#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: extend_gtf.sh <input.gtf> <output.gtf> <extend_bp>" >&2
  exit 2
fi

input_gtf=$1
output_gtf=$2
extend_bp=$3

[[ -f "$input_gtf" ]] || { echo "input GTF not found: $input_gtf" >&2; exit 2; }
[[ "$extend_bp" =~ ^[0-9]+$ ]] || { echo "extend_bp must be a non-negative integer" >&2; exit 2; }

output_dir=$(dirname -- "$output_gtf")
mkdir -p "$output_dir"
temp_output="${output_gtf}.tmp.$$"
trap 'rm -f "$temp_output"' EXIT

# Pass 1 records every transcript 5-prime boundary. Pass 2 extends every
# feature from the same chromosome, strand and gene that shares that original
# boundary. Each input record is evaluated once, preventing cascading extension
# when different isoform TSSs are separated by exactly extend_bp.
awk -v FS='\t' -v OFS='\t' -v extend_bp="$extend_bp" '
function gene_id_from_attributes(attributes, value) {
  if (match(attributes, /gene_id[[:space:]]+"[^"]+"/)) {
    value = substr(attributes, RSTART, RLENGTH)
    sub(/^gene_id[[:space:]]+"/, "", value)
    sub(/"$/, "", value)
    return value
  }
  return ""
}

NR == FNR {
  if ($0 ~ /^#/ || NF == 0) next
  if (NF != 9) {
    print "invalid GTF record with " NF " columns at input line " FNR > "/dev/stderr"
    exit 2
  }
  if ($3 == "transcript" && ($7 == "+" || $7 == "-")) {
    gene_id = gene_id_from_attributes($9)
    if (gene_id == "") {
      print "transcript lacks gene_id at input line " FNR > "/dev/stderr"
      exit 2
    }
    tss = ($7 == "+" ? $4 : $5)
    transcript_tss[$1 SUBSEP $7 SUBSEP gene_id SUBSEP tss] = 1
  }
  next
}

$0 ~ /^#/ || NF == 0 {
  print
  next
}

NF != 9 {
  print "invalid GTF record with " NF " columns at input line " FNR > "/dev/stderr"
  exit 2
}

{
  gene_id = gene_id_from_attributes($9)
  if ($7 == "+" && (($1 SUBSEP $7 SUBSEP gene_id SUBSEP $4) in transcript_tss)) {
    new_start = $4 - extend_bp
    $4 = (new_start < 1 ? 1 : new_start)
  } else if ($7 == "-" && (($1 SUBSEP $7 SUBSEP gene_id SUBSEP $5) in transcript_tss)) {
    $5 += extend_bp
  }
  print
}
' "$input_gtf" "$input_gtf" > "$temp_output"

mv -f "$temp_output" "$output_gtf"
trap - EXIT
