#!/usr/bin/env bash
# Rebuild a genome-specific TSO-like BED. This is an optional reference utility,
# not part of the normal tsscatcher mkref or run workflow.
set -euo pipefail

fasta= bowtie2_index= output= work_root= threads=10 keep_work=no
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) fasta=$2; shift 2;;
    --bowtie2-index) bowtie2_index=$2; shift 2;;
    --output) output=$2; shift 2;;
    --work-dir) work_root=$2; shift 2;;
    --threads) threads=$2; shift 2;;
    --keep-work) keep_work=yes; shift;;
    --help)
      echo "Usage: $0 --fasta genome.fa --bowtie2-index INDEX_PREFIX --output TSO_mismatch_2_sort.bed [--threads 10] [--work-dir DIR] [--keep-work]"
      exit 0;;
    *) echo "ERROR: unknown option: $1" >&2; exit 1;;
  esac
done
[[ -n "$fasta" && -n "$bowtie2_index" && -n "$output" ]] || {
  echo "ERROR: --fasta, --bowtie2-index and --output are required" >&2
  exit 1
}
[[ -s "$fasta" ]] || { echo "ERROR: FASTA not found or empty: $fasta" >&2; exit 1; }
[[ -e "${bowtie2_index}.1.bt2" || -e "${bowtie2_index}.1.bt2l" ]] || {
  echo "ERROR: Bowtie2 index prefix not found: $bowtie2_index" >&2
  exit 1
}
for command in python3 bowtie2 samtools bedtools sort awk; do
  command -v "$command" >/dev/null || { echo "ERROR: missing command: $command" >&2; exit 1; }
done

fasta=$(realpath "$fasta")
output=$(realpath -m "$output")
mkdir -p "$(dirname "$output")"
[[ -n "$work_root" ]] || work_root="$(dirname "$output")/tso_build_work"
mkdir -p "$work_root"
work=$(mktemp -d "$work_root/TSOlike.XXXXXX")
if [[ "$keep_work" == yes ]]; then
  trap 'echo "Work directory retained: $work" >&2' EXIT
else
  trap 'status=$?; if [[ $status -eq 0 ]]; then rm -rf -- "$work"; rmdir "$work_root" 2>/dev/null || true; else echo "Build failed; work directory retained: $work" >&2; fi' EXIT
fi

# Generate the original 13-nt TSO-derived candidates. Only substitutions are
# allowed: <=3 in the first 10 nt and <=1 in the final 3 nt.
python3 - "$work/TSO.similar.fa" "$work/TSO.similar.fastq" <<'PY'
import itertools
import sys

sequence = "TTTCTTATATGGG"
bases = "ATCG"

def variants(seq, maximum):
    answer = set()
    for distance in range(maximum + 1):
        for positions in itertools.combinations(range(len(seq)), distance):
            choices = [[base for base in bases if base != seq[pos]] for pos in positions]
            for replacements in itertools.product(*choices):
                value = list(seq)
                for pos, base in zip(positions, replacements):
                    value[pos] = base
                answer.add("".join(value))
    return answer

full = sorted(a + b for a in variants(sequence[:10], 3) for b in variants(sequence[10:], 1))
with open(sys.argv[1], "w") as fasta, open(sys.argv[2], "w") as fastq:
    for index, seq in enumerate(full):
        name = f"variant_{index}"
        fasta.write(f">{name}\n{seq}\n")
        fastq.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
print(f"Generated {len(full)} unique candidate sequences", file=sys.stderr)
PY

[[ -s "${fasta}.fai" ]] || samtools faidx "$fasta"
cut -f1 "${fasta}.fai" > "$work/fasta.chromosomes.txt"

bowtie2 -k 999 -p "$threads" -L 13 --end-to-end --score-min C,0,0 \
  --no-1mm-upfront -x "$bowtie2_index" -U "$work/TSO.similar.fastq" \
  2> "$work/bowtie2.log" \
  | samtools view -b -o "$work/TSO.similar_raw.bam" -

bedtools bamtobed -i "$work/TSO.similar_raw.bam" \
  | awk 'BEGIN{FS=OFS="\t"} NR==FNR{keep[$1]=1;next} $1 in keep' \
      "$work/fasta.chromosomes.txt" - \
  > "$work/TSO.similar.bed"

# bedtools returns each 13-mer in the alignment strand. Keep genomic sites with
# <=2 mismatches in the first 10 nt and <=1 in the final 3 nt. BED columns 7-8
# record those two mismatch counts.
bedtools getfasta -fi "$fasta" -bed "$work/TSO.similar.bed" -s -bedOut \
  | awk 'BEGIN{FS=OFS="\t"; target="TTTCTTATATGGG"}
    {
      seq=toupper($7); if(length(seq)!=13) next
      first=0; last=0
      for(i=1;i<=10;i++) if(substr(seq,i,1)!=substr(target,i,1)) first++
      for(i=11;i<=13;i++) if(substr(seq,i,1)!=substr(target,i,1)) last++
      if(first<=2 && last<=1) print $1,$2,$3,".",".",$6,first,last
    }' \
  | LC_ALL=C sort -u -k1,1V -k2,2n -k3,3n -k6,6 \
  > "$output"

[[ -s "$output" ]] || { echo "ERROR: no TSO-like sites were generated" >&2; exit 1; }
echo "TSO-like BED created: $output ($(wc -l < "$output") unique sites)" >&2
