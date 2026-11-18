# TSScatcher

TSScatcher calls and quantifies strand-aware transcription start sites (TSSs) from paired-end 5′ single-cell RNA-seq. In these libraries, the cell barcode and unique molecular identifier (UMI) are carried by R1, but the usable cDNA sequence in R1 depends on the library design. TSScatcher supports two read structures:

- `long`: R1 extends beyond the cell barcode and UMI into cDNA and is long enough to be aligned. TSScatcher uses the R1 alignment to identify transcript 5-prime ends. In the 10x Genomics 5′ workflow, this corresponds to equal-length R1/R2 data processed with Cell Ranger chemistry `SC5P-PE`.
- `short`: R1 contains the cell barcode and UMI but little or no alignable cDNA sequence (commonly 26 bp), while R2 contains the cDNA insert. TSScatcher obtains the signal from R2 and applies transcriptome-based shift estimation to recover the TSS coordinate. In the 10x Genomics 5′ workflow, this corresponds to Cell Ranger chemistry `SC5P-R2`.

The `long` and `short` labels therefore describe the R1 read structure, not different biological classes of TSS.

This README is arranged as a start-to-finish tutorial. Complete Steps 1-5 in
order for the first run. The parameter tables and algorithm notes are reference
material and are not required to start.

## Step 1: install TSScatcher and its dependencies

Linux is required. Clone the repository, enter it, and create the supplied micromamba environment:

```bash
git clone https://github.com/yswutan/TSScatcher.git
cd TSScatcher
micromamba create -f env/environment.yml
micromamba activate tsscatcher
```

Cell Ranger is distributed separately by 10x Genomics and is not installed by micromamba. TSScatcher was validated with Cell Ranger 7.2.0. Install it as follows:

1. Open the 10x Genomics Cell Ranger download page, accept the license, and download the Linux tarball.
2. Extract it into a user-controlled software directory:

   ```bash
   mkdir -p /path/to/software
   tar -xzf cellranger-7.2.0.tar.gz -C /path/to/software
   ```

3. Locate the executable and verify it:

   ```bash
   find /path/to/software/cellranger-7.2.0 -type f -name cellranger -perm -111
   /path/to/software/cellranger-7.2.0/cellranger --version
   ```

4. Either add its directory to `PATH`:

   ```bash
   export PATH=/path/to/software/cellranger-7.2.0:$PATH
   ```

   or provide the full executable path in every `tsscatcher run` and `tsscatcher mkref` command:

   ```bash
   --cellranger /path/to/software/cellranger-7.2.0/cellranger
   ```

The runtime never installs packages. If a dependency is absent, the command stops
and reports it before the affected step.

Add TSScatcher to `PATH`:

```bash
export PATH=/path/to/TSScatcher/bin:$PATH
```

## Step 2: build the reference

Long- and short-R1 libraries require different reference bundles. Long R1 uses only the genome, reference TSS and Cell Ranger index:

```bash
tsscatcher mkref \
  --type long \
  --fasta /reference/genome.fa \
  --gtf /reference/genes.gtf \
  --genome-name GRCh38 \
  --threads 20 \
  --output /reference/tsscatcher_GRCh38_long
```

Short R1 additionally builds the extended GTF and RSEM/STAR reference:

```bash
tsscatcher mkref \
  --type short \
  --fasta /reference/genome.fa \
  --gtf /reference/genes.gtf \
  --genome-name GRCh38 \
  --extend-bp 300 \
  --threads 20 \
  --output /reference/tsscatcher_GRCh38_short
```

Use the short reference bundle for a mixed long/short project because it contains the shared files plus the short-R1-only files. Genome-scale indexes are not stored in this repository.

### What does `mkref` create?

`mkref` creates a new, self-contained directory at `--output`; it does not modify the input FASTA or GTF. Keep the entire directory together. Pass the top-level directory, not its nested Cell Ranger directory, to `run --reference`.

```text
<long-reference>/
├── reference_type.txt                  # contains the word: long
├── genome/
│   ├── genome.fa                        # copy of --fasta
│   ├── genome.fa.fai                    # samtools FASTA index
│   ├── genes.gtf                        # copy of --gtf
│   ├── chrom.sizes                     # chromosome and length
│   └── reference_tss.bed               # annotation TSS, BED6
└── cellranger/<genome-name>/
    ├── fasta/
    ├── genes/
    ├── star/
    └── reference.json

<short-reference>/
├── reference_type.txt                  # contains the word: short
├── genome/
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── genes.gtf
│   ├── chrom.sizes
│   ├── reference_tss.bed
│   └── genes.extend<extend-bp>.gtf    # e.g. genes.extend300.gtf
├── cellranger/<genome-name>/            # same structure as above
└── rsem/                               # RSEM and STAR indexes
    ├── reference.grp
    ├── reference.ti
    ├── reference.seq
    ├── reference.transcripts.fa
    ├── reference.idx.fa
    ├── reference.n2g.idx.fa
    ├── Genome
    ├── SA
    ├── SAindex
    └── other STAR index tables
```

| Output | Purpose |
|---|---|
| `reference_type.txt` | Records whether `mkref` built a `long` or `short` bundle; it is also a simple human-readable check of the bundle type. |
| `genome/genome.fa` and `genes.gtf` | Copies of the supplied inputs, kept so one bundle contains the exact genome and annotation used to build it. |
| `genome/genome.fa.fai` | FASTA random-access index made by `samtools faidx`. |
| `genome/chrom.sizes` | Two columns, chromosome and length; used for genomic boundaries and bigWig creation. |
| `genome/reference_tss.bed` | Deduplicated one-base TSS sites derived from GTF `transcript` records. BED6 columns are chromosome, 0-based start, end, generated site ID, `.`, and strand. It is used for reference rescue after model filtering. |
| `genome/genes.extend<extend-bp>.gtf` | Short only. Extends the transcript 5-prime boundary upstream on `+` and downstream on `-`; used for transcript-to-genome conversion during shift estimation. |
| `cellranger/<genome-name>/` | Complete Cell Ranger reference built by `cellranger mkref`; used for read alignment. Do not edit its internal files. |
| `rsem/reference.*` | Short only. RSEM transcript reference and metadata. |
| `rsem/Genome`, `SA`, `SAindex`, and related tables | Short only. STAR index used for transcript assignment and library-specific shift estimation. |

A long bundle can run long samples only. A short bundle contains all shared files and can run long-only, short-only, or mixed projects. If `mkref --extend-bp` is not 300, use the same value in `run --extend-bp`.

### How to confirm that `mkref` finished correctly

The final terminal line should report `long reference bundle created at ...` or `short reference bundle created at ...`. Then check the important files:

```bash
REF=/reference/tsscatcher_GRCh38_short   # change to your output directory
cat "$REF/reference_type.txt"
test -s "$REF/genome/reference_tss.bed" && echo "reference TSS: OK"
test -d "$REF/cellranger/GRCh38/star" && echo "Cell Ranger index: OK"
test -s "$REF/rsem/SA" && echo "short-R1 RSEM/STAR index: OK"
```

For a long bundle, omit the last RSEM/STAR check. The `rsem/` directory and extended GTF are intentionally absent.

**Expected time and disk space.** Reference construction is a one-time, genome-scale operation. In the hg38 validation on this server with 20 threads and network storage, the long bundle took approximately 1 hour 52 minutes and occupied 19 GB; the short bundle took approximately 2 hours 48 minutes and occupied 48 GB. As a practical estimate, reserve 2-4 hours and at least 25 GB for long or 60 GB for short. CPU speed, available memory, filesystem performance, FASTA/GTF size, and Cell Ranger version can change these times substantially. Cell Ranger indexing dominates long builds; short builds additionally run GTF extension and RSEM/STAR indexing.

### TSO-like sites and strand-invasion filtering

The template-switching oligo (TSO) is required by 10x 5-prime chemistry, but a TSO-derived sequence can sometimes anneal to a complementary internal sequence instead of marking a genuine transcript 5-prime end. This **TSO strand invasion** can produce an artificial pile-up of read starts that resembles a TSS. TSScatcher therefore uses a strand-aware BED of genomic TSO-like sequences to remove likely strand-invasion artifacts during peak calling.

TSScatcher provides `ref/TSOlike/TSO_mismatch_2_sort.bed` for hg38 primary chromosomes (`chr1`-`chr22`, `chrX`, and `chrY`). Columns 1-6 are BED6 and are used by TSScatcher; columns 7 and 8 contain the mismatch counts in the first 10 nt and final 3 nt.

The TSO-like BED is genome-build specific and is not generated by `mkref`. Do not use the bundled hg38 BED with hg19, a different chromosome naming convention, or another species. For another genome, build a Bowtie2 index from exactly the same FASTA supplied to `mkref`, and then run the optional reproduction script. Bowtie2 is needed only for this one-time auxiliary procedure and is deliberately not installed by the TSScatcher micromamba environment.

Install Bowtie2 separately, make an index in a user-selected reference directory, and verify the index files:

```bash
# The exact installation command depends on the user's software setup.
bowtie2 --version

mkdir -p /reference/custom/bowtie2
bowtie2-build --threads 20 \
  /reference/custom/genome.fa \
  /reference/custom/bowtie2/genome

ls /reference/custom/bowtie2/genome*.bt2*
```

Then use the standalone helper script:

```bash
scripts/build_tsolike_bed.sh \
  --fasta /reference/custom/genome.fa \
  --bowtie2-index /reference/custom/bowtie2/genome \
  --threads 20 \
  --work-dir /reference/custom/tso_work \
  --output /reference/custom/TSO_mismatch_2_sort.bed
```

The helper reproduces the original construction method for the TSO-derived sequence `TTTCTTATATGGG`:

1. generate substitution-only candidates with up to three mismatches in the first 10 nt and up to one mismatch in the final 3 nt (36,760 unique sequences);
2. align them end-to-end with Bowtie2, retaining up to 999 placements per sequence;
3. extract each strand-oriented genomic 13-mer from the target FASTA;
4. retain sites with no more than two mismatches in the first 10 nt and no more than one mismatch in the final 3 nt;
5. deduplicate and coordinate-sort the BED. The final first-part cutoff gives the output its `mismatch_2` name.

This script is a reference-building utility, not part of `tsscatcher mkref` or the normal analysis workflow. It does not require additional Python packages: candidate generation uses the Python standard library and mismatch counting uses `awk`. It does require external `bowtie2`, plus `samtools` and `bedtools` already present in the TSScatcher environment. Its run-specific files stay under `--work-dir` and are removed after success; use `--keep-work` to retain them. A failed build retains its work directory for diagnosis.

Use the generated BED in the analysis:

```bash
tsscatcher run \
  ... \
  --reference /reference/custom_short_or_long \
  --tsolike-bed /reference/custom/TSO_mismatch_2_sort.bed \
  --output /analysis/custom_project
```


## Step 3: prepare the FASTQs and sample sheet

The tab-separated file has exactly these first four columns:

```text
sample_id    fastq_dir       fastq_sample    type
sample01     /data/run_01    sample01        long
sample02     /data/run_02    sample02        short
```

Each sample occurs once and has one type. A sample is never counted twice toward
the replicate threshold.

## Step 4: run TSScatcher

```bash
tsscatcher run \
  --samples samples.tsv \
  --reference /reference/tsscatcher_GRCh38_short \
  --output /analysis/project01 \
  --threads 20 \
  --memory 128 \
  --deep-learning yes \
  --dl-threshold 0.5 \
  --rescue-distance 100 \
  --min-replicates 2 \
  --merge-distance 100 \
  --shift-max-distance 1000
```

Single-sample direct input is also supported:

```bash
tsscatcher run \
  --id sample01 \
  --fastqs /data/run_01 \
  --sample sample01 \
  --type long \
  --reference /reference/tsscatcher_GRCh38_long \
  --output /analysis/sample01
```

Before Cell Ranger runs, TSScatcher samples reads from both R1 and R2. Equal
median lengths are expected for `long`; unequal lengths with R1 shorter than R2
and no more than 35 bp are expected for `short`. A mismatch stops the run. After
manual inspection, an intentional mismatch can be accepted with:

```bash
--confirm-type-mismatch
```

A successful run ends with `TSScatcher run completed: <output>`. Confirm the two principal results before continuing:

```bash
OUT=/analysis/project01   # change to your --output directory
test -s "$OUT/5_merged_tss/final_merged_tss.bed.gz" && echo "merged TSS: OK"
test -s "$OUT/6_quantification/all_samples/matrix.mtx.gz" && echo "count matrix: OK"
```

If the run stops, read the last error first. Temporary files are created only under `<output>/tmp/`; the run-specific temporary directory is removed after a successful run.

## Step 5: inspect the results

All files are written below the user-supplied `--output` directory.

```text
<output>/
├── 1_run_information/
│   ├── samples.tsv
│   └── merge_manifest.tsv
├── 2_cellranger/<sample_id>/outs/
├── 3_preprocessed_bam/<sample_id>/
│   ├── <sample_id>.filtered_reads.bam
│   └── <sample_id>.filtered_reads.bam.bai
├── 4_sample_tss/<sample_id>/
│   ├── 1_read_sites/
│   │   ├── <sample_id>.read_sites.bed.gz
│   │   ├── <sample_id>.read_sites.before_shift.bed.gz  # short only
│   │   ├── <sample_id>.read_sites.after_shift.bed.gz   # short only
│   │   ├── <sample_id>.shift_distance.txt              # short only
│   │   └── <sample_id>.shift_distance.pdf              # short only
│   ├── 2_coverage_tracks/
│   │   ├── <sample_id>.read_coverage.plus.bw            # long
│   │   ├── <sample_id>.read_coverage.minus.bw           # long
│   │   └── <sample_id>.read_coverage.{before_shift,after_shift}.{plus,minus}.bw  # short
│   ├── 3_peak_calling/
│   │   ├── <sample_id>.paraclu_raw_summits.bed
│   │   └── <sample_id>.tso_filtered_peaks.before_dl.bed
│   └── 4_final_sample_tss/
│       ├── <sample_id>.final_sample_tss.bed
│       └── qc/<sample_id>.dl_filter_status.tsv.gz
├── 5_merged_tss/
│   ├── final_merged_tss.bed.gz
│   └── final_merged_tss.supporting_samples.tsv.gz
├── 6_quantification/
│   ├── per_sample/<sample_id>/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
│   └── all_samples/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
└── tmp/TSScatcher.XXXXXX/  # deleted after success; retained after failure
```

### Text-file columns

`1_run_information/samples.tsv`:

| Column | Meaning |
|---|---|
| `sample_id` | Unique output identifier |
| `fastq_dir` | Input FASTQ directory |
| `fastq_sample` | FASTQ filename sample prefix passed to Cell Ranger |
| `type` | `long` or `short` |

`1_run_information/merge_manifest.tsv` has two columns: `sample_id` and the corresponding sample-level `peak_file` used for cross-sample merging.

All `*.read_sites*.bed.gz` files contain seven tab-separated columns:

| Column | Meaning |
|---:|---|
| 1 | Chromosome |
| 2 | 0-based start coordinate |
| 3 | End coordinate |
| 4 | Read name |
| 5 | BED score (the input alignment score, often `.`) |
| 6 | TSS strand (`+` or `-`) |
| 7 | Cell barcode |

`*.paraclu_raw_summits.bed` contains eight columns: chromosome, 0-based start, end, peak ID, total peak count, strand, maximum summit count, and summit coordinate(s). `*.tso_filtered_peaks.before_dl.bed` and `*.final_sample_tss.bed` are BED6: chromosome, 0-based start, end, peak ID, total peak count, and strand.

`final_merged_tss.bed.gz` is BED6: chromosome, 0-based start, end, final TSS ID, number of supporting samples, and strand. `final_merged_tss.supporting_samples.tsv.gz` contains `final_tss_id`, `support_samples` (count), and `sample_ids` (comma-separated identifiers).

`*.dl_filter_status.tsv.gz` contains `peak_id`, model `probability`, and `status`: `dl_kept`, `rescued`, or `dl_rejected`.

`features.tsv.gz` contains final TSS ID, genomic label (`chrom:start-end:strand`), and feature type (`TSS`). `barcodes.tsv.gz` has one barcode per row; the combined matrix prefixes it with `sample_id#`. In `matrix.mtx.gz`, the dimension line is number of TSS features, number of barcodes, and number of nonzero entries; each subsequent row is feature index, barcode index, and integer count.

`shift_distance.txt` contains one unrounded median shift value. BAM/BAI, bigWig, PDF, and Cell Ranger files are binary or tool-standard formats and do not have TSScatcher-defined text columns. Minus-strand bigWigs store negative coverage for mirrored browser display.

See [docs/OUTPUTS.md](docs/OUTPUTS.md) for the same layout as a standalone reference.

## Parameter reference (read after the first successful run)

### `tsscatcher mkref`

| Parameter | Required/default | Meaning |
|---|---|---|
| `--type` | required: `long` or `short` | Selects the reference bundle. `long` builds genome/TSS/Cell Ranger files only; `short` additionally builds the extended GTF and RSEM/STAR reference. |
| `--fasta` | required | Input genomic FASTA. Sequence names must match the GTF chromosome field. |
| `--gtf` | required | Input gene annotation GTF used for Cell Ranger, reference TSS extraction, and short-R1 transcript conversion. |
| `--output` | required | New reference-bundle directory. |
| `--genome-name` | `TSScatcher` | Name passed to `cellranger mkref`; it becomes the Cell Ranger reference subdirectory name. |
| `--threads` | `10` | Threads for RSEM/STAR preparation and Cell Ranger reference construction. |
| `--extend-bp` | `300` | Short reference only: bases added upstream of transcript 5-prime ends in the extended GTF. Use the same value for `run --extend-bp`. |
| `--cellranger` | `cellranger` from `PATH` | Full Cell Ranger executable path when it is not on `PATH`. |

### `tsscatcher run`

Input can be supplied either with `--samples` or with all four direct-input parameters `--id`, `--fastqs`, `--sample`, and `--type`. Do not mix the two modes.

| Parameter | Required/default | Meaning |
|---|---|---|
| `--samples` | conditional | Tab-separated multi-sample sheet whose first four columns are exactly `sample_id`, `fastq_dir`, `fastq_sample`, and `type`. |
| `--id` | direct mode | Unique TSScatcher output identifier. This need not equal the FASTQ prefix. |
| `--fastqs` | direct mode | Directory containing paired R1/R2 FASTQ files. |
| `--sample` | direct mode | Prefix used to locate FASTQs and passed to Cell Ranger `--sample`; for `test_long_S1_L001_R1_001.fastq.gz`, use `test_long`. |
| `--type` | direct mode | `long` for equal-length R1/R2 or `short` for short R1 and longer R2. Determines Cell Ranger chemistry and downstream algorithm. |
| `--reference` | required | TSScatcher reference bundle. A project containing any short sample requires a short bundle. |
| `--output` | required | Project output directory. All results and temporary data are created below it. Use a new directory for a new run. |
| `--threads` | `10` | Local CPU cores passed to Cell Ranger and used by downstream tools. |
| `--memory` | `64` | Cell Ranger local memory limit in GB. |
| `--cellranger` | `cellranger` from `PATH` | Full path to the Cell Ranger executable. |
| `--tsolike-bed` | bundled hg38 BED | Custom genome-specific TSO-like BED. Required when the reference genome is not hg38. |
| `--confirm-type-mismatch` | off | Continue despite disagreement between declared type and sampled R1/R2 median lengths. Use only after inspecting the reported lengths. |
| `--peak-count` | `5` | Paraclu minimum total count used to call candidate clusters. |
| `--peak-summit` | `3` | Minimum maximum per-position read count required at a candidate summit. |
| `--deep-learning` | `yes` | `yes` applies sequence-model filtering plus reference rescue; `no` merges the TSO-filtered pre-DL peaks directly. |
| `--dl-threshold` | `0.5` | Model probability must be strictly greater than this value to be DL-kept. |
| `--rescue-distance` | `100` | A DL-rejected peak is restored when its same-strand distance to a reference TSS is between 0 and this many bases, inclusive. |
| `--extend-bp` | `300` | Short R1 only. Must match the extended GTF generated by `mkref`; participates in genomic shift conversion. |
| `--shift-max-distance` | `1000` | Short R1 only. Only transcript-assignment distances below this value enter the median shift calculation. The median is not rounded. |
| `--min-replicates` | `2` | Minimum number of distinct samples simultaneously supporting a strand-aware interval. Automatically becomes 1 for a one-sample run. |
| `--merge-distance` | `100` | After sample-support filtering, adjacent retained partitions on the same chromosome and strand are joined when their gap is at most this value. |

Before mapping, `run` examines up to 10,000 reads per mate. `long` requires equal R1/R2 median lengths. `short` requires R1 shorter than R2 and R1 median length no greater than 35 bp. These two internal defaults are exposed by `check-fastq`/`map`, not by `run`.

### Advanced subcommands

Most users should use `run`; the commands below expose individual stages.

#### `tsscatcher check-fastq`

| Parameter | Required/default | Meaning |
|---|---|---|
| `--fastqs` | required | FASTQ directory. |
| `--sample` | required | FASTQ filename prefix. |
| `--type` | required | Declared `long` or `short` type. |
| `--reads` | `10000` | Maximum reads sampled from each mate for median-length inference. |
| `--short-r1-max-length` | `35` | Largest R1 median accepted as `short`. |
| `--confirm-type-mismatch` | off | Report the mismatch but return success. |

#### `tsscatcher map`

| Parameter | Required/default | Meaning |
|---|---|---|
| `--id` | required | Cell Ranger run/output ID. |
| `--fastqs`, `--sample`, `--type` | required | Same meanings as in direct `run` mode. |
| `--transcriptome` | required | Cell Ranger reference directory itself, not the TSScatcher bundle parent. |
| `--output` | required | Parent directory for the Cell Ranger run. |
| `--threads` | `10` | Cell Ranger `--localcores`. |
| `--memory` | `64` | Cell Ranger `--localmem` in GB. |
| `--cellranger` | `cellranger` | Cell Ranger executable. |
| `--type-check-reads` | `10000` | Reads per mate checked before mapping. |
| `--short-r1-max-length` | `35` | Short-R1 cutoff. |
| `--confirm-type-mismatch` | off | Allow an intentional length/type mismatch. |

#### `tsscatcher preprocess`

| Parameter | Required/default | Meaning |
|---|---|---|
| `--bam` | required | Cell Ranger position-sorted BAM. |
| `--id` | required | Sample output ID. |
| `--barcodes` | required | Filtered Cell Ranger barcode list (`barcodes.tsv` or `.gz`). |
| `--type` | required | `long` or `short`; selects the read flag and long-R1 TSO preprocessing. |
| `--output` | required | Directory for filtered BAM/BAI. |
| `--threads` | `5` | Samtools threads. |
| `--tmp-dir` | `<output>/tmp` | Temporary working directory. `run` overrides this to its unique project-local temp directory. |

#### `tsscatcher call-peaks`

Common required parameters are `--type`, `--bam`, `--id`, `--output`, `--chrom-sizes`, `--tsolike`, `--fasta`, and `--reference-tss`. `--threads` defaults to 5; `--peak-count`, `--peak-summit`, `--deep-learning`, `--dl-threshold`, and `--rescue-distance` have the same meanings/defaults as `run`; `--tmp-dir` defaults to `<output>/tmp`. Short type additionally requires `--rsem-ref` and `--extended-gtf`, and accepts `--extend-bp 300` and `--shift-max-distance 1000`.

#### `tsscatcher merge`

| Parameter | Required/default | Meaning |
|---|---|---|
| `--manifest` | required | TSV with header `sample_id` and `peak_file`; each sample ID may occur once. |
| `--output-bed` | required | Gzipped final BED6 path. |
| `--output-support` | required | Gzipped support-detail TSV path. |
| `--min-replicates` | `2` | Distinct-sample support threshold; automatically 1 for one manifest sample. |
| `--merge-distance` | `100` | Maximum gap used to join retained partitions. |

#### `tsscatcher quantify`

| Parameter | Required/default | Meaning |
|---|---|---|
| `--id` | required | Sample ID used in status messages. |
| `--reads` | required | Seven-column read-site BED or BED.gz with cell barcode in column 7. |
| `--final-sites` | required | Final BED6/BED6.gz whose feature ID is column 4. |
| `--output` | required | Output matrix directory. |
| `--counts-output` | optional | Also copy the intermediate three-column feature/barcode/count table here. |
| `--tmp-dir` | `<output>/tmp` | Temporary working directory. |


## Deep-learning filtering is optional

The packaged model is used only when `--deep-learning yes` (default). Set `--deep-learning no` to skip model scoring and reference rescue entirely; the strand-aware TSO-filtered peaks in `*.tso_filtered_peaks.before_dl.bed` then enter multi-sample merging directly. No `final_sample_tss.bed` or `dl_filter_status.tsv.gz` is produced in that mode.

## Algorithm notes

The exact confirmed workflow is documented in
[docs/ALGORITHM.md](docs/ALGORITHM.md). Every genomic comparison, rescue,
replicate calculation, merge and quantification is strand-aware.

## Test with the included real FASTQs

The repository includes anonymized real long-R1 and short-R1 FASTQ subsets under `tests/data/`. Follow the installation, reference preparation, and test commands below; both commands execute the complete workflow from Cell Ranger through final TSS quantification.

On the validation server using 10 threads, the 100,000-read-pair long-R1 test took approximately 4 minutes and the short-R1 test approximately 7 minutes. Short R1 takes longer because it additionally performs RSEM/STAR assignment and shift estimation. Runtime depends on CPU, storage speed, Cell Ranger version and system load.

```bash
# Long R1
tsscatcher run \
  --id test_long \
  --fastqs "$PWD/tests/data/long" \
  --sample test_long \
  --type long \
  --reference /reference/tsscatcher_GRCh38_long \
  --output "$PWD/test_output/long" \
  --threads 10 \
  --memory 64 \
  --cellranger /path/to/cellranger

# Short R1
tsscatcher run \
  --id test_short \
  --fastqs "$PWD/tests/data/short" \
  --sample test_short \
  --type short \
  --reference /reference/tsscatcher_GRCh38_short \
  --output "$PWD/test_output/short" \
  --threads 10 \
  --memory 64 \
  --cellranger /path/to/cellranger
```

The included subsets contain 100,000 paired reads each and are intended to verify installation and file flow, not to reproduce biological results from the complete library.

## License

TSScatcher is released under the GNU General Public License v3.0 (GPL-3.0); see [LICENSE](LICENSE). Third-party software, including Cell Ranger, remains subject to its own license and is not redistributed by TSScatcher.
