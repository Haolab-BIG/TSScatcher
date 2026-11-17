# Output layout

All project-generated files are placed below the directory supplied with `--output`. Each sample is processed as one library type (`long` or `short`).

```text
<output>/
├── 1_run_information/
│   ├── samples.tsv
│   └── merge_manifest.tsv
├── 2_cellranger/
│   └── <sample_id>/outs/
├── 3_preprocessed_bam/
│   └── <sample_id>/
│       ├── <sample_id>.filtered_reads.bam
│       └── <sample_id>.filtered_reads.bam.bai
├── 4_sample_tss/
│   └── <sample_id>/
│       ├── 1_read_sites/
│       │   ├── <sample_id>.read_sites.bed.gz
│       │   ├── <sample_id>.read_sites.before_shift.bed.gz  # short only
│       │   ├── <sample_id>.read_sites.after_shift.bed.gz   # short only
│       │   ├── <sample_id>.shift_distance.txt              # short only
│       │   └── <sample_id>.shift_distance.pdf              # short only
│       ├── 2_coverage_tracks/
│       │   ├── <sample_id>.read_coverage.plus.bw            # long
│       │   ├── <sample_id>.read_coverage.minus.bw           # long
│       │   └── <sample_id>.read_coverage.{before_shift,after_shift}.{plus,minus}.bw  # short
│       ├── 3_peak_calling/
│       │   ├── <sample_id>.paraclu_raw_summits.bed
│       │   └── <sample_id>.tso_filtered_peaks.before_dl.bed
│       └── 4_final_sample_tss/
│           ├── <sample_id>.final_sample_tss.bed
│           └── qc/<sample_id>.dl_filter_status.tsv.gz
├── 5_merged_tss/
│   ├── final_merged_tss.bed.gz
│   └── final_merged_tss.supporting_samples.tsv.gz
├── 6_quantification/
│   ├── per_sample/<sample_id>/
│   │   ├── matrix.mtx.gz
│   │   ├── barcodes.tsv.gz
│   │   └── features.tsv.gz
│   └── all_samples/
│       ├── matrix.mtx.gz
│       ├── barcodes.tsv.gz
│       └── features.tsv.gz
└── tmp/TSScatcher.XXXXXX/   # removed after success; retained after failure
```

`final_sample_tss.bed` is the final BED6 for one sample after DL filtering and reference rescue. `final_merged_tss.bed.gz` is the cross-sample BED6 used for quantification. Its score column is the number of supporting samples. `final_merged_tss.supporting_samples.tsv.gz` records their exact sample IDs.

The QC table contains `peak_id`, DL `probability`, and `status` (`dl_kept`, `rescued`, or `dl_rejected`). It is not used as the final TSS input.

Minus-strand bigWigs contain negative coverage for mirrored browser display. Read names and cell barcodes are retained in the seven-column read-site BED files. Temporary bedGraph files are deleted only after a non-empty bigWig is created.

When `--deep-learning no` is used, `tso_filtered_peaks.before_dl.bed` enters multi-sample merging directly, so no `final_sample_tss.bed` or DL status table is created.
