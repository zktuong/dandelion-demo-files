# dandelion-remap

Code used to generate the files processed by `Dandelion` preprocessing pipeline.

## abTCR
```bash
cd /path/to/dandelion-demo-files/dandelion_manuscript/data/dandelion-remap
singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess \
	--meta="abTCR_metadata.csv" --file_prefix="all" --chain="TR" \
	--keep_trailing_hyphen_number --skip_format_header --filter_to_high_confidence
```

## BCR 
```bash
singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess \
	--meta="BCR_metadata.csv" --sep="-" --file_prefix="all" \
	--keep_trailing_hyphen_number --filter_to_high_confidence
```

Please note that for the BCR folders, the GEX barcode was appended as a prefix to the `sequence_id` and `cell_id` whereas the TCR folders have their original barcodes (due to `--skip_format_header`).

## gdTCR

### original

This refers to output from cellranger 6.1.2 ran with the 10X GRCh38 5.0.0 V(D)J reference, with the contigs identified by cellranger as high confidence subsequently re-annotated with `Dandelion`.
```bash
cd /path/to/dandelion-demo-files/dandelion_manuscript/data/dandelion-remap/gdT/original
singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess \
	--meta="gdTCR_metadata.csv" --sep="-" --file_prefix="all" --chain="TR" \
	--keep_trailing_hyphen_number --skip_format_header --filter_to_high_confidence
```

### gdhack

This refers to output from cellranger 6.1.2 ran with 5.0.0 reference modified to obtain annotated γδTCR contigs as per 10X Genomics' instructions.

## Failed samples

Please note that there are 4 `abTCR` samples (in `abTCR_metadata_failed.csv`) where the libraries failed, producing little to no contigs. These were not further processed and ignored for the downstream analyses.
