# cycloheximide

Code used to generate the files processed by `Dandelion` preprocessing pipeline.

```bash
# BCR
# cd /path/to/dandelion-demo-files/dandelion_manuscript/data/cycloheximide
singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess \
	--meta="cycloheximideBCR_metadata.csv" --sep="-" --file_prefix="all" \
	--filter_to_high_confidence

# TCR
singularity run -B $PWD /path/to/sc-dandelion_latest.sif dandelion-preprocess \
	--meta="cycloheximideTCR_metadata.csv" --sep="-" --file_prefix="all" --chain="TR" \
	--filter_to_high_confidence
```
