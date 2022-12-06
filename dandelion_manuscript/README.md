# dandelion_manuscript

This folder contains all the notebooks and relevant files required to generata analysis results and plots for the `Dandelion` manuscript.

## Contents

- `notebook`: contains scripts and notebooks used for data processing and analysis.

- `metadata`: contains metadata relevant for samples.

- `gene_list`: lists of genes required in running some codes within the notebooks.

- `utils`: utility functions required within the notebooks.

- `data`:
  - `dandelion-remap`: remapped TCR/BCR files from 'Suo et al 2022 Science' used in the manuscript.
  - `cycloheximide`: anndata and TCR/BCR files used in the manuscript post treatment.
    Due to large file size limitations on github, the anndata is now hosted on a permanent ftp site:
    ftp://ftp.sanger.ac.uk/pub/users/kp9/PBMC_cycloheximide_predictions_compressed.h5ad