# Training-only references for the cutoff experiment

`folds.csv` assigns each of the 130 original training samples to one held-out fold.
`full_train_de.tsv.gz` contains the complete training DESeq2 result. The five
`fold_*_de.tsv.gz` tables repeat DESeq2 using only the other four folds. Each table
contains every input gene, including genes with unavailable p values.

`reference.json` records the exact input hashes, software versions, fit sample
identities, and blind parametric VST coefficients. Freeze each fit's geometric
means and VST coefficients before transforming its held-out samples. The test and
external-validation sets are never read by this exporter.

`padj` uses the original repository's default DESeq2 independent filtering at
`alpha=0.05`. Thresholding this column varies the selection cutoff while keeping
that original independent filter fixed. `padj_unfiltered` is also provided: it
uses BH adjustment without independent filtering and can support a separately
labelled sensitivity analysis. Both columns retain DESeq2's default outlier rules.

To regenerate from the repository root:

```bash
.venv/bin/python Files_For_Manuscript/gene_cutoff_reference/generate_folds.py
Rscript Files_For_Manuscript/gene_cutoff_reference/export_reference.R "$PWD"
```

Python requires the existing notebook environment. R requires DESeq2 and jsonlite.
The exporter writes compact reference tables and JSON, with no expression matrices.
