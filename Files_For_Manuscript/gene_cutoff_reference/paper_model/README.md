# Reconstructed published six-gene model

This benchmark follows `Scripts_For_Manuscript/Figure_2.ipynb`, cells 7, 11 and 13
(zero-based): GBP5, BNIP3L, KLF6, DYSF, LASP1 and PCBP1; binomial-logit GLM;
no intercept, no penalty and no feature scaling. The original code independently
fits blind DESeq2 VST in each split. Its Youden threshold is selected on training
samples only.

`model.json` records the training-refit coefficients, full-precision threshold,
VST parameters, software versions and file hashes. Original serialized model
weights were not available. `paper_features.tsv` contains only the six normalized
features, with sample identities and split labels. `r_predictions.csv` and
`r_metrics.csv` allow the Python notebook to verify its own scoring:

```python
probability = expit(features[genes].to_numpy() @ coefficients)
predicted = probability > threshold
```

The exporter requires the three ROC-AUCs and threshold to match published rounded
values before writing these files. This historical benchmark is separate from
the new cross-validation search. Applying its coefficients to frozen training
normalization is a different preprocessing protocol and must be labelled as such.

Regenerate from the repository root:

```bash
Rscript Files_For_Manuscript/gene_cutoff_reference/paper_model/export_reference.R "$PWD"
```

[Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC11164910/).
