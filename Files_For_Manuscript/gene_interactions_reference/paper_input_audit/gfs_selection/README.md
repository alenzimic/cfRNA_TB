# Audit of the original greedy forward selection

Running the actual selection loop from `Scripts_For_Manuscript/Figure_2.ipynb`
(cell 9, zero-based) over all **142 reconstructed candidates** recovers the six
published genes in the same order: **GBP5, BNIP3L, KLF6, DYSF, LASP1, PCBP1**.
The seventh proposed gene, JMJD1C, increases the original composite score by
0.000721, below its strictly-greater-than-0.01 acceptance threshold. All 994
candidate fits are retained. This was a selection rerun, not a fixed-six-gene fit.

The algorithm preserves the original unpenalized binomial-logit GLM without an
intercept or additional scaling, pROC defaults, first maximum Youden threshold,
score of AUC + sensitivity + specificity, and correlation tie breaker. The
correlation field is NA at the first step because there are no selected genes.
The full candidate loop includes already-selected genes, as the source does;
R's formula handling removes duplicated terms.

The input list in `../training_candidates.tsv` was independently reconstructed
with the article's **R 4.1.3 / DESeq2 1.34.0**. The original training VST reference
was fitted with **DESeq2 1.30.1**. The first audit used a cached full VST matrix;
the portable rerun below instead computes its training size factors from all
training genes and applies the saved parametric trend from `../../reference.json`.
No large VST file is included. The portable rerun was checked against all 994
cached-matrix candidate fits; selection, scores, and stopping decisions agree
within floating-point precision. The gene list is not truncated to 113.

## Rerun

Install `pROC`, `dplyr`, and `jsonlite` in the R environment, then run:

```bash
Rscript Files_For_Manuscript/gene_interactions_reference/paper_input_audit/gfs_selection/export_audit.R /path/to/cfRNA_TB
```

An optional second argument writes results to another directory. The saved run
uses R 4.1.3, pROC 1.18.0, and dplyr 1.0.10; `audit.json` records every version.
The script verifies the candidate and training-input hashes before fitting.
Neither test nor validation samples are used by selection. After selection,
the known six-gene reference matrices are scored only to verify predictions.

## Files

- `export_audit.R`: portable selection and verification script.
- `selected_genes.tsv`: selected IDs, names, DE statistics, and selection order.
- `selection_steps.csv`: seven proposed additions, including the rejected one.
- `all_candidate_scores.csv`: all 142 candidates at each of seven steps.
- `coefficients.csv`: the selected final model's unpenalized coefficients.
- `audit.json`: source, settings, versions, input/output SHA256 hashes, selected
  genes, match checks, stop decision, and rescored reference AUCs.
- `cached_vs_portable_verification.json`: comparison with the first cached-VST
  run, including numerical differences and its original matrix hash.

This establishes that the repository's 142-candidate pipeline recovers the final
published GFS model. It **does not recover or certify the reported 113-gene
candidate pool**, or establish identical inputs for the historical RF,
ExtraTrees, and other reported Figure 1F classifiers.
