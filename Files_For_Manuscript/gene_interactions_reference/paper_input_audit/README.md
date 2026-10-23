# Audit of the reported 113-gene input

The original `Figure_2.ipynb` training filter selects **142 genes** with both
DESeq2 1.30.1 and the article's stated **R 4.1.3 / DESeq2 1.34.0** versions.
`training_candidates.tsv` records the exact-version result. It is not a recovered
113-gene paper list.

`audit.json` records the versions, selection rule, counts, checksums and diagnostic
checks. All 130 training samples pass QC. The selected genes have unique IDs and
symbols, no missing annotation and no overlap with the repository blacklist.
Selecting only protein-coding genes would give 138, not 113. No additional filter
or truncation was applied. The 312 outlier replacements and six reported marker
statistics agree with the original notebook output.

Reproduce from the repository root:

```bash
/programs/R-4.1.3-r9/bin/Rscript Files_For_Manuscript/gene_interactions_reference/paper_input_audit/export_audit.R "$PWD"
```

On another machine, use an R 4.1.3 environment with DESeq2 1.34.0 and jsonlite.
The exporter verifies these versions. An optional cached DDS argument avoids
refitting but verifies its original counts, sample labels and DESeq2 version.

[Article and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11164910/).
