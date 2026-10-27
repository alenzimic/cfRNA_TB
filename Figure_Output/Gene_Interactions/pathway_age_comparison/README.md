# Pathway analysis and continuous-age comparison

This is an archived **aggregate result snapshot** from the executed local
`cfRNA_TB_overview/notebooks/02_pathway_enrichment_by_age.ipynb`, together with a
comparison against `Scripts_For_Manuscript/gene_interactions_pathways.ipynb`.
The snapshot contains model estimates, tests, coverage, aggregate counts,
figures and provenance. It does not include participant metadata, individual
pathway scores, raw counts or sample-level QC records.

## Findings and interpretation

- The earlier manuscript analysis supports general TB-associated expression,
  most consistently for interferon-gamma after excluding the selected genes.
  The focused comparison table retains all five previously specified programs,
  including uncertain and untested comparisons. The exhaustive table retains
  all 49 shared eligible Hallmark sets.
- The age analysis found **0 of
  250 age-interaction tests at FDR <0.05**
  (smallest q=0.985744). Its separate omnibus
  test for any TB difference identified 6
  pathway–study associations. General TB association does not demonstrate an
  age-dependent TB difference.
- All **251 earlier participants overlap** the
  812-person age analysis. GHL, R2D2 and TBSQUARE are
  exactly this overlapping group; the 561 additional people
  are ADTB and PROTECT. Of the 6 general-TB
  discoveries, 6
  occur in overlapping adult studies and
  0 occur in ADTB or PROTECT.
  Agreement is **same-participant analytical
  concordance**, not independent replication or confirmed transfer.
- Children and adults were sampled in different studies. The age curves are
  fitted separately by study and do not infer a continuous biological transition
  across unsampled adolescence. Sparse ages are masked, including internal gaps.

## Different scores and questions

| Feature | Earlier manuscript pathways | Continuous-age analysis |
|---|---|---|
| Participants | 251, divided into training/test/validation | 812, divided into five studies |
| Expression | Training-reference DESeq2 VST | Within-person ranks of length-adjusted counts |
| Participant pathway score | Average standardized expression; six selected genes excluded in focused comparison | Directional single-sample rank enrichment; all measured set members included |
| SD scale | Training-reference SD | Pooled SD across 812 people |
| Primary comparison | TB–nonTB level difference | Whether TB–nonTB difference changes with age |
| Adjustment in the comparison | Recruiting study + estimated RNA sources | Smooth age, sex, HIV, country and count-read length as estimable within study |
| Multiple testing | 30 focused pathway-score comparisons; GSEA FDR separately per run | 250 pathway–study tests for each primary question |

The old GSEA normalized enrichment score, the old participant mean-expression
score and the new single-sample rank score are different quantities. Their effect
sizes and adjusted p-values are not interchangeable. New curve extrema in the
comparison CSV are descriptive summaries of supported ages, not additional tests.
A GSEA output of zero reflects finite-permutation resolution, not a literal zero
probability. Neither analysis establishes biochemical activity or causation.

## Contents

- `focused_pathway_comparison.csv`: five prespecified older programs, their old
  effects/intervals/FDR and the corresponding new study tests/curve ranges.
- `all_shared_hallmark_comparison.csv`: all shared eligible Hallmarks under the
  old specification adjusting for study and estimated RNA sources and excluding
  the six XGBoost genes. This prevents selecting only concordant pathways.
- `participant_overlap_summary.csv`: aggregate verified overlap by old split and
  study. Participant identities were checked locally but are not exported here.
- `age_analysis_any_TB_hits.csv`: all new general-TB discoveries, with a flag
  indicating whether each set was eligible in the old GSEA.
- `primary_tests.tsv` and `primary_curves.tsv`: complete new test family and
  supported age-grid estimates. Missing curves are unsupported, not zero.
- `sensitivity_*`, `model_diagnostics.json`: robustness variants and explicit
  reasons for any skipped model fits.
- `pathway_timelines_all.pdf`: all eligible pathways, one pathway per page, with
  all studies on the same age axis. Color identifies study; solid lines are TB
  and dashed lines are nonTB. Pale ribbons are pointwise 95% confidence intervals.
  The support strip below each plot shows observed ages and supported comparisons
  separately for each study, preserving internal gaps and sparse tails. These are
  separate study fits, not a pooled lifespan model; study differences are not
  isolated age effects.
- `pathway_timeline_index.tsv`: pathway-to-page lookup, with a flag identifying
  the focus pathways shown in the notebook. Their `timeline_*.png` and `.pdf`
  figures are included individually.
- `age_pathway_analysis.html`: archived HTML report with embedded figures and
  aggregate outputs. It loads MathJax/RequireJS from a CDN for rendering and is
  not an independently rerunnable analysis package.
- `age_analysis_provenance.json`: original overview source/software/output hashes.
- `comparison_audit.json`: local input hashes, matching checks and comparison scope.
- `snapshot_manifest.json`: exact SHA-256 values of every published snapshot file.

## Rebuild

From the `cfRNA_TB` repository root, after executing the original overview and
manuscript notebooks:

```bash
.venv/bin/python Scripts_For_Manuscript/export_pathway_age_comparison.py   --overview ../cfRNA_TB_overview
```

The exporter reads `../cfRNA_TB_overview/results/pathway_age/`, the manuscript
`Figure_Output/Gene_Interactions/focused/pathway*.csv` results, and old/new
metadata locally to verify overlap. These source paths and hashes are recorded
without copying participant records. Use the original overview notebook and its
local data/resources to refit the age models; this exporter only copies and
compares completed results.
