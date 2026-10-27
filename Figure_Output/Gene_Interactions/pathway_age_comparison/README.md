# Pathway analysis and continuous-age comparison

This is an archived **aggregate result snapshot** from the executed local
`cfRNA_TB_overview/notebooks/02_pathway_enrichment_by_age.ipynb`, together with a
comparison against `Scripts_For_Manuscript/gene_interactions_pathways.ipynb`.
The snapshot contains model estimates, tests, coverage, aggregate counts,
figures and provenance. Faint patient dots show observed distributions in the
figures; the interactive explorer embeds them as raster images. It does not
include participant identifiers, individual score tables, raw counts or
sample-level QC records.

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
- Children and adults were sampled in different studies. The primary age curves are
  fitted separately by study and do not infer a continuous biological transition
  across unsampled adolescence. Sparse ages are masked, including internal gaps.
- The declared-age follow-up estimates TB-minus-nonTB differences at reporting
  ages 5, 10, 20, 30, 40, 50, 60 using the original study models and exact spline
  contrasts. **0/750 supported comparisons pass age-contrast
  FDR <0.05** (minimum q=0.151223). One BH family
  covers all supported pathways, studies and reporting ages; unsupported planned
  rows retain reasons/counts but have no estimate or p-value. These are evaluation
  ages for an exploratory follow-up, not fitting bins or preregistered cutoffs.
  The age-specific tests, global age interaction and any-TB omnibus tests answer
  different questions. No continuous significant age interval is claimed.
- The supplementary age-bin analysis uses declared ranges **0–<5, 5–<10, 10–<18, 18–<30, 30–<45, 45–<60, 60+ years**
  and fits each study/bin independently, adjusting for linear age and clinical/
  protocol covariates. **7/550 tested comparisons pass bin FDR
  <0.05** (minimum q=0.024212). One separate BH family includes
  all estimable pathway/study/bin comparisons. Sparse or unreliable fits remain
  unavailable. These are age-range associations under a constant TB difference
  within each bin, not exact-age effects or evidence of an age interaction.
  The bin model, population and number of tests differ from the exact-age
  analysis; discoveries do not establish that binning is preferable. Cutoffs
  were declared before these tests and were not optimized for significance.

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
- `pathway_timelines_all.pdf`: all eligible pathways, one per page, with two
  group curves on one age axis. Blue solid shows TB; orange dashed shows nonTB.
  Children retain the original PROTECT model. Adults use an exploratory summary
  with fixed study/clinical composition, study-specific baseline and TB levels,
  and shared smooth age shapes. These are new summary estimates, not averages of
  whichever studies happen to cover an age. The adolescent gap stays unconnected.
  Whiskers show pointwise model intervals conditional on the model assumptions;
  study-omission sensitivity is reported. At ages covered by fewer studies, the
  summary relies on the shared-shape assumption to represent the other studies.
  Faint dots show actual patient ages and observed standardized scores, before
  covariate adjustment. They remain visible where sparse support withholds curves.
  The closer default vertical scale retains all model intervals; counts of patient
  dots outside the displayed limits are printed explicitly. All observations remain
  in the models. `pathway_patient_distributions_all.pdf` shows the full patient
  score range, using the same pathway-to-page index.
- `summary_*`: separate summary-model tests, complete predictions, diagnostics,
  and leave-one-adult-study-out checks. Its FDR families differ from the original
  study-specific primary analysis, which remains the basis of the manuscript
  comparison tables. The shared-age-shape assumption is not validated: both
  full-range and supported-age formal heterogeneity checks fail the leverage
  safeguards and are explicitly unavailable. Omission checks include changes in
  the reference population and do not constitute independent validation.
- `pathway_age_explorer.html`: self-contained interactive view with one age axis,
  pathway/source selectors, optional confidence bands and a patient-dot toggle.
  Its vertical-scale selector switches between Focus on curves and All patient
  scores. Original-study views include exact-age effect/CI/q tables; these tests
  are not assigned to the exploratory adult-summary or omission curves.
  A separate age-bin evidence strip and table show each original study's
  independent within-bin tests, with counts and unavailable reasons.
  Dots follow the study/omission selection and are embedded as transparent raster
  images, without patient identifiers or machine-readable patient rows. Numerical
  model estimates are aggregate results. The file requires no external scripts.
  The primary figure
  design follows two-group lifespan displays such as
  [Bethlehem et al., 2022](https://www.nature.com/articles/s41586-022-04554-y),
  while preserving the limitations of the much smaller, disjoint-age TB cohorts.
- `pathway_timeline_index.tsv`: pathway-to-page lookup, with a flag identifying
  the focus pathways shown in the notebook. Their `timeline_*.png` and `.pdf`
  figures are included individually.
- `age_contrasts.tsv`: all planned study/pathway/reporting-age combinations,
  with direct HC3/t effects, pointwise 95% CIs, p-values, global BH q-values and
  support information. Original singscore-unit effects re-express the same SD
  contrast; they are not fold changes, GSEA NES or independent evidence.
  `age_contrast_diagnostics.json` documents the family and unchanged interaction
  tests; `age_contrast_summary.tsv` summarizes results by study.
- `age_contrast_heatmap.*`: all pathways at supported reporting ages; color is
  effect size and a black dot would mark age-contrast q <0.05. Confidence intervals
  are pointwise, so a CI can exclude zero without an FDR discovery.
  `age_contrast_examples.*` compares group levels and direct difference curves
  for four pairs selected by minimum age-specific q; the TSV records that choice.
- `age_bins.tsv`: all planned study/pathway/bin comparisons, adjusted effects,
  pointwise HC3/t intervals, raw p-values and globally corrected bin q-values.
  `age_bin_counts.tsv`, `age_bin_summary.tsv` and `age_bin_diagnostics.json`
  retain group counts, complete availability reasons and design assumptions.
  Tests require at least 10 patients per diagnosis, overlapping observed ages,
  at least 15 residual degrees of freedom, leverage below 0.95 and an identifiable
  TB effect after nuisance adjustment. Confidence intervals are not adjusted
  for multiplicity; a CI excluding zero is not equivalent to a bin FDR discovery.
- `pathway_age_bins_all.pdf`: original two-curve timelines with separate study
  evidence rows underneath, using the same 50-pathway page order. Horizontal
  spans identify bins, not significant continuous-age windows. Filled markers
  identify bin q<0.05 (blue higher TB, orange higher nonTB), open markers mean
  no FDR discovery and crosses mean unavailable. `binned_timeline_*.png` and
  `.pdf` are the same focus pathways shown in the original timeline figures.
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
