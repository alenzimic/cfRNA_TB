# Circulating Cell-Free RNA in Blood as a Host Response Biomarker for the Detection of Tuberculosis
## This repository contains the pipeline and analysis scripts used to evaluate TB markers in plasma cfRNA.

The alignment pipeline for cfRNA has an individual Readme with instructions on how to run the pipeline. Full sequencing files will be made available on dbGAP.

All scripts and files to recreate the paper figures, including metadata and count matrices, can be found within the "Files_For_Manuscript" and "Scripts_For_Manuscript" folders.

R version 4.2.0 (2022-04-22)
Platform: x86_64-conda-linux-gnu

## Compare the biology of the two six-gene panels

Open the executed [gene-panel biology notebook](Scripts_For_Manuscript/gene_panel_biology.ipynb) for gene functions, curated pathway membership and observed adult TB associations. It reuses saved results without refitting models.

The notebook exports the [function diagram](Figure_Output/Gene_Interactions/focused/two_panel_functions.png), [pathway-membership map](Figure_Output/Gene_Interactions/focused/two_panel_pathway_membership.png) and [adult-association plot](Figure_Output/Gene_Interactions/focused/two_panel_tb_associations.png). Each figure is also available as a PDF in the same folder.

## Inspect model scores and HIV-positive samples

The active [model comparison](Scripts_For_Manuscript/gene_interactions.ipynb),
[cutoff interpretation](Scripts_For_Manuscript/gene_cutoff_interpretation.ipynb), and
[XGBoost interpretation](Scripts_For_Manuscript/gene_interactions_interpretation.ipynb)
notebooks place patient-level panels after the relevant evaluation cells.
Each panel combines its models, shows train/test/validation score distributions,
and marks HIV-positive participants. `[HIV+]` means recorded positive;
`[HIV?]` means unknown status, not negative.

The original violin/dot plots and every ranked patient bar are retained. Wider
ranked panels place errors and HIV status in separate rows below the bars, in the
same patient order: black marks an incorrect prediction; purple marks HIV positive
and gold marks unknown status. Blank cells mean correct prediction or HIV negative
in their respective rows. Score units, patient order and decision boundaries are
unchanged; the main legend explains both point markers and annotation-row colors.
The participants in these panels all have recorded ages of at least 18 years,
so no pediatric row is shown here. The shared renderer adds that row when a
plotted cohort contains participants recorded as pediatric; missing age stays unknown.

Open the [paper-GLM panel](Figure_Output/Model_Evaluation/paper.pdf),
[all six-gene classifiers](Figure_Output/Model_Evaluation/six_genes.pdf), or
[two gene-panel comparison](Figure_Output/Model_Evaluation/six_gene_panels.pdf).
These are local generated outputs; `Figure_Output/Model_Evaluation/` is excluded
from Git. The notebooks retain their embedded summary figures. Rerun the plotting
cells to regenerate the local files. Every panel has `_samples.pdf` pages with readable IDs and a searchable
`_samples.html` table; `_metrics.tsv` provides HIV subgroup counts and metrics.
One stored model fit gives one score per patient, not a repeated-test success rate.
See [decisions and validation scope](docs/DECISIONS.md).
