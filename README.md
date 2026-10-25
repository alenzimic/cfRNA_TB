# Circulating Cell-Free RNA in Blood as a Host Response Biomarker for the Detection of Tuberculosis
## This repository contains the pipeline and analysis scripts used to evaluate TB markers in plasma cfRNA.

The alignment pipeline for cfRNA has an individual Readme with instructions on how to run the pipeline. Full sequencing files will be made available on dbGAP.

All scripts and files to recreate the paper figures, including metadata and count matrices, can be found within the "Files_For_Manuscript" and "Scripts_For_Manuscript" folders.

R version 4.2.0 (2022-04-22)
Platform: x86_64-conda-linux-gnu

## Compare the biology of the two six-gene panels

Open the executed [gene-panel biology notebook](Scripts_For_Manuscript/gene_panel_biology.ipynb) for gene functions, curated pathway membership and observed adult TB associations. It reuses saved results without refitting models.

The notebook exports the [function diagram](Figure_Output/Gene_Interactions/focused/two_panel_functions.png), [pathway-membership map](Figure_Output/Gene_Interactions/focused/two_panel_pathway_membership.png) and [adult-association plot](Figure_Output/Gene_Interactions/focused/two_panel_tb_associations.png). Each figure is also available as a PDF in the same folder.
