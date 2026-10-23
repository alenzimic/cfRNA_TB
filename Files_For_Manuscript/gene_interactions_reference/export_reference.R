# Regenerate the compact original-DESeq2 reference used by gene_interactions.ipynb.
# Usage (from any working directory):
# Rscript /path/to/cfRNA_TB/Files_For_Manuscript/gene_interactions_reference/export_reference.R /path/to/cfRNA_TB
# Requires DESeq2, jsonlite, and the sha256sum command. The notebook itself is Python.
# This script writes only training_de_genes.tsv and reference.json; no large matrices.
suppressPackageStartupMessages(library(DESeq2))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install the R package jsonlite.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript export_reference.R /path/to/cfRNA_TB")
repo_root <- normalizePath(args[1], mustWork = TRUE)
input_dir <- file.path(repo_root, "Files_For_Manuscript")
output_dir <- file.path(input_dir, "gene_interactions_reference")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_counts <- function(split) {
  as.matrix(read.delim(file.path(input_dir, paste0("counts.", split, ".tsv")),
                       row.names = 1, check.names = FALSE))
}
counts_train <- read_counts("train")
metadata <- read.delim(file.path(input_dir, "metadata.train.tsv"))
rownames(metadata) <- metadata$sample_id
stopifnot(identical(colnames(counts_train), rownames(metadata)))

# Figure_2.ipynb: training-only DESeq2, contrast positive versus negative.
dds <- DESeqDataSetFromMatrix(counts_train,
                              colData = metadata[, c("sample_id", "tb")], design = ~tb)
dds <- DESeq(dds)
results_train <- as.data.frame(results(dds, alpha = 0.05,
                                       contrast = c("tb", "positive", "negative")))
results_train$gene_id <- rownames(results_train)
annotation <- read.delim(file.path(input_dir, "gencode.biotype.name.key.tsv"))
results_train <- merge(results_train, annotation, by = "gene_id")
candidate_genes <- results_train[!is.na(results_train$padj) &
                                 results_train$padj < 0.05 &
                                 results_train$baseMean > 100, ]
candidate_genes <- candidate_genes[order(candidate_genes$padj, candidate_genes$gene_id),
                                    c("gene_id", "gene_name", "baseMean", "log2FoldChange", "padj")]
write.table(candidate_genes, file.path(output_dir, "training_de_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(dds)

# These are the exact dispersion-fit steps used by the original notebooks'
# varianceStabilizingTransformation(as.matrix(counts)) with blind=TRUE.
# Only the two fitted coefficients are retained; Python computes the transform.
blind_vst_coefficients <- function(counts) {
  object <- DESeqDataSetFromMatrix(counts, DataFrame(row.names = colnames(counts)), ~1)
  object <- estimateSizeFactors(object)
  object <- estimateDispersionsGeneEst(object, quiet = TRUE)
  object <- estimateDispersionsFit(object, quiet = TRUE, fitType = "parametric")
  if (attr(dispersionFunction(object), "fitType") != "parametric") {
    stop("The fitted VST is not parametric; the Python closed-form transform is invalid.")
  }
  as.list(attr(dispersionFunction(object), "coefficients"))
}
training_vst <- blind_vst_coefficients(counts_train)
counts_test <- read_counts("test")
counts_val <- read_counts("val")
stopifnot(identical(rownames(counts_train), rownames(counts_test)),
          identical(rownames(counts_train), rownames(counts_val)))
pooled_vst <- blind_vst_coefficients(cbind(counts_train, counts_test, counts_val))

input_paths <- file.path("Files_For_Manuscript",
                        c(paste0("counts.", c("train", "test", "val"), ".tsv"),
                          paste0("metadata.", c("train", "test", "val"), ".tsv"),
                          "gencode.biotype.name.key.tsv"))
sha256 <- function(relative_path) {
  result <- system2("sha256sum", shQuote(file.path(repo_root, relative_path)), stdout = TRUE)
  hash <- substr(result[1], 1, 64)
  if (!grepl("^[a-f0-9]{64}$", hash)) stop("Could not hash input: ", relative_path)
  hash
}
input_hashes <- as.list(setNames(vapply(input_paths, sha256, character(1)), input_paths))
reference <- list(
  r_version = as.character(getRversion()),
  deseq2_version = as.character(packageVersion("DESeq2")),
  published_deseq2_version = "1.34.0",
  input_sha256 = input_hashes,
  training_vst = training_vst,
  pooled_vst = pooled_vst,
  paper_reported_candidates = 113L,
  reproduced_candidates = nrow(candidate_genes),
  candidate_filter = "Training only: padj < 0.05 and baseMean > 100, as in Figure_2.ipynb.",
  discrepancy_note = paste(
    "The 2024 article reports 113 candidates, but its repository's exact training filter",
    "reproduces 142 with the audited inputs; the six saved marker statistics match.",
    "An exact 113-gene list was not found in the inspected article supplements.",
    "This reference reproduces the repository filter and does not certify the reported 113-gene set."
  ),
  inputsource = list(
    article = "https://www.nature.com/articles/s41467-024-49245-6",
    training_analysis = "Scripts_For_Manuscript/Figure_2.ipynb: training DESeq2 and candidate_genes cells",
    pooled_analysis = "Scripts_For_Manuscript/Figure_1.ipynb and Files_For_Manuscript/DAA_scripts.R: blind VST",
    regeneration = "Files_For_Manuscript/gene_interactions_reference/export_reference.R"
  ),
  size_factor_method = paste(
    "Default DESeq2 ratio: exp(median(log(count) - rowMean(log(reference counts))))",
    "over genes positive in every reference sample; do not geometrically center size factors."
  ),
  normalization_scope = paste(
    "Pooled coefficients are for descriptive Figure 1C only.",
    "New ML comparisons freeze training geometric means and training coefficients for test/validation.",
    "The original Figure 2 notebook instead independently fitted VST for each split."
  )
)
jsonlite::write_json(reference, file.path(output_dir, "reference.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("Wrote compact reference; training candidates: ", nrow(candidate_genes))
