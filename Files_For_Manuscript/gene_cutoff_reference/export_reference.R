# Export training-only DESeq2 statistics and fold-specific VST coefficients.
# Usage: Rscript export_reference.R /path/to/cfRNA_TB
# Requires DESeq2, jsonlite, sha256sum, and the supplied folds.csv.
# No test or external-validation counts or labels are read by this script.
suppressPackageStartupMessages(library(DESeq2))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript export_reference.R /path/to/cfRNA_TB")
repo_root <- normalizePath(args[1], mustWork = TRUE)
input_dir <- file.path(repo_root, "Files_For_Manuscript")
output_dir <- file.path(input_dir, "gene_cutoff_reference")

counts <- as.matrix(read.delim(file.path(input_dir, "counts.train.tsv"),
                               row.names = 1, check.names = FALSE))
metadata <- read.delim(file.path(input_dir, "metadata.train.tsv"), check.names = FALSE)
rownames(metadata) <- metadata$sample_id
annotation <- read.delim(file.path(input_dir, "gencode.biotype.name.key.tsv"))
folds <- read.csv(file.path(output_dir, "folds.csv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts), metadata$sample_id),
          identical(folds$sample_id, metadata$sample_id),
          identical(folds$tb_binary, as.integer(metadata$tb_binary)),
          identical(sort(unique(folds$validation_fold)), 0:4),
          !anyDuplicated(annotation$gene_id), !anyDuplicated(folds$sample_id))

sha256 <- function(relative_path) {
  result <- system2("sha256sum", shQuote(file.path(repo_root, relative_path)), stdout = TRUE)
  hash <- substr(result[1], 1, 64)
  if (!grepl("^[a-f0-9]{64}$", hash)) stop("Could not hash input: ", relative_path)
  hash
}

blind_vst_coefficients <- function(fit_counts) {
  object <- DESeqDataSetFromMatrix(fit_counts, DataFrame(row.names = colnames(fit_counts)), ~1)
  object <- estimateSizeFactors(object)
  object <- estimateDispersionsGeneEst(object, quiet = TRUE)
  object <- estimateDispersionsFit(object, quiet = TRUE, fitType = "parametric")
  if (attr(dispersionFunction(object), "fitType") != "parametric") {
    stop("VST fit is not parametric; the Python closed-form transform is invalid.")
  }
  as.list(attr(dispersionFunction(object), "coefficients"))
}

export_fit <- function(name, training_samples, held_out_samples = character()) {
  started <- Sys.time()
  message("Starting ", name, " with ", length(training_samples), " training samples")
  fit_counts <- counts[, training_samples, drop = FALSE]
  fit_metadata <- metadata[training_samples, c("sample_id", "tb"), drop = FALSE]
  fit_metadata$tb <- factor(fit_metadata$tb, levels = c("negative", "positive"))
  dds <- DESeqDataSetFromMatrix(fit_counts, fit_metadata, design = ~tb)
  dds <- DESeq(dds, quiet = TRUE)
  filtered <- results(dds, alpha = 0.05, contrast = c("tb", "positive", "negative"))
  unfiltered <- results(dds, independentFiltering = FALSE,
                        contrast = c("tb", "positive", "negative"))
  table <- as.data.frame(filtered)
  table$gene_id <- rownames(table)
  table$gene_name <- annotation$gene_name[match(table$gene_id, annotation$gene_id)]
  table$padj_unfiltered <- unfiltered$padj
  table <- table[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE",
                     "stat", "pvalue", "padj", "padj_unfiltered")]
  filename <- paste0(name, "_de.tsv.gz")
  connection <- gzfile(file.path(output_dir, filename), "wt")
  write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  close(connection)
  rm(dds)
  coefficients <- blind_vst_coefficients(fit_counts)
  reference <- list(
    de_file = filename,
    de_sha256 = sha256(file.path("Files_For_Manuscript/gene_cutoff_reference", filename)),
    training_samples = training_samples,
    held_out_samples = held_out_samples,
    n_training_samples = length(training_samples),
    n_positive = sum(fit_metadata$tb == "positive"),
    n_negative = sum(fit_metadata$tb == "negative"),
    n_genes = nrow(table),
    n_repository_de_genes = sum(!is.na(table$padj) & table$padj < 0.05 & table$baseMean > 100),
    independent_filter_threshold = as.numeric(S4Vectors::metadata(filtered)$filterThreshold),
    vst = coefficients,
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  message("Finished ", name, " in ", round(reference$seconds, 1), " seconds")
  reference
}

input_paths <- file.path("Files_For_Manuscript", c("counts.train.tsv", "metadata.train.tsv",
                           "gencode.biotype.name.key.tsv", "gene_cutoff_reference/folds.csv"))
reference <- list(
  r_version = as.character(getRversion()),
  deseq2_version = as.character(packageVersion("DESeq2")),
  input_sha256 = as.list(setNames(vapply(input_paths, sha256, character(1)), input_paths)),
  folds = list(method = "sklearn.model_selection.StratifiedKFold", sklearn_version = "1.6.1",
                n_splits = 5L, shuffle = TRUE, random_state = 87L,
                column = "validation_fold", index_base = 0L),
  design = "~tb", contrast = c("tb", "positive", "negative"),
  padj = "DESeq2 default independent filtering with alpha=0.05, reproducing the repository filter.",
  padj_unfiltered = "DESeq2 BH-adjusted p values with independentFiltering=FALSE; provided explicitly for alternative sensitivity analyses.",
  normalization = paste("Each fit uses only its training samples for blind parametric VST coefficients.",
                        "Python must compute geometric means from those same training samples, then",
                        "freeze geometric means and these coefficients to transform held-out samples."),
  scope = "Only counts.train.tsv and metadata.train.tsv are read. No test or external-validation fit or selection.",
  fits = list()
)
fit_samples <- c(list(full_train = metadata$sample_id),
                 setNames(lapply(0:4, function(fold) folds$sample_id[folds$validation_fold != fold]),
                          paste0("fold_", 0:4)))
for (name in names(fit_samples)) {
  training_samples <- fit_samples[[name]]
  held_out <- metadata$sample_id[!metadata$sample_id %in% training_samples]
  reference$fits[[name]] <- export_fit(name, training_samples, held_out)
  jsonlite::write_json(reference, file.path(output_dir, "reference.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = NA)
}
message("Wrote all six training-only references to ", output_dir)
