# Reconstruct the paper's final six-gene model from Figure_2.ipynb.
# Usage: Rscript export_reference.R /path/to/cfRNA_TB
# R fits the original model once; the Python notebook scores the exported coefficients.
suppressPackageStartupMessages(library(DESeq2))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript export_reference.R /path/to/cfRNA_TB")
repo_root <- normalizePath(args[1], mustWork = TRUE)
input_dir <- file.path(repo_root, "Files_For_Manuscript")
output_dir <- file.path(input_dir, "gene_cutoff_reference/paper_model")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(87)

genes <- c("ENSG00000154451.14", "ENSG00000104765.16", "ENSG00000067082.15",
           "ENSG00000135636.15", "ENSG00000002834.18", "ENSG00000169564.7")
gene_names <- c("GBP5", "BNIP3L", "KLF6", "DYSF", "LASP1", "PCBP1")

prepare_split <- function(split) {
  message("Preparing original independent VST for ", split)
  counts <- as.matrix(read.delim(file.path(input_dir, paste0("counts.", split, ".tsv")),
                                 row.names = 1, check.names = FALSE))
  metadata <- read.delim(file.path(input_dir, paste0("metadata.", split, ".tsv")))
  stopifnot(identical(colnames(counts), metadata$sample_id), all(genes %in% rownames(counts)))
  # These are the exact steps inside varianceStabilizingTransformation(matrix),
  # with its defaults blind=TRUE and fitType="parametric". Keep the fit coefficients.
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(row.names = colnames(counts)), ~1)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
  dds <- estimateDispersionsFit(dds, quiet = TRUE, fitType = "parametric")
  stopifnot(attr(dispersionFunction(dds), "fitType") == "parametric")
  vst <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  features <- as.data.frame(t(vst[genes, , drop = FALSE]))
  features$sample_id <- rownames(features)
  features$split <- split
  features$tb_binary <- as.integer(metadata$tb == "positive")
  features <- features[order(features$sample_id), c("sample_id", "split", "tb_binary", genes)]
  list(features = features,
       vst = as.list(attr(dispersionFunction(dds), "coefficients")))
}

splits <- c("train", "test", "val")
prepared <- setNames(lapply(splits, prepare_split), splits)
training <- prepared$train$features
training$tb <- factor(training$tb_binary, levels = c(0, 1))
formula <- as.formula(paste("tb ~", paste(genes, collapse = " + "), "- 1"))
model <- glm(formula, data = training, family = binomial())
stopifnot(model$converged, length(coef(model)) == 6L, all(is.finite(coef(model))))

features <- do.call(rbind, lapply(prepared, function(value) value$features))
rownames(features) <- NULL
probabilities <- predict(model, newdata = features, type = "response")
train_probabilities <- probabilities[features$split == "train"]
ordered_probabilities <- sort(unique(train_probabilities))
# pROC's '<' direction uses ascending midpoint thresholds. Its first maximum
# Youden index matches which.max(coords(roc, "all")$youden) in the original code.
thresholds <- c(-Inf, (head(ordered_probabilities, -1L) + tail(ordered_probabilities, -1L)) / 2, Inf)
youden <- vapply(thresholds, function(threshold) {
  predicted <- train_probabilities > threshold
  mean(predicted[training$tb_binary == 1]) + mean(!predicted[training$tb_binary == 0]) - 1
}, numeric(1))
threshold <- thresholds[which.max(youden)]

roc_auc <- function(labels, scores) {
  positive <- labels == 1
  n_positive <- sum(positive)
  (sum(rank(scores)[positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * sum(!positive))
}
metrics <- do.call(rbind, lapply(splits, function(split) {
  keep <- features$split == split
  labels <- features$tb_binary[keep]
  scores <- probabilities[keep]
  predicted <- as.integer(scores > threshold)
  data.frame(split = split, n_samples = sum(keep), roc_auc = roc_auc(labels, scores),
             accuracy = mean(predicted == labels),
             sensitivity = mean(predicted[labels == 1] == 1),
             specificity = mean(predicted[labels == 0] == 0))
}))
stopifnot(identical(round(metrics$roc_auc, 3), c(0.921, 0.947, 0.946)),
          round(threshold, 3) == 0.498)

write.table(features, file.path(output_dir, "paper_features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
predictions <- features[, c("sample_id", "split", "tb_binary")]
predictions$probability <- probabilities
predictions$predicted <- as.integer(probabilities > threshold)
write.csv(predictions, file.path(output_dir, "r_predictions.csv"), row.names = FALSE)
write.csv(metrics, file.path(output_dir, "r_metrics.csv"), row.names = FALSE)

sha256 <- function(relative_path) {
  output <- system2("sha256sum", shQuote(file.path(repo_root, relative_path)), stdout = TRUE)
  hash <- substr(output[1], 1, 64)
  if (!grepl("^[a-f0-9]{64}$", hash)) stop("Could not hash ", relative_path)
  hash
}
inputs <- c(file.path("Files_For_Manuscript", c(paste0("counts.", splits, ".tsv"),
                                                paste0("metadata.", splits, ".tsv"))),
            "Scripts_For_Manuscript/Figure_2.ipynb")
outputs <- file.path("Files_For_Manuscript/gene_cutoff_reference/paper_model",
                     c("paper_features.tsv", "r_predictions.csv", "r_metrics.csv"))
reference <- list(
  label = "Paper six-gene GLM: reconstructed original configuration",
  provenance = paste("No serialized fitted model or full-precision coefficients were found in the",
                     "repository or inspected supplements. These coefficients are a training-only",
                     "refit using the exact final six genes and GLM specification in Figure_2.ipynb.",
                     "The saved paper AUCs and threshold are reproduced to published precision."),
  article = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11164910/",
  source_notebook = "Scripts_For_Manuscript/Figure_2.ipynb",
  source_cells_zero_based = list(normalization = 7L, selection = 9L, gene_order = 11L, model = 13L),
  r_version = as.character(getRversion()),
  deseq2_version = as.character(packageVersion("DESeq2")),
  published_r_version = "4.1.3", published_deseq2_version = "1.34.0",
  input_sha256 = as.list(setNames(vapply(inputs, sha256, character(1)), inputs)),
  output_sha256 = as.list(setNames(vapply(outputs, sha256, character(1)), outputs)),
  genes = genes, gene_names = gene_names,
  coefficients = as.list(coef(model)), intercept = 0,
  formula = paste(deparse(formula), collapse = " "),
  family = "binomial", link = "logit", penalty = "none", scaling = "none",
  fit_intercept = FALSE, method = "glm.fit", solver = "IRLS",
  glm_control = glm.control(), iterations = model$iter, converged = model$converged,
  threshold = threshold, threshold_rule = "P(TB positive) > training Youden threshold",
  published_threshold_rounded = 0.498,
  published_auc_rounded = list(train = 0.921, test = 0.947, val = 0.946),
  split_vst = lapply(prepared, function(value) value$vst),
  preprocessing = paste("Original Figure_2.ipynb: independent blind DESeq2 VST in each split,",
                        "including separate size-factor geometric means and dispersion fits.",
                        "This reproduces the paper's batch-dependent preprocessing and differs",
                        "from the cutoff study's frozen training preprocessing."),
  selection_note = paste("The six genes were selected by GFS on the full original training set.",
                         "The paper selected GFS among fifteen families using test ROC-AUC.",
                         "Treat this fixed benchmark as a historical comparator, not a candidate",
                         "in the new cross-validation selection procedure."),
  gfs_settings = list(candidate_filter = "padj < 0.05 and baseMean > 100",
                       score = "training ROC-AUC + sensitivity + specificity at Youden threshold",
                       minimum_score_gain = 0.01,
                       tie_breaker = "lowest absolute mean correlation with previously selected genes")
)
jsonlite::write_json(reference, file.path(output_dir, "model.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = NA)
print(metrics)
message("Training Youden threshold: ", format(threshold, digits = 17))
message("Exported the reconstructed paper model and six-gene features.")
