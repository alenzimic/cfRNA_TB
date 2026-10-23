# Audit the reported 113-gene input using the paper's stated R/DESeq2 versions.
# Usage: /programs/R-4.1.3-r9/bin/Rscript export_audit.R /path/to/cfRNA_TB
# An optional second argument reuses a DDS fitted by this exact version. Its raw
# counts, sample identities, labels and recorded DESeq2 version are checked.
suppressPackageStartupMessages(library(DESeq2))
stopifnot(as.character(getRversion()) == "4.1.3",
          as.character(packageVersion("DESeq2")) == "1.34.0")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) stop("Usage: Rscript export_audit.R REPOSITORY [CACHED_DDS]")
repo_root <- normalizePath(args[1], mustWork = TRUE)
input_dir <- file.path(repo_root, "Files_For_Manuscript")
output_dir <- file.path(input_dir, "gene_interactions_reference/paper_input_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_counts <- as.matrix(read.delim(file.path(input_dir, "counts.train.tsv"),
                                   row.names = 1, check.names = FALSE))
samples <- read.delim(file.path(input_dir, "metadata.train.tsv"))
rownames(samples) <- samples$sample_id
annotation <- read.delim(file.path(input_dir, "gencode.biotype.name.key.tsv"))
blacklist <- read.delim(file.path(input_dir, "genes_to_exlude.tsv"), header = FALSE)
previous <- read.delim(file.path(input_dir, "gene_interactions_reference/training_de_genes.tsv"))
stopifnot(identical(colnames(raw_counts), samples$sample_id),
          !anyDuplicated(annotation$gene_id))

if (length(args) == 2L) {
  dds <- readRDS(args[2])
  stopifnot(identical(dimnames(raw_counts), dimnames(counts(dds))),
            all(raw_counts == counts(dds)),
            identical(as.character(dds$tb), samples$tb),
            as.character(S4Vectors::metadata(dds)$version) == "1.34.0")
} else {
  dds <- DESeqDataSetFromMatrix(raw_counts, samples[, c("sample_id", "tb")], ~tb)
  dds <- DESeq(dds)
}

# Figure_2.ipynb cell 5: the first result is unused; the second supplies res_df.
unused_default <- results(dds)
result <- as.data.frame(results(dds, alpha = 0.05,
                                contrast = c("tb", "positive", "negative")))
result$gene_id <- rownames(result)
result <- merge(result, annotation, by = "gene_id")
# Figure_2.ipynb cell 9: this is the complete, unmodified candidate filter.
selected <- result[!is.na(result$padj) & result$padj < 0.05 & result$baseMean > 100, ]
selected <- selected[order(selected$padj), ]
candidate_file <- "training_candidates.tsv"
write.table(selected, file.path(output_dir, candidate_file),
            sep = "\t", quote = FALSE, row.names = FALSE)

sha256 <- function(relative_path) {
  value <- system2("sha256sum", shQuote(file.path(repo_root, relative_path)), stdout = TRUE)
  hash <- substr(value[1], 1, 64)
  if (!grepl("^[a-f0-9]{64}$", hash)) stop("Could not hash ", relative_path)
  hash
}
input_paths <- c(file.path("Files_For_Manuscript", c("counts.train.tsv", "metadata.train.tsv",
                     "gencode.biotype.name.key.tsv", "genes_to_exlude.tsv",
                     "gene_interactions_reference/training_de_genes.tsv")),
                 "Scripts_For_Manuscript/Figure_2.ipynb")
matched <- selected[match(previous$gene_id, selected$gene_id), ]
statistic_difference <- vapply(c("baseMean", "log2FoldChange", "padj"), function(column) {
  max(abs(matched[[column]] - previous[[column]]), na.rm = TRUE)
}, numeric(1))
audit <- list(
  status = "published_113_gene_set_unresolved",
  r_version = as.character(getRversion()),
  deseq2_version = as.character(packageVersion("DESeq2")),
  selected_gene_count = nrow(selected), paper_reported_gene_count = 113L,
  candidate_file = candidate_file,
  candidate_sha256 = sha256(file.path("Files_For_Manuscript/gene_interactions_reference/paper_input_audit",
                                      candidate_file)),
  input_sha256 = as.list(setNames(vapply(input_paths, sha256, character(1)), input_paths)),
  article = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11164910/",
  source = list(notebook = "Scripts_For_Manuscript/Figure_2.ipynb",
                 de_cell_zero_based = 5L, candidate_cell_zero_based = 9L),
  filter = "Training-only DESeq2 ~tb; positive versus negative; results(alpha=0.05); padj < 0.05 and baseMean > 100; order by padj.",
  n_training_samples = ncol(raw_counts), n_input_genes = nrow(raw_counts),
  n_positive = sum(samples$tb == "positive"), n_negative = sum(samples$tb == "negative"),
  outlier_replaced_genes = sum(mcols(dds)$replace, na.rm = TRUE),
  all_training_samples_pass_qc = all(samples$passQC),
  selected_genes_missing_annotation = sum(is.na(selected$gene_name)),
  duplicate_selected_ids = anyDuplicated(selected$gene_id),
  duplicate_selected_symbols = anyDuplicated(selected$gene_name),
  selected_gene_types = as.list(table(selected$gene_type)),
  selected_blacklist_overlap = sum(sub("\\..*$", "", selected$gene_id) %in% blacklist$V2),
  padj_exactly_cutoff = sum(result$padj == 0.05, na.rm = TRUE),
  base_mean_exactly_cutoff = sum(result$baseMean == 100, na.rm = TRUE),
  same_ids_as_existing_142_reference = setequal(selected$gene_id, previous$gene_id),
  max_absolute_statistic_difference_from_existing_reference = as.list(statistic_difference),
  unused_default_results_candidate_count = sum(!is.na(unused_default$padj) &
                                                unused_default$padj < 0.05 & unused_default$baseMean > 100),
  conclusion = paste("The source filter reproduces 142 genes even with the article's stated R and DESeq2 versions.",
                     "This is not an identified 113-gene paper input list. No arbitrary truncation or extra filter",
                     "has been applied. The exact 113-gene benchmark arm remains unavailable."),
  protocol_note = paste("DE is fitted to raw training counts before VST in the original notebook.",
                        "Changing prediction normalization therefore cannot resolve this candidate-count discrepancy.",
                        "The unused default results(dds) call is reported only as a provenance check.")
)
jsonlite::write_json(audit, file.path(output_dir, "audit.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = NA)
message("Audited ", nrow(selected), " candidates; paper-reported count: 113.")
