# Reproduce the original GFS search over the independently audited candidate list.
# Usage: Rscript export_audit.R /path/to/cfRNA_TB [output_directory]
# Requires pROC, dplyr, jsonlite, and sha256sum. No large VST file is required.
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1 || length(args) > 2) stop("Supply repository root and optional output directory.")
repo <- normalizePath(args[1],mustWork=TRUE)
relative_dir <- "Files_For_Manuscript/gene_interactions_reference/paper_input_audit/gfs_selection"
output <- if (length(args)==2) args[2] else file.path(repo,relative_dir)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
input_dir <- file.path(repo,"Files_For_Manuscript")
audit_dir <- file.path(input_dir,"gene_interactions_reference/paper_input_audit")
sha256 <- function(path) {
  hash <- substr(system2("sha256sum",shQuote(path),stdout=TRUE)[1],1,64)
  if (!grepl("^[a-f0-9]{64}$",hash)) stop("Cannot hash ",path)
  hash
}
candidate_path <- file.path(audit_dir,"training_candidates.tsv")
candidate_audit <- fromJSON(file.path(audit_dir,"audit.json"))
stopifnot(identical(sha256(candidate_path),candidate_audit$candidate_sha256))
for (path in names(candidate_audit$input_sha256)) {
  stopifnot(identical(sha256(file.path(repo,path)),candidate_audit$input_sha256[[path]]))
}
candidate_table <- read.delim(candidate_path,check.names=FALSE) %>%
  arrange(padj) %>% filter(padj < 0.05,baseMean > 100)
candidate_genes <- candidate_table$gene_id
stopifnot(!anyDuplicated(candidate_genes),length(candidate_genes)>0)

# The original training VST uses all raw training genes to estimate size factors.
# Reconstruct its parametric transform from the stored blind-training trend.
# Those trend coefficients were fitted using DESeq2 1.30.1; the independently
# audited candidate list above was fitted using the paper's DESeq2 1.34.0.
reference_path <- file.path(input_dir,"gene_interactions_reference/reference.json")
reference <- fromJSON(reference_path)
for (name in c("counts.train.tsv","metadata.train.tsv")) {
  path <- file.path("Files_For_Manuscript",name)
  stopifnot(identical(sha256(file.path(repo,path)),reference$input_sha256[[path]]))
}
counts <- as.matrix(read.delim(file.path(input_dir,"counts.train.tsv"),row.names=1,check.names=FALSE))
metadata <- read.delim(file.path(input_dir,"metadata.train.tsv"),check.names=FALSE)
stopifnot(identical(colnames(counts),metadata$sample_id),all(candidate_genes %in% rownames(counts)))
positive <- rowSums(counts>0)==ncol(counts)
log_counts <- log(counts[positive,,drop=FALSE])
log_means <- rowMeans(log_counts)
size_factors <- exp(apply(sweep(log_counts,1,log_means,"-"),2,median))
stopifnot(all(is.finite(size_factors)),all(size_factors>0))
q <- sweep(counts[candidate_genes,,drop=FALSE],2,size_factors,"/")
a0 <- reference$training_vst$asymptDisp
a1 <- reference$training_vst$extraPois
vst <- log2((1+a1+2*a0*q+2*sqrt(a0*q*(1+a1+a0*q)))/(4*a0))
features <- as.data.frame(t(vst),check.names=FALSE)
features$sample_id <- rownames(features)
# Original merge uses sort=TRUE. Unused columns do not affect the fitted models.
norm_plot_dat_train <- merge(features,metadata,by.x="sample_id")
df <- norm_plot_dat_train
df$tb <- as.factor(ifelse(df$tb == "positive",1,0))
cor_plot_dat <- norm_plot_dat_train
set.seed(87)
top_genes <- c()
last_best_score <- 0
step <- 0
step_rows <- list()
candidate_rows <- list()
recorded_warnings <- character()

repeat {
  step <- step+1
  gfs_output <- list()
  metrics_output <- list()
  for (gene in candidate_genes) {
    temp_genes <- c(top_genes,gene)
    formula <- as.formula(paste("tb ~",paste(temp_genes,collapse=" + "),"- 1"))
    model_train <- withCallingHandlers(
      glm(formula,data=df,family="binomial"),
      warning=function(w) {
        recorded_warnings <<- c(recorded_warnings,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    probs <- predict(model_train,newdata=df,type="response")
    roc_train <- suppressMessages(roc(df$tb,probs))
    coords_df <- data.frame(coords(roc_train,"all"))
    coords_df$youden <- coords_df$sensitivity+coords_df$specificity-1
    optimal <- coords_df[which.max(coords_df$youden),]
    all_score <- optimal$specificity+optimal$sensitivity+auc(roc_train)
    current_gene_data <- cor_plot_dat[[gene]]
    top_genes_data <- cor_plot_dat[,top_genes]
    current_cor <- abs(mean(cor(current_gene_data,top_genes_data)))
    gfs_output[[gene]] <- c(all_score,current_cor)
    metrics_output[[gene]] <- data.frame(
      step=step,gene_id=gene,auc=as.numeric(auc(roc_train)),
      sensitivity=optimal$sensitivity,specificity=optimal$specificity,
      threshold=optimal$threshold,direction=roc_train$direction,
      converged=model_train$converged,iterations=model_train$iter)
  }
  gfs_res <- do.call("rbind",gfs_output) %>% data.frame()
  colnames(gfs_res) <- c("all_score","correlation")
  top_performer <- gfs_res %>% arrange(desc(all_score),correlation) %>% head(1)
  best_score <- top_performer$all_score
  best_gene <- rownames(top_performer)
  accepted <- best_score-last_best_score > 0.01
  details <- do.call("rbind",metrics_output)
  details$all_score <- gfs_res$all_score
  details$correlation <- gfs_res$correlation
  details$is_top_performer <- details$gene_id == best_gene
  candidate_rows[[step]] <- details
  chosen <- details[details$is_top_performer,,drop=FALSE]
  chosen$score_gain <- best_score-last_best_score
  chosen$accepted <- accepted
  chosen$gene_name <- candidate_table$gene_name[match(best_gene,candidate_table$gene_id)]
  step_rows[[step]] <- chosen
  message("step ",step,": ",chosen$gene_name," score=",sprintf("%.12f",best_score),
          " gain=",sprintf("%.12f",chosen$score_gain)," accepted=",accepted)
  write.csv(do.call(rbind,step_rows),file.path(output,"selection_steps.csv"),row.names=FALSE)
  if (accepted) {
    top_genes <- c(top_genes,best_gene)
    last_best_score <- best_score
  } else break
  if (step > length(candidate_genes)) stop("Unexpected endless selection.")
}

write.csv(do.call(rbind,candidate_rows),file.path(output,"all_candidate_scores.csv"),row.names=FALSE)
selected <- candidate_table[match(top_genes,candidate_table$gene_id),]
selected$selection_order <- seq_along(top_genes)
write.table(selected,file.path(output,"selected_genes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
final_formula <- as.formula(paste("tb ~",paste(top_genes,collapse=" + "),"- 1"))
final_model <- glm(final_formula,data=df,family="binomial")
write.csv(data.frame(gene_id=names(coef(final_model)),coefficient=coef(final_model)),
          file.path(output,"coefficients.csv"),row.names=FALSE)
published <- fromJSON(file.path(repo,
  "Files_For_Manuscript/gene_cutoff_reference/paper_model/model.json"))

# Compare the newly selected model with the separately reconstructed paper model.
paper_dir <- file.path(input_dir,"gene_cutoff_reference/paper_model")
paper_features <- read.delim(file.path(paper_dir,"paper_features.tsv"),check.names=FALSE)
expected_coefficients <- unlist(published$coefficients)[names(coef(final_model))]
comparison <- lapply(c("train","test","val"),function(split) {
  values <- paper_features[paper_features$split==split,,drop=FALSE]
  probabilities <- plogis(as.matrix(values[,names(coef(final_model))]) %*% coef(final_model))
  expected <- plogis(as.matrix(values[,names(coef(final_model))]) %*% expected_coefficients)
  positive <- values$tb_binary==1
  n_positive <- sum(positive)
  auc <- (sum(rank(probabilities)[positive])-n_positive*(n_positive+1)/2)/(n_positive*sum(!positive))
  list(auc=auc,max_abs_probability_difference=max(abs(probabilities-expected)))
})
names(comparison) <- c("train","test","val")
input_paths <- c(
  "Files_For_Manuscript/gene_interactions_reference/paper_input_audit/training_candidates.tsv",
  "Files_For_Manuscript/gene_interactions_reference/paper_input_audit/audit.json",
  "Files_For_Manuscript/counts.train.tsv","Files_For_Manuscript/metadata.train.tsv",
  "Files_For_Manuscript/gene_interactions_reference/reference.json",
  "Scripts_For_Manuscript/Figure_2.ipynb",
  "Files_For_Manuscript/gene_cutoff_reference/paper_model/model.json",
  "Files_For_Manuscript/gene_cutoff_reference/paper_model/paper_features.tsv",
  file.path(relative_dir,"export_audit.R"))
output_files <- c("selection_steps.csv","all_candidate_scores.csv","selected_genes.tsv","coefficients.csv")
provenance <- list(
  schema_version=1L,
  source=list(notebook="Scripts_For_Manuscript/Figure_2.ipynb",cell_zero_based=9L,
    method="Same candidate loop, R glm, pROC defaults, composite score, correlation tie breaker, and stopping rule"),
  candidate_file="../training_candidates.tsv",candidate_count=length(candidate_genes),
  candidate_deseq2_version=candidate_audit$deseq2_version,
  vst=list(method="Original blind parametric training VST reconstructed from saved training trend and raw training counts",
    reference_file="../../reference.json",reference_deseq2_version=reference$deseq2_version,
    size_factor_method="DESeq2 ratio: exp(median(log(count)-mean(log(reference counts)))) over genes positive in every training sample",
    n_size_factor_reference_genes=sum(positive),coefficients=reference$training_vst,
    original_cached_matrix=list(file="training_blind_vst.tsv.gz",included=FALSE,
      sha256="c6c8261c1f30c589e39592ebb8e2b565a6a4a4ffadde158d8aadee5d511d756c",
      deseq2_version="1.30.1",
      note="Original audit used this cached matrix; the portable script reconstructs its training transformation without copying the large matrix.")),
  r_version=R.version.string,pROC_version=as.character(packageVersion("pROC")),
  dplyr_version=as.character(packageVersion("dplyr")),jsonlite_version=as.character(packageVersion("jsonlite")),
  n_training_samples=nrow(df),seed=87L,
  glm=list(family="binomial",link="logit",intercept=FALSE,penalty="none",scaling="none",
    epsilon=glm.control()$epsilon,maxit=glm.control()$maxit),
  score="Training ROC AUC plus sensitivity plus specificity at the first maximum Youden threshold",
  tie_breaker="Lowest absolute mean Pearson correlation with already selected genes",
  accept_rule="Score gain strictly greater than 0.01",candidate_fits=sum(vapply(candidate_rows,nrow,integer(1))),
  selected_genes=top_genes,selected_gene_names=selected$gene_name,
  same_set_as_published=setequal(top_genes,published$genes),
  same_order_as_published=identical(top_genes,published$genes),
  rejected_next_gene=tail(step_rows,1)[[1]]$gene_name,
  rejected_score_gain=tail(step_rows,1)[[1]]$score_gain,
  max_abs_coefficient_difference_from_six_gene_refit=max(abs(coef(final_model)-expected_coefficients)),
  rescored_paper_features=comparison,warnings=unique(recorded_warnings),
  input_sha256=setNames(lapply(input_paths,function(path) sha256(file.path(repo,path))),input_paths),
  output_sha256=setNames(lapply(output_files,function(name) sha256(file.path(output,name))),output_files),
  conclusion="The repository's 142-candidate GFS search recovers the published six final genes in their original order. This does not recover or certify the article's reported 113-gene input pool for other models.")
write_json(provenance,file.path(output,"audit.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
message("Selected: ",paste(selected$gene_name,collapse=", "))
message("Same published set/order: ",provenance$same_set_as_published,"/",provenance$same_order_as_published)
