#!/usr/bin/env Rscript

# Differential expression per study_id using Salmon quant.sf -> tximport -> DESeq2.
# The input metadata keeps the existing columns: sample_id, condition, study_id.
# Each row must represent one independent biological replicate (one SRA run per sample).

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else getwd()
metadata_file <- Sys.getenv("PROMETA_METADATA", unset = file.path(project_dir, "metadata_geral.csv"))
quant_dir <- file.path(project_dir, "quant")
tx2gene_file <- file.path(project_dir, "reference", "tx2gene.csv")
output_dir <- file.path(project_dir, "differential_expression")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(metadata_file)) stop("Não encontrei metadata_geral.csv em: ", metadata_file)
if (!file.exists(tx2gene_file)) stop("Não encontrei reference/tx2gene.csv. Arquivo deve conter colunas TXNAME e GENEID.")
meta <- read_csv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
required <- c("sample_id", "condition", "study_id")
if (!all(required %in% names(meta))) stop("Metadados precisam conter: ", paste(required, collapse = ", "))
meta <- meta %>%
  mutate(condition = tolower(trimws(condition)), study_id = trimws(study_id), sample_id = trimws(sample_id))
if (anyDuplicated(meta$sample_id)) stop("sample_id deve ser único (uma linha por amostra/run).")
if (any(!meta$condition %in% c("control", "treatment"))) stop("condition deve ser control ou treatment.")
tx2gene <- read_csv(tx2gene_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
if (!all(c("TXNAME", "GENEID") %in% names(tx2gene))) stop("tx2gene.csv deve ter colunas TXNAME e GENEID.")

all_results <- list()
for (study in unique(meta$study_id)) {
  sub <- meta %>% filter(study_id == study)
  sub$condition <- factor(sub$condition, levels = c("control", "treatment"))
  n_by_condition <- table(sub$condition)
  if (length(n_by_condition) != 2 || any(n_by_condition < 2)) {
    warning("Pulando ", study, ": são necessárias pelo menos 2 amostras biológicas em cada condição.")
    next
  }
  files <- file.path(quant_dir, sub$sample_id, "quant.sf")
  names(files) <- sub$sample_id
  if (any(!file.exists(files))) {
    warning("Pulando ", study, ": quant.sf ausente para ", paste(sub$sample_id[!file.exists(files)], collapse = ", "))
    next
  }
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")],
                  countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
  rownames(sub) <- sub$sample_id
  dds <- DESeqDataSetFromTximport(txi, colData = sub, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "treatment", "control"), alpha = 0.05)
  result_df <- as.data.frame(res)
  result_df$gene_id <- rownames(result_df)
  result_df$study_id <- study
  result_df <- result_df %>% arrange(padj)
  safe_study <- gsub("[^A-Za-z0-9._-]+", "_", study)
  write_csv(result_df, file.path(output_dir, paste0(safe_study, "_DESeq2_treatment_vs_control.csv")))
  all_results[[study]] <- result_df
  message("Concluído: ", study, " (", nrow(result_df), " genes testados)")
}
if (!length(all_results)) stop("Nenhuma comparação foi analisada. Confira réplicas, quant.sf e metadados.")
write_csv(bind_rows(all_results), file.path(output_dir, "DESeq2_todos_estudos.csv"))
message("Resultados salvos em ", output_dir)
