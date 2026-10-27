############################## IMPORT LIBRARIES ############################### 
library(readr)
library(biomaRt)
library(tibble)
library(DESeq2)
library(survival)
library(survminer)
library(timeROC)
library(glmnet)
library(mclust)
library(writexl)
library(openxlsx)
library(ggVennDiagram)
library(ggplot2)
library(grid)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(msigdbr)
library(purrr)
library(cowplot)
library(ggplotify)
library(ggrepel)
library(patchwork)
library(tidyr)
library(dplyr)
select <- dplyr::select

################################## LOAD DATA ###################################

df_phenotype <- read_tsv("Phenotype_GBM.tsv")
df_gene_count <- read_tsv("Gene_GBM.tsv")

################################# FILTER DATA ##################################

# keep original data frames
df_phenotype_original <- df_phenotype 
df_count_original <- df_gene_count 

# filter and arrange to specifications 
df_phenotype <- df_phenotype %>% 
  filter(tumor_descriptor.samples == "Primary", tissue_or_organ_of_origin.diagnoses == "Brain, NOS", tissue_type.samples == "Tumor", prior_treatment.diagnoses == "No", vital_status.demographic == "Dead") %>%
  arrange(age_at_index.demographic) 

# filter phenotype to only keep samples which also have a gene_count 
samples_in_df_gene_count <- colnames(df_gene_count)[-1]  # remove if first column is not sample names
df_phenotype <- df_phenotype[df_phenotype$sample %in% samples_in_df_gene_count, ]

# keep one sample per patient (one patient has two vials: keep the first, 01A)
df_phenotype <- df_phenotype %>%
  arrange(sample) %>%
  distinct(case_id, .keep_all = TRUE) %>%
  arrange(age_at_index.demographic)

stopifnot(!anyDuplicated(df_phenotype$case_id))
nrow(df_phenotype)   # should be 122

########################## CATEGORIZE DATA #########################
# categorize 1: Age 
    # <55 = Young
    # >= 55 = Old
# further categorize into 2 categories:
    # =< 1 year = SS
    # >1 = LS

df_phenotype <- df_phenotype %>%
  mutate(
    age_groups = if_else(age_at_index.demographic < 55, "Young", "Old"),
    survival_groups = if_else(days_to_death.demographic <= 365, "SS", "LS")
  )

# test table 
print(
  test_table <- table(df_phenotype$age_groups, df_phenotype$survival_groups)
)

############### MATCH SAMPLES KEPT IN PHENOTYPE TO GENE_COUNT###############

samples_to_keep <- df_phenotype$sample
df_gene_count <- df_gene_count[, c(TRUE, colnames(df_gene_count)[-1] %in% samples_to_keep)]

######### EXPRESSION VALUE CONVERSION FROM log2(count+1) TO (count) ############
# keep gene ID col separate
gene_ids <- df_gene_count[, 1]

# convert only numeric val
expr <- df_gene_count[, -1] 
expr_conversion <- (2^expr) -1

# combine gene IDs to converted vals
df_gene_count <- cbind(gene_ids, expr_conversion)

############## USE Ensembl TO CONVERT Ensembl_ID TO GENE ID ####################
################## AND FILTER FOR PROTEIN CODING GENES #########################

# 1. Connect to Ensembl for human genome
# 1. Connect to Ensembl for human genome (Ensembl 115 archive - US East mirror retired)
mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host    = "https://e115.ensembl.org"
)

# 2. Remove version nr (i.e .15)
df_gene_count$ensembl_gene_id <- sub("\\..*$", "", df_gene_count$Ensembl_ID)

# 3. Get gene biotype and gene symbol from BioMart package
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = unique(df_gene_count$ensembl_gene_id),
  mart = mart
)

# 4. Add gene symbol column
df_gene_count$gene_symbol <- annot$external_gene_name[
  match(df_gene_count$ensembl_gene_id, annot$ensembl_gene_id)
]

# 5. Add gene biotype column
df_gene_count$gene_biotype <- annot$gene_biotype[
  match(df_gene_count$ensembl_gene_id, annot$ensembl_gene_id)
]

# 6. Keep only protein-coding genes
df_gene_count <- df_gene_count[df_gene_count$gene_biotype == "protein_coding", ]

# 7. Reorder column
df_gene_count <- df_gene_count[, c(
  "Ensembl_ID",
  "ensembl_gene_id",
  "gene_symbol",
  "gene_biotype",
  setdiff(colnames(df_gene_count),
          c("Ensembl_ID", "ensembl_gene_id", "gene_symbol", "gene_biotype"))
)]

# 8. Remove row where gene symbol is na/empty & duplicates from gene_symbol col 
# (note: keep copy with largest count)
df_gene_count <- df_gene_count %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = rowSums(across(where(is.numeric)), na.rm = TRUE), n = 1, with_ties = FALSE) %>%
  ungroup()

# 9. Save object
saveRDS(df_gene_count, "df_gene_count.rds")

# 10. Load ENSEMBL 115 Annotation (from saved object)
# Annotation (symbol, biotype), protein-coding filter and duplicate removal were
# done once against Ensembl 115 and saved, so BioMart is not queried again.
df_gene_count <- readRDS("df_gene_count.rds")

# drop the removed vial (01B) from the cached table
df_gene_count <- df_gene_count[, c(colnames(df_gene_count)[1:4], samples_to_keep)]

stopifnot(
  nrow(df_gene_count) == 19414,
  ncol(df_gene_count) == 4 + length(samples_to_keep),
  all(df_gene_count$gene_biotype == "protein_coding"),
  !anyDuplicated(df_gene_count$gene_symbol)
)

#################### RUNNING DESeq2 (FOR GENE EXPRESSION ANALYSIS) ######################

# 1. Update df_gene_count to only keep samples
gene_matrix <- df_gene_count %>%
  select(-Ensembl_ID, -ensembl_gene_id, -gene_biotype)

gene_matrix <- gene_matrix %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

gene_matrix <- round(gene_matrix)
mode(gene_matrix) <- "integer"

# 2. set refrence groups Old and LS
phenotype_matrix <- df_phenotype %>%
  select(sample, age_groups, survival_groups) # only keep necessary columns

phenotype_matrix <- phenotype_matrix %>%
  mutate(
    age_groups = factor(age_groups, levels = c("Old", "Young")),
    survival_groups = factor(survival_groups, levels = c("LS", "SS"))
  )
rownames(phenotype_matrix) <- phenotype_matrix$sample

# match order of column and row in gene_matrix and phenotype_matrix
common_samples <- intersect(colnames(gene_matrix), rownames(phenotype_matrix))
gene_matrix <- gene_matrix[, common_samples]
phenotype_matrix <- phenotype_matrix[common_samples, ]

  
# 3. Setup DESeq2 model
dds <- DESeqDataSetFromMatrix(
  countData = gene_matrix,
  colData = phenotype_matrix,
  design = ~ age_groups * survival_groups
)

# 4. Run DESeq2 model
dds <- DESeq(dds)
print(
  resultsNames(dds)
)

# 5. Run analysis for 4 different groups
# 5.1: Old: SS vs LS (survival effect in Old group)
# 5.2: Young: SS vs LS (survival effect in Young group)
# 5.3: SS: Young vs Old (age effect in SS)
# 5.4: LS: Young vs Old (age effect in LS)

# 5.1: Old: SS vs LS
# Does survival affect gene expression in Old patients
results_Old_survival <- results(dds, 
  name = "survival_groups_SS_vs_LS",
  independentFiltering = FALSE
  )

# 5.2: Young: SS vs LS (survival effect in Young group)
# Does survival affect gene expression in Young patients?
results_Young_survival <- results(dds,
  contrast = list(
    c("survival_groups_SS_vs_LS",
      "age_groupsYoung.survival_groupsSS"),
    character(0)
  ),
  independentFiltering = FALSE
)

# 5.3: SS: Young vs Old (age effect in SS group)
# Does age affect gene expression in SS patients?
results_SS_age <- results(dds,
  contrast = list(
    c("age_groups_Young_vs_Old", "age_groupsYoung.survival_groupsSS"),
    character(0)
  ),
  independentFiltering = FALSE
)

# 5.4: LS: Young vs Old (age effect in LS group)
# Does age affect gene expression in LS patients?
results_LS_age <- results(dds, 
  name = "age_groups_Young_vs_Old",
  independentFiltering = FALSE
)

################ CREATE EXCEL SHEET WITH DESeq2 GENE RESULTS ###################
# Make function to create excel book with 3 sheets for significant genes:
# Use padj <= 0.05 and |log2fc| > 1
# 1. Top_genes
# 2. Upregulated_genes
# 3. Downregulated_genes

dseq_data_excel <- function(results_Old_survival, 
                            results_Young_survival, 
                            results_SS_age,       
                            results_LS_age,
                            file_name,
                            padj_cutoff = 0.05,
                            log2fc_cutoff = 1) {
  
  process_results <- function(results, label) { 
    df <- as.data.frame(results)
    df$gene_symbol <- rownames(df)
    df <- df[, c("gene_symbol", setdiff(names(df), "gene_symbol"))]
    df <- df[!is.na(df$padj), ]
    
    # Add name to clarify where results came from
    df$comparison <- label
    
    return(df)
  }
  
  # Process each result
  df_Old <- process_results(results_Old_survival, "Old_SS_vs_LS")
  df_Young <- process_results(results_Young_survival, "Young_SS_vs_LS")
  df_SS <- process_results(results_SS_age, "SS_Young_vs_Old")
  df_LS <- process_results(results_LS_age, "LS_Young_vs_Old")
  
  
  # Combine all
  combined <- rbind(df_Old, df_Young, df_SS, df_LS)
  
  # Top Genes
  top_genes <- combined[
    combined$padj <= padj_cutoff &
      abs(combined$log2FoldChange) >= log2fc_cutoff,
  ]
  top_genes <- top_genes[order(top_genes$padj), ]
  
  # Upregulated Genes
  upregulated_genes <- combined[
    combined$padj <= padj_cutoff &
      combined$log2FoldChange >= log2fc_cutoff,
  ]
  upregulated_genes <- upregulated_genes[order(upregulated_genes$padj), ]
  
  # Downregulated Genes
  downregulated_genes <- combined[
    combined$padj <= padj_cutoff &
      combined$log2FoldChange <= -log2fc_cutoff,
  ]
  downregulated_genes <- downregulated_genes[order(downregulated_genes$padj), ]
  
  # Create Excel Workbook
  excel_workbook <- createWorkbook()
  
  addWorksheet(excel_workbook, "Top_genes")
  addWorksheet(excel_workbook, "Upregulated_genes")
  addWorksheet(excel_workbook, "Downregulated_genes")
  
  writeData(excel_workbook, "Top_genes", top_genes)
  writeData(excel_workbook, "Upregulated_genes", upregulated_genes)
  writeData(excel_workbook, "Downregulated_genes", downregulated_genes)
  
  saveWorkbook(excel_workbook, file_name, overwrite = TRUE)
}

# Run function to create workbook
dseq_data_excel(
  results_Old_survival,
  results_Young_survival,
  results_SS_age,       
  results_LS_age,
  "Significant_DEG.xlsx"
)

############## DESeq2 RESULTS COUNT FOR DIFFERENT log2FC VALUES ################

# list of results
results_list <- list(
  "Old SS vs Old LS"       = results_Old_survival,
  "Young SS vs Young LS"   = results_Young_survival,
  "Young SS vs Old SS"     = results_SS_age,
  "Young LS vs Old LS"     = results_LS_age
)

# Create function to count significant genes
count_signif_genes <- function(results, padj_cutoff = 0.05, log2fc_cutoff) {
  df <- as.data.frame(results)
  df <- df[!is.na(df$padj), ]
  
  signif <- df[df$padj <= padj_cutoff & abs(df$log2FoldChange) >= log2fc_cutoff, ]
  up  <- df[df$padj <= padj_cutoff & df$log2FoldChange >= log2fc_cutoff, ]
  down <- df[df$padj <= padj_cutoff & df$log2FoldChange <= -log2fc_cutoff, ]
  
  return(c(
    Total = nrow(signif),
    Up    = nrow(up),
    Down  = nrow(down)
  ))
}

# Build summary table for 4 log2FC cutoffs (0, 0.5, 1, 1,5)
log2fc_cutoffs <- c(0, 0.5, 1, 1.5)

summary_rows <- list()

for (cutoff in log2fc_cutoffs) {
  for (comparison in names(results_list)) {
    counts <- count_signif_genes(results_list[[comparison]], log2fc_cutoff = cutoff)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Comparison    = comparison,
      padj_cutoff   = 0.05,
      log2FC_cutoff = cutoff,
      Total_Signif     = counts["Total"],
      Upregulated   = counts["Up"],
      Downregulated = counts["Down"]
    )
  }
}

summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
print(summary_table)

# create excel workbook
excel_workbook <- createWorkbook()
addWorksheet(excel_workbook, "Summary")
writeData(excel_workbook, "Summary", summary_table)
saveWorkbook(excel_workbook, "DESeq2_DIFFERENT_LOG2FC _Table.xlsx", overwrite = TRUE)

################### VISUAL REPRESENTATION OF DESeq2 RESULTS ####################
# 1.1 Volcano Plot Old SS vs Old LS
# 1.2 Volcano Plot Young SS vs Young LS
# 1.3 Volcano Plot Young SS vs Old SS 
# 1.4 Volcano Plot Young LS vs Old LS 
# 2. Venn Diagram
# 2.1 Significant genes in age groups
# 2.2 Significant genes across survival groups
# 2.3 Comparison of LS, SS, Old, Young

# create folder for plots
output_dir_volcano <- "volcano_plots"
if (!dir.exists(output_dir_volcano)) dir.create(output_dir_volcano)

# 1. create function to draw volcano plot with the same parameters
plot_volcano <- function(res, title, n_top = 20) {
  
  res <- res[!is.na(res$pvalue), ]
  
  # top n_top significant genes for this comparison, ranked by padj
  sig_genes <- rownames(res)[res$padj < 0.05 & abs(res$log2FoldChange) >= 1]
  top_genes <- head(sig_genes[order(res[sig_genes, "padj"])], n_top)
  
  EnhancedVolcano(res,
                  lab       = rownames(res),
                  x         = "log2FoldChange",
                  y         = "padj",
                  selectLab = unique(c(top_genes)),
                  pCutoff   = 0.05,
                  FCcutoff  = 1,
                  title     = title,
                  ylab = bquote(-log[10]~italic(padj)),
                  subtitle  = NULL,
                  caption   = paste0("padj < 0.05  |  |log2FC| > 1  |  n = ", nrow(res)),
                  drawConnectors  = TRUE,
                  maxoverlapsConnectors = 20
  )
}

# Generate all 5 plots
VP1 <- plot_volcano(results_Old_survival,   "Old: SS vs LS")
VP2 <- plot_volcano(results_Young_survival, "Young: SS vs LS")
VP3 <- plot_volcano(results_SS_age,         "SS: Young vs Old")
VP4 <- plot_volcano(results_LS_age,         "LS: Young vs Old")

print(VP1)
print(VP2)
print(VP3)
print(VP4)

# Save into output_dir_volcano, same naming pattern as before
ggsave(filename = file.path(output_dir_volcano, "volcano_Old_survival.pdf"),   plot = VP1, width = 9, height = 8, dpi = 300)
ggsave(filename = file.path(output_dir_volcano, "volcano_Young_survival.pdf"), plot = VP2, width = 9, height = 8, dpi = 300)
ggsave(filename = file.path(output_dir_volcano, "volcano_SS_age.pdf"),         plot = VP3, width = 9, height = 8, dpi = 300)
ggsave(filename = file.path(output_dir_volcano, "volcano_LS_age.pdf"),         plot = VP4, width = 9, height = 8, dpi = 300)

cat("volcano_Old_survival.pdf written:",   file.exists(file.path(output_dir_volcano, "volcano_Old_survival.pdf")), "\n")
cat("volcano_Young_survival.pdf written:", file.exists(file.path(output_dir_volcano, "volcano_Young_survival.pdf")), "\n")
cat("volcano_SS_age.pdf written:",         file.exists(file.path(output_dir_volcano, "volcano_SS_age.pdf")), "\n")
cat("volcano_LS_age.pdf written:",         file.exists(file.path(output_dir_volcano, "volcano_LS_age.pdf")), "\n")

# 2. Venn Diagrams
# create folder for diagrams 
output_dir_venn <- "venn_diagrams"
if (!dir.exists(output_dir_venn)) dir.create(output_dir_venn)

# Create function to extract significant genes
get_signif_genes <- function(results, padj_cutoff = 0.05, log2fc_cutoff = 1) {
  df <- as.data.frame(results)
  df$gene_symbol <- rownames(df)
  df <- df[!is.na(df$padj), ]
  df$gene_symbol[df$padj <= padj_cutoff & abs(df$log2FoldChange) >= log2fc_cutoff]
}

# Extract significant genes for 5.1 - 5.4 in DESeq2 analysis
signif_Old_survival   <- get_signif_genes(results_Old_survival)
signif_Young_survival <- get_signif_genes(results_Young_survival)
signif_SS_age         <- get_signif_genes(results_SS_age)
signif_LS_age         <- get_signif_genes(results_LS_age)

# Draw set names above/below the shapes instead of beside them
add_set_labels <- function(p, sets, size = 4, pad = 0.10, center = length(sets) <= 2) {
  d    <- process_data(Venn(sets))
  lab  <- as.data.frame(venn_setlabel(d))
  edge <- as.data.frame(venn_setedge(d))
  
  ymin <- min(edge$Y); ymax <- max(edge$Y); h <- ymax - ymin
  
  ord    <- order(-lab$Y, seq_len(nrow(lab)))
  n_top  <- ceiling(nrow(lab) / 2)
  is_top <- seq_len(nrow(lab)) %in% ord[seq_len(n_top)]
  lab$Y2 <- ifelse(is_top, ymax + pad * h, ymin - pad * h)
  
  if (center) lab$X <- mean(range(edge$X))   # centre two-set labels
  
  p + geom_text(data = lab, aes(x = X, y = Y2, label = name),
                size = size, fontface = "bold", inherit.aes = FALSE)
}

shrink_venn <- function(p, mult = 0.15) {
  p +
    scale_x_continuous(expand = expansion(mult = mult)) +
    scale_y_continuous(expand = expansion(mult = mult))
}

# 2.1
sets_1 <- list("Old: SS vs LS" = signif_Old_survival,
               "Young: SS vs LS" = signif_Young_survival)
venn_1 <- ggVennDiagram(sets_1, set_size = 0, label_size = 4) +
  scale_fill_gradient(low = "lightgreen", high = "steelblue") +
  ggtitle("Survival effect within each age group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
venn_1 <- shrink_venn(add_set_labels(venn_1, sets_1), mult = 0.2)

# 2.2
sets_2 <- list("SS: Young vs Old" = signif_SS_age,
               "LS: Young vs Old" = signif_LS_age)
venn_2 <- ggVennDiagram(sets_2, set_size = 0, label_size = 4) +
  scale_fill_gradient(low = "pink", high = "coral") +
  ggtitle("Age effect within each survival group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
venn_2 <- shrink_venn(add_set_labels(venn_2, sets_2), mult = 0.2)

# 2.3
sets_3 <- list("Old: SS vs LS"    = signif_Old_survival,
               "Young: SS vs LS"  = signif_Young_survival,
               "SS: Young vs Old" = signif_SS_age,
               "LS: Young vs Old" = signif_LS_age)
venn_3 <- ggVennDiagram(sets_3, set_size = 0, label_size = 3.5) +
  scale_fill_gradient(low = "blue", high = "purple") +
  ggtitle("All Comparisons") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
venn_3 <- shrink_venn(add_set_labels(venn_3, sets_3, size = 3.5), mult = 0.2)

print(venn_1); print(venn_2); print(venn_3)

# alter layout of each plot to make them readable
# venn_1 and venn_2 (each about 3 in wide in your layout)
ggsave(file.path(output_dir_venn, "venn_survival_effect_young_vs_old.pdf"), venn_1, width = 4.5, height = 4.5)
ggsave(file.path(output_dir_venn, "venn_age_effect_ss_vs_ls.pdf"),venn_2, width = 4.5, height = 4.5)
# venn_3 (wider, about 5 in on the page)
ggsave(file.path(output_dir_venn, "venn_all_comparisons.pdf"), venn_3, width = 6.5, height = 5.5)

cat("venn_survival_effect_young_vs_old.pdf written:", file.exists(file.path(output_dir_venn, "venn_survival_effect_young_vs_old.pdf")), "\n")
cat("venn_age_effect_ss_vs_ls.pdf written:",          file.exists(file.path(output_dir_venn, "venn_age_effect_ss_vs_ls.pdf")), "\n")
cat("venn_all_comparisons.pdf written:",               file.exists(file.path(output_dir_venn, "venn_all_comparisons.pdf")), "\n")

########################## ENRICHMENT ANALYSIS - ORA (GO) ###########################
# For enrichment analysis using all significant genes: create list to be used as input argument
# All significant genes per comparison (full significant DESeq2 result)
all_gene_sets <- list(
  Old_survival   = signif_Old_survival,
  Young_survival = signif_Young_survival,
  SS_age         = signif_SS_age,
  LS_age         = signif_LS_age
)


# Create function to run GO analysis
run_go_ora <- function(gene_sets, set_type = "all", output_file = NULL) {
  
  if (is.null(output_file)) {
    output_file <- paste0("GO_ORA_", set_type, "_results.xlsx")
  }
  
  # Set background universe
  background_genes <- rownames(as.data.frame(results_Old_survival))
  background_entrez <- mapIds(org.Hs.eg.db,
                              keys = background_genes,
                              column = "ENTREZID",
                              keytype = "SYMBOL",
                              multiVals = "first") %>% na.omit()
  
  wb <- createWorkbook()
  all_results <- list()
  go_objects <- list()
  
  for (comparison_name in names(gene_sets)) {
    genes <- gene_sets[[comparison_name]]
    
    if (length(genes) < 2) {
      message("Skipping ", comparison_name, ": fewer than 2 ", set_type, " genes")
      next
    }
      
    entrez <- mapIds(org.Hs.eg.db,
                     keys = genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first") %>% na.omit()
      
    go <- enrichGO(gene = entrez,
                    universe = background_entrez,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
    
    go_objects[[comparison_name]] <- go  # object for dotplots 
    go_df <- as.data.frame(go) # dataframe for excel
      
    sheet_name <- paste(comparison_name, set_type, "GO", sep = "_")
    sheet_name <- substr(sheet_name, 1, 31)
      
    all_results[[sheet_name]] <- go_df
      
    addWorksheet(wb, sheet_name)
    if (nrow(go_df) > 0) {
      writeData(wb, sheet_name, go_df)
    } else {
      writeData(wb, sheet_name, "No significant GO terms")
    }
  }
  saveWorkbook(wb, output_file, overwrite = TRUE)
  message("Saved: ", output_file)
  
  return(list(tables = all_results, objects = go_objects)) 
}

# Run GO for all_significant
GO_ORA_all_fc1_results <- run_go_ora(all_gene_sets, set_type = "all_fc1")

###################### ENRICHMENT ANALYSIS - ORA (KEGG) ########################

# Create function to run KEGG analysis
run_kegg_ora <- function(gene_sets, set_type = "all", output_file = NULL) {
  
  if (is.null(output_file)) {
    output_file <- paste0("KEGG_ORA_", set_type, "_results.xlsx")
  }
  
  # Set background universe
  background_genes <- rownames(as.data.frame(results_Old_survival))
  background_entrez <- mapIds(org.Hs.eg.db,
                              keys = background_genes,
                              column = "ENTREZID",
                              keytype = "SYMBOL",
                              multiVals = "first") %>% na.omit()
  
  wb <- createWorkbook()
  all_results <- list()
  kegg_objects <- list()
  
  for (comparison_name in names(gene_sets)) {
    genes <- gene_sets[[comparison_name]]
    
    if (length(genes) < 2) {
      message("Skipping ", comparison_name, ": fewer than 2 ", set_type, " genes")
      next
    }
    
    entrez <- mapIds(org.Hs.eg.db,
                     keys = genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first") %>% na.omit()
    
    kegg <- enrichKEGG(gene = entrez,
                       universe = background_entrez,
                       organism = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
    
    kegg_objects[[comparison_name]] <- kegg # object for dotplots
    kegg_df <- as.data.frame(kegg) # dataframe for excel
    
    sheet_name <- paste(comparison_name,set_type, "KEGG", sep = "_")
    sheet_name <- substr(sheet_name, 1, 31)
    
    all_results[[sheet_name]] <- kegg_df
    
    addWorksheet(wb, sheet_name)
    if (nrow(kegg_df) > 0) {
      writeData(wb, sheet_name, kegg_df)
    } else {
      writeData(wb, sheet_name, "No significant KEGG pathways")
    }
  }
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
  message("Saved: ", output_file)
  
  return(list(tables = all_results, objects = kegg_objects)) 
}

# Run KEGG for all_significant
KEGG_ORA_all_fc1_results <- run_kegg_ora(all_gene_sets, set_type = "all_fc1")

########################### VISUALIZATION OF ORA (GO) ##########################

# Create function for GO dotplots (showing top 20 results based on significance per ontology)
plot_GO_dotplots <- function(go_objects, output_dir = "dotplots_GO", showCategory = 20) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  for (comparison_name in names(go_objects)) {
    
    go <- go_objects[[comparison_name]]
    
    if (is.null(go) || nrow(as.data.frame(go)) == 0) {
      message("Skipping ", comparison_name, ": no significant GO terms")
      next
    }
    
    # Split dotplots into ontology categories
    p <- dotplot(go, showCategory = showCategory, split = "ONTOLOGY") +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
      ggtitle(paste(comparison_name, "- GO (BP / MF / CC)"))
    
    ggsave(filename = file.path(output_dir, paste0("dotplot_", comparison_name, ".pdf")),
           plot = p, width = 9, height = 10, dpi = 300)
    
    print(p)
  }
}

# Dot plot for all_significant
plot_GO_dotplots(GO_ORA_all_fc1_results$objects, output_dir = "dotplots_GO_all_fc1")

########################### VISULIZATION OF ORA (KEGG) ##########################

# Create function for KEGG dotplots (showing top 20 results based on padj)
plot_KEGG_dotplots <- function(kegg_objects, output_dir = "dotplots_KEGG", showCategory = 20) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  for (comparison_name in names(kegg_objects)) {
    
    kegg <- kegg_objects[[comparison_name]]
    
    if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) {
      message("Skipping ", comparison_name, ": no significant KEGG pathways")
      next
    }
    
    p <- dotplot(kegg, showCategory = showCategory) +
      ggtitle(paste(comparison_name, "- KEGG Pathways"))
    
    ggsave(filename = file.path(output_dir, paste0("dotplot_", comparison_name, ".pdf")),
           plot = p, width = 9, height = 8, dpi = 300)
    
    print(p)
  }
}

# Dot plot for all_significant
plot_KEGG_dotplots(KEGG_ORA_all_fc1_results$objects, output_dir = "dotplots_KEGG_all_fc1")

################ ENRICHMENT ANALYSIS AGAINST HALLMARK - GSEA ###################

# Step 1: Get Hallmark gene sets from MSigDB (category "H" = Hallmark, 50 curated gene sets)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

hallmark_term2gene <- hallmark_sets %>%
  rename(term = gs_name, gene = gene_symbol)

# Step 2 Build ranked gene lists for GSEA
# Ranked gene list uses the stat column (not log2FC) - keeps every gene, no significance filtering
prepare_gsea_ranks_symbol <- function(results_df) {
  results_df <- as.data.frame(results_df)
  results_df <- results_df %>%
    filter(!is.na(stat)) %>% # use stat coloumn instead of log2FC
    tibble::rownames_to_column("gene") %>%
    distinct(gene, .keep_all = TRUE)
  
  ranks <- results_df$stat
  names(ranks) <- results_df$gene
  ranks <- sort(ranks, decreasing = TRUE)
  return(ranks)
}

# Run function - one ranked list per comparison
ranks_Old_survival    <- prepare_gsea_ranks_symbol(results_Old_survival)
ranks_Young_survival  <- prepare_gsea_ranks_symbol(results_Young_survival)
ranks_SS_age          <- prepare_gsea_ranks_symbol(results_SS_age)
ranks_LS_age          <- prepare_gsea_ranks_symbol(results_LS_age)


# Step 3:  GSEA function: run single comparison against Hallmark 
run_hallmark_gsea_single <- function(ranks, term2gene) {
  gsea_result <- GSEA(geneList = ranks,
                      TERM2GENE = term2gene,
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,  # only keep significant pathways 
                      pAdjustMethod = "BH",
                      seed = TRUE)
  return(gsea_result)
}

# Step 4: Run GSEA across all comparisons and collect results
run_hallmark_gsea_all <- function(term2gene) {
  results <- list(
    Old_survival   = run_hallmark_gsea_single(ranks_Old_survival, term2gene),
    Young_survival = run_hallmark_gsea_single(ranks_Young_survival, term2gene),
    SS_age         = run_hallmark_gsea_single(ranks_SS_age, term2gene),
    LS_age         = run_hallmark_gsea_single(ranks_LS_age, term2gene)
  )
  return(results)
}

hallmark_gsea_all <- run_hallmark_gsea_all(hallmark_term2gene)

# Step 5: Combine all comparisons into one table and write to excel
hallmark_gsea_combined <- purrr::imap_dfr(hallmark_gsea_all, function(res, comparison_name) {
  as.data.frame(res) %>%
    mutate(comparison = comparison_name)
})

export_hallmark_gsea_excel <- function(gsea_combined, file_path) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "Hallmark_GSEA_All")
  writeData(wb, "Hallmark_GSEA_All", gsea_combined)
  
  for (comp in unique(gsea_combined$comparison)) {
    sheet_name <- comp
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, gsea_combined %>% filter(comparison == comp))
  }
  
  saveWorkbook(wb, file_path, overwrite = TRUE)
}

export_hallmark_gsea_excel(hallmark_gsea_combined, "hallmark_gsea_results.xlsx")

######### NES BAR PLOT SUMMARY - HALLMARK GSEA (ONE PLOT PER COMPARISON) ############

plot_hallmark_nes_barplot <- function(gsea_objects, output_dir = "barplots_hallmark",
                                      showCategory = 20) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  for (comparison_name in names(gsea_objects)) {
    
    gsea_result <- gsea_objects[[comparison_name]]
    gsea_df <- as.data.frame(gsea_result)
    
    if (is.null(gsea_result) || nrow(gsea_df) == 0) {
      message("Skipping ", comparison_name, ": no significant Hallmark pathways")
      next
    }
    
    # Keep top N by p.adjust (in case there are many significant pathways)
    plot_df <- gsea_df %>%
      arrange(p.adjust) %>%
      slice_head(n = showCategory) %>%
      mutate(
        pathway = gsub("HALLMARK_", "", ID),
        pathway = gsub("_", " ", pathway),
        direction = ifelse(NES > 0, "Positively Enriched", "Negatively Enriched"),
        pathway = factor(pathway, levels = pathway[order(NES)])  # order by NES for clean bar layout
      )
    
    p <- ggplot(plot_df, aes(x = NES, y = pathway, fill = direction)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
      scale_fill_manual(values = c("Positively Enriched" = "red", "Negatively Enriched" = "blue")) +
      labs(title = paste(comparison_name, "- Hallmark GSEA"),
           x = "Normalized Enrichment Score (NES)", y = NULL, fill = NULL) +
      theme_classic(base_size = 12) +
      theme(legend.position = "top")
    
    ggsave(filename = file.path(output_dir, paste0("barplot_", comparison_name, ".pdf")),
           plot = p, width = 9, height = max(4, 0.35 * nrow(plot_df) + 1.5), dpi = 300)
    
    print(p)
  }
}

plot_hallmark_nes_barplot(hallmark_gsea_all)

############################ SURVIVAL ANALYSIS #################################

####################### DATA PREP FOR SURVIVAL ANALYSIS ########################X

# Step 1: Extract patient-level clinical data for survival analysis from df_phenotype
df_survival_stats <- df_phenotype %>%
  select(patient_id = sample,
         days_to_death = days_to_death.demographic,
         vital_status = vital_status.demographic)
# 
# # Step 2: Report any patients missing a survival time
# # do not drop or impute
n_missing_days_to_death <- sum(is.na(df_survival_stats$days_to_death))
if (n_missing_days_to_death > 0) {
  warning(n_missing_days_to_death, " patient(s) missing days_to_death: ",
          paste(df_survival_stats$patient_id[is.na(df_survival_stats$days_to_death)], collapse = ", "))
}

# Step 3: All patients are deceased (df_phenotype is already filtered to vital_status.demographic == "Dead")
df_survival_stats <- df_survival_stats %>%
  mutate(vital_status = 1) # Dead = 1

# Step 4: Keep the full TCGA barcode and add a patient-level ID (first 3 barcode segments) for
# matching against the expression matrix at the patient level;
df_survival_stats <- df_survival_stats %>%
  mutate(patient_id_full = patient_id,
         patient_id = sub("^([^-]+-[^-]+-[^-]+)-.*$", "\\1", patient_id))

# Step 5: Check duplicate samples
stopifnot(!anyDuplicated(df_survival_stats$patient_id))   # one row per patient

df_survival_stats <- df_survival_stats %>%
  select(patient_id, days_to_death, vital_status)

cat("Patients in survival table:", nrow(df_survival_stats), "\n")   # expect 122

# Step 6: VST-normalize the DESeq2 counts for use in the Cox model (dds is left untouched)
vst_matrix <- assay(vst(dds, blind = FALSE))

cat("vst_matrix dimensions:", dim(vst_matrix), "\n")

# Step 7: renamed from full barcode to patient-level ID
vst_patient_id <- sub("^([^-]+-[^-]+-[^-]+)-.*$", "\\1", colnames(vst_matrix))

vst_matrix_patient <- vst_matrix
colnames(vst_matrix_patient) <- sub("^([^-]+-[^-]+-[^-]+)-.*$", "\\1", colnames(vst_matrix))

stopifnot(!anyDuplicated(colnames(vst_matrix_patient)))

cat("vst_matrix_patient dimensions:", dim(vst_matrix_patient), "\n")   # expect 122 columns
cat("vst_matrix_patient column names match df_survival_stats patient_id set:",
    setequal(colnames(vst_matrix_patient), df_survival_stats$patient_id), "\n")

# Step 8: Filter each of the 4 DESeq2 results to significant DEGs (padj < 0.05, |log2FC| >= 1 -
# same fc1 threshold, using the existing get_signif_genes() function
deg_Old_survival   <- get_signif_genes(results_Old_survival)
deg_Young_survival <- get_signif_genes(results_Young_survival)
deg_SS_age         <- get_signif_genes(results_SS_age)
deg_LS_age         <- get_signif_genes(results_LS_age)

cat("Old_survival DEGs:", length(deg_Old_survival), "\n")
cat("Young_survival DEGs:", length(deg_Young_survival), "\n")
cat("SS_age DEGs:", length(deg_SS_age), "\n")
cat("LS_age DEGs:", length(deg_LS_age), "\n")

# Step 9: Union across all 4 comparisons, all unique significantly expressed genes
deg_union_genes <- unique(c(deg_Old_survival, deg_Young_survival, deg_SS_age, deg_LS_age))
cat("DEG union (deduplicated):", length(deg_union_genes), "\n")

# Step 10: Confirm gene ID format matches between DESeq2 results and vst_matrix_patient row
# names before subsetting - Note: a format mismatch (e.g. gene symbols vs Ensembl IDs) would return zero matches
cat("Example DEG union gene names:", head(deg_union_genes, 5), "\n")
cat("Example vst_matrix_patient row names:", head(rownames(vst_matrix_patient), 5), "\n")
cat("All DEG union genes found in vst_matrix_patient row names:",
    all(deg_union_genes %in% rownames(vst_matrix_patient)), "\n")

# Step 11: Subset vst_matrix_patient to the DEG union genes, then transpose so rows =
# patients and columns = genes (needed to merge with df_survival_stats by patient_id)
vst_deg_matrix <- vst_matrix_patient[deg_union_genes, ]

vst_deg_t <- as.data.frame(t(vst_deg_matrix)) %>%
  rownames_to_column("patient_id")

# Step 12: Merge patient-level expression data with survival data by patient_id
df_cox_input <- df_survival_stats %>%
  inner_join(vst_deg_t, by = "patient_id")

cat("df_cox_input rows:", nrow(df_cox_input), "\n")
cat("df_cox_input columns:", ncol(df_cox_input),
    "(expect 3 metadata +", length(deg_union_genes), "genes =", 3 + length(deg_union_genes), ")\n")
cat("any missing values in gene expression columns:",
    any(is.na(df_cox_input[, 4:ncol(df_cox_input)])), "\n")
cat("first 6 column names:", head(colnames(df_cox_input), 6), "\n")

####################### UNIVARIANT COX REGRESSION ANALYSIS #####################

# Step 1: Run a univariate Cox regression for each of the 237 DEG union genes
# (columns 4 onward in df_cox_input), extracting HR, 95% CI, and padj-value for each
gene_cols <- colnames(df_cox_input)[4:ncol(df_cox_input)]

cox_results <- lapply(gene_cols, function(gene) {
  formula_gene <- as.formula(paste0("Surv(days_to_death, vital_status) ~ `", gene, "`"))
  fit <- coxph(formula_gene, data = df_cox_input)
  summary_fit <- summary(fit)

  data.frame(
    gene = gene,
    HR = summary_fit$coefficients[1, "exp(coef)"],
    CI_lower = summary_fit$conf.int[1, "lower .95"],
    CI_upper = summary_fit$conf.int[1, "upper .95"],
    p_value = summary_fit$coefficients[1, "Pr(>|z|)"]
  )
})

univariate_cox_results <- do.call(rbind, cox_results)

# Step 2: Benjamini-Hochberg FDR correction across all 237 p-values
univariate_cox_results$padj <- p.adjust(univariate_cox_results$p_value, method = "BH")

# Step 3: Sort by padj ascending
univariate_cox_results <- univariate_cox_results[order(univariate_cox_results$padj), ]
rownames(univariate_cox_results) <- NULL

cat("Genes passing padj < 0.05:", sum(univariate_cox_results$padj < 0.05), "\n")
print(head(univariate_cox_results, 10))

# Step 4: Export to Excel (openxlsx, matching the export style used throughout the project)
wb <- createWorkbook()
addWorksheet(wb, "Univariate_Cox")
writeData(wb, "Univariate_Cox", univariate_cox_results)
saveWorkbook(wb, "univariate_cox_results.xlsx", overwrite = TRUE)

cat("univariate_cox_results.xlsx written:", file.exists("univariate_cox_results.xlsx"), "\n")

############################### LASSO-COX ANALYSIS #############################

# Step 1: Build the predictor matrix (x) and survival object (y) from df_cox_input.
# gene_cols (237 DEG union genes) was already defined in the univariate Cox section above.
x <- as.matrix(df_cox_input[, gene_cols])
y <- Surv(df_cox_input$days_to_death, df_cox_input$vital_status)

cat("x dimensions:", dim(x), "\n")
cat("y length:", length(y), "\n")

# Step 2: Set the project-wide standard seed before the cross-validation step
set.seed(123)

# Step 3: Run cv.glmnet (Cox model, 10-fold CV) to select lambda
cv_fit <- cv.glmnet(x, y, family = "cox", nfolds = 10)

# Step 4: Report lambda.min and lambda.1se
cat("lambda.min:", cv_fit$lambda.min, "\n")
cat("lambda.1se:", cv_fit$lambda.1se, "\n")

# Step 5: Extract non-zero coefficient genes at lambda.1se and lambda.min
coef_1se <- coef(cv_fit, s = "lambda.1se")
coef_min <- coef(cv_fit, s = "lambda.min")

# Use which() on the plain coefficient vector (not Matrix row-subsetting) so gene
# names stay character(0) rather than NULL if no genes survive at a given lambda 
# NULL would silently drop the "gene" column from data.frame() below
idx_1se <- which(coef_1se[, 1] != 0)
idx_min <- which(coef_min[, 1] != 0)

cat("Non-zero genes at lambda.1se:", length(idx_1se), "\n")
cat("Non-zero genes at lambda.min:", length(idx_min), "\n")

# lambda.1se was intended as the primary (more conservative) signature, but it
# shrinks all 237 gene coefficients to zero given n = 122 patients (kept here for
# the record only, not exported or used as the reported signature)
lasso_cox_signature_1se <- data.frame(
  gene = rownames(coef_1se)[idx_1se],
  coefficient = as.numeric(coef_1se[idx_1se, 1])
)
lasso_cox_signature_1se <- lasso_cox_signature_1se[order(-abs(lasso_cox_signature_1se$coefficient)), ]
rownames(lasso_cox_signature_1se) <- NULL

cat("NOTE: lambda.1se yielded 0 non-zero genes - using lambda.min signature as primary going forward\n")
print(lasso_cox_signature_1se)

# Step 6: Save the lambda.min signature (gene names + coefficients) as the primary
# reported signature, sorted by effect size (largest |coefficient| first)
lasso_cox_signature <- data.frame(
  gene = rownames(coef_min)[idx_min],
  coefficient = as.numeric(coef_min[idx_min, 1])
)
lasso_cox_signature <- lasso_cox_signature[order(-abs(lasso_cox_signature$coefficient)), ]
rownames(lasso_cox_signature) <- NULL

print(lasso_cox_signature)

# Step 7: Export the primary (lambda.min) signature to Excel (openxlsx, same style
# as univariate_cox_results.xlsx)
wb <- createWorkbook()
addWorksheet(wb, "Lasso_Cox_Signature")
writeData(wb, "Lasso_Cox_Signature", lasso_cox_signature)
saveWorkbook(wb, "lasso_cox_signature.xlsx", overwrite = TRUE)

cat("lasso_cox_signature.xlsx written:", file.exists("lasso_cox_signature.xlsx"), "\n")

################################ RISK SCORE MODEL ##############################

# Step 1: Subset df_cox_input to the 6 lasso_cox_signature genes, in signature order
signature_genes <- lasso_cox_signature$gene
cat("Signature genes found in df_cox_input:", all(signature_genes %in% colnames(df_cox_input)), "\n")

sig_expr <- as.matrix(df_cox_input[, signature_genes])

# Step 2: Compute risk_score = sum(coefficient_i * gene_expression_i) per patient
risk_score <- as.numeric(sig_expr %*% lasso_cox_signature$coefficient)

# Step 3: Save patient_id, days_to_death, vital_status, and risk_score into df_risk_score
df_risk_score <- data.frame(
  patient_id = df_cox_input$patient_id,
  days_to_death = df_cox_input$days_to_death,
  vital_status = df_cox_input$vital_status,
  risk_score = risk_score
)

cat("df_risk_score rows:", nrow(df_risk_score), "\n")

# Step 4: Split patients into high/low risk groups using the median risk_score as
# cutoff (at or above median = high risk, below = low risk)
median_risk_score <- median(df_risk_score$risk_score)
df_risk_score$risk_group <- ifelse(df_risk_score$risk_score >= median_risk_score, "high", "low")

cat("Median risk score (cutoff):", median_risk_score, "\n")
cat("Patients per risk group:\n")
print(table(df_risk_score$risk_group))

cat("risk_score summary - mean:", mean(df_risk_score$risk_score),
    "min:", min(df_risk_score$risk_score),
    "max:", max(df_risk_score$risk_score), "\n")

# Step 5: Export df_risk_score to Excel 
wb <- createWorkbook()
addWorksheet(wb, "Risk_Score")
writeData(wb, "Risk_Score", df_risk_score)
saveWorkbook(wb, "risk_score_results.xlsx", overwrite = TRUE)

cat("risk_score_results.xlsx written:", file.exists("risk_score_results.xlsx"), "\n")

# C-index of the original risk score (reverse = TRUE: higher score = higher risk)
c_original <- concordance(
  Surv(days_to_death, vital_status) ~ risk_score,
  data = df_risk_score, reverse = TRUE
)
cat("Original C-index:", c_original$concordance, " SE:", sqrt(c_original$var), "\n")

############################### KAPLAN-MEIER CURVE #############################

# Note: all 122 patients in this cohort are deceased (no censoring), so both KM
# curves will reach 0% survival probability rather than leveling off with censored
# patients still at risk this is expected given the cohort design and differs
# from typical partially-censored TCGA survival plots

# Step 1: Fit survival curves for the high vs low risk groups
km_fit <- survfit(Surv(days_to_death, vital_status) ~ risk_group, data = df_risk_score)

# Step 2: Log-rank test comparing the two risk groups
km_logrank <- survdiff(Surv(days_to_death, vital_status) ~ risk_group, data = df_risk_score)

# Step 3: Extract the log-rank p-value (survdiff gives a chi-square stat on 1 df
# for a 2-group comparison)
km_logrank_pvalue <- pchisq(km_logrank$chisq, df = 1, lower.tail = FALSE)

cat("Log-rank test p-value:", km_logrank_pvalue, "\n")
cat("Significant at p < 0.05:", km_logrank_pvalue < 0.05, "\n")

# Step 4: Plot with survminer::ggsurvplot - high/low curves in distinct colors,
# log-rank p-value on the plot, and a risk table below it. risk_group levels default
# to alphabetical order ("high", "low"), so palette/legend.labs follow that order
km_plot <- ggsurvplot(
  km_fit,
  data = df_risk_score,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("red", "blue"),
  legend.title = "Risk group",
  legend.labs = c("High risk", "Low risk"),
  xlab = "Time (days)",
  ylab = "Survival probability"
)

print(km_plot)

# Step 5: Save to file (pdf + print, so the curve and risk table are both included)
pdf("kaplan_meier_curve.pdf", width = 8, height = 8)
print(km_plot)
dev.off()

cat("kaplan_meier_curve.pdf written:", file.exists("kaplan_meier_curve.pdf"), "\n")

########################### TIME-DEPENDENT ROC CURVES ##########################

# Note: since all patients are deceased (no censoring), the time-dependent ROC
# calculation doesn't need to handle censoring weights the way it typically would
# in a partially-censored cohort - simplifying the estimation compared to the
# standard case in most published examples

# Step 1: Compute time-dependent ROC/AUC at 1 year, 3 year, and 5 year time points,
# using risk_score as the continuous predictor
roc_timepoints <- c(365, 1095, 1825)

roc_fit <- timeROC(
  T = df_risk_score$days_to_death,
  delta = df_risk_score$vital_status,
  marker = df_risk_score$risk_score,
  cause = 1,
  weighting = "marginal",
  times = roc_timepoints,
  ROC = TRUE,
  iid = FALSE
)

# Step 2: Report AUC at each timepoint
cat("AUC at 1-year (365 days):", roc_fit$AUC[1], "\n")
cat("AUC at 3-year (1095 days):", roc_fit$AUC[2], "\n")
cat("AUC at 5-year (1825 days):", roc_fit$AUC[3], "\n")

# Step 3: Build a long data frame of FP/TP (1 - specificity / sensitivity) across
# all three time points, with AUC values folded into the time point labels so they
# show up in the plot legend
roc_labels <- sprintf("%s (AUC = %.3f)",
                       c("1-year", "3-year", "5-year"),
                       roc_fit$AUC)

roc_df <- data.frame(
  FP = as.vector(roc_fit$FP),
  TP = as.vector(roc_fit$TP),
  timepoint = factor(rep(roc_labels, each = nrow(roc_fit$FP)), levels = roc_labels)
)

# Step 4: Plot all three ROC curves on one figure - distinct colors/line types per
# timepoint, AUC in the legend (via roc_labels), and a diagonal reference line
roc_plot <- ggplot(roc_df, aes(x = FP, y = TP, color = timepoint, linetype = timepoint)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted") +
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       color = "Timepoint", linetype = "Timepoint",
       title = "Time-dependent ROC curves (risk score)") +
  theme_minimal()

print(roc_plot)

# Step 5: Save to file
ggsave(filename = "time_dependent_roc.pdf", plot = roc_plot, width = 7, height = 6)

cat("time_dependent_roc.pdf written:", file.exists("time_dependent_roc.pdf"), "\n")

############################### RISK SCORE PLOT #############################

# Step 1: Order patients by risk_score (ascending) - all 3 panels share this order
order_idx <- order(df_risk_score$risk_score)
plot_df <- df_risk_score[order_idx, ]
plot_df$rank <- seq_len(nrow(plot_df))

# Step 2: Panel 1 - risk score distribution, colored by risk_group, with the
# median cutoff shown as a dashed horizontal line
p1 <- ggplot(plot_df, aes(x = rank, y = risk_score, color = risk_group)) +
  geom_point() +
  geom_hline(yintercept = median_risk_score, linetype = "dashed", color = "black") +
  scale_color_manual(values = c(high = "red", low = "blue")) +
  labs(y = "Risk score", color = "Risk group") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Step 3: Panel 2 - survival time (days_to_death) in the same patient order. Colored
# by risk_group rather than vital_status, since vital_status is constant (all 122
# patients are deceased) and wouldn't show any contrast - risk_group is what actually
# reveals whether shorter survival times cluster on the high-risk side
p2 <- ggplot(plot_df, aes(x = rank, y = days_to_death, color = risk_group)) +
  geom_point() +
  scale_color_manual(values = c(high = "red", low = "blue")) +
  labs(y = "Days to death", color = "Risk group") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Step 4: Panel 3 - signature gene expression heatmap, same patient order. Each gene
# is z-scored across patients (per-gene relative expression) so the diverging
# blue-white-red scale is centered meaningfully at 0
heatmap_z <- as.data.frame(scale(sig_expr[order_idx, ]))
heatmap_z$rank <- plot_df$rank

heatmap_long <- pivot_longer(heatmap_z, cols = -rank, names_to = "gene", values_to = "z_expression")
heatmap_long$gene <- factor(heatmap_long$gene, levels = rev(signature_genes))

p3 <- ggplot(heatmap_long, aes(x = rank, y = gene, fill = z_expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Patient rank (sorted by risk score, ascending)", y = NULL, fill = "Z-score\n(VST expr)") +
  theme_classic(base_size = 12)

# Step 5: Stack all 3 panels vertically, aligned on patient rank (panel 1 on top)
risk_plot <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))

print(risk_plot)

# Step 6: Save to file
ggsave(filename = "risk_score_plot.pdf", plot = risk_plot, width = 8, height = 9)

cat("risk_score_plot.pdf written:", file.exists("risk_score_plot.pdf"), "\n")

##################### INDIVIDUAL SIGNATURE GENE KM CURVES ######################

# Step 1: Diagnostic - check each gene's median split BEFORE running the full loop.
# OR10Q1 is expected to fail the median split (majority of samples at the VST floor);
# it is handled separately in Step 2.
for (gene in signature_genes) {
  vals <- df_cox_input[[gene]]
  med <- median(vals)
  grp <- ifelse(vals >= med, "high", "low")
  cat(gene, "- median:", med, "- group counts: ")
  print(table(grp))
}

# Step 2: For each of the 6 signature genes, split patients into high/low
# expression (median split; OR10Q1 split at the VST floor)
km_gene_plots <- list()
km_fits <- list()
km_data <- list()
km_labs <- list()
km_gene_results <- data.frame()


for (gene in signature_genes) {
  
  gene_data <- df_cox_input %>%
    select(patient_id, days_to_death, vital_status, all_of(gene))
  
  median_expr <- median(gene_data[[gene]])
  
  # Grouping: median split, except OR10Q1 (majority of samples at the VST floor)
  if (gene == "OR10Q1") {
    floor_val <- min(gene_data[[gene]])
    cat("OR10Q1 at floor:", sum(gene_data[[gene]] == floor_val), "of", nrow(gene_data), "\n")
    gene_data$expr_group <- ifelse(gene_data[[gene]] > floor_val, "high", "low")
    legend_labs <- c("Above floor", "At floor")
  } else {
    gene_data$expr_group <- ifelse(gene_data[[gene]] >= median_expr, "high", "low")
    legend_labs <- c("High", "Low")
  }
  
  cat(gene, "- group counts:\n")
  print(table(gene_data$expr_group))
  
  if (length(unique(gene_data$expr_group)) < 2) {
    warning(paste0(gene, ": only one group - skipped."))
    next
  }
  
# Step 3: Fit KM curve and log-rank test for this gene
  fit <- survfit(Surv(days_to_death, vital_status) ~ expr_group, data = gene_data)
  logrank <- survdiff(Surv(days_to_death, vital_status) ~ expr_group, data = gene_data)
  pval <- pchisq(logrank$chisq, df = 1, lower.tail = FALSE)
  
  # Record the cutoff actually used, plus group sizes
  cutoff <- if (gene == "OR10Q1") floor_val else median_expr
  
  km_gene_results <- rbind(km_gene_results,
                           data.frame(gene = gene,
                                      cutoff = cutoff,
                                      n_high = sum(gene_data$expr_group == "high"),
                                      n_low  = sum(gene_data$expr_group == "low"),
                                      logrank_p = pval))
  
# Step 4: Store fit, data and legend labels
  km_fits[[gene]] <- fit
  km_data[[gene]] <- gene_data
  km_labs[[gene]] <- legend_labs
} 

# Step 5: Adjust the 6 log-rank p-values (BH) across the 6 genes
km_gene_results$logrank_p_adj <- p.adjust(km_gene_results$logrank_p, method = "BH")

# Step 6: Plot one KM curve per gene, showing raw and adjusted p-values
for (g in km_gene_results$gene) {
  r <- km_gene_results[km_gene_results$gene == g, ]
  km_gene_plots[[g]] <- ggsurvplot(
    km_fits[[g]], data = km_data[[g]],
    pval = paste0("p-value = ", signif(r$logrank_p, 2),
                  "\npadj = ", signif(r$logrank_p_adj, 2)),
    pval.coord = c(1300, 0.9),
    palette = c("red", "blue"),
    legend.title = paste0(g, " expression"),
    legend.labs = km_labs[[g]],
    xlab = "Time (days)", ylab = "Survival probability", title = g)
}

# Step 7: Report p-values for all 6 genes
cat("Log-rank p-values by signature gene (median split; OR10Q1 split at VST floor):\n")
print(km_gene_results)

writexl::write_xlsx(km_gene_results, "signature_gene_km_results.xlsx")

# Step 8: Arrange all 6 curves into one figure (2 rows x 3 cols)
km_gene_arranged <- arrange_ggsurvplots(km_gene_plots, ncol = 3, nrow = 2, print = FALSE)

ggsave(filename = "signature_gene_km_curves.pdf", plot = km_gene_arranged, width = 15, height = 9)

cat("signature_gene_km_curves.pdf written:", file.exists("signature_gene_km_curves.pdf"), "\n")

##########################SIGNATURE GENES BOXPLOTS #############################

# Step 1: Build a patient-level lookup of the 4 combined groups (age_groups x
# survival_groups) from df_phenotype, trimmed to patient-level ID to match df_cox_input
df_groups_lookup <- df_phenotype %>%
  select(patient_id = sample, age_groups, survival_groups) %>%
  mutate(patient_id = sub("^([^-]+-[^-]+-[^-]+)-.*$", "\\1", patient_id)) %>%
  distinct(patient_id, age_groups, survival_groups) %>%
  mutate(combined_group = paste(age_groups, survival_groups, sep = "_"))

cat("df_groups_lookup rows:", nrow(df_groups_lookup), "\n")
cat("combined_group counts:\n")
print(table(df_groups_lookup$combined_group))

# Step 2: Merge the 6 signature genes' expression (from df_cox_input) with
# combined_group, into a new object - df_cox_input itself is left untouched
df_signature_boxplot_data <- df_cox_input %>%
  select(patient_id, all_of(signature_genes)) %>%
  inner_join(df_groups_lookup, by = "patient_id")

cat("df_signature_boxplot_data rows:", nrow(df_signature_boxplot_data), "\n")
cat("any missing combined_group:", any(is.na(df_signature_boxplot_data$combined_group)), "\n")

# Step 3: Reshape to long format - one row per patient per gene, so all 6
# genes can be faceted in one figure and calculate BH-padj
df_signature_long <- df_signature_boxplot_data %>%
  pivot_longer(cols = all_of(signature_genes), names_to = "gene", values_to = "expression")

# Per-gene test + BH adjustment across the 6 signature genes
adj_p_labels <- function(df_long, group_col, test = c("wilcox", "kruskal")) {
  test <- match.arg(test)
  df_long %>%
    group_by(gene) %>%
    summarise(p = if (test == "wilcox")
      wilcox.test(expression ~ .data[[group_col]], exact = FALSE)$p.value
      else
        kruskal.test(expression ~ factor(.data[[group_col]]))$p.value,
      .groups = "drop") %>%
    mutate(p_adj = p.adjust(p, method = "BH"),
           label = paste0("padj = ", formatC(p_adj, format = "g", digits = 2)))
}

# Step 4: Boxplots by combined_group (Old_LS, Old_SS, Young_LS, Young_SS) - one
# facet per gene, with an overall Kruskal-Wallis test p-value shown per facet
# (4-group comparison, so a single omnibus test rather than one pairwise p-value)
lab_combined <- adj_p_labels(df_signature_long, "combined_group", "kruskal")

signature_boxplots_combined <- ggplot(df_signature_long, aes(x = combined_group, y = expression, fill = combined_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  geom_text(data = lab_combined, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "VST expression", fill = "Group") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "signature_boxplots_combined.pdf", plot = signature_boxplots_combined, width = 11, height = 7)

cat("signature_boxplots_combined.pdf written:", file.exists("signature_boxplots_combined.pdf"), "\n")

# Step 5: Boxplots by risk_group (high vs low) - one facet per gene, 
# NOTE: risk_group is derived directly from these same 6 genes' weighted expression
# (risk_score = sum of coefficient_i * gene_i expression), so a strong separation
# here is expected by construction,-
# this is a descriptive/illustrative companion to the risk score plot's heatmap panel

df_signature_riskgroup_data <- df_cox_input %>%
  select(patient_id, all_of(signature_genes)) %>%
  inner_join(df_risk_score %>% select(patient_id, risk_group), by = "patient_id")

cat("df_signature_riskgroup_data rows:", nrow(df_signature_riskgroup_data), "\n")

df_signature_riskgroup_long <- df_signature_riskgroup_data %>%
  pivot_longer(cols = all_of(signature_genes), names_to = "gene", values_to = "expression")

signature_boxplots_riskgroup <- ggplot(df_signature_riskgroup_long, aes(x = risk_group, y = expression, fill = risk_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(high = "red", low = "blue")) +
  labs(x = NULL, y = "VST expression", fill = "Risk group",
       caption = "Descriptive only: risk groups are defined by these genes.") +
  theme_classic(base_size = 12)

#print(signature_boxplots_riskgroup)

ggsave(filename = "signature_boxplots_riskgroup.pdf", plot = signature_boxplots_riskgroup, width = 10, height = 7)

cat("signature_boxplots_riskgroup.pdf written:", file.exists("signature_boxplots_riskgroup.pdf"), "\n")

# Step 6: Boxplots by survival_groups (SS vs LS) - one facet per gene, with a
# Wilcoxon rank-sum test p-value shown per facet.
# NOTE: these 6 genes were selected via survival-based methods (DESeq2 survival
# contrasts and LASSO-Cox on survival time), so a strong SS/LS separation is
# expected by construction. This plot is descriptive/illustrative, consistent with the high/low risk circularity

df_signature_survival_data <- df_cox_input %>%
  select(patient_id, all_of(signature_genes)) %>%
  inner_join(df_groups_lookup %>% select(patient_id, survival_groups), by = "patient_id")

cat("df_signature_survival_data rows:", nrow(df_signature_survival_data), "\n")

df_signature_survival_long <- df_signature_survival_data %>%
  pivot_longer(cols = all_of(signature_genes), names_to = "gene", values_to = "expression")

lab_survival <- adj_p_labels(df_signature_survival_long, "survival_groups", "wilcox")

signature_boxplots_survival <- ggplot(df_signature_survival_long, aes(x = survival_groups, y = expression, fill = survival_groups)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  geom_text(data = lab_survival, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "VST expression", fill = "Survival group",
       caption = "Descriptive only: genes were selected using survival.") +
  theme_classic(base_size = 12)

print(signature_boxplots_survival)

writexl::write_xlsx(list(combined_groups = lab_combined, survival_groups = lab_survival),
                    "signature_boxplot_pvalues.xlsx")

ggsave(filename = "signature_boxplots_survival.pdf", plot = signature_boxplots_survival, width = 10, height = 7)

cat("signature_boxplots_survival.pdf written:", file.exists("signature_boxplots_survival.pdf"), "\n")

############################### INDEPENDENT-PREDICTOR CHECK #############################

# Goal: check whether risk_score is still a significant predictor of survival once
# age is adjusted for, i.e. that its prognostic value isn't just standing in for age

# Step 1: Build an age lookup table from df_phenotype, using the continuous age at
# diagnosis (not the Old/Young split, which was built for the DESeq2 design)
#  truncate the full TCGA barcode to the same patient-level ID used elsewhere
df_age_lookup <- df_phenotype %>%
  select(patient_id = sample, age = age_at_index.demographic) %>%
  mutate(patient_id = sub("^([^-]+-[^-]+-[^-]+)-.*$", "\\1", patient_id))

# Step 2: Confirm duplicate vials agree on age before collapsing to one row per patient
age_dup_check <- df_age_lookup %>%
  group_by(patient_id) %>%
  summarise(n_distinct_age = n_distinct(age), .groups = "drop") %>%
  filter(n_distinct_age > 1)

if (nrow(age_dup_check) > 0) {
  stop("Duplicate patient_id rows disagree on age - check before merging")
}

df_age_lookup <- df_age_lookup %>% distinct(patient_id, age)

cat("df_age_lookup rows:", nrow(df_age_lookup), "\n")

# Step 3: Merge age into df_risk_score by patient_id, into a new object 
# Note: df_risk_score itself is left untouched
df_multivariate_input <- df_risk_score %>%
  inner_join(df_age_lookup, by = "patient_id")

cat("df_multivariate_input rows:", nrow(df_multivariate_input), "\n")
cat("any missing age values:", any(is.na(df_multivariate_input$age)), "\n")

# Step 4: Fit the multivariate Cox model - risk_score adjusted for age
# Note: SS/LS survival group is NOT included as a covariate it's derived from the
# outcome itself, so including it would be circular 
multivariate_cox_model <- coxph(Surv(days_to_death, vital_status) ~ risk_score + age,
                                 data = df_multivariate_input)

print(summary(multivariate_cox_model))

# Step 5: Extract HR, 95% CI, and p-value for each variable
cox_summary <- summary(multivariate_cox_model)

multivariate_cox_results <- data.frame(
  variable = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[, "exp(coef)"],
  CI_lower = cox_summary$conf.int[, "lower .95"],
  CI_upper = cox_summary$conf.int[, "upper .95"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]
)

cat("Multivariate Cox model results:\n")
print(multivariate_cox_results)

cat("risk_score remains significant after adjusting for age (p < 0.05):",
    multivariate_cox_results$p_value[multivariate_cox_results$variable == "risk_score"] < 0.05, "\n")

# Step 6: Forest plot showing both variables' HR and CI
forest_plot <- ggforest(multivariate_cox_model, data = df_multivariate_input)

print(forest_plot)

ggsave(filename = "forest_plot.pdf", plot = forest_plot, width = 7, height = 5)

cat("forest_plot.pdf written:", file.exists("forest_plot.pdf"), "\n")

# Step 7: Export HR/CI/p-value table to Excel
wb <- createWorkbook()
addWorksheet(wb, "Multivariate_Cox")
writeData(wb, "Multivariate_Cox", multivariate_cox_results)
saveWorkbook(wb, "multivariate_cox_results.xlsx", overwrite = TRUE)

cat("multivariate_cox_results.xlsx written:", file.exists("multivariate_cox_results.xlsx"), "\n")

################### DIAGNOSTICS - PROPORTIONAL HAZARDS CHECK ###################

# Goal: check whether multivariate_cox_model satisfies the proportional hazards (PH)
# assumption that each variable's effect on the hazard is constant over time.
# A significant p-value (p < 0.05) here flags a violation of the assumption.

# Step 1: Test PH assumption per variable and globally, based on scaled Schoenfeld residuals
ph_test <- cox.zph(multivariate_cox_model)

print(ph_test)

# Step 2: Plot scaled Schoenfeld residuals over time for each variable
ph_plot <- ggcoxzph(ph_test)

print(ph_plot)

# Step 3: Save to file (pdf + print, same pattern as the KM curve, since ggcoxzph
# returns a multi-panel object rather than a single ggplot)
pdf("ph_diagnostics_plot.pdf", width = 8, height = 8)
print(ph_plot)
dev.off()

cat("ph_diagnostics_plot.pdf written:", file.exists("ph_diagnostics_plot.pdf"), "\n")

# Step 4: Plain-language readout of whether the PH assumption holds (p > 0.05 = holds)
ph_pvalues <- ph_test$table[, "p"]

cat("PH assumption holds for risk_score (p > 0.05):", ph_pvalues["risk_score"] > 0.05, "\n")
cat("PH assumption holds for age (p > 0.05):", ph_pvalues["age"] > 0.05, "\n")
cat("PH assumption holds globally (p > 0.05):", ph_pvalues["GLOBAL"] > 0.05, "\n")

######################### INTERNAL BOOTSTRAP VALIDATION #######################

# Step 1: Set the standard seed before the bootstrap loop (same seed as before)
set.seed(123)

# Step 2: Set up storage for the 1,000 bootstrap iterations - NA marks iterations
# where the resampled model selects zero genes (C-index is undefined there)
n_boot <- 1000
n_patients <- nrow(df_cox_input)

# original data objects and original C-index
x_orig <- as.matrix(df_cox_input[, gene_cols])
y_orig <- Surv(df_cox_input$days_to_death, df_cox_input$vital_status)
c_orig <- c_original$concordance #from the risk score block

bootstrap_cindex_results <- rep(NA_real_, n_boot)   # apparent (on resample)
test_cindex_results      <- rep(NA_real_, n_boot)   # refitted model on original data
optimism                 <- rep(NA_real_, n_boot)
lambda_min_results       <- rep(NA_real_, n_boot)
n_genes_results          <- rep(NA_integer_, n_boot)
zero_gene_iterations <- 0

# Step 3: Bootstrap loop - resample the 122 patients with replacement, refit the
# LASSO-Cox model on the resampled data (same approach as the original: cv.glmnet,
# family = "cox", 10-fold CV, lambda.min), and compute the C-index
for (i in seq_len(n_boot)) {
  boot_idx  <- sample(seq_len(n_patients), size = n_patients, replace = TRUE)
  boot_data <- df_cox_input[boot_idx, ]
  x_boot <- as.matrix(boot_data[, gene_cols])
  y_boot <- Surv(boot_data$days_to_death, boot_data$vital_status)
  
  cv_boot  <- cv.glmnet(x_boot, y_boot, family = "cox", nfolds = 10)
  coef_boot <- coef(cv_boot, s = "lambda.min")
  idx_boot  <- which(coef_boot[, 1] != 0)
  
  lambda_min_results[i] <- cv_boot$lambda.min
  n_genes_results[i]    <- length(idx_boot)
  
  if (length(idx_boot) == 0) {
    zero_gene_iterations <- zero_gene_iterations + 1
    next
  }
  
  lp_boot <- as.numeric(predict(cv_boot, newx = x_boot, s = "lambda.min", type = "link"))
  lp_orig <- as.numeric(predict(cv_boot, newx = x_orig, s = "lambda.min", type = "link"))
  
  c_boot <- concordance(y_boot ~ lp_boot, reverse = TRUE)$concordance
  c_test <- concordance(y_orig ~ lp_orig, reverse = TRUE)$concordance
  
  bootstrap_cindex_results[i] <- c_boot
  test_cindex_results[i]      <- c_test
  optimism[i]                 <- c_boot - c_test
}

mean_optimism <- mean(optimism, na.rm = TRUE)
c_corrected   <- c_orig - mean_optimism

cat("Apparent (original) C-index:", c_orig, "\n")
cat("Mean optimism:", mean_optimism, "\n")
cat("Optimism-corrected C-index:", c_corrected, "\n")
cat("Median lambda.min:", median(lambda_min_results), " (original: 0.1887)\n")
cat("Median number of genes selected:", median(n_genes_results), "\n")

# Step 4: Summarize the valid (non-NA) C-index values
valid_cindex <- bootstrap_cindex_results[!is.na(bootstrap_cindex_results)]
n_valid <- length(valid_cindex)
mean_cindex <- mean(valid_cindex)
ci_lower <- quantile(valid_cindex, 0.025)
ci_upper <- quantile(valid_cindex, 0.975)

cat("Valid (non-empty model) iterations:", n_valid, "\n")
cat("Mean C-index:", mean_cindex, "\n")
cat("95% CI:", ci_lower, "-", ci_upper, "\n")

# Step 5: Export summary statistics to Excel
bootstrap_summary <- data.frame(
  metric = c("n_valid_iterations", "mean_cindex", "ci_lower_2.5", "ci_upper_97.5",
             "zero_gene_iterations", "original_cindex", "mean_optimism",
             "optimism_corrected_cindex", "median_lambda_min", "median_n_genes"),
  value = c(n_valid, mean_cindex, ci_lower, ci_upper, zero_gene_iterations,
            c_orig, mean_optimism, c_corrected,
            median(lambda_min_results), median(n_genes_results))
)

wb <- createWorkbook()
addWorksheet(wb, "Bootstrap_Summary")
writeData(wb, "Bootstrap_Summary", bootstrap_summary)
saveWorkbook(wb, "bootstrap_validation_results.xlsx", overwrite = TRUE)

cat("bootstrap_validation_results.xlsx written:", file.exists("bootstrap_validation_results.xlsx"), "\n")
