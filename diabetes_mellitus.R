library(tidyverse)
library(DESeq2)
library(readxl)
library(ggrepel)
library(VennDiagram)


genesets_dir <- file.path('.', 'genesets')

pval_cutoff <- 0.05
lfc_cutoff <- 0.585


# --------------------
# Explore counts data 
# --------------------

meta_data <- read_csv('PRJNA748404_meta.txt') %>% 
  mutate(treatment = recode_factor(treatment, 'counterpart' = 'control', 'post treatment 16 weeks' = 'treated'))

data <- read_tsv('counts_mRNA_2_deduped_UCSC.txt', skip = 1)

sample_names <- sub('.+(SRR\\d+)_(?:sorted|deduped)\\.bam', '\\1', names(data)) %>% gsub(pattern = '\\.', replacement = '_')

names(data) <- sample_names

genes_data <- data[, c(1:6)]
counts_data <- data %>% select(starts_with('SRR'))

length(counts_data)  # 6 samples
nrow(counts_data)  # 43432 genes


# ---------------------------
# Create DESeqDataSet object
# ---------------------------

get_DESeqDataSet_obj <- function(cnts, design_formula) {
  
  col_data <- meta_data %>%
    filter(Run %in% names(cnts)) %>%
    select(Run, treatment)
  
  col_data <- col_data[match(names(cnts), col_data$Run),] %>%
    column_to_rownames(var = 'Run')
  
  print(all(rownames(col_data) %in% colnames(cnts)))
  print(all(rownames(col_data) == colnames(cnts)))
  
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(cnts)),
                                colData = col_data,
                                design = design_formula)
  
  print(summary(dds))
  mcols(dds)$basepairs <- genes_data$Length
  
  # Filter rows with low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  print(summary(dds))
  
  dds
}


dds <- get_DESeqDataSet_obj(data.frame(counts_data, row.names = data$Geneid), ~ treatment)
colData(dds)

counts_fpkm <- fpkm(dds) %>%
  data.frame() %>% 
  rownames_to_column(var = 'ensembl_gene_id') %>% 
  as_tibble()
counts_fpkm$meanFPKM <- rowMeans(counts_fpkm[, c(2:7)])
counts_fpkm$highFPKM <- if_else(counts_fpkm$SRR15202006 >= 0.5 & counts_fpkm$SRR15202007 >= 0.5 & counts_fpkm$SRR15202008 >= 0.5 & counts_fpkm$SRR15202009 >= 0.5 & counts_fpkm$SRR15202010 >= 0.5 & counts_fpkm$SRR15202011 >= 0.5, 1, 0)


# --------
# Results
# --------

dds_full <- DESeq(dds)
colData(dds_full)

res <- results(dds_full, contrast = c('treatment', 'treated', 'control'), alpha = pval_cutoff)
res

mcols(res)
summary(res)  # 220 up, 166 down, 5 outliers, 2328 low counts

res[which(is.na(res$pvalue)), ]  # 5 rows where pvalue is NA and padj is NA
res[which(is.na(res$padj) & !is.na(res$pvalue)), ]  # 2328 rows where only padj is NA

res_df <- res %>%  # dataframe version
  data.frame() %>% 
  rownames_to_column(var = 'ensembl_gene_id') %>% 
  as_tibble()

res_sig_df <- res_df %>%  # filter only 386 DE genes (278 if apply log2FoldChange threshold, 267 if add FPKM/236 if other FPKM)
  filter(
    padj < pval_cutoff,
    (log2FoldChange <= -lfc_cutoff | log2FoldChange >= lfc_cutoff),
    ensembl_gene_id %in% (counts_fpkm %>% filter(meanFPKM >= 0.5))$ensembl_gene_id  # highFPKM == 1
  ) %>% 
  arrange(padj)


# -----------
# Comparison
# -----------

mm_ensembl_symbol <- readRDS(file.path(genesets_dir, 'mm_ensembl_symbol.GRCm38.p6.RDS'))

# Paper identified 186 DE genes
DE_genes_df <- read_excel('DE_RNA.xlsx', sheet = 'mRNA')

# Gene symbols we don't have conversion to Ensembl ID for
DE_genes_df %>% 
  anti_join(mm_ensembl_symbol, by = c('gene_symbol' = 'mgi_symbol')) %>% 
  View()

# Gene symbols whose corresponding Ensembl ID is not in results
DE_genes_df %>% 
  inner_join(mm_ensembl_symbol, by = c('gene_symbol' = 'mgi_symbol')) %>% 
  filter(!ensembl_gene_id %in% res_df$ensembl_gene_id) %>% 
  View()

# Merge gene symbols w/ Ensembl ID
res_sig_merge_df <- res_sig_df %>% 
  left_join(mm_ensembl_symbol, by = 'ensembl_gene_id')

# Upregulated genes
dev.off()
grid::grid.draw(
  venn.diagram(
    list(
      paper = (DE_genes_df %>% filter(log2FoldChange > 0))$gene_symbol,
      project = (res_sig_merge_df %>% filter(log2FoldChange > 0))$mgi_symbol
    ),
    fill = c('#F2F0F7', '#DADAEB'),
    filename = NULL
  )
)

# Downregulated genes
dev.off()
grid::grid.draw(
  venn.diagram(
    list(
      paper = (DE_genes_df %>% filter(log2FoldChange < 0))$gene_symbol,
      project = (res_sig_merge_df %>% filter(log2FoldChange < 0))$mgi_symbol
    ),
    fill = c('#F2F0F7', '#DADAEB'),
    filename = NULL
  )
)


# --------------
# Volcano plots
# --------------

# p-adjusted value
res_df %>% 
  mutate(
    sig_threshold = if_else(
      padj < pval_cutoff & abs(log2FoldChange) >= lfc_cutoff,
      if_else(log2FoldChange > 0, 'DE-up', 'DE-down'),
      'non-DE'
    ),
    gene_symbol = case_when(
      ensembl_gene_id == 'ENSMUSG00000072941' ~ 'Sod3',
      ensembl_gene_id == 'ENSMUSG00000048583' ~ 'Igf2',
      ensembl_gene_id == 'ENSMUSG00000019997' ~ 'Ctgf',  # Ccn2
      ensembl_gene_id == 'ENSMUSG00000056758' ~ 'Hmga2',
      ensembl_gene_id == 'ENSMUSG00000029919' ~ 'Hpgds',
      ensembl_gene_id == 'ENSMUSG00000063011' ~ 'Msln'
    )
  ) %>% 
  filter(!is.na(sig_threshold)) %>% 
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = sig_threshold)) +
  geom_text_repel(aes(label = gene_symbol, x = log2FoldChange, y = -log10(padj)), box.padding = 1.2, segment.size = 0.3, segment.color = 'black') +
  scale_color_manual(values = c('green', 'red', 'gray')) +
  xlab('log2 fold change') + 
  ylab('-log10 adjusted p-value') +
  theme_minimal()


# p value
res_df %>% 
  mutate(
    sig_threshold = if_else(
      pvalue < pval_cutoff & abs(log2FoldChange) >= lfc_cutoff,
      if_else(log2FoldChange > 0, 'DE-up', 'DE-down'),
      'non-DE'
    ),
    gene_symbol = case_when(
      ensembl_gene_id == 'ENSMUSG00000072941' ~ 'Sod3',
      ensembl_gene_id == 'ENSMUSG00000048583' ~ 'Igf2',
      ensembl_gene_id == 'ENSMUSG00000019997' ~ 'Ctgf',  # Ccn2
      ensembl_gene_id == 'ENSMUSG00000056758' ~ 'Hmga2',
      ensembl_gene_id == 'ENSMUSG00000029919' ~ 'Hpgds',
      ensembl_gene_id == 'ENSMUSG00000063011' ~ 'Msln'
    )
  ) %>% 
  filter(!is.na(sig_threshold)) %>% 
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = sig_threshold)) +
  geom_text_repel(aes(label = gene_symbol, x = log2FoldChange, y = -log10(pvalue)), box.padding = 1.2, segment.size = 0.3, segment.color = 'black') +
  scale_color_manual(values = c('green', 'red', 'gray')) +
  xlab('log2 fold change') + 
  ylab('-log10 p-value') +
  theme_minimal()
