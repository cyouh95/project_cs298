library(tidyverse)


data_dir <- file.path('.', 'data')
genesets_dir <- file.path('.', 'genesets')


# ----------------------------
# Load counts data (paper #1)
# ----------------------------

meta_data_1 <- read_csv('PRJNA694971_meta.txt') %>%
  arrange(experimental_class_type, Run)


# Paper data
counts_files <- list.files(data_dir)

data_paper_1 <- read_tsv(file.path(data_dir, counts_files[1])) %>% select(-EU181)

for (counts_file in counts_files) {
  m <- str_match(counts_file, '(GSM\\d+)_(EU\\d+)')
  run_id <- meta_data_1[meta_data_1$`Sample Name` == m[, 2], ]$Run
  
  tmp_df <- read_tsv(file.path(data_dir, counts_file)) %>%
    dplyr::rename(!!run_id := m[, 3]) %>%
    select(Geneid, !!run_id)
  
  data_paper_1 <- data_paper_1 %>% left_join(tmp_df, by = 'Geneid')
}


# Project data
data_project_1 <- read_tsv('counts_mRNA_1_deduped.txt', skip = 1)
sample_names <- sub('.+(SRR\\d+)_(?:deduped|sorted)\\.bam', '\\1', names(data_project_1)) %>% gsub(pattern = '\\.', replacement = '_')
names(data_project_1) <- sample_names


# ----------------------------
# Load counts data (paper #2)
# ----------------------------

meta_data_2 <- read_csv('PRJNA748404_meta.txt') %>% 
  mutate(treatment = recode(treatment, 'counterpart' = 'control', 'post treatment 16 weeks' = 'diabetes'))


# Project data
data_project_2 <- read_tsv('counts_mRNA_2_deduped_UCSC.txt', skip = 1)
sample_names <- sub('.+(SRR\\d+)_(?:sorted|deduped)\\.bam', '\\1', names(data_project_2)) %>% gsub(pattern = '\\.', replacement = '_')
names(data_project_2) <- sample_names


# ----------------
# Load other data
# ----------------

meta_data_miRNA <- read_csv('PRJNA748413_meta.txt') %>% 
  mutate(treatment = recode(treatment, 'counterpart' = 'control', 'post treatment 16 weeks' = 'diabetes'))


# miRNA data
data_miRNA_1 <- read_csv('counts_miRNA_1.csv') %>%
  dplyr::rename('miRNA' = 'X1')

data_miRNA_2 <- read_tsv('counts_miRNA_2.csv') %>% 
  dplyr::rename('miRNA' = '#miRNA') %>% 
  # mutate(miRNA = str_replace(`#miRNA`, 'mmu-', '')) %>%
  select(miRNA, `000`, `001`, `002`, `003`, `004`, `005`) %>% 
  group_by(miRNA) %>% 
  summarise_all('mean')
names(data_miRNA_2) <- c('miRNA', meta_data_miRNA$Run)


# Genesets data
mm_ensembl_entrez <- readRDS(file.path(genesets_dir, 'mm_ensembl_entrez.GRCm38.p6.RDS'))
mm_ensembl_symbol <- readRDS(file.path(genesets_dir, 'mm_ensembl_symbol.GRCm38.p6.RDS'))

mm_h <- readRDS(file.path(genesets_dir, 'Mm.h.all.v7.1.entrez.rds'))
mm_c5_bp <- readRDS(file.path(genesets_dir, 'Mm.c5.bp.v7.1.entrez.rds'))
mm_c5_cc <- readRDS(file.path(genesets_dir, 'Mm.c5.cc.v7.1.entrez.rds'))
mm_c5_mf <- readRDS(file.path(genesets_dir, 'Mm.c5.mf.v7.1.entrez.rds'))


# ----------------
# Prep input data
# ----------------

# Meta data
col_data_1 <- meta_data_1 %>% select(-treatment, Run, experimental_class_type, regime) %>% dplyr::rename('treatment' = 'experimental_class_type', 'label' = 'regime')
col_data_1$treatment <- factor(col_data_1$treatment, levels = c('E', 'C', 'A'))
col_data_1$label <- factor(col_data_1$label, levels = c('on land', 'in space with gravity', 'in space without gravity'))

col_data_2 <- meta_data_2 %>% select(Run, treatment) %>% mutate(label = treatment)
col_data_2$treatment <- factor(col_data_2$treatment, levels = c('control', 'diabetes'))


# Counts data
counts_paper_1_AC <- data.frame(data_paper_1 %>% select((meta_data_1 %>% filter(experimental_class_type %in% c('A', 'C'), `Assay Type` == 'RNA-Seq'))$Run), row.names = data_paper_1$Geneid)
counts_paper_1_AE <- data.frame(data_paper_1 %>% select((meta_data_1 %>% filter(experimental_class_type %in% c('A', 'E'), `Assay Type` == 'RNA-Seq'))$Run), row.names = data_paper_1$Geneid)

counts_project_1_AC <- data.frame(data_project_1 %>% select((meta_data_1 %>% filter(experimental_class_type %in% c('A', 'C'), `Assay Type` == 'RNA-Seq'))$Run), row.names = data_project_1$Geneid)
counts_project_1_AE <- data.frame(data_project_1 %>% select((meta_data_1 %>% filter(experimental_class_type %in% c('A', 'E'), `Assay Type` == 'RNA-Seq'))$Run), row.names = data_project_1$Geneid)

counts_miRNA_1_AC <- data.frame(data_miRNA_1 %>% select((meta_data_1 %>% filter(experimental_class_type %in% c('A', 'C'), `Assay Type` == 'miRNA-Seq'))$Run), row.names = data_miRNA_1$miRNA)
counts_miRNA_1_AE <- data.frame(data_miRNA_1 %>% select((meta_data_1 %>% filter(experimental_class_type %in% c('A', 'E'), `Assay Type` == 'miRNA-Seq'))$Run), row.names = data_miRNA_1$miRNA)

counts_project_2 <- data.frame(data_project_2 %>% select(contains('SRR')), row.names = data_project_2$Geneid)
counts_project_2 <- data.frame(data_project_2 %>% select(contains('SRR')), row.names = data_project_2$Geneid)

counts_miRNA_2 <- data.frame(data_miRNA_2 %>% select(contains('SRR')), row.names = data_miRNA_2$miRNA)
counts_miRNA_2 <- data.frame(data_miRNA_2 %>% select(contains('SRR')), row.names = data_miRNA_2$miRNA)


# -----------------
# Generate reports
# -----------------

params_default <- list(
  mm_ensembl_entrez = mm_ensembl_entrez,
  mm_ensembl_symbol = mm_ensembl_symbol,
  mm_h = mm_h,
  mm_c5_bp = mm_c5_bp,
  mm_c5_cc = mm_c5_cc,
  mm_c5_mf = mm_c5_mf
)


# Paper 1 data
rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_paper_A_vs_C.html',
  params = c(
    params_default,
    col_data = list(col_data_1),
    counts_data = list(counts_paper_1_AC),
    condition = 'A',
    control = 'C',
    lfc_cutoff = 1
  )
)

rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_paper_A_vs_E.html',
  params = c(
    params_default,
    col_data = list(col_data_1),
    counts_data = list(counts_paper_1_AE),
    condition = 'A',
    control = 'E',
    lfc_cutoff = 1
  )
)


# Project 1 data
rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_project_A_vs_C.html',
  params = c(
    params_default,
    col_data = list(col_data_1),
    counts_data = list(counts_project_1_AC),
    condition = 'A',
    control = 'C',
    lfc_cutoff = 1
  )
)

rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_project_A_vs_E.html',
  params = c(
    params_default,
    col_data = list(col_data_1),
    counts_data = list(counts_project_1_AE),
    condition = 'A',
    control = 'E',
    lfc_cutoff = 1
  )
)


# Comparison
rmarkdown::render(
  input = 'comparison_report.Rmd',
  output_file = 'comparison_report_A_vs_C_miRNA.html',
  params = c(
    params_default,
    col_data = list(col_data_1),
    counts_paper = list(counts_paper_1_AC),
    counts_project = list(counts_project_1_AC),
    counts_miRNA = list(counts_miRNA_1_AC),
    condition = 'A',
    control = 'C',
    lfc_cutoff = 1
  )
)

rmarkdown::render(
  input = 'comparison_report.Rmd',
  output_file = 'comparison_report_A_vs_E_miRNA.html',
  params = c(
    params_default,
    col_data = list(col_data_1),
    counts_paper = list(counts_paper_1_AE),
    counts_project = list(counts_project_1_AE),
    counts_miRNA = list(counts_miRNA_1_AE),
    condition = 'A',
    control = 'E',
    lfc_cutoff = 1
  )
)


# Project 2 data
rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_project_diabetes_vs_control_UCSC.html',
  params = c(
    params_default,
    col_data = list(col_data_2),
    counts_data = list(counts_project_2),
    condition = 'diabetes',
    control = 'control',
    lfc_cutoff = 0.585,
    genes_length = list(data_project_2$Length),
    use_fpkm = T
  )
)
