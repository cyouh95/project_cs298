library(tidyverse)


data_dir <- file.path('.', 'data')
genesets_dir <- file.path('.', 'genesets')


# -----------------
# Load counts data
# -----------------

# Paper data
meta_data <- read_csv('PRJNA694971_meta.txt')
meta_data$regime <- factor(meta_data$regime, levels = c('on land', 'in space with gravity', 'in space without gravity'))
meta_data$experimental_class_type <- factor(meta_data$experimental_class_type, levels = c('F', 'E', 'D', 'C', 'B', 'A'))

counts_files <- list.files(data_dir)

data_paper <- read_tsv(file.path(data_dir, counts_files[1])) %>% select(-EU181)

for (counts_file in counts_files) {
  m <- str_match(counts_file, '(GSM\\d+)_(EU\\d+)')
  run_id <- meta_data[meta_data$`Sample Name` == m[, 2], ]$Run
  
  tmp_df <- read_tsv(file.path(data_dir, counts_file)) %>%
    dplyr::rename(!!run_id := m[, 3]) %>%
    select(Geneid, !!run_id)
  
  data_paper <- data_paper %>% left_join(tmp_df, by = 'Geneid')
}


# Project data
data_project <- read_tsv('counts.txt', skip = 1)

sample_names <- sub('.+(SRR\\d+)_(?:deduped|sorted)\\.bam', '\\1', names(data_project)) %>% gsub(pattern = '\\.', replacement = '_')

names(data_project) <- sample_names


# miRNA data
data_miRNA <- read_csv('counts_data_miRNA.csv')


# Genesets data
conversion_df <- readRDS(file.path(genesets_dir, 'mm_gene_conversions.RDS')) %>%
  dplyr::select(-ensembl_transcript_id) %>% unique()

mm_h <- readRDS(file.path(genesets_dir, 'Mm.h.all.v7.1.entrez.rds'))
mm_c5_bp <- readRDS(file.path(genesets_dir, 'Mm.c5.bp.v7.1.entrez.rds'))
mm_c5_cc <- readRDS(file.path(genesets_dir, 'Mm.c5.cc.v7.1.entrez.rds'))
mm_c5_mf <- readRDS(file.path(genesets_dir, 'Mm.c5.mf.v7.1.entrez.rds'))


# -----------------
# Generate reports
# -----------------

params_default <- list(
  meta_data = meta_data,
  conversion_df = conversion_df,
  mm_h = mm_h,
  mm_c5_bp = mm_c5_bp,
  mm_c5_cc = mm_c5_cc,
  mm_c5_mf = mm_c5_mf
)


# Paper data
rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_paper_A_vs_C.html',
  params = c(params_default, counts_data = list(data_paper), condition = 'A', control = 'C')
)

rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_paper_A_vs_E.html',
  params = c(params_default, counts_data = list(data_paper), condition = 'A', control = 'E')
)


# Project data
rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_project_A_vs_C.html',
  params = c(params_default, counts_data = list(data_project), condition = 'A', control = 'C')
)

rmarkdown::render(
  input = 'analysis_report.Rmd',
  output_file = 'analysis_report_project_A_vs_E.html',
  params = c(params_default, counts_data = list(data_project), condition = 'A', control = 'E')
)


# Comparison
rmarkdown::render(
  input = 'comparison_report.Rmd',
  output_file = 'comparison_report_A_vs_C_miRNA.html',
  params = c(params_default, counts_paper = list(data_paper), counts_project = list(data_project), counts_miRNA = list(data_miRNA), condition = 'A', control = 'C')
)

rmarkdown::render(
  input = 'comparison_report.Rmd',
  output_file = 'comparison_report_A_vs_E_miRNA.html',
  params = c(params_default, counts_paper = list(data_paper), counts_project = list(data_project), counts_miRNA = list(data_miRNA), condition = 'A', control = 'E')
)
