library(tidyverse)
library(biomaRt)


data_dir <- file.path('.', 'genesets')


# ----------------------
# Current version (104)
# ----------------------

# Choosing database
listEnsembl()
listEnsemblArchives()

ensembl <- useEnsembl(biomart = 'ensembl', host = 'http://may2021.archive.ensembl.org')


# Choosing dataset
listDatasets(ensembl) %>% filter(dataset == 'mmusculus_gene_ensembl')  # GRCm39

ensembl <- useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl)


# View available fields
View(listAttributes(ensembl))


# Save ensembl_gene_id to entrezgene_id conversion
mm_ensembl_entrez <- getBM(
  attributes = c('ensembl_gene_id', 'entrezgene_id'),
  mart = ensembl
)
mm_ensembl_entrez %>% count(is.na(entrezgene_id))
saveRDS(mm_ensembl_entrez, file = file.path(data_dir, 'mm_ensembl_entrez.GRCm39.RDS'))


# Save ensembl_gene_id to mgi_symbol conversion
mm_ensembl_symbol <- getBM(
  attributes = c('ensembl_gene_id', 'mgi_symbol'),
  mart = ensembl
)
mm_ensembl_symbol[mm_ensembl_symbol == ''] <- NA
mm_ensembl_symbol %>% count(is.na(mgi_symbol))
saveRDS(mm_ensembl_symbol, file = file.path(data_dir, 'mm_ensembl_symbol.GRCm39.RDS'))


# -----------------------
# Archived version (102)
# -----------------------

# Choosing database
listEnsembl()
listEnsemblArchives()

ensembl <- useEnsembl(biomart = 'ensembl', host = 'http://nov2020.archive.ensembl.org')


# Choosing dataset
listDatasets(ensembl) %>% filter(dataset == 'mmusculus_gene_ensembl')  # GRCm38.p6

ensembl <- useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl)


# View available fields
View(listAttributes(ensembl))


# Save ensembl_gene_id to entrezgene_id conversion
mm_ensembl_entrez <- getBM(
  attributes = c('ensembl_gene_id', 'entrezgene_id'),
  mart = ensembl
)
mm_ensembl_entrez %>% count(is.na(entrezgene_id))
saveRDS(mm_ensembl_entrez, file = file.path(data_dir, 'mm_ensembl_entrez.GRCm38.p6.RDS'))


# Save ensembl_gene_id to mgi_symbol conversion
mm_ensembl_symbol <- getBM(
  attributes = c('ensembl_gene_id', 'mgi_symbol'),
  mart = ensembl
)
mm_ensembl_symbol[mm_ensembl_symbol == ''] <- NA
mm_ensembl_symbol %>% count(is.na(mgi_symbol))
saveRDS(mm_ensembl_symbol, file = file.path(data_dir, 'mm_ensembl_symbol.GRCm38.p6.RDS'))
