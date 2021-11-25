library(tidyverse)
library(scales)
library(RColorBrewer)


# Trimmomatic logs
get_trimmomatic_logs <- function(dir_name) {
  trimmomatic_dir <- file.path(dir_name, 'trimmomatic')
  
  trimmomatic_logs <- list.files(trimmomatic_dir)
  
  trimmomatic_df <- data.frame(
    sample_id = character(),
    raw_cnt = numeric(0),
    trimmed_survived_cnt = numeric(0),
    trimmed_survived_pct = numeric(0),
    trimmed_forward_cnt = numeric(0),
    trimmed_forward_pct = numeric(0),
    trimmed_reverse_cnt = numeric(0),
    trimmed_reverse_pct = numeric(0),
    trimmed_dropped_cnt = numeric(0),
    trimmed_dropped_pct = numeric(0),
    stringsAsFactors = F
  )
  
  for (idx in seq_along(trimmomatic_logs)) {
    log_text <- read_file(file.path(trimmomatic_dir, trimmomatic_logs[[idx]]))
    
    trimmomatic_df[idx, 1] <- str_extract(trimmomatic_logs[[idx]], 'SRR\\d+')
    trimmomatic_df[idx, 2:10] <- str_match(log_text, 'Input Read Pairs: (\\d+) Both Surviving: (\\d+) \\(([\\d.]+)%\\) Forward Only Surviving: (\\d+) \\(([\\d.]+)%\\) Reverse Only Surviving: (\\d+) \\(([\\d.]+)%\\) Dropped: (\\d+) \\(([\\d.]+)%\\)')[1, 2:10] %>% as.numeric()
  }
  
  trimmomatic_df
}


# STAR logs
get_star_logs <- function(dir_name) {
  star_dir <- file.path(dir_name, 'star')
  
  star_logs <- list.files(star_dir)
  
  star_df <- data.frame(
    sample_id = character(),
    mapped_input_cnt = numeric(0),
    mapped_unique_cnt = numeric(0),
    mapped_unique_pct = numeric(0),
    mapped_multi_cnt = numeric(0),
    mapped_multi_pct = numeric(0),
    mapped_many_cnt = numeric(0),
    mapped_many_pct = numeric(0),
    stringsAsFactors = F
  )
  
  for (idx in seq_along(star_logs)) {
    log_text <- read_file(file.path(star_dir, star_logs[[idx]]))
    
    star_df[idx, 1] <- str_extract(star_logs[[idx]], 'SRR\\d+')
    star_df[idx, 2:8] <- str_match(log_text, 'Number of input reads \\|\\s+(\\d+)[\\s\\S]+Uniquely mapped reads number \\|\\s+(\\d+)\\s+Uniquely mapped reads % \\|\\s+([\\d.]+)[\\s\\S]+Number of reads mapped to multiple loci \\|\\s+(\\d+)\\s+% of reads mapped to multiple loci \\|\\s+([\\d.]+)[\\s\\S]+Number of reads mapped to too many loci \\|\\s+(\\d+)\\s+% of reads mapped to too many loci \\|\\s+([\\d.]+)[\\s\\S]+')[1, 2:8] %>% as.numeric()
  }
  
  star_df
}


# Picard logs
get_picard_logs <- function(dir_name) {
  picard_dir <- file.path(dir_name, 'picard')
  
  picard_logs <- list.files(picard_dir)
  
  picard_df <- data.frame(
    sample_id = character(),
    deduped_dropped_pct = numeric(0),
    stringsAsFactors = F
  )
  
  for (idx in seq_along(picard_logs)) {
    log_text <- read_file(file.path(picard_dir, picard_logs[[idx]]))
    
    picard_df[idx, 1] <- str_extract(picard_logs[[idx]], 'SRR\\d+')
    picard_df[idx, 2] <- round(str_match(log_text, 'Unknown Library[\\s\\d]+\\t([\\d.]+)')[1, 2] %>% as.numeric() * 100, 2)
  }
  
  picard_df %>% mutate(
    deduped_survived_pct = 100 - deduped_dropped_pct
  )
}


# featureCounts logs
get_featurecounts_logs <- function(dir_name) {

  log_text <- read_file(file.path(dir_name, 'featurecounts', 'counts_deduped.txt.summary'))
  
  sample_id <- str_extract_all(log_text, 'SRR\\d+', simplify = T)[1, ]
  
  assigned_text <- str_extract(log_text, 'Assigned(.+)')
  assigned_cnt <- str_extract_all(assigned_text, '\\d+', simplify = T)[1, ]  %>% as.numeric()
  
  data.frame(
    sample_id = sample_id,
    assigned_cnt = assigned_cnt,
    stringsAsFactors = F
  )
}


# Combine dataframes
get_samples_df <- function(meta_data, dir_name) {
  trimmomatic_df <- get_trimmomatic_logs(dir_name)
  star_df <- get_star_logs(dir_name)
  picard_df <- get_picard_logs(dir_name)
  featurecounts_df <- get_featurecounts_logs(dir_name)
  
  samples_df <- meta_data %>%
    dplyr::select(Run, treatment) %>%
    right_join(trimmomatic_df, by = c('Run' = 'sample_id')) %>% 
    mutate(
      trimmed_removed_cnt = raw_cnt - trimmed_survived_cnt,
      trimmed_removed_pct = round(trimmed_removed_cnt / raw_cnt * 100, 2)
    ) %>% 
    left_join(star_df, by = c('Run' = 'sample_id')) %>% 
    mutate(
      mapped_nonunique_cnt = mapped_multi_cnt + mapped_many_cnt,
      mapped_nonunique_pct = mapped_multi_pct + mapped_many_pct,
      mapped_mapped_cnt = mapped_unique_cnt + mapped_nonunique_cnt,
      mapped_mapped_pct = mapped_unique_pct + mapped_nonunique_pct,
      mapped_unmapped_cnt = mapped_input_cnt - (mapped_unique_cnt + mapped_multi_cnt + mapped_many_cnt),
      mapped_unmapped_pct = round(mapped_unmapped_cnt / mapped_input_cnt * 100, 2)
    ) %>% 
    left_join(picard_df, by = c('Run' = 'sample_id')) %>% 
    left_join(featurecounts_df, by = c('Run' = 'sample_id')) %>% 
    mutate(
      assigned_pct = round(assigned_cnt / raw_cnt * 100, 2)
    ) %>% 
    dplyr::rename('sample_id' = 'Run') 
  
  samples_df$sample_id <- factor(samples_df$sample_id, levels = rev(samples_df$sample_id))
  samples_df$trimmed_survived_cnt == samples_df$mapped_input_cnt
  
  samples_df
}


# Raw counts
get_raw_counts <- function(samples_df, dir_name, height = 3) {
  samples_df %>%
    select(sample_id, treatment, raw_cnt) %>%
    View()
  
  samples_df %>% 
    ggplot(aes(y = sample_id, x = raw_cnt, fill = treatment)) +
    geom_col(width = 0.8) +
    geom_text(aes(x = raw_cnt * 0.99, label = format(raw_cnt, big.mark = ',')), color = '#444444', size = 2.5, hjust = 1, vjust = 0.5, show.legend = F) +
    xlab('') + ylab('') + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(suffix = 'M', scale = 1e-6)) +
    scale_fill_manual(values = brewer.pal(3, 'YlGnBu'), name = 'Experimental conditions') +
    theme(
      text = element_text(size = 8, family = 'Times'),
      legend.title = element_text(size = 7, face = 'bold'),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(face = 'bold')
    )
  
  ggsave(file.path(dir_name, 'figures', 'figure_1.png'), width = 8, height = height, units = 'in')
}


# Trimmed counts
get_trimmed_counts <- function(samples_df, dir_name, height = 3) {
  samples_df %>%
    select(sample_id, raw_cnt, trimmed_survived_cnt, trimmed_survived_pct, trimmed_removed_cnt, trimmed_removed_pct) %>%
    View()
  
  samples_df %>%
    select(sample_id, treatment, trimmed_survived_cnt, trimmed_survived_pct, trimmed_removed_cnt, trimmed_removed_pct) %>%
    pivot_longer(
      -c(sample_id, treatment),
      names_to = c('trim_res', '.value'),
      names_pattern = 'trimmed_([a-z]+)_([a-z]+)'
    ) %>% 
    mutate(
      trim_res = recode(trim_res, 'survived' = 'Survived', 'removed' = 'Dropped')
    ) %>% 
    ggplot(aes(x = cnt, y = sample_id, alpha = trim_res, fill = treatment)) + 
    geom_bar(position = 'stack', stat = 'identity', width = 0.8) +
    geom_text(aes(x = cnt * 0.99, label = if_else(trim_res == 'Survived', str_c(sprintf('%.1f', pct), '%'), '')), color = '#444444', size = 2.5, hjust = 1, vjust = 0.5, show.legend = F) +
    xlab('') + ylab('') + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(suffix = 'M', scale = 1e-6)) +
    scale_fill_manual(values = brewer.pal(3, 'YlGnBu'), name = 'Experimental conditions') +
    scale_alpha_manual(values = c(0.5, 1), name = 'Trimming results') +
    theme(
      text = element_text(size = 8, family = 'Times'),
      legend.title = element_text(size = 7, face = 'bold'),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(face = 'bold')
    ) +
    guides(alpha = guide_legend(reverse = T, order = 2), fill = guide_legend(order = 1))
  
  ggsave(file.path(dir_name, 'figures', 'figure_2.png'), width = 8, height = height, units = 'in')
}


# Mapped counts
get_mapped_counts <- function(samples_df, dir_name, height = 3) {
  samples_df %>%
    select(sample_id, mapped_input_cnt, mapped_unmapped_pct, mapped_unique_pct, mapped_nonunique_pct, mapped_multi_pct, mapped_many_pct) %>%
    View()
  
  samples_df %>%
    select(sample_id, treatment, trimmed_removed_cnt, trimmed_removed_pct, mapped_unique_cnt, mapped_unique_pct, mapped_nonunique_cnt, mapped_nonunique_pct, mapped_unmapped_cnt, mapped_unmapped_pct) %>%
    pivot_longer(
      -c(sample_id, treatment),
      names_to = c('map_res', '.value'),
      names_pattern = '_([a-z]+)_([a-z]+)'
    ) %>% 
    mutate(
      map_res = recode_factor(map_res, 'removed' = 'Dropped', 'unmapped' = 'Unmapped', 'nonunique' = 'Multi-mapped', 'unique' = 'Uniquely mapped')
    ) %>% 
    ggplot(aes(x = cnt, y = sample_id, alpha = map_res, fill = treatment)) + 
    geom_bar(position = 'stack', stat = 'identity', width = 0.8) +
    geom_text(aes(x = cnt, label = if_else(map_res != 'Dropped', str_c(sprintf('%.1f', pct), '%'), '')), color = '#444444', size = 2.5, hjust = 0.5, position = position_stack(vjust = 0.5, reverse = F), show.legend = F) +
    xlab('') + ylab('') + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(suffix = 'M', scale = 1e-6)) +
    scale_fill_manual(values = brewer.pal(3, 'YlGnBu'), name = 'Experimental conditions') +
    scale_alpha_manual(values = c(0.5, 0.8, 0.9, 1), name = 'Mapping results') +
    theme(
      text = element_text(size = 8, family = 'Times'),
      legend.title = element_text(size = 7, face = 'bold'),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(face = 'bold')
    ) +
    guides(alpha = guide_legend(reverse = T, order = 2), fill = guide_legend(order = 1))
  
  ggsave(file.path(dir_name, 'figures', 'figure_3.png'), width = 8, height = height, units = 'in')
}


# Deduped counts
get_deduped_counts <- function(samples_df, dir_name, height = 3) {
  samples_df %>%
    select(sample_id, mapped_mapped_cnt, deduped_survived_pct, deduped_dropped_pct) %>%
    View()
  
  samples_df %>%
    select(sample_id, treatment, trimmed_removed_cnt, trimmed_removed_pct, mapped_unmapped_cnt, mapped_unmapped_pct, mapped_mapped_cnt, deduped_survived_pct, deduped_dropped_pct) %>%
    mutate(
      deduped_survived_cnt = deduped_survived_pct / 100 * mapped_mapped_cnt,
      deduped_dropped_cnt = deduped_dropped_pct / 100 * mapped_mapped_cnt
    ) %>% 
    select(-mapped_mapped_cnt) %>% 
    pivot_longer(
      -c(sample_id, treatment),
      names_to = c('dedupe_res', '.value'),
      names_pattern = '_([a-z]+)_([a-z]+)'
    ) %>% 
    mutate(
      dedupe_res = recode_factor(dedupe_res, 'removed' = 'Dropped', 'unmapped' = 'Unmapped', 'dropped' = 'Duplicates', 'survived' = 'Unique')
    ) %>% 
    ggplot(aes(x = cnt, y = sample_id, alpha = dedupe_res, fill = treatment)) + 
    geom_bar(position = 'stack', stat = 'identity', width = 0.8) +
    geom_text(aes(x = cnt, label = if_else(!dedupe_res %in% c('Dropped', 'Unmapped'), str_c(sprintf('%.1f', pct), '%'), '')), color = '#444444', size = 2.5, hjust = 0.5, position = position_stack(vjust = 0.5, reverse = F), show.legend = F) +
    xlab('') + ylab('') + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(suffix = 'M', scale = 1e-6)) +
    scale_fill_manual(values = brewer.pal(3, 'YlGnBu'), name = 'Experimental conditions') +
    scale_alpha_manual(values = c(0.5, 0.6, 0.9, 1), name = 'Deduping results') +
    theme(
      text = element_text(size = 8, family = 'Times'),
      legend.title = element_text(size = 7, face = 'bold'),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(face = 'bold')
    ) +
    guides(alpha = guide_legend(reverse = T, order = 2), fill = guide_legend(order = 1))
  
  ggsave(file.path(dir_name, 'figures', 'figure_4.png'), width = 8, height = height, units = 'in')
}


# Assigned counts
get_assigned_counts <- function(samples_df, dir_name, height = 3) {
  samples_df %>%
    select(sample_id, raw_cnt, assigned_cnt, assigned_pct) %>%
    View()
  
  samples_df %>%
    select(sample_id, treatment, raw_cnt, assigned_cnt, assigned_pct) %>%
    mutate(
      rest_cnt = raw_cnt - assigned_cnt,
      rest_pct = 100 - assigned_pct
    ) %>% 
    select(-raw_cnt) %>% 
    pivot_longer(
      -c(sample_id, treatment),
      names_to = c('assign_res', '.value'),
      names_pattern = '([a-z]+)_([a-z]+)'
    ) %>% 
    mutate(
      assign_res = recode_factor(assign_res, 'rest' = 'Total', 'assigned' = 'Assigned')
    ) %>% 
    ggplot(aes(x = cnt, y = sample_id, alpha = assign_res, fill = treatment)) + 
    geom_bar(position = 'stack', stat = 'identity', width = 0.8) +
    geom_text(aes(x = cnt * 1.02, label = if_else(assign_res == 'Assigned', str_c(sprintf('%.1f', pct), '%'), '')), color = '#444444', size = 2.5, hjust = 0, vjust = 0.5, show.legend = F) +
    xlab('') + ylab('') + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.01)), labels = label_number(suffix = 'M', scale = 1e-6)) +
    scale_fill_manual(values = brewer.pal(3, 'YlGnBu'), name = 'Experimental conditions') +
    scale_alpha_manual(values = c(0.5, 1), name = 'Quantification results') +
    theme(
      text = element_text(size = 8, family = 'Times'),
      legend.title = element_text(size = 7, face = 'bold'),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(face = 'bold')
    ) +
    guides(alpha = guide_legend(reverse = T, order = 2), fill = guide_legend(order = 1))
  
  ggsave(file.path(dir_name, 'figures', 'figure_5.png'), width = 8, height = height, units = 'in')
}


# Generate figures
meta_data_1 <- read_csv('PRJNA694971_meta.txt') %>%
  arrange(Run) %>% 
  select(-treatment) %>% 
  dplyr::rename('treatment' = 'experimental_class_type')
meta_data_2 <- read_csv('PRJNA748404_meta.txt') %>% 
  mutate(treatment = recode(treatment, 'counterpart' = 'Control', 'post treatment 16 weeks' = 'Diabetic'))

# scp chan@spartan01.sjsu.edu:/home/chan/diabetes_mellitus/outputs/star/*.Log.final.out diabetes_mellitus/star/
# scp chan@spartan01.sjsu.edu:/home/chan/diabetes_mellitus/outputs/slurm/trimmomatic* diabetes_mellitus/trimmomatic/
# scp chan@spartan01.sjsu.edu:/home/chan/diabetes_mellitus/outputs/picard/*_deduped.log diabetes_mellitus/picard/
# scp chan@spartan01.sjsu.edu:/home/chan/diabetes_mellitus/outputs/featurecounts/counts_deduped.txt.summary diabetes_mellitus/featurecounts/
# scp chan@spartan01.sjsu.edu:/home/chan/mice_microgravity/outputs_primary/star/*.Log.final.out mice_microgravity/star/
# scp chan@spartan01.sjsu.edu:/home/chan/mice_microgravity/outputs_primary/slurm/trimmomatic* mice_microgravity/trimmomatic/
# scp chan@spartan01.sjsu.edu:/home/chan/mice_microgravity/outputs_primary/picard/*_deduped.log mice_microgravity/picard/
# scp chan@spartan01.sjsu.edu:/home/chan/mice_microgravity/outputs_primary/featurecounts/counts_deduped.txt.summary mice_microgravity/featurecounts/

dir_name_1 <- 'mice_microgravity'
dir_name_2 <- 'diabetes_mellitus'


samples_df_1 <- get_samples_df(meta_data_1, dir_name_1)

get_raw_counts(samples_df_1, dir_name_1)
get_trimmed_counts(samples_df_1, dir_name_1)
get_mapped_counts(samples_df_1, dir_name_1)
get_deduped_counts(samples_df_1, dir_name_1)
get_assigned_counts(samples_df_1, dir_name_1)


samples_df_2 <- get_samples_df(meta_data_2, dir_name_2)

get_raw_counts(samples_df_2, dir_name_2, height = 2.5)
get_trimmed_counts(samples_df_2, dir_name_2, height = 2.5)
get_mapped_counts(samples_df_2, dir_name_2, height = 2.5)
get_deduped_counts(samples_df_2, dir_name_2, height = 2.5)
get_assigned_counts(samples_df_2, dir_name_2, height = 2.5)
