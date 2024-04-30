library(tidyverse)
library(janitor)

devel_file <- list.files(pattern = 'sample_percent_development_by_species.csv', recursive = TRUE)
devel_df <- read_csv(devel_file)

devel_df <- devel_df %>% 
  clean_names() %>% 
  rename(Cochliomyia = c_macellaria) %>% 
  rename(Chrysomya = c_rufifacies) %>% 
  rename(Lucilia = l_cuprina) %>% 
  rename(Plinthopyga = b_plinthopyga)

devel_df <- devel_df %>% 
  mutate(across(-stage, ~na_if(., "-"))) %>% 
  mutate(across(-stage, ~gsub("%", "", .))) %>% 
  mutate(across(-stage, ~as.numeric(.)))

str(devel_df)


norm_data_cleaner <- function(filename, signames=NULL) {
  # filename should be a complete, unique path string.
  # signames is an optional chr vector of miRNA names to analyze.
  
  # read in data and standardize 1st column's name
  normdata <- read_csv(filename)
  normdata <- normdata %>% 
    rename(miR = 1) 
  
  # remove all non-significant miR data, if sig names are provided
  if(!is.null(signames)) {
    normdata <- normdata %>% 
      filter(miR %in% signames)
  }
  
  # flip the x & y axes and name them appropriately
  normdata <- normdata %>% 
    mutate(miR = paste0('miR_', miR)) %>% 
    t() %>%
    as.data.frame() %>% 
    row_to_names(1) %>% 
    rownames_to_column(var = 'stage') %>% 
    mutate(stage = gsub('\\.\\d+$', '', stage))
  
  # make development stage naming consistent
  normdata <- normdata %>% 
    mutate(stage = case_when(grepl("Early", stage) | grepl("early", stage) ~ "Early",
                             grepl("Mid", stage) | grepl("mid", stage) ~ "Mid",
                             grepl("Late", stage) | grepl("late", stage) ~ "Late",
                             TRUE ~ "unrecognized stage"))
  
  # make beginning < middle < end (avoids alphabetical ordering)
  ordlevels <- unique(normdata$stage)
  normdata <- normdata %>%
    mutate(stage = factor(stage, levels=ordlevels, ordered=FALSE)) %>% 
    mutate_if(is.character, as.numeric)
  
  # add species
  normdata$species = gsub('/.*?.csv', '', filename)
  
  return(normdata)
}


infile <- list.files(pattern = 'chry_allrufcounts_py_consol_normalized.csv', recursive = TRUE)
chry_dat <- norm_data_cleaner(infile)
chry_dat$development <- devel_df$Chrysomya[!is.na(devel_df$Chrysomya)]

infile <- list.files(pattern = 'coch_cmacwtcombined_nolarvae_py_consol_normalized.csv', recursive = TRUE)
coch_dat <- norm_data_cleaner(infile)
coch_dat$development <- devel_df$Cochliomyia[!is.na(devel_df$Cochliomyia)]

infile <- list.files(pattern = 'plinthocountsnosex_py_consol_normalized.csv', recursive = TRUE)
plin_dat <- norm_data_cleaner(infile)
plin_dat$development <- devel_df$Plinthopyga[!is.na(devel_df$Plinthopyga)]

infile <- list.files(pattern = 'luciliawt_py_consol_normalized.csv', recursive = TRUE)
luci_dat <- norm_data_cleaner(infile)
luci_dat$development <- devel_df$Lucilia[!is.na(devel_df$Lucilia)]


chry_coch <- full_join(chry_dat, coch_dat)
chry_coch_plin <- full_join(chry_coch, plin_dat)
all_species_df <- full_join(chry_coch_plin, luci_dat)



all_species_df %>% 
  select(c(miR_100, miR_125, `miR_let-7`, miR_bft, species, development)) %>% 
  pivot_longer(starts_with('miR_'), 
               names_to = 'miR', 
               values_to = 'norm_read') %>% 
  ggplot(aes(x=development, y=norm_read, color=species)) +
  geom_smooth(size=2, se = FALSE) +
  facet_wrap(~miR, scales = 'free')


minimums <- all_species_df %>%
  group_by(species, stage) %>% 
  filter(stage=='Early') %>% 
  summarize(min_devel=min(development))
minimums <- minimums %>% 
  select(-stage)

all_species_df <- full_join(all_species_df, minimums)

all_species_df <- all_species_df %>% 
  mutate(adjusted_devel = (development - min_devel)/(100 - min_devel))



alldat_file <- list.files(pattern = 'all_norm_wt_data.csv', recursive = TRUE)
alldat <- read_csv(alldat_file)

#adj_plot <- all_species_df %>% 
#deseq_mirs_expression <- alldat %>% 
rf_mirs_expression <- alldat %>% 
  clean_names() %>% 
  #select(c(mir_100, mir_125, mir_let_7, mir_bft, species, adjusted_devel)) %>% 
  select(c(mir_315a, mir_375, mir_305, mir_9a, species, adjusted_devel)) %>% 
  pivot_longer(starts_with('miR_'), 
               names_to = 'miR', 
               values_to = 'norm_read') %>% 
  ggplot(aes(x=adjusted_devel, y=norm_read, color=species)) +
  geom_smooth(linewidth=2, se = FALSE) +
  facet_wrap(~miR, scales = 'free') +
  labs(x='Adjusted Pupal Development (0-100%)', 
       y='Normalized miR Count', 
       title='Random Forest miR Expression across Pupal Development, by Species') +
  scale_color_viridis_d(end=0.85)

#ggsave('norm_data_by_adj_devel_percent.pdf', deseq_mirs_expression, scale = 2.5)
ggsave('norm_data_by_adj_devel_percent_rf.pdf', rf_mirs_expression, scale = 2.5)

