library(forcats)
library(ggpmisc)
library(tidyverse)
library(janitor)
library(ranger)
library(vip)
library(UpSetR)
set.seed(2024)

#### DATA & FUNCTION SETUP ####

get_devel_data <- function() {
  # This function imports fly development % data and formats it for use with miR count data
  
  # This file lists the % development of each fly sample at each pupal stage
  devel_file <- list.files(pattern = 'sample_percent_development_by_species.csv', recursive = TRUE)
  devel_data <- read_csv(devel_file)
  
  # make the names match the data folders
  devel_data <- devel_data %>% 
    clean_names() %>% 
    rename(Cochliomyia = c_macellaria) %>% 
    rename(Chrysomya = c_rufifacies) %>% 
    rename(Lucilia = l_cuprina) %>% 
    rename(Plinthopyga = b_plinthopyga)
  
  # convert chr % values to NAs/decimals, as appropriate
  devel_data <- devel_data %>% 
    mutate(across(-stage, ~na_if(., "-"))) %>% 
    mutate(across(-stage, ~gsub("%", "", .))) %>% 
    mutate(across(-stage, ~as.numeric(.))) %>% 
    mutate(across(-stage, ~./100))
}

devel_df <- get_devel_data()
str(devel_df)


norm_data_cleaner <- function(filename, signames=NULL) {
  # filename should be a complete, unique path string of a normalized miRNA count file.
  # signames is an optional chr vector of miRNA names to analyze.
  
  # read in data and standardize 1st column's name
  normdata <- read_csv(filename)
  normdata <- normdata %>% 
    rename(mir = 1) 
  
  # remove all non-significant miR data, if sig names are provided
  if(!is.null(signames)) {
    normdata <- normdata %>% 
      filter(mir %in% signames)
  }
  
  # flip the x & y axes and name them appropriately
  normdata <- normdata %>% 
    mutate(mir = paste0('mir_', mir)) %>% 
    t() %>%
    as.data.frame() %>% 
    row_to_names(1) %>% 
    rownames_to_column(var = 'stage') %>% 
    mutate(stage = gsub('\\.\\d+$', '', stage)) %>% 
    clean_names()
  
  # make development stage naming consistent
  normdata <- normdata %>% 
    mutate(stage = case_when(grepl("Early", stage) | grepl("early", stage) ~ "Early",
                             grepl("Mid", stage) | grepl("mid", stage) ~ "Mid",
                             grepl("Late", stage) | grepl("late", stage) ~ "Late",
                             TRUE ~ "unrecognized stage")) %>% 
    mutate(stage = factor(stage)) %>% 
    mutate_if(is.character, as.numeric) # convert miR reads into numbers
  
  # add species
  species_name = gsub('/.*?.csv', '', filename)
  normdata$species = species_name
  
  # add development percentage
  normdata$development <- devel_df[[species_name]][!is.na(devel_df[[species_name]])]
  # WARNING: devel_df must be in the global code environment for this to work!
  
  # convert lifetime development % to puparial development %
  min_devel = min(normdata$development)
  normdata <- normdata %>% 
    mutate(adjusted_devel = (development - min_devel)/(1 - min_devel))
  
  return(normdata)
}


#test_file <- list.files(pattern = 'chry_allrufcounts_py_consol_normalized.csv', recursive = TRUE)
#test_df <- norm_data_cleaner(test_file)


# These files house all normalized wildtype fly miRNA counts, including "low-quality" samples
wt_files <- c('chry_allrufcounts_py_consol_normalized.csv',
              'coch_cmacwtcombined_nolarvae_py_consol_normalized.csv',
              'plinthocountsnosex_py_consol_normalized.csv',
              'luciliawt_py_consol_normalized.csv')

species_names <- c('Chrysomya', 'Cochliomyia', 'Plinthopyga', 'Lucilia')
# WARNING: name order here^ must match order of wt_files list


# make a list of full file paths, based on the above file names
wt_filepaths <- 
  wt_files %>%
  map(\(x) list.files(pattern=x, recursive=TRUE))

# make a data frame from each species' file and collect the frames in a list
df_list <- 
  wt_filepaths %>% 
  map(norm_data_cleaner)
names(df_list) <- species_names

# make that list of individual data frames into one big frame
alldat <- reduce(df_list, full_join)

result_file <- paste0(getwd(), '/all_norm_wt_data.csv')
write_csv(alldat, result_file)


# grow a random forest that models % development across all species
#all_devel_mod <- 
#  alldat %>% 
#  select(-c(stage, adjusted_devel)) %>% 
#  select_if(~ !any(is.na(.))) %>% # eliminate columns with NAs
#  ranger(formula = development ~ ., importance = 'permutation')
#vars <- vi(all_devel_mod)
#vip(vars, num_features = 25) +
#  ggtitle(paste0('Most Predictive Variables Across Species'))


forest_grower <- function(dataframe, species_name, devel=NULL) {
  # grows a random forest that models species_name's % development
  # dataframe must have been through norm_data_cleaner() to work here
  # devel is optional; if 'adjusted', analysis will use adjusted_devel

  devel_mod <- 
    dataframe %>% 
    filter(species == species_name) %>% 
    select_if(~ !any(is.na(.))) # eliminate columns with NAs
  
  # p represents the number of independent variables that the trees can branch on
  p <- devel_mod %>% 
    select(starts_with('mir_')) %>% 
    ncol()

  if (is.null(devel)) {
    devel_mod <- devel_mod %>% 
      select(-c(stage, species, adjusted_devel)) %>% 
      ranger(formula = development ~ ., 
             importance = 'permutation', # we're modeling regression, not classification
             min.node.size = 4, # fits our sample groupings better than default 5
             mtry = floor(p/3), # recommended m for regression is p/3, not default p^0.5
             scale.permutation.importance = TRUE) # seemed like a good idea?
  } else if (devel == 'adjusted') {
    devel_mod <- devel_mod %>% 
      select(-c(stage, species, development)) %>% 
      ranger(formula = adjusted_devel ~ ., 
             importance = 'permutation', 
             min.node.size = 4,
             mtry = floor(p/3),
             scale.permutation.importance = TRUE)
  }

  return(devel_mod)
}


#### SINGLE FOREST RESULTS ####

# save graphs showing top 10 predictive miRs by species' total development %
# from a single forest
pdf(paste0(getwd(), '/randomforest_devel_predictors_by_species.pdf'))
for (name in species_names) {
  model <- forest_grower(alldat, name)
  graph <- vip(model, num_features = 25) +
    ggtitle(paste0('Most Predictive miRNAs of ', name, ' Total Development'))
  plot(graph)
}
dev.off()

# save graphs showing top 10 predictive miRs by species' adjusted development %
# from a single forest
pdf(paste0(getwd(), '/randomforest_adj_devel_predictors_by_species.pdf'))
for (name in species_names) {
  model <- forest_grower(alldat, name, devel='adjusted')
  graph <- vip(model, num_features = 25) +
    ggtitle(paste0('Most Predictive miRNAs of ', name, ' Puparial Development'))
  plot(graph)
}
dev.off()

# save graphs showing how well random forest model predicts development %
# from a single forest
pdf(paste0(getwd(), '/randomforest_devel_predictions_by_species.pdf'))
for (name in species_names) {
  model <- forest_grower(df_list[[name]], name)
  preds <- predict(model, df_list[[name]])
  df_list[[name]]$predicted_devel <- preds$predictions
  graph <- df_list[[name]] %>% 
    ggplot(aes(x=development, y=predicted_devel)) +
    ggtitle(paste0('Predicted versus Actual Total Development of ', name)) +
    labs(x='Actual Development', y='Predicted Development') +
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "R2"))) +
    geom_point()
  plot(graph)
  df_list[[name]]$predicted_devel <- NULL
}
dev.off()

# save graphs showing how well random forest model predicts adjusted development %
# from a single forest
pdf(paste0(getwd(), '/randomforest_adj_devel_predictions_by_species.pdf'))
for (name in species_names) {
  model <- forest_grower(df_list[[name]], name, devel='adjusted')
  preds <- predict(model, df_list[[name]])
  df_list[[name]]$predicted_devel <- preds$predictions
  graph <- df_list[[name]] %>% 
    ggplot(aes(x=adjusted_devel, y=predicted_devel)) +
    ggtitle(paste0('Predicted versus Actual Puparial Development of ', name)) +
    labs(x='Actual Development', y='Predicted Development') +
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "R2"))) +
    geom_point() +
    lims(x=c(0,1.05), y=c(0,1.05))
  plot(graph)
  df_list[[name]]$predicted_devel <- NULL
}
dev.off()


# determine which miRNAs are important across species
# from a single forest
signame_list <- list()
for (name in species_names) {
  model <- forest_grower(alldat, name, devel='adjusted')
  print(model)
  name_import <- vi(model)
  name_import <- head(name_import, n=25)
  signame_list[[name]] <- name_import$Variable
}
upset(fromList(signame_list), nsets=length(signame_list))
mir_allspecies <- Reduce(intersect, signame_list)


#### MULTI FOREST RESULTS ####

# running randomForest a bunch of times
#chry_forests <- list()
#for(i in 1:10){
#chry_forests[[i]] <- forest_grower(alldat,species_name = species_names[1])
#}
#vi_list_chry <- purrr::map(chry_forests,function(z){vi(z)[['Variable']] %>% head(25)})


# getting a bunch of randomForest results in one big table
vip_mirs <- list()
vip_mir_counts <- list()

for(species in species_names){
  vip_mirs[[species]] <- list()
  
  # run randomForest a bunch of times, grabbing the top 25 miRs each time
  for(i in 1:100){
    forest <- forest_grower(alldat, species_name=species, devel='adjusted')
    vip_mirs[[species]][[i]] <- 
      vi(forest)[['Variable']] %>% 
      head(25)
  }
  
  # consolidate the top miR results into a frequency table
  vip_mir_counts[[species]] <- 
    vip_mirs[[species]] %>% 
    unlist() %>% 
    table() %>% 
    data.frame() %>% 
    clean_names() %>% 
    rename(mir = x) %>% 
    mutate(species = species) %>% 
    arrange(desc(freq))
}

# combine randomForest results from all species into one table
vip_mir_allcounts <- 
  reduce(vip_mir_counts, full_join)

# sum how many times each miR is among the top variables across all species
total_mir_counts <- 
  vip_mir_allcounts %>% 
    group_by(mir) %>% 
    summarise(total_freq = sum(freq)) %>% 
    arrange(desc(total_freq))

# make the total miR sums a new column in the table of individual counts
vip_mir_allcounts <- 
  vip_mir_allcounts %>% 
    full_join(total_mir_counts) %>% 
  arrange(desc(total_freq))

result_file <- paste0(getwd(), '/randomforest_adj_devel_400_trees.csv')
write_csv(vip_mir_allcounts, result_file)


# visualize how many times each miR appeared in each species
# but only the miRs that appeared at least 100 summed times across all species
vip_mir_allcounts %>% 
  filter(total_freq >= 100) %>% 
  #filter(freq > 95) %>% 
  ggplot(aes(x=fct_reorder(mir, total_freq), y=freq)) +
  geom_col() +
  facet_wrap(~species, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

# visualize how many species each miR appeared in
# but only the miRs that appeared at least 35 times in a species
vip_mir_allcounts %>% 
  filter(freq >= 35) %>% 
  ggplot(aes(x=fct_reorder(mir, total_freq))) +
  geom_bar() +
  coord_flip()


# visualize how many miRs appeared in each grouping of species
# then print the miRs that appeared in all species AND in at least 300/400 forests
forest_name_results <- list()
for (species in species_names){
  forest_name_results[[species]] <- as.character(vip_mir_counts[[species]][["mir"]])
}

pdf(paste0(getwd(), '/randomforest_adj_devel_upset.pdf'))
forest_upset <- upset(fromList(forest_name_results), nsets=length(forest_name_results))
forest_upset
dev.off()

upset_all_species <- reduce(forest_name_results, intersect)
frequent_mirs <- 
  total_mir_counts %>% 
  filter(total_freq >= 300) %>% 
  pull(mir)
freq_upset_adj_devel_mirs <- intersect(upset_all_species, frequent_mirs)
#freq_upset_devel_mirs


# visualize miRNA expression patterns for those intersecting frequent miRNAs
adj_plot <- alldat %>% 
  select(c(mir_315a, mir_375, mir_305, mir_9a, species, adjusted_devel)) %>% 
  pivot_longer(starts_with('mir_'), 
               names_to = 'mir', 
               values_to = 'norm_read') %>% 
  ggplot(aes(x=adjusted_devel, y=norm_read, color=species)) +
  geom_smooth(size=2, se = FALSE) +
  facet_wrap(~mir, scales = 'free')
ggsave(paste0(getwd(), '/norm_data_by_adj_devel_percent_rf.pdf'), adj_plot, scale = 3)


