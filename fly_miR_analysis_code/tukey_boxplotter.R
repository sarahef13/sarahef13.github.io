library(tidyverse)
library(janitor)
library(ggplot2)
library(multcompView)

#### OLD CODE ####

#infile <- list.files(pattern = 'chry_allrufcounts_py_consol_normalized', recursive = TRUE)
#chry1 <- read_csv(infile)

# Get normalized data into a usable format for this analysis:
chry_cleaned <- chry1 %>% 
  rename(miR = 1) %>% 
  filter(miR %in% c('274', '100', '285', '125')) %>% 
  mutate(miR = paste0('miR_', miR)) %>% 
  t() %>%
  as.data.frame() %>% 
  row_to_names(1) %>% 
  rownames_to_column(var = 'stage') %>% 
  mutate(stage = gsub('\\.\\d+$', '', stage)) %>% 
  mutate(stage = factor(stage, levels=c("Early", "Mid", "Late"), ordered=TRUE)) %>% 
  mutate_if(is.character, as.numeric)

colnames(coch1[1])
?rename


model1=lm(chry_cleaned$miR_100 ~ chry_cleaned$stage)
ANOVA1=aov(model1)

TUKEY <- TukeyHSD(x=ANOVA1, 'chry_cleaned$stage', conf.level=0.95) # Tukey test to study each pair of treatment
labels <- generate_label_df(TUKEY, 'chry_cleaned$stage') # generate labels using function
names(labels) <- c('Letters','stage') #rename columns for merging
yvalue <- aggregate(.~stage, data=chry_cleaned, mean) # obtain letter position for y axis using means
final1 <- merge(labels, yvalue) # merge dataframes

ggplot(chry_cleaned, aes(x = stage, y = miR_100)) +
  labs(x = 'Stage', y = 'Normalized Expression') +
  ggtitle("Chrysomya miR_100 expression")+
  theme(plot.title = element_text(face='bold')) +
  geom_boxplot(fill = 'slateblue', alpha=0.2, stat = "boxplot") +
  geom_text(data = final1, aes(x = stage, y = miR_100, label = Letters), vjust=-8, hjust=-0.5) +
  stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme(plot.title = element_text(vjust=-0.6))

str(final1)

#### NEW CODE ####

# process _py_consol_normalized.csv files to be usable for Tukey:
norm_data_cleaner <- function(filename, signames) {
  # filename should be a complete, unique path string
  # signames should be a chr vector of miR to analyze
  normdata <- read_csv(filename)
  normdata <- normdata %>% 
    rename(miR = 1) %>% 
    filter(miR %in% signames) %>% 
    mutate(miR = paste0('miR_', miR)) %>% 
    t() %>%
    as.data.frame() %>% 
    row_to_names(1) %>% 
    rownames_to_column(var = 'stage') %>% 
    mutate(stage = gsub('\\.\\d+$', '', stage))
  ordlevels <- unique(normdata$stage)
  normdata <- normdata %>%
    mutate(stage = factor(stage, levels=ordlevels, ordered=TRUE)) %>% 
    mutate_if(is.character, as.numeric)
  return(normdata)
}


generate_label_df <- function(TUKEY, variable) {
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  #TODO: put the labels in the same order as in the boxplot
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}


single_mir_analyzer <- function(mircol, indf) {
  # indf is a cleaned df of normalized mir data
  # mircol is the column of the cleaned data to analyze
  model1 = lm(indf[[mircol]] ~ indf$stage)
  ANOVA1 = aov(model1)
  
  model1=lm(indf[[mircol]] ~ indf$stage)
  ANOVA1=aov(model1)
  
  TUKEY <- TukeyHSD(x=ANOVA1, 'indf$stage', conf.level=0.95) # Tukey test to study each pair of treatment
  labels <- generate_label_df(TUKEY, 'indf$stage') # generate labels using function
  names(labels) <- c('Letters','stage') #rename columns for merging
  yvalue <- aggregate(.~stage, data=indf, max) # obtain letter position for y axis using means
  final1 <- merge(labels, yvalue) # merge dataframes
  return(final1)
}


single_mir_grapher <- function(mircol, final1, prefinal, speciesname) {
  col_to_graph <- colnames(final1)[mircol]
  tukey_plot <- ggplot(prefinal, aes(x = stage, y = !!sym(col_to_graph))) +
    geom_boxplot(fill = 'slateblue', alpha=0.2, stat = "boxplot") +
    labs(x = 'Stage', y = 'Normalized Expression') +
    ggtitle(paste0(speciesname, " ", colnames(final1)[mircol], " expression"))+
    theme(plot.title = element_text(face='bold')) +
    geom_text(data = final1, aes(x = stage, y = !!sym(col_to_graph)+0.5, label = Letters), hjust=-0.5) +
    stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme(plot.title = element_text(vjust=-0.6))
  return(tukey_plot)
}


#### MAKING PLOTS FROM DATA FILES ####

#infile <- list.files(pattern = 'plinthocountsnosex_py_consol_normalized.csv', recursive = TRUE)
#infile <- list.files(pattern = 'chry_allrufcounts2_py_consol_normalized.csv', recursive = TRUE)
#infile <- list.files(pattern = 'coch_cmacwtcombined_py_consol_normalized.csv', recursive = TRUE)
#infile <- list.files(pattern = 'coch_cmacwtcombined_nolarvae_py_consol_normalized.csv', recursive = TRUE)
#infile <- list.files(pattern = 'chry_allrufcounts_py_consol_normalized.csv', recursive = TRUE)
infile <- list.files(pattern = 'luciliawt_py_consol_normalized.csv', recursive = TRUE)
#innames <- c('274', '100', '285', '125')
#innames <- c("let-7", "125", "100", "bft") #for Plinthopyga
#innames <- c("let-7", "125", "100", "bft", "315a") #for Chrysomya & Cochliomyia
innames <- c("bft", "92b") #for low quality Lucilia
outdf <- norm_data_cleaner(infile, innames)
species <- gsub('/.*?.csv', '', infile)


pdf_name <- gsub('\\.csv', '_tukey.pdf', infile)
pdf(pdf_name)
#for (colname in names(outdf[-1])) {
for (i in 2:length(outdf)) {
  analysis <- single_mir_analyzer(i, outdf)
  graph <- single_mir_grapher(i+1, analysis, outdf, species)
  plot(graph)
}
dev.off()


#model1=lm(outdf[["miR_100"]] ~ outdf$stage)
#ANOVA1=aov(outdf[["miR_100"]] ~ outdf$stage)
#summary(ANOVA1)


#### MAKING PLOTS OF ALL SIGNIFICANT MIRS BY SPECIES ####

name_adder <- function(filename) {
  species <- gsub('\\_py_consol_sig.csv', '', filename)
  sig_full <- read_csv(filename)
  colnames(sig_full)[1] <- 'miR'
  sig_names <- sig_full$miR
  return(sig_names)
}

#signames_file <- list.files(pattern = 'plinthocountsnosex_py_consol_sig.csv', recursive = TRUE)
#infile <- list.files(pattern = 'plinthocountsnosex_py_consol_normalized.csv', recursive = TRUE)

#signames_file <- list.files(pattern = 'chry_allrufcounts2_py_consol_sig.csv', recursive = TRUE)
#infile <- list.files(pattern = 'chry_allrufcounts2_py_consol_normalized.csv', recursive = TRUE)

#signames_file <- list.files(pattern = 'coch_cmacwtcombined_nolarvae_py_consol_sig.csv', recursive = TRUE)
#infile <- list.files(pattern = 'coch_cmacwtcombined_nolarvae_py_consol_normalized.csv', recursive = TRUE)

#signames_file <- list.files(pattern = 'luciliawt_py_consol_sig.csv', recursive = TRUE)
#infile <- list.files(pattern = 'luciliawt_py_consol_normalized.csv', recursive = TRUE)

signames_file <- list.files(pattern = 'chry_allrufcounts_py_consol_sig.csv', recursive = TRUE)
infile <- list.files(pattern = 'chry_allrufcounts_py_consol_normalized.csv', recursive = TRUE)

innames <- name_adder(signames_file)
outdf <- norm_data_cleaner(infile, innames)
species <- gsub('/.*?.csv', '', infile)


pdf_name <- gsub('\\.csv', '_tukey_allsigs.pdf', infile)
pdf(pdf_name)
for (i in 2:length(outdf)) {
  analysis <- single_mir_analyzer(i, outdf)
  graph <- single_mir_grapher(i+1, analysis, outdf, species)
  plot(graph)
}
dev.off()


#### MAKING PLOTS BASED ON RANDOM FOREST RESULTS ####

innames <- c("315a", "375", "305", "9a")
#infile <- list.files(pattern = 'chry_allrufcounts_py_consol_normalized.csv', recursive = TRUE)
#infile <- list.files(pattern = 'coch_cmacwtcombined_nolarvae_py_consol_normalized.csv', recursive = TRUE)
#infile <- list.files(pattern = 'plinthocountsnosex_py_consol_normalized.csv', recursive = TRUE)
infile <- list.files(pattern = 'luciliawt_py_consol_normalized.csv', recursive = TRUE)

outdf <- norm_data_cleaner(infile, innames)
species <- gsub('/.*?.csv', '', infile)

pdf_name <- gsub('\\.csv', '_tukey_rf.pdf', infile)
pdf(pdf_name)
for (i in 2:length(outdf)) {
  analysis <- single_mir_analyzer(i, outdf)
  graph <- single_mir_grapher(i+1, analysis, outdf, species)
  plot(graph)
}
dev.off()

