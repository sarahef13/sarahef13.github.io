library(tidyverse)
library(UpSetR)

####Functions that set up data for upset()####

name_adder <- function(filename) {
  species <- gsub('\\_py_consol_sig.csv', '', filename)
  sig_full <- read_csv(filename)
  #sig_names <- as.data.frame(sig_full$...1)
  colnames(sig_full)[1] <- 'miR'
  sig_names <- sig_full$miR
  #colnames(sig_names) <- species
  #assign(species, sig_names)
  return(sig_names)
}

#test case for name_adder():
#filename <- list.files(pattern = 'chry_allrufcounts_py_consol_sig.csv', recursive = TRUE)

get_names <- function(filename) {
  species <- gsub('\\_py_consol_sig.csv', '', filename)
  return (species)
}

name_compiler <- function(infile) {
  files_to_analyze <- list.files(pattern = infile, recursive = TRUE)
  sig_names <- lapply(files_to_analyze, name_adder)
  file_names <- lapply(files_to_analyze, get_names)
  names_together <- setNames(sig_names, file_names)
  return(names_together)
}


####Upset plot for all species' significant consolidated miR names####

all_names_together <- name_compiler('py_consol_sig.csv')
#all_overlap_names <- Reduce(intersect, all_names_together)

pdf(paste0(getwd(), '/all_consol_upset.pdf'))
upset(fromList(all_names_together), nsets=length(all_names_together))
dev.off()


####Upset plot for Chrysomya's significant consolidated miR names####

chry_names_together <- name_compiler('chry_.*?py_consol_sig.csv')
chry_overlap_names <- Reduce(intersect, chry_names_together)

pdf(paste0(getwd(), '/Chrysomya/chry_consol_upset.pdf'))
upset(fromList(chry_names_together), nsets=length(chry_names_together))
dev.off()


####Upset plot for Cochliomyia's significant consolidated miR names####

coch_names_together <- name_compiler('coch_.*?py_consol_sig.csv')
coch_overlap_names <- Reduce(intersect, coch_names_together)
#coch_mf_overlap <- Reduce(intersect, coch_names_together[2:3])

pdf(paste0(getwd(), '/Cochliomyia/coch_consol_upset.pdf'))
upset(fromList(coch_names_together[-2]), nsets=length(coch_names_together[-2]))
dev.off()


####Upset plot for Lucilia's significant consolidated miR names####

luci_names_together <- name_compiler('luci.*?py_consol_sig.csv')
luci_overlap_names <- Reduce(intersect, luci_names_together)
#luci_wt_overlap <- Reduce(intersect, luci_names_together[2:3])

pdf(paste0(getwd(), '/Lucilia/luci_consol_upset.pdf'))
upset(fromList(luci_names_together), nsets=length(luci_names_together))
dev.off()


####Upset plot for Plinthopyga's significant consolidated miR names####

plin_names_together <- name_compiler('plin.*?py_consol_sig.csv')
plin_overlap_names <- Reduce(intersect, plin_names_together)

pdf(paste0(getwd(), '/Plinthopyga/plin_consol_upset.pdf'))
upset(fromList(plin_names_together), nsets=length(plin_names_together))
dev.off()


####Upset plot for found overlapping miR names by species####
multi_names_together <- list("Chrysomya" = chry_overlap_names, 
                             "Cochliomyia" = coch_mf_overlap,
                             "Lucilia" = luci_wt_overlap,
                             "Plinthopyga" = plin_overlap_names)
chpl_overlap_names <- intersect(chry_overlap_names, plin_overlap_names)

pdf(paste0(getwd(), '/species_consol_upset.pdf'))
upset(fromList(multi_names_together), nsets=length(multi_names_together))
dev.off()


####Upset plot for overlapping wildtype miR names####
wt_names_together <- list("Chrysomya" = chry_names_together[[2]], 
                             "Cochliomyia" = coch_names_together[[1]],
                             "Lucilia" = luci_names_together[[3]],
                             "Plinthopyga" = plin_names_together[[2]])

upset(fromList(wt_names_together), nsets=length(wt_names_together))
wt_overlap_names1 <- intersect(chry_names_together[[2]], 
                              coch_names_together[[1]])
wt_overlap_names2 <- intersect(wt_overlap_names1, 
                               plin_names_together[[2]])


####Upset plot for overlapping wildtype miR names, using high+low quality data####
wt_any_qual_names_together <- list("Chrysomya" = chry_names_together[[1]], 
                          "Cochliomyia" = coch_names_together[[1]],
                          "Lucilia" = luci_names_together[[2]],
                          "Plinthopyga" = plin_names_together[[2]])

upset(fromList(wt_any_qual_names_together), nsets=length(wt_any_qual_names_together))
wt_any_qual_chry_coch <- intersect(chry_names_together[[1]], 
                               coch_names_together[[1]])
wt_any_qual_luci_plin <- intersect(luci_names_together[[2]], 
                               plin_names_together[[2]])
wt_any_qual_chry_coch_plin <- intersect(wt_any_qual_chry_coch, 
                                        plin_names_together[[2]])
#wt_any_qual_all_species <- intersect(wt_any_qual_chry_coch,
#                                        wt_any_qual_luci_plin)
wt_any_qual_coch_luci_plin <- intersect(coch_names_together[[1]],
                                         wt_any_qual_luci_plin)
wt_any_qual_luci_chry <- intersect(luci_names_together[[2]],
                                   chry_names_together[[1]])



