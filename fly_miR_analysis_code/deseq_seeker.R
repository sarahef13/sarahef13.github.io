# This code will run DESeq on all files that have names containing 
# the pattern in files_to_analyze (line 73).
#
# Each input file will result in:
#   - 'file_resOrdered.csv', all DESeq results in order of adjusted p value
#   - 'file_sig.csv', only DESeq results with low adjusted p & high log fold change
#   - 'file_normalized.csv', normalized versions of all counts
#   - 'file_libSize_MA.pdf', a bar plot of sample library sizes and a dot plot of counts vs fold change
#
# If at least 2 significantly differential miRNAs are found,
# this code also generates:
#   - 'file_heatmap.pdf', a heatmap of significant miRNA counts per sample


library(DESeq2)
library(RColorBrewer)
library(gplots)


DESeeker <- function(csvFile) {
  
  #### Getting Data From File ####
  countData <- read.csv(csvFile, header=T, row.names=1)
  column_names <- gsub('\\.\\d+$', '', names(countData))
  
  
  #### Running DESeq On Data ####
  colData <- DataFrame(condition=factor(column_names))
  dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
  dds <- DESeq(dds)
  
  res <- results(dds)
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.10 & abs(resOrdered$log2FoldChange)>=1,]
  head(sig)
  selected <- rownames(sig)
  # TODO: make final big function return selected miR names
  
  
  #### Saving Results To CSVs ####
  resOrdered_name <- gsub('\\.csv', '_resOrdered.csv', csvFile)
  write.csv(resOrdered, resOrdered_name)
  
  sig_name <- gsub('\\.csv', '_sig.csv', csvFile)
  write.csv(sig, sig_name)
  
  normalized_name<- gsub('\\.csv', '_normalized.csv', csvFile)
  write.csv(counts(dds,normalized=TRUE), normalized_name)
  
  
  #### Plotting Data & DESeq Results ####
  pdf_name1 <- gsub('\\.csv', '_libSize_MA.pdf', csvFile)
  pdf(pdf_name1)
  barplot(colSums(countData)*1e-6, names=column_names, ylab="Library size (millions)")
  plotMA(dds)
  dev.off()
  
  if (nrow(sig) > 1) {
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
    pdf_name2 <- gsub('\\.csv', '_heatmap.pdf', csvFile)
    pdf(pdf_name2)
    heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]+1), 
               col = hmcol, scale="row", trace="none", margin=c(8,8), cexRow=.9, cexCol=.9, keysize=1.2 )
    dev.off()
  }
  
  
  return(paste('# significant miRNAs found in', csvFile, '=', nrow(sig), sep=' '))
}


files_to_analyze <- list.files(pattern = 'py_consol.csv', recursive = TRUE)
length(files_to_analyze)
lapply(files_to_analyze, DESeeker)

#file_to_analyze <- list.files(pattern = 'coch_cmacwtcombined_nolarvae_py_consol', recursive = TRUE)
#DESeeker(file_to_analyze)


