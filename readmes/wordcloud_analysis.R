getwd()
library("dplyr")
setwd("../sophiep/Desktop/Extended Research Project in Applied Bioinformatics /GO_analysis/")

################################1. Loading in the upreg vs non_upreg GO analysis for each files  
files <- list.files("GO_results_per_celltype/", full.names = TRUE)

GO_upreg_vs_non_upreg <- lapply(files, function(file) {
  cell_type <- gsub("^GO_terms_(.*)_genes.csv", "\\1", basename(file))
  df <- read.csv(file)
  df$New_cellType <- cell_type
  return(df)
}) %>% 
  bind_rows()


############################## 2. Loading in the macro_upreg vs upreg GO analysis for each files 
files_2 <- list.files("GO_results_per_celltype_combined/", full.names = TRUE)

GO_combined_upreg_vs_upreg <- lapply(files_2, function(file) {
  cell_type <- gsub("^GO_terms_(.*)_genes.csv", "\\1", basename(file))
  df <- read.csv(file)
  df$New_cellType <- cell_type
  return(df)
}) %>%
  bind_rows()

### Finding the main terms for each cell type 
#Using word cloud 

if(!require("wordcloud")){
  install.packages("wordcloud")
}
library("wordcloud")

if(!require("wordcloud2")){
  install.packages("wordcloud2")
}
library("wordcloud2")

install.packages(c("tm", "slam"))
library("tm")
library("slam")

if(!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}
library("RColorBrewer")

################################## Trying just with Apical Progenitors first
png("wordcloud_upreg_vs_non_upreg_AP.png", width = 900, height = 750)
apical_GO <- subset(GO_upreg_vs_non_upreg, New_cellType == "Apical_progenitors")
wordcloud(apical_GO$term_name) #Shows all the terms in apical progenitors for upreg vs non-upreg 
dev.off()



################################# Creating a For Loop to run through all cell types 
celltypes <- unique(GO_upreg_vs_non_upreg$New_cellType)
output_dir <- "wordcloud_upreg_vs_non_upreg_genes"
dir.create(output_dir, showWarnings = FALSE)

for (celltype in celltypes) {
  
  celltype_subset <- subset(GO_upreg_vs_non_upreg_grouped, New_cellType == celltype)
  
  output_file <- file.path(output_dir, paste0(celltype, ".png"))
  png(output_file, width = 900, height = 750)
  wordcloud <- wordcloud(words = celltype_subset$term_name, main = celltype, random.order = FALSE, scale = c(5, 0.9))
  
  title(main = celltype, cex.main = 3, line = -3)
  
  dev.off()
}



################################## Identifying which terms are shared in both 
apical_GO_combined <- subset(GO_combined_upreg_vs_upreg, New_cellType == "Apical_progenitors")


png("wordcloud_shared_big.png", width = 1800, height = 1500)
words <- as.character(apical_GO$term_name)
highlight_terms <- as.character(apical_GO_combined$term_name)
freqs <- rep(1, length(words))  # Assuming equal frequency

colors <- ifelse(words %in% highlight_terms, "deeppink3", "black")

wordcloud(
  words = words,
  freq = freqs,
  colors = colors,
  ordered.colors = TRUE,
  random.order = FALSE,
  scale = c(0.9, 0.8)
)

dev.off()


#Cross-checking 
intersect(apical_GO$term_name, apical_GO_combined$term_name)

