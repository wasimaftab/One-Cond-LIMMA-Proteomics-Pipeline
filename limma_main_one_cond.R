#This is a pipeline to analyze proteiomic data in a proteinGroups.txt (MaxQuant output) file with one condition
#Author:Wasim Aftab

cat('\014')
rm(list = ls())

## Installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
list.of.packages <- c("limma", "qvalue")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  BiocManager::install(new.packages)

## Installing CRAN packages
list.of.packages <-
  c(
    "dplyr",
    "stringr",
    "MASS",
    "matlab",
    "plotly",
    "htmlwidgets",
    "rstudioapi",
    "webshot",
    "matrixStats"
  )
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)
if (is.null(webshot:::find_phantom())) {
  webshot::install_phantomjs()
}

library(dplyr)
library(stringr)
library(MASS)
library(matlab)
library(plotly)
library(limma)
library(qvalue)
library(htmlwidgets)
library(rstudioapi)
library(sva)

## Chdir to source dir
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))
cur_dir <- getwd()

source("limma_helper_functions.R")
## Load the proteingroups file
# myFilePath <- file.choose()
# temp <- unlist(strsplit(myFilePath, "\\", fixed = TRUE))
# myFilePath <- "/home/wasim/Downloads/One-cond-limma/Example/proteinGroups.txt"
# myFilePath <- "prodat_SILAC.xlsx"
myFilePath <- "Parabolic_volcano/1110_data_before.xlsx"

data <-
  as.data.frame(readxl::read_xlsx(myFilePath, col_names = TRUE))

colnames_data <- colnames(data)

## Temporary change in column names
colnames(data)[2:ncol(data)] <- c("Log2R_1", "Log2R_2", "Log2R_3", "Log2R_4")

## Display data to faciliate choice of treatment and control
print(colnames(data))
# browser()

## Extract gene/protein names
# proteins <- data$Proteins

#Extract required data for Limma
condition <-
  readline('Enter condition name (case insensitive) as it appeared in the column names printed above= ')
# browser()
## check condition
if(is.null(condition)) stop("Need a condition.")
if(length(condition) != 1) stop("Need exactly one condition.")

## select condition
sel <- grep(condition, colnames(data), ignore.case = TRUE)
if(length(sel) < 2) stop("Need at least two replicates for one-sample differential expression")

## removing blank rows
print(paste("data rows before removing blank rows = ", NROW(data)))
temp <-
  as.matrix(rowSums(apply(data[,2:ncol(data)], 2, as.numeric), na.rm = TRUE))
idx <- which(is.na(temp))
if (length(idx)) {
  data <- data[-idx,]
  # proteins <- proteins[-idx]
}
print(paste("data rows after removing blank rows = ", NROW(data)))

print(head(data))

FC_Cutoff <- readfloat("Enter the log fold change cut off=")


## Extract only those proteins/peptides that has intensity values in
## K out of N replicates in each group
dat <- data[, sel]
rep_treats <- ncol(dat)
k_out_N_treatment <- read_k_out_of_N(rep_treats, 'treatment')
idx_nz_treatment <-
  which(rowSums(!is.na(dat)) >= k_out_N_treatment)

idx_nz_both <- idx_nz_treatment
if (length(idx_nz_both)) {
  dat <- data[idx_nz_both, c(1,sel)]
  # proteins <- proteins[idx_nz_both]
}
boxplot(dat[,2:ncol(dat)], names = c("rep1", "rep2", "rep3", "rep4"), main = "data before batch correction and median normalization")

## Batch correction and Normalization
# For comBat to work, you need to impute missing values, so the steps are
# 1. Record the indexes of missing values
data_matrix <- as.matrix(dat[,2:ncol(dat)])
data_matrix[is.infinite(data_matrix)] <- NA
nan_idx <- which(is.na(data_matrix))

# 2. Impute missing values using a small normal distribution at the left tail of the original data distribution
# Note, you can replace this algorithm with the one that is most applicable to you context
data_matrix <- impute_from_normal_dist(data_matrix, nan_idx)

# 3. Perform batch correction using comBat on the full data (no missing values)
# define the batch vector to indicate which batch each sample belongs to
batch <- c("Batch1", "Batch1", "Batch2", "Batch2")
data_matrix <- ComBat(dat = data_matrix, batch = batch, par.prior = TRUE, prior.plots = FALSE)

# 4. Replace the values at indexes where earlier missing values were with NA
data_matrix[nan_idx] <- NA
boxplot(data_matrix, names = c("rep1", "rep2", "rep3", "rep4"), main = "data after batch correction but before median normalization")
norm_data <- median_normalization(data_matrix)
dat[,2:ncol(dat)] <- norm_data
boxplot(dat[,2:ncol(dat)], names = c("rep1", "rep2", "rep3", "rep4"), main = "data after batch correction and median normalization")
# browser()

## Limma main code
p <- readinteger_binary(
  cat(
    'Enter 1: if you want to use limma moderated pvalues after BH adjustment in the volcano plot\n',
    '\bEnter 0: if you want use limma moderated pvalues without adjustment in the volcano plot ='
  )
)
res.eb <- one_cond_limma_SILAC(tab = dat, condition = condition, sig.level = 0.05, FC_Cutoff, p)

## Save the data file
final_data <- res.eb
filename_final_data <-
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), '_final_data')

##Create plotly object and save plot as html
filename_mod <-
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), '_limma_plot')

if (p) {
  title <- paste0(toupper(condition), ": volcano plot of moderated p-values after BH adjustment")
}else {
  title <- paste0(toupper(condition), ": volcano plot of moderated p-values")
}

fullpath_subDir <- display_plotly_fig_one_cond(final_data, FC_Cutoff, filename_mod, title = title)

write.table(
  final_data,
  paste0(fullpath_subDir, "/", filename_final_data, ".tsv"),
  sep = '\t',
  row.names = FALSE,
  col.names = TRUE
)
