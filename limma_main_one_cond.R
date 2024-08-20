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
myFilePath <- "TestData/data_quad.xlsx"

ratios <-
  as.data.frame(readxl::read_xlsx(myFilePath, col_names = TRUE))

data <- ratios

## Display data to faciliate choice of treatment and control
print(colnames(data))

## Extract gene/protein names
# proteins <- ratios$Proteins

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

# ## Extract only those proteins/peptides that has intensity values in
# ## K out of N replicates in each group
# rep_treats <- ncol(dat)
# k_out_N_treatment <- read_k_out_of_N(rep_treats, 'treatment')
# 
# idx_nz_treatment <-
#   which(rowSums(data[, 1:rep_treats] != 0) >= k_out_N_treatment)
# 
# idx_nz_both <- idx_nz_treatment
# if (length(idx_nz_both)) {
#   dat <- dat[idx_nz_both, ]
#   proteins <- proteins[idx_nz_both]
# }

FC_Cutoff <- readfloat("Enter the log fold change cut off=")
# temp <- log2(as.matrix(data[,2:ncol(data)]))
# temp[is.infinite(temp)] <- NA
# data[,2:ncol(data)] <- temp

## Impute data
# data_limma <- log2(as.matrix(data))
# data_limma[is.infinite(data_limma)] <- NA
# nan_idx <- which(is.na(data_limma))
# temp <- reshape(temp, nrow(data_limma)*ncol(data_limma), 1)
# hist(temp, na.rm = TRUE, xlab = "log2(intensity)", ylab = "Frequency",
#      main =  "All data: before imputation")
# fit <- fitdistr(c(na.exclude(data_limma)), "normal")
# mu <- as.double(fit$estimate[1])
# sigma <- as.double(fit$estimate[2])
# sigma_cutoff <- 6
# new_width_cutoff <- 0.3
# downshift <- 1.8
# width <- sigma_cutoff * sigma
# new_width <- width * new_width_cutoff
# new_sigma <- new_width / sigma_cutoff
# new_mean <- mu - downshift * sigma
# imputed_vals_my = rnorm(length(nan_idx), new_mean, new_sigma)
# scaling_factor <- readfloat_0_1("Enter a number > 0 and <=1 to scale imputed values = ")
# data_limma[nan_idx] <- imputed_vals_my*scaling_factor
# data_limma[nan_idx] <- imputed_vals_my

# ## Median Normalization Module
# want_normalization <- as.integer(readline(
#   cat(
#     'Enter 1: if you want to normalize the protein intensities in each experiemnt by substrating the median of the corresponding experiment\n',
#     '\bEnter 2: if you want to perform column wise median normalization of the data matrix\n',
#     '\b(for definition of "column wise median normalization" see README)= '
#   )
# ))
# data_limma <- data
# if (want_normalization == 1) {
#   boxplot(data_limma, main = "data before median substract normalization")
#   col_med <- matrixStats::colMedians(data_limma)
#   med_mat <- matlab::repmat(col_med, nrow(data_limma), 1)
#   data_limma <- data_limma - med_mat
#   boxplot(data_limma, main = "data after median substract normalization")
# } else if (want_normalization == 2) {
#   boxplot(data_limma, main = "data before column wise median normalization")
#   data_limma <- median_normalization(data_limma)
#   boxplot(data_limma, main = "data after column wise median normalization")
# }

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

## Limma main code
p <- readinteger_binary(
  cat(
    'Enter 1: if you want to use limma moderated pvalues after BH adjustment in the volcano plot\n',
    '\bEnter 0: if you want use limma moderated pvalues without adjustment in the volcano plot ='
  )
)
res.eb <- one_cond_limma_SILAC(tab = dat, condition = condition, sig.level = 0.05, FC_Cutoff, p)
# res.eb <- one_cond_limma_SILAC(tab = dat) ## calling Stephan's version

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
