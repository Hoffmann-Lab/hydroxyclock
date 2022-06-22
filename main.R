library(ggplot2)
library(glmnet)
library(doParallel)
library(caret)
library(data.table)
library(dplyr) 
library(gridExtra)
library(pals)
library(rtracklayer)


source("bin/testHyperparameters.R")
source("bin/loocv.R")
source("bin/lotocv.R")


# Load Data


import_matrix <- fread("hg38_2000bp_bins_woBl_wPeakOl_Cui_He_featureMatrix_normalized.txt", data.table = F)

sample_info_df <- import_matrix[,c("tissue", "gender", "dataset")]
Y <-  import_matrix[,"age"]
X <- as.matrix(import_matrix[,!colnames(import_matrix) %in% c("tissue", "gender", "dataset", "age")])
# remove features with bins on X- and Y-chromosome 
X <- X[,!grepl("chrY|chrX",colnames(X))]

# Hyperparameter Optimization

df_testAlphas <- testAlphas(X = X, Y = Y, sample_info = sample_info_df, k_outer = 10, k_inner = 10, family = "poisson",
                            intercept_bool = T, threads = 10)

testAlphas_boxplot(df_testAlphas)


# LOOCV

sample_df_loocv <- loocv(X = X, Y = Y, sample_info = sample_info_df, k_inner = 10, family = "poisson",
                         intercept_bool = T, alpha = 0.5, lambda = "min", threads = 10)
sample_df_loocv <- readRDS("temp.rds")

LOOCV_scatterplot(sample_df_loocv)

# LOTOCV

sample_df_lotocv <- lotocv(X = temp, Y = Y, sample_info = sample_info_df, k_inner = 10, family = "poisson",
                           intercept_bool = T, alpha = 0.5, lambda = "min", threads = 10)

LOTOCV_boxplot(sample_df_lotocv)
