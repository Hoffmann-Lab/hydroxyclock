

innerFold <- function(X_train, X_test, Y_train, k, family, intercept_bool, alpha, lambda, threads){ 
  # X_train: input matrix of samples in train set with independent variables in columns (features: log-CPM in 2kb bins)
  # X_test: input matrix of samples in test set with independent variables in columns (features: log-CPM in 2kb bins)
  # Y_train: vector of samples in train set containing dependent variable (age in years)
  # Y_test: vector of samples in test set containing dependent variable (age in years)
  # k: integer stating the number of folds in inner cross-validation
  # family: string stating the family for the glmnet model (e.g., "poisson" or "gaussian")
  # intercept_bool: boolean stating whether the model should fit an intercept
  # alpha: numeric 0<= alpha <= 1, defines the mixing percentage between ridge and lasso regularization for elastic net (alpha=0 --> ridge, alpha=1 --> lasso)
  # lambda: "lambda.min" or "lambda.1se", where min uses lambda with optimal fit and 1se lambda for which error is within 1 standard error 
  # threads: number of cores that are used for parallel execution (doParallel)
  
  registerDoParallel(threads)
  
  message(paste("   Inner ", k, "-fold CV for model training and lambda optimization ...", sep = ""))  
  cvfit <- cv.glmnet(x = X_train, y = Y_train , family = family, type.measure = "mse",  intercept = intercept_bool, 
                     nfolds = k, alpha = alpha, trace.it = T, parallel = T, standardize = F) 
  
  message(paste("   Applying model to test data ... "))  
  prediction_test <- predict(cvfit, type = "response", newx = X_test, s = ifelse(lambda=="min", cvfit$lambda.min, ifelse(lambda=="1se", cvfit$lambda.1se, stop("lambda has to be a string: either 'min' or '1se'"))))
  
  return(prediction_test)
}


loocv <- function(X, Y, sample_info, k_inner, family, intercept_bool, alpha, lambda, threads){
  # X: input matrix with independent variables in columns (features: log-CPM in 2kb bins)
  # Y: vector containing dependent variable (age in years)
  # sample_info: matrix with information about samples in X (tissue, gender, dataset)
  # k_inner: integer stating the number of folds in inner cross-validation
  # family: string stating the family for the glmnet model (e.g., "poisson" or "gaussian")
  # intercept_bool: boolean stating whether the model should fit an intercept
  # alpha: numeric 0<= alpha <= 1, defines the mixing percentage between ridge and lasso regularization for elastic net (alpha=0 --> ridge, alpha=1 --> lasso)
  # lambda: "min" or "1se", where min uses lambda with optimal fit and 1se lambda for which error is within 1 standard error 
  # threads: number of cores that are used for parallel execution (doParallel)
  
  message("Performing Leave-One-Out Cross-Validation")  
  
  predictions <- sapply(1:nrow(X), function(test_sample) {
    message(paste(" Training a model on all samples except ", test_sample, "/", nrow(X), sep = ""))

    prediction_test <- innerFold(X_train = X[-test_sample,], X_test = X[test_sample,], Y_train = Y[-test_sample], 
                                 k = k_inner, family, intercept_bool, alpha, lambda, threads)
    return(prediction_test)
  })
  

  sample_df <- data.frame(tissue = sample_info$tissue,
                          gender =  sample_info$gender,
                          dataset = sample_info$dataset,
                          age = Y,
                          prediction = unlist(predictions))
  sample_df$err <- sample_df$prediction -  sample_df$age
  sample_df$abs_err <- abs(sample_df$err)
  
  return(sample_df)
}



LOOCV_scatterplot <- function(df_loocv){
  # df_loocv: output of loocv
  ggscatter(df_loocv, x = "age", y = "prediction", color = "tissue",  shape = "dataset" , alpha = .75, size = 3)  + 
    xlim(min(df_loocv[,c("age", "prediction")]), max(df_loocv[,c("age", "prediction")])) + ylim(min(df_loocv[,c("age", "prediction")]), max(df_loocv[,c("age", "prediction")])) + 
    annotate("text", x = min(df_loocv[,c("age", "prediction")]), y = max(df_loocv[,c("age", "prediction")]), 
             label = paste("cor = ", signif(cor$estimate, 2), " (p = ", signif(cor[["p.value"]],2), ")\nMAE = ", signif(median(absolute_errors),2), sep = ""), hjust = 0, vjust = 0.75, size = 5) + 
    geom_abline(intercept = 0, color = "black", linetype="dashed")  + xlab("chronological age \n[years]") + ylab("hydroxymethylation age \n[years]") +   
    scale_color_manual(values=as.vector(pals::alphabet(length(unique(df_loocv$tissue)))), name = "tissue  ") + 
    theme_bw(base_size = 14, base_line_size = 12/30) + 
    theme(legend.position="bottom", legend.spacing.x = unit(0.2, 'cm'), legend.spacing.y = unit(0.2, 'cm'),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(),
          axis.title = element_text(face="bold"), legend.text=element_text(color="grey20"),  legend.title=element_text(color="grey20"), 
          legend.box="vertical", legend.box.just = "left", legend.margin = margin(0,1.5,0,0, unit="cm")) + 
    guides(color = guide_legend(override.aes=list(shape = 15, size = 5), title.position="left", nrow = 5), shape = guide_legend(title.position="left"))
}







