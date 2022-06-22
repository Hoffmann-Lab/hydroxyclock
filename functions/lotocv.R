


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

lotocv <- function(X, Y, sample_info, k_inner, family, intercept_bool, alpha, lambda, threads){
  # X: input matrix with independent variables in columns (features: log-CPM in 2kb bins)
  # Y: vector containing dependent variable (age in years)
  # sample_info: matrix with information about samples in X (tissue, gender, dataset)
  # k_inner: integer stating the number of folds in inner cross-validation
  # family: string stating the family for the glmnet model (e.g., "poisson" or "gaussian")
  # intercept_bool: boolean stating whether the model should fit an intercept
  # alpha: numeric 0<= alpha <= 1, defines the mixing percentage between ridge and lasso regularization for elastic net (alpha=0 --> ridge, alpha=1 --> lasso)
  # lambda: "min" or "1se", where min uses lambda with optimal fit and 1se lambda for which error is within 1 standard error
  # threads: number of cores that are used for parallel execution (doParallel)

  message("Performing Leave-One-Tissue-Out Cross-Validation")


  predictions <- sapply(unique(sample_info$tissue), function(test_tissue) {
    message(paste(" Training a model on all tissues except:  ", test_tissue))
    test_tissue_ids <- sample_info$tissue %in% test_tissue
    prediction_test <- innerFold(X_train = X[!test_tissue_ids,], X_test = X[test_tissue_ids,], Y_train = Y[!test_tissue_ids],
                                 k = k_inner, family, intercept_bool, alpha, lambda, threads)
    return(prediction_test)
    })


  ids  <- unlist(sapply(unique(sample_info$tissue), function(x) which(sample_info$tissue %in% x)))
  sample_df <- data.frame(tissue = sample_info[ids,]$tissue,
                          gender =  sample_info[ids,]$gender,
                          dataset = sample_info[ids,]$dataset,
                          age = Y[ids],
                          prediction = unlist(predictions))
  sample_df$err <- sample_df$prediction -  sample_df$age
  sample_df$abs_err <- abs(sample_df$err)

  return(sample_df)
}


stat_box_data <- function(y, upper_limit = max(sample_df$err) * 1.2){
  return(data.frame(y =  1 * upper_limit, label = length(y)))
}


LOTOCV_boxplot <- function(df_lotocv){
  # df_lotocv: output of lotocv
  ggplot(df_lotocv, aes(x = tissue, y = err)) + geom_hline(yintercept = 0, alpha = 0.5) +
    geom_boxplot(alpha=0.5, aes(fill = tissue, color = tissue), show.legend = FALSE, outlier.shape = NA) +  theme_bw(base_size = 14.4, base_line_size = 12/30) +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), axis.title = element_text(face="bold"),
          legend.position = "bottom",  legend.text=element_text(color="grey20", size=12),  legend.title=element_text(color="grey20")) +
    ylab("error\n[years]") + xlab(" ") +
    scale_fill_manual(values=as.vector(pals::alphabet(26))) + scale_color_manual(values=as.vector(pals::alphabet(26)))  +
    stat_summary(fun.data = "mean_cl_normal",aes(shape = "age acceleration\n(mean error)"), colour = "black", size = 3, geom = "point") +
    scale_shape_manual("", values = c("age acceleration\n(mean error)" = 18)) +
    stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.2, size = 4, fontface = "bold") +
    geom_jitter(size = 2, alpha = 0.35, shape = 16, aes(color = tissue), show.legend = FALSE, width = 0.3)
}
