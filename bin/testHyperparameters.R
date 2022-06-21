

testAlphas <- function(X, Y, sample_info, k_outer, k_inner, family, intercept_bool, threads){
  # X: input matrix with independent variables in columns (features: log-CPM in 2kb bins)
  # Y: vector containing dependent variable (age in years)
  # sample_info: matrix with information about samples in X (tissue, gender, dataset)
  # k_outer: integer stating the number of folds in outer cross-validation
  # k_inner: integer stating the number of folds in inner cross-validation
  # family: string stating the family for the glmnet model (e.g., "poisson" or "gaussian")
  # intercept_bool: boolean stating whether the model should fit an intercept
  # threads: number of cores that are used for parallel execution (doParallel)

  message("Performing Nested Cross-Validation with Grid Search over alpha")

  registerDoParallel(threads)


  ### NESTED CROSS VALIDATION
  ### OUTER AND INNER CROSS VALIDATION
  ### OUTER: to validate the performance in k-folds with the respective optimal lambda values for each fold
  ### INNER: using cv.glmnet cross-validation for k-fold cross-validation to find the best lambda

  ## OUTER
  # Create equally sized folds with balanced sample sized
  folds <- createFolds(Y, k = k_outer) # indices for each fold

  perf_alphas <- (matrix(ncol = 5, nrow = length(seq(0,1,0.1)*10)*10*2))
  id <- 1

  for(alpha in seq(0,1,0.1)){
    message(paste(" ", k_outer,"-CV for alpha=", alpha, sep = ""))
    #Perform k fold cross validation
    for(i in 1:k_outer){
      message(paste("  Start fold", i, "of ", k_outer))
      # split data in training and testing set according to fold
      fold <- folds[[i]]

      # Use the train data partitions for the specific fold to perform cv for lambda selection and model fitting ...
      ### INNER
      cvfit <- cv.glmnet(x = X[-fold,], y = Y[-fold], family = family,  intercept = intercept_bool,
                         nfolds = k_inner, alpha = alpha,  standardize = F, trace.it = T, parallel = T)
      # cvfit <- cv.glmnet(x = X[-fold,], y = Y[-fold], family = "poisson",  intercept = T,
                         # nfolds = 10, alpha = alpha,  standardize = F, trace.it = T, parallel = T)
      # Apply model to test dataset for the specific fold
      prediction_test_min <- predict(cvfit, type = "response", newx = testData[,2:ncol(testData)], s = cvfit$lambda.min)
      prediction_test_1se <- predict(cvfit, type = "response", newx = testData[,2:ncol(testData)], s = cvfit$lambda.1se)
      perf_alphas[id,] <-  c(alpha,"min", i,  median(abs(prediction_test_min - testData[,1])), cor.test(prediction_test_min,testData[,1])$estimate)
      perf_alphas[id+1,] <-  c(alpha,"1se", i,  median(abs(prediction_test_1se - testData[,1])), cor.test(prediction_test_1se,testData[,1])$estimate)
      id <- id + 1
    }
  }

  perf_alphas <- perf_alphas %>% data.frame() %>% setNames( c("alpha", "lambda", "fold" , "err" , "cor"))
  perf_alphas <- perf_alphas %>%  mutate_at(vars("alpha", "lambda", "fold"), ~as.factor(.))
  perf_alphas <- perf_alphas %>%  mutate_at(vars("err" , "cor"), ~as.numeric(.))

  return(perf_alphas)
}




testAlphas_boxplot <- function(df_testAlphas){
  # df_loocv: output of testAlphas
  plot1 <- ggplot(df_testAlphas, aes(x = alpha, y = err)) +
    geom_boxplot(alpha=0.5, aes(fill = lambda, color = lambda), show.legend = F) + theme_bw(base_size = 14, base_line_size = 12/30) +
    ylab("median absolute error \n[years]")  +  xlab('\u03b1 \n Ridge <-----> Lasso') +
    scale_fill_manual(values=c("aquamarine4", "forestgreen")) + scale_color_manual(values=c("aquamarine4", "forestgreen")) +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(),
          axis.title = element_text(face="bold"), axis.text.x = element_text(angle = 45, hjust=1))

  plot2 <- ggplot(df_testAlphas, aes(x = alpha, y = cor)) +
    geom_boxplot(alpha=0.5, aes(fill = lambda, color = lambda), show.legend = T) + theme_bw(base_size = 14, base_line_size = 12/30) +
    ylab("Pearson correlation")  + xlab('\u03b1 \n Ridge <-----> Lasso') +
    scale_fill_manual(values=c("aquamarine4", "forestgreen")) + scale_color_manual(values=c("aquamarine4", "forestgreen")) +
    theme(legend.text=element_text(color="grey20"),  legend.title=element_text(color="grey20"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(),
          axis.title = element_text(face="bold"), legend.spacing.y = unit(0.2, 'cm'), axis.text.x = element_text(angle = 45, hjust=1))

  grid.arrange(plot1, plot2, ncol=2, widths=c(5.7,7))
}
