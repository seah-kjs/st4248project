#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### HELPER FUNCTIONS ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MSE <- function(true_value, pred_value) {
  # true_value: observed y value
  # pred_value: yhat
  non_na_index <- !is.na(pred_value)
  mean((true_value[non_na_index] - pred_value[non_na_index])^2)
}

random_split_data <- function(dataset, train_ratio, seed) {
  # dataset: entire dataframe
  # train_ratio: value between 0 and 1
  # seed: for reproducibility
  
  if(!missing(seed)) set.seed(seed)
  
  n <- nrow(dataset)
  train_index <- sample(n, size = train_ratio * n)
  train_set <- dataset[train_index, ]
  test_set <- dataset[-train_index, ]
  list(train = train_set, test = test_set)
}


predict.regsubsets <- function(object, newdata, id,...) {
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- coef(object, id = id)
  xvars <- names(coefi)
  mat[, xvars]%*%coefi
}


# run, then save the images, then run generate_thumbnails
# histogram for single numerical vars
# scatterplots for numerical vars
# boxplot for numerical and factor 
# barchart for numerical and factor
# factor and factor summarise into table


lazysummary <- function(dataset, columns_to_look = "all", verbose = FALSE) {
  ## FOR cleaned dataset
  
  ## convert to data.frame
  if(!("data.frame" %in% class(dataset))) {
    dataset <- data.frame(dataset)
  }
  
  size <- dim(dataset)
  var_names <- names(dataset)
  
  ## peek at data
  if(size[2] > 10) {
    small_look <- dataset %>% 
      dplyr::select(seq(1,10,1)) %>%
      head()
  } else {
    small_look <- dataset %>% head()
  }
  
  print(small_look %>%
    kableExtra::kable() %>%
    kableExtra::kable_styling())
  
  
  #dataset %>% 
  #  head() %>%
  #  kableExtra::kable() %>%
  #  kableExtra::kable_styling()
  
  ## summary of data
  dplyr::glimpse(dataset)
  
  #class(dataset$variable)
  # integer, factor, logical, numeric, Date
  
  # group all variable classes tgt
  #
  
  ## graphical summaries
  ggplot(dataset, aes_string())
  
}











kfoldCV <- function(k, x, y, df){
  n <- length(y)
  shuffled_sample <- sample(n)
  MSE_p <- numeric(k)
  for(i in seq(k)) {
    top_index <- (n %/% k)*(i - 1) + 1
    if(i == k) btm_index <- n
    else btm_index <- (n %/% k)*i
    
    val_index <- shuffled_sample[top_index:btm_index]
    test_y <- y[val_index]
    train_y <- y[-val_index]
    train_x <- x[-val_index]
    test_x <- x[val_index]
    lm_fit <- lm(train_y ~ bs(train_x, df = df))
    MSE_p[i] <- mean((test_y - predict(lm_fit, newdata = list(train_x = test_x)))^2)
  }
  mean(MSE_p)
}

gaussian_kernel <- function(x, y, sigma = 1) {
  exp(-(sum((x - y)^2)) / sigma^2)
}


gram_matrix <- function(X, n, kernel_fn, kernel_args = NULL) {
  Kx <- matrix(0, nrow = n, ncol = n)
  for(i in seq(n)) {
    for(j in seq(n)) {
      if(is.null(dim(X))) {   ## X is 1 dimension
        if(is.null(kernel_args)) {
          Kx[i, j] <- do.call(kernel_fn, args = list(X[i], X[j]))
        } else {
          Kx[i, j] <- do.call(kernel_fn, args = list(X[i], X[j], kernel_args))
        }
      } else {    ## X is more than 1 dimension
        if(is.null(kernel_args)) {
          Kx[i, j] <- do.call(kernel_fn, args = list(X[i, ], X[j, ]))
        } else {
          Kx[i, j] <- do.call(kernel_fn, args = list(X[i, ], X[j, ], kernel_args))
        }
      }
    }
  }
  Kx
}

rkhs_estimator <- function(xgrid, xdata, ydata, kernel_fn, lambda, kernel_args) {
  n <- length(ydata)
  diag_matrix <- diag(1, nrow = n, ncol = n)
  
  Kx <- gram_matrix(xdata, n, kernel_fn = kernel_fn, kernel_args = kernel_args)
  front_part <- t(ydata) %*% solve(Kx + lambda * diag_matrix, tol = 1e-60)
  
  if(is.null(dim(xgrid)) | is.null(dim(xdata))){
    f_hat <- sapply(
      xgrid, function(xx) front_part %*% 
        (Map(function(yy) kernel_fn(xx, yy, kernel_args), xdata) %>% unlist(.)))
    f_hat
  } else{
    f_hat <- apply(
      xgrid, 1, function(xx) front_part %*% 
        (apply(xdata, 1, function(yy) kernel_fn(xx, yy, kernel_args)) ))
    f_hat
  }
}






rkhs_gCV <- function(xdata, ydata, kernel_fn, kernel_args, lambda) {
  n <- length(ydata)
  Kx <- gram_matrix(xdata, n, kernel_fn, kernel_args)
  kmatrix <- solve(Kx + diag(n) * lambda, tol = 1e-60)
  hat_matrix <- kmatrix %*% Kx
  
  if(nrow(hat_matrix) != ncol(hat_matrix)){
    stop("Error in calculating hat matrix. It is not a square matrix.")
  }
  
  yhat <- hat_matrix %*% ydata
  
  ((1/n) * sum((ydata - yhat)^2)) / ((1/n) * sum(diag(diag(1, n) - hat_matrix)))^2
}






