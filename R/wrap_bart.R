# Getting the BART wrapped function
#' @export
rbart <- function(x_train,
                 y,
                 x_test,
                 n_tree = 200,
                 n_mcmc = 2000,
                 n_burn = 500,
                 alpha = 0.95,
                 beta = 2,
                 df_splines = 10,
                 df = 3,
                 sigquant = 0.9,
                 kappa = 2,
                 scale_bool = TRUE) {

     # Verifying if x_train and x_test are matrices
     if(!is.matrix(x_train) || !is.matrix(x_test)){
          x_train <- as.matrix(x_train)
          x_test <- as.matrix(x_test)
     }

     if(is.null(colnames(x_train)) || is.null(colnames(x_test))) {
          stop("Insert valid NAMED matrix for x_train or x_test")
     }

     # Scaling x
     x_min <- apply(x_train,2,min)
     x_max <- apply(x_train,2,max)

     # Storing the original
     x_train_original <- x_train
     x_test_original <- x_test

     # Creating the scaled vesion
     x_train_scale <- x_train
     x_test_scale <- x_test

     # Normalising all the columns
     for(i in 1:ncol(x_train)){
             x_train_scale[,i] <- normalize_covariates_bart(y = x_train[,i],a = x_min[i], b = x_max[i])
             x_test_scale[,i] <- normalize_covariates_bart(y = x_test[,i],a = x_min[i], b = x_max[i])
     }


     # Scaling the y
     min_y <- min(y)
     max_y <- max(y)

     # Creating the B spline
     B_train <- as.matrix(splines::bs(x = x_train_scale,df = df_splines,intercept = TRUE))
     # B_test <- as.matrix(splines::bs(x = x_test_scale,df = df_splines,Boundary.knots = range(x_train_scale),intercept = TRUE))
     B_test <- as.matrix(predict(B_train,newx = x_test_scale))

     # Transforming back to matrix
     # x_train_scale <- as.matrix(B_train)
     # x_test_scale <- as.matrix(B_test)

     # Scaling "y"
     if(scale_bool){
        y_scale <- normalize_bart(y = y,a = min_y,b = max_y)
     } else {
        y_scale <- y
     }

     # Calculating \tau_{mu}
     tau_b <- tau_mu <- (4*n_tree*(kappa^2))
     # tau_b <- n_tree

     # Getting the naive sigma value
     nsigma <- naive_sigma(x = x_train_scale[,,drop = FALSE],y = y_scale)

     # Calculating tau hyperparam
     a_tau <- df/2

     # Calculating lambda
     qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
     lambda <- (nsigma*nsigma*qchi)/df
     d_tau <- (lambda*df)/2

     # Defining a_tau_b and d_tau_b
     a_tau_b <- n_tree
     d_tau_b <- 1

     # Call the bart function
     tau_init <- nsigma^(-2)
     mu_init <- mean(y_scale)

     # Creating the vector that stores all trees
     all_tree_post <- vector("list",length = round(n_mcmc-n_burn))

     # Generating the BART obj
     bart_obj <- sbart(x_train_scale,
          y_scale,
          x_test_scale,
          B_train = B_train,
          B_test = B_test,
          n_tree,
          n_mcmc,
          n_burn,
          tau_init,
          mu_init,
          tau_mu,
          tau_b,
          tau_b, # Same initial value as tau_b
          alpha,
          beta,
          a_tau,d_tau,
          a_tau_b,d_tau_b)


     if(scale_bool){
             # Tidying up the posterior elements
             y_train_post <- unnormalize_bart(z = bart_obj[[1]],a = min_y,b = max_y)
             y_test_post <- unnormalize_bart(z = bart_obj[[2]],a = min_y,b = max_y)
             for(i in 1:round(n_mcmc-n_burn)){
                     all_tree_post[[i]] <-  unnormalize_bart(z = bart_obj[[4]][,,i],a = min_y,b = max_y)
             }
             tau_post <- bart_obj[[3]]/((max_y-min_y)^2)
             tau_b_post <-  bart_obj[[5]]/((max_y-min_y)^2)
             tau_b_post_intercept <-  bart_obj[[6]]/((max_y-min_y)^2)

     } else {
             y_train_post <- bart_obj[[1]]
             y_test_post <- bart_obj[[2]]
             tau_post <- bart_obj[[3]]
             for(i in 1:round(n_mcmc-n_burn)){
                     all_tree_post[[i]] <-  bart_obj[[4]][,,i]
             }
             tau_b_post <-  bart_obj[[5]]
             tau_b_post_intercept <-  bart_obj[[6]]

     }

     # Return the list with all objects and parameters
     return(list(y_hat = y_train_post,
                 y_hat_test = y_test_post,
                 tau_post = tau_post,
                 tau_b_post = tau_b_post,
                 tau_b_post_intercept = tau_b_post_intercept,
                 all_tree_post = all_tree_post,
                 prior = list(n_tree = n_tree,
                              alpha = alpha,
                              beta = beta,
                              tau_mu = tau_mu,
                              a_tau = a_tau,
                              d_tau = d_tau),
                 mcmc = list(n_mcmc = n_mcmc,
                             n_burn = n_burn),
                 data = list(x_train = x_train,
                             y = y,
                             x_test = x_test)))
}


#
