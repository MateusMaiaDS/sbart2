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
                  scale_bool = TRUE,
                  # Hyperparam for tau_b and tau_b_0
                  df_tau_b = 5,
                  prob_tau_b = 0.9) {

     # Verifying if x_train and x_test are matrices
     if(!is.data.frame(x_train) || !is.data.frame(x_test)){
          stop("Insert valid data.frame for both data and xnew.")
     }

     # Creating the col
     original_p <- 1:ncol(x_train)-1


     # Getting the valid
     dummy_x <- caret::dummyVars(~.,data = x_train)

     col_names <- attr(dummy_x$terms,"term.labels")
     dummy_x_train_m <- predict(object = dummy_x,newdata = x_train)
     dummy_x_test_m <- predict(object = dummy_x,newdata = x_test)

     recode_names <- recode_vars(x_train = x_train,dummy_obj = dummy_x)
     n_levels <- c(table(recode_names)-1) # It also counts the zero

     x_train_scale <- dummy_x_train_m
     x_test_scale <- dummy_x_test_m

     # Scaling x
     x_min <- apply(as.matrix(x_train_scale),2,min)
     x_max <- apply(as.matrix(x_train_scale),2,max)

     # Storing the original
     x_train_original <- x_train
     x_test_original <- x_test


     # Normalising all the columns
     for(i in 1:ncol(x_train)){
             x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
             x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
     }


     # Scaling the y
     min_y <- min(y)
     max_y <- max(y)

     # Creating the B spline
     B_train <- as.matrix(splines::bs(x = x_train_scale[,col_names[!(col_names %in% dummy_x$facVars)], drop = FALSE],df = df_splines,intercept = TRUE))
     B_test <- as.matrix(predict(B_train,newx = x_test_scale[,col_names[!(col_names %in% dummy_x$facVars)], drop = FALSE]))

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
     tau_b_0 <- tau_b <- tau_mu <- 0.1*(4*n_tree*(kappa^2))
     # tau_b <- n_tree

     # Getting the naive sigma value
     nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

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

     rate_d_tau <- optim(par = 1,d_tau_b_rate,method = "L-BFGS-B",lower = 0,
                         df = df_tau_b, prob = prob_tau_b,kappa = kappa,
                         n_tree = n_tree)

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
          tau_b_0, # Same initial value as tau_b
          alpha,
          beta,
          a_tau,d_tau,
          a_tau_b = df_tau_b*0.5,d_tau_b = rate_d_tau$par, # Hypeparameters from tau_b and tau_b_0
          original_p, # Getting the p available variables
          n_levels # Getting the sample levels
          )


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
