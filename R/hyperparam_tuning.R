# Setting the parameters and other comparisons
rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/sbart.cpp")
source("R/wrap_bart.R")
source("R/other_functions.R")
source("R/cv.R")

n_ <- 100
df_splines_ <- 8
n_tree
set.seed(42)
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_*100))
colnames(x) <- "x"
colnames(x_new) <- "x"
y <- sin(x) + rnorm(n = n_,sd = 0.1)
y[x<0] <- y[x<0] + 2
y[x>0] <- y[x>0] - 2

# Doing the same for the motor data
library(boot)
data("motor")
x <- motor$times %>% as.matrix
y <- motor$accel %>% as.matrix()
colnames(x) <- "x"

x <- as.data.frame(x)
data <- cbind(x,y)
colnames(data) <- c("x","y")

kfold_obj <- k_fold(data = data,dependent_variable = "y",k_partitions = 10,seed = 42,as_data_frame = TRUE)

# Metrics data.frame
metrics_df <- data.frame(n_rep = NA, rmse = NA, crps = NA, pi = NA, model = NA)
metrics_df <- metrics_df[-1,]

for(i in 1:10){

     # Doing for sbart
     sbartmod <- rbart(x_train = kfold_obj[[i]]$x_train,
                       y = unlist(c(kfold_obj[[i]]$y_train)),
                       x_test = kfold_obj[[i]]$x_test,
                       n_tree = 10,n_mcmc = 2500,n_burn = 500,
                       alpha = 0.95,beta = 2,df_splines = df_splines_)

     sd_post_vec <- c((sbartmod$tau_post^(-1/2))[-2000],
                      mean((sbartmod$tau_post[-2000]^(-1/2))))

     rmse_sbart <- rmse(x = unlist(c(kfold_obj[[i]]$y_test)),
                        y = rowMeans(sbartmod$y_hat_test))

     crps_sbart <- crps(y = unlist(c(kfold_obj[[i]]$y_test)),
                        means = rowMeans(sbartmod$y_hat_test),
                        sds = rep(mean(sd_post_vec),length(rowMeans(sbartmod$y_hat_test))))$CRPS

     pi_c_sbart <- pi_coverage(y = unlist(c(kfold_obj[[i]]$y_test)),
                               y_hat_post = t(sbartmod$y_hat_test),
                               sd_post = sd_post_vec,prob = 0.5)

     metrics_df <- rbind(metrics_df, data.frame(n_rep = i,
                                                rmse = rmse_sbart,
                                                crps = crps_sbart,
                                                pi = pi_c_sbart,
                                                model = "SBART"))

     # Doing the same for BART
     bartmod <- dbarts::bart(x.train = kfold_obj[[i]]$x_train,
                             y.train = unlist(c(kfold_obj[[i]]$y_train)),
                             x.test = kfold_obj[[i]]$x_test)

     rmse_bart <- rmse(x = unlist(c(kfold_obj[[i]]$y_test)),
                        y = bartmod$yhat.test.mean)

     crps_bart <- crps(y = unlist(c(kfold_obj[[i]]$y_test)),
                        means = bartmod$yhat.test.mean,
                        sds = rep(mean(bartmod$sigma),length(bartmod$yhat.test.mean)))$CRPS

     pi_c_bart <- pi_coverage(y = unlist(c(kfold_obj[[i]]$y_test)),
                               y_hat_post = bartmod$yhat.test,
                               sd_post = bartmod$sigma,prob = 0.5)

     metrics_df <- rbind(metrics_df, data.frame(n_rep = i,
                                                rmse = rmse_bart,
                                                crps = crps_bart,
                                                pi = pi_c_bart,
                                                model = "BART"))

}

plot_rmse <- metrics_df %>%
     ggplot(mapping = aes(x = model, y = rmse))+
     geom_boxplot()+
     theme_bw()

plot_crps <- metrics_df %>%
     ggplot(mapping = aes(x = model, y = crps))+
     geom_boxplot()+
     theme_bw()

plot_pi <- metrics_df %>%
     ggplot()+
     geom_hline(yintercept = 0.5, col = "red", lty = "dashed")+
     geom_boxplot(mapping = aes(x = model, y = pi))+
     theme_bw()

cowplot::plot_grid(plot_rmse,plot_crps,plot_pi,nrow = 1)
