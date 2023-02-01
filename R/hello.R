rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/sbart.cpp")
source("R/wrap_bart.R")
source("R/other_functions.R")
n_ <- 100
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_*10))
y <- sin(x) + rnorm(n = n_,sd = 0.1)
y[x<0] <- y[x<0] + 2
y[x>0] <- y[x>0] - 2

colnames(x) <- "x"
colnames(x_new) <- "x"

# Testing over the motorbike data
library(boot)
data("motor")
x <- motor$times %>% as.matrix
y <- motor$accel %>% as.matrix()
x_new <- seq(min(x),max(x),length.out = 500) %>% as.matrix()
colnames(x) <- "x"
colnames(x_new) <- "x"

# Testing the GP-BART
bart_test <- rbart(x_train = x,y = y,x_test = x_new,n_tree = 10,n_mcmc = 2000,
                   alpha = 0.95,beta = 2,
                   n_burn = 500,scale_bool = TRUE)



bartmod <- dbarts::bart(x.train = x,y.train = y,ntree = 20,x.test = x_new)

# All ll trees prediction plot
all_tree_posterior_mean <- Reduce("+",bart_test$all_tree_post)/length(bart_test$all_tree_post)
if(!(ncol(all_tree_posterior_mean) %>% is_null())){
        colnames(all_tree_posterior_mean) <- paste0("tree.",1:ncol(all_tree_posterior_mean))
} else {
        all_tree_posterior_mean <- all_tree_posterior_mean %>% as.matrix()
        colnames(all_tree_posterior_mean) <- "tree.1"
}
all_tree_posterior_mean <- all_tree_posterior_mean %>% as.data.frame %>% add_column(x) %>% pivot_longer(starts_with("tree"))


# Getting quantiles for ribon
pi_ <- bart_test$y_hat_test %>% apply(1,function(x){quantile(x,probs = c(0.025,0.975))}) %>% t %>% cbind(x_new)
colnames(pi_) <- c("lower","upper","x")
pi_ <- pi_ %>% as.data.frame()
# Replicating the same plot on ggplot
ggplot()+
     geom_ribbon(data = pi_,mapping = aes(x = x,ymin = lower, ymax = upper), col = NA, alpha = 0.1, lty = "dashed",fill = "blue")+
     geom_point(data = data.frame(x = x, y= y ), mapping = aes(x = x, y =y ))+
     geom_point(data = data.frame(x = x, y = apply(bart_test[[1]],1,mean)), mapping = aes(x = x,  y = y), col = "darkblue", alpha = 0.7, pch= 3)+
     geom_line(data = data.frame(x = x_new, y = apply(bart_test[[2]],1,mean)), mapping = aes(x = x, y = y) , col = "blue") +
     geom_line(data = data.frame(x = x_new, y = bartmod$yhat.test.mean), mapping = aes(x = x, y = y), col ="red")+
     geom_line(data = all_tree_posterior_mean,
               mapping = aes(x = x, y = value, col = name), alpha = 0.5,show.legend = FALSE)+
     ylim(y = range(y)*1.4)+
     theme_bw()

# plot(bart_test$tau_post, type = "l")
# microbenchmark::microbenchmark(bart(x_train = x,y_train = y,x_test = x,n_tree = 200,n_mcmc = 5000,
#                                     n_burn = 0,tau = 1,mu = 1,
#                                     tau_mu = 4*4*200,naive_sigma = 1,alpha = 0.95,
#                                     beta = 2,a_tau = 1,d_tau = 1,nsigma = 1),
#                                dbarts::bart(x.train = x,y.train = y,x.test = ,ntree = 200,ndpost = 5000,nskip = 0),times = 10)
plot(bart_test$tau_b_post[-length(bart_test$tau_b_post)], type = "l")


# Getting new B-splines function
n_basis <- 10
x<- normalize_covariates_bart(y = x,a = min(x),b = max(x))
B <- bspline(x = x,x_obs = seq(from = sort(x)[1], to = sort(x)[length(x)],length.out = n_basis))
B <- B[,-11]
B_new <- bspline(x = x_new,x_obs = seq(from = sort(x)[1], to = sort(x)[length(x)-1],length.out = n_basis))
# B_new <- bspline(x = x_new,x_obs = x)


# BS from package
B <- splines::bs(x,knots = seq(from = sort(x)[2], to = sort(x)[length(x)-1],length.out = n_basis))
B <- splines::bs(x,knots = seq(from = sort(x)[2], to = sort(x)[length(x)-1],length.out = n_basis))
B <- splines
# B <- mgcv::s(x,k = 10)

B_new <- splines::bs(x_new,knots = seq(from = sort(x)[1], to = sort(x)[length(x)-1],length.out = n_basis))
B_new <- splines::bs(x_new,df = 10)
# Getting the values of B
colnames(B) <- paste0("B.",1:ncol(B))
plot_basis <- B %>% cbind(x) %>% as.data.frame() %>% pivot_longer(starts_with("B"))

colnames(B_new) <- paste0("B.",1:ncol(B_new))
plot_basis_new <- B_new %>% cbind(x_new) %>% as.data.frame() %>% pivot_longer(starts_with("B"))

ggplot(plot_basis)+
        geom_line(mapping = aes(x = x, y = value, col = name))+
        theme_bw()


# B <- B[,-n_basis]
# B_new <- B_new[,-n_basis]
# Getting y_hat
H <- solve(crossprod(B) )%*%crossprod(B,y)
y_train_hat <- B%*%H
y_new_hat <- B_new%*%H
par(mfrow = c(1,2))
plot(x,y, ylim = range(y))
plot(x_new,y_new_hat,col = "red", ylim = range(y))


