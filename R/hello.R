rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/sbart.cpp")
source("R/other_functions.R")
source("R/wrap_bart.R")
n_ <- 100
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_*1))
colnames(x) <- "x"
colnames(x_new) <- "x"
# x <- as.data.frame(x)
# x_test <- as.data.frame(x_new)
y <- sin(2*x) + rnorm(n = n_,sd = 0.1)
y[x<0] <- y[x<0] + 2
y[x>0] <- y[x>0] - 2
# add_class <- rnorm(n = n_)
# add_class <- factor(ifelse(add_class>0,"A","B"))

# x <- x_new <- cbind(x,add_class)
# y[x$add_class=="A",1] <- y[x$add_class=="A",1] + 2
# y[x$add_class=="B",1] <- y[x$add_class=="B",1] - 2

# x <- cbind(1,x)
# colnames(x) <- c("x.0","x.1")
# x_new <- cbind(1,x_new)
# colnames(x_new) <- c("x.0","x.1")

# Testing over the motorbike data
# library(boot)
# data("motor")
# x <- motor$times %>% as.matrix
# y <- motor$accel %>% as.matrix()
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
# x_new <- x


# Iris test
# x <- dplyr::select(iris,Sepal.Length,Species)
# x_test <- x
# y <- iris$Sepal.Width

x <- as.data.frame(x)
x_test <- as.data.frame(x_new)


# Testing the GP-BART
bart_test <- rbart(x_train = x,y = unlist(c(y)),x_test = x_test,
                   n_tree = 10,n_mcmc = 2000,
                   alpha = 0.95,beta = 2,nIknots = 10,
                   df_tau_b = 10,prob_tau_b = 0.9,
                   n_burn = 500,scale_bool = TRUE)

# Convergence plots
par(mfrow = c(3,1))
plot(bart_test$tau_post[-2500],type = "l", main = expression(tau),ylab=  "")
plot(bart_test$tau_b_post[-2500],type = "l", main = expression(tau[b]),ylab=  "")
plot(bart_test$tau_b_post_intercept[-2500],type = "l", main = expression(tau[b[0]]),ylab=  "")

bartmod <- dbarts::bart(x.train = x,y.train = unlist(c(y)),ntree = 20,x.test = x_test)

library(MOTRbart)
motrbartmod <- MOTRbart::motr_bart(x = cbind(1,x),y = c(y))
pred_motrbart = predict_motr_bart(motrbartmod,cbind(1,x_test),type = "all")

# par(mfrow=c(2,1))
# plot(y$x,bart_test$y_hat %>% rowMeans())
# plot(y$x,bartmod$yhat.train.mean)

plot(x$x,bart_test$y_hat %>% rowMeans(), col = "blue")
plot(x$x,bartmod$yhat.train.mean, col = "red")

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
     geom_point(data = data.frame(x = x, y = y ), mapping = aes(x = x, y =y ))+
     geom_point(data = data.frame(x = x, y = apply(bart_test[[1]],1,mean)), mapping = aes(x = x,  y = y), col = "darkblue", alpha = 0.7, pch= 3)+
     geom_line(data = data.frame(x = x_new, y = apply(bart_test[[2]],1,mean)), mapping = aes(x = x, y = y) , col = "blue") +
     # geom_line(data = data.frame(x = x_new, y = bartmod$yhat.test.mean), mapping = aes(x = x, y = y), col ="red")+
     geom_line(data = data.frame(x = x_new, y = pred_motrbart %>% colMeans()), mapping = aes(x = x, y = y), col ="darkgreen")+

     geom_line(data = all_tree_posterior_mean,
               mapping = aes(x = x, y = value, col = name), alpha = 0.5,show.legend = FALSE)+
     ylim(y = range(y)*2)+
     theme_bw()

# plot(bart_test$tau_post, type = "l")
# microbenchmark::microbenchmark(bart(x_train = x,y_train = y,x_test = x,n_tree = 200,n_mcmc = 5000,
#                                     n_burn = 0,tau = 1,mu = 1,
#                                     tau_mu = 4*4*200,naive_sigma = 1,alpha = 0.95,
#                                     beta = 2,a_tau = 1,d_tau = 1,nsigma = 1),
#                                dbarts::bart(x.train = x,y.train = y,x.test = ,ntree = 200,ndpost = 5000,nskip = 0),times = 10)
# plot(bart_test$tau_b_post[-length(bart_test$tau_b_post)], type = "l")



# b_no_intercept <- bs(x = x)
# b_with_intercept <- bs(x = x,intercept = TRUE)
# Getting the values of B
# B <- B_train
# B_new <- B_test
# colnames(B) <- paste0("B.",1:ncol(B))
# plot_basis <- B %>% cbind(x) %>% as.data.frame() %>% pivot_longer(starts_with("B"))
#
# colnames(B_new) <- paste0("B.",1:ncol(B_new))
# plot_basis_new <- B_new %>% cbind(x_new) %>% as.data.frame() %>% pivot_longer(starts_with("B"))
#
# ggplot(plot_basis)+
#         geom_line(mapping = aes(x = x, y = value, col = name))+
#         theme_bw()

par(mfrow=c(1,1))
