
# UCT ASSIGNMENT
# Author: Julian Albert
# Date:   01/03/2019

# Clean Environment
rm(list = ls()); dev.off()
# Assignment prelim script, useful packages and custom ggplot theme
source("../../../../PrelimScript/assignment_prelim_script.R")
## 0. Load all the packages for the assignment question
if (!require("pacman")) install.packages("pacman")
p_load(splines, RColorBrewer, psych)

# Question 1 ----
# P-Splines

## 1.1 Simulate the data ----
set.seed(123) # for reproducibility

n <- 100 # number of data points
x_dat <- runif(n, 0, 1) # x from U(0, 1)
x_dat_sorted <- sort(x_dat)

func.mu <- function(x)
{
  mu <-  5 + sin(3*pi*(x - 0.6))
  return(mu)
}

mu_dat <- func.mu(x_dat_sorted) # mean for y
sigma_dat <- 0.5 # sd for y
y_dat <- rnorm(n, mean = mu_dat, sd = sigma_dat) # y ~ N(mu, 0.5^2)

# plot(y_dat~x_dat_sorted) # quick check of data

## 1.2 Start with a B-Spline Basis ----

## 1.2.1 FROM eilers2010splines

tpower <- function(x, t, p){
  (x-t)^p * (x > t)
}

newBase <- function(x, xl, xr, ndx, deg){
  dx <- (xr - xl)/ndx
  knots <- seq(xl - deg*dx, xr + deg*dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1)/(gamma(deg+1)*dx^deg)
  B <- (-1)^(deg + 1)* P%*%t(D)
  return(B)
}

x_dat_plot <- seq(0, 1, length.out = 50000) # increase n for smooth basis plot

newBS_plot <- newBase(x_dat_plot, x_dat_plot[1], 
                      x_dat_plot[length(x_dat_plot)], 21, 3)

color <- colorRampPalette(c('red', 'orange', 'springgreen', 'royalblue'))(NCOL(newBS_plot))

## 1.2.2 B-Spline from bs(splines)

n.knots <- 20
# we need 20 evenly spaced interior knots
## our boundary knots are on 0 and 1 find even step from 0 to 1
even_space_size <- diff(seq(0, 1, length.out = 22))[1]
## our interior knots are spaced equally from 0 +step:1-step
range.knots <- c(even_space_size, 1-even_space_size)
## construct basis
B_basis <- bs(x_dat_plot,
              knots = seq(range.knots[1],
                          range.knots[2],
                          length.out = 20),
              intercept = TRUE,
              Boundary.knots = c(0, 1))

## Compare th Figures for each method

#setwd('../Figs')
#pdf('newBS.#pdf') # eilers2010splines
plot(newBS_plot[,1]~x_dat_plot , typ = 'l', 
     ylim = c(0, 1), xlab = 'x', ylab = '', 
     main = 'Cubic B-Spline Basis Functions')
for(i in 1:NCOL(newBS_plot))
  lines(newBS_plot[,i]~x_dat_plot, col = color[i], lwd = 2)
dev.off()

#pdf('OldBS.#pdf') # bs()
plot(B_basis[,1]~x_dat_plot , typ = 'l',
     ylim = c(0, 1), xlab = 'x', ylab = '', 
     main = 'Cubic B-Spline Basis Functions')
for(i in 1:NCOL(B_basis))
  lines(B_basis[,i]~x_dat_plot, col = color[i], lwd = 2)
dev.off()
#setwd('../Code')

## 1.2.3 FOR ANALYSIS WE WILL USE 100 POINTS DATA

newBS <- newBase(x_dat_sorted, x_dat_sorted[1], 
                 x_dat_sorted[length(x_dat_sorted)], 21, 3)

### For a B-Spline Basis -> no penalty returns same so see later
beta_hat <- solve(crossprod(newBS)) %*% t(newBS) %*% y_dat # betas
f_hat <- newBS %*% beta_hat # predictions

## 1.3 P-spline is B-spline + difference penalty ----

## need to specify a difference operation matrix and lambda
## the P_matrix is constructed P = sum(b_{i+1} - b_i)^2
k <- NROW(beta_hat)
P_matrix <- diff(diag(k), differences = 1)
S_matrix <- t(P_matrix) %*% P_matrix # penalty matrix

## 1.3.1 Calc. GCV for range of lambdas

lambda <- seq(0, 100, length.out = 5000) # range of lambda values
## need a smoothing matrix from S_mat in notes
Omega_mat <- S_matrix # for large n more efficient to use P_spline as penalty matrix

## initialise GCV and define N
GCV <- numeric()
N <- length(y_dat)

## 1.3.2 Get GCV values for different Lambda's and find the minimum
for(i in 1:length(lambda))
{
  tmp.lambda <- lambda[i]
  tmp.S_mat <- newBS %*% solve(crossprod(newBS) + tmp.lambda*Omega_mat) %*% t(newBS)
  tmp.f_hat <- tmp.S_mat %*% y_dat
  GCV[i] <- (1/N) * sum( ((y_dat - tmp.f_hat)/(1 - tr(tmp.S_mat)/N))^2 )
}

## find the minimum point
df_GCV <- tibble(lambda = lambda, 
                 GCV = GCV, 
                 MinGCV.lambda = lambda[which(GCV == min(GCV))],
                 MinGCV.GCV = min(GCV)) 

df_GCV$MinGCV.lambda[1]; df_GCV$MinGCV.GCV[1]

## 1.3.3 plot the GCV vs Lambda
#setwd("../Figs")
#pdf("GCVvsLambda.#pdf")
ggplot(df_GCV, aes(x = lambda)) + 
  geom_line(aes(y = GCV), col = "dodgerblue3") +
  theme_university() +
  geom_point(aes(x = MinGCV.lambda[1], y = MinGCV.GCV[1]), 
             col = "darkorange1", pch = 13, size = 3) +
  ggtitle( "GCV Score for Different Lambda Values") +
  xlab(expression(paste(lambda))) +
  ylab('GCV') + 
  ylim(c(min(GCV)*0.95, max(GCV)))
dev.off()
#setwd("../Code")

## 1.4 Final plot for Lambda = 0, optimal, Inf ----

## define different lambda values
lambda_vec <- c(0, df_GCV$MinGCV.lambda[1], 100000)
## functions for betas and predictions
func.beta <- function(lambda) newBS %*% solve(crossprod(newBS) + lambda*S_matrix) %*% t(newBS)
func.pred <- function(beta_hat) beta_hat %*% y_dat

penalised_beta <- lapply(lambda_vec, func.beta)
penalised_fhat <- lapply(penalised_beta, func.pred)

final_df_P_basis <- penalised_fhat %>%
  bind_cols() %>%
  rename(Lambda0 = names(.)[1],
         LambdaOpt = names(.)[2],
         LambdaInf = names(.)[3]) %>%
  gather(Lambda_value, f_hat)  %>%
  mutate(x_data =  rep(x_dat_sorted, 3), 
         y_data = rep(y_dat, 3))

#setwd("../Figs")
#pdf("P_Spline_diff_Lambdas.#pdf")
ggplot(final_df_P_basis, aes(x =x_data)) + 
  geom_point(aes(y = y_data), col = "black") +
  geom_line(aes(y = f_hat, col = factor(Lambda_value)), size = 0.75) +
  theme_university() +
  ggtitle( "P-Spline Fit on Simulated Data") +
  xlab("x") +
  ylab('f(x)') + 
  scale_colour_discrete(name="", 
                        breaks=c("Lambda0", 
                                 "LambdaOpt", 
                                 "LambdaInf"),
                        labels=c(expression(paste(lambda, ' = ', 0)), 
                                 expression(paste(lambda, ' = ', 2.14)), 
                                 expression(paste(lambda, ' = ', infinity)))) +
  guides(colour = guide_legend(override.aes = list(size = c(2,2,2))))
dev.off()
#setwd("../Code")

# --- END --- #
