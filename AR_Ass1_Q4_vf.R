
# UCT ASSIGNMENT
# Author: Julian Albert
# Date:   20/03/2019

# Clean Environment
rm(list = ls()); dev.off()
# Assignment prelim script, useful packages and custom ggplot theme
source("../../../../PrelimScript/assignment_prelim_script.R")
## 0. Load all the packages for the assignment question
if (!require("pacman")) install.packages("pacman")
p_load(wavethresh, wmtsa, coop, rwavelet, RMThreshold)

# Question 4 ---- 
# Wavelets

## 4.1 Simulate Data ----

set.seed(123) # reproducible

n <- 1024
x_dat <- seq(0, 0.99999, length.out = n) # x from U(0, 1)

func.mu <- function(x)
{
  mu <-  5 + sin(3*pi*(x - 0.6))
  return(mu)
}

mu_dat <- func.mu(x_dat) # mean for y
sigma_dat <- 0.5 # sd for y
y_dat <- rnorm(n, mean = mu_dat, sd = sigma_dat) # y ~ N(mu, 0.5^2)

## 4.2 construct a HAAR Basis function ----

## 4.2.1 mother wavelet
func.mother_haar <- function(x)
{
  ifelse((x >= 0) & (x < 1/2), 1, 
         ifelse((x >= 1/2) & (x < 1), -1, 0))
}

## 4.2.2 daughter wavelet
func.daughter_haar <- function(j, k, x)
{
  return( 2^(j/2) * func.mother_haar((2^j)*x - k))
}

## 4.2.3 father wavelet
func.father_haar <- function(x)
{
  ifelse(x >= 0 & x <= 1, 1, 0)
}

## 4.2.4 So constructing a HAAR Basis
func.haar_basis <- function(J, data)
{
  
haar_basis <- matrix(0, length(data), 2^J) # size of basis matrix

c00 <- func.father_haar(data)
haar_basis[, 1] <- c00 # specify first column as father

# update the matrix using the double sum boi
for(j in 0:(J-1)){
  for(k in 0:((2^j) - 1)){
    index <- 2^(j) + k + 1
    haar_basis[, index] <- func.daughter_haar(j, k, data)
  }
}

return(haar_basis)
}

## 4.3 Predicting Using our HAAR basis for different J ----

## 4.3.1 Prediction function
func.wave_pred <- function(x_data, y_data, basis.func, J)
{
  func.haar_basis <- match.fun(basis.func)
  
  haar_basis <- func.haar_basis(J, x_data)
  beta_hat <- solve(t(haar_basis) %*% haar_basis) %*% t(haar_basis) %*% y_data
  f_hat <- haar_basis %*% beta_hat
  
  return(list(pred = f_hat,
              basis = haar_basis))
}

## 4.3.2 Results

resolutions_vec <- c(2, 4, 6) # different resolutions

haar_bases <- lapply(resolutions_vec, 
                     function(x) func.wave_pred(x_dat, y_dat, func.haar_basis, x)$basis)

fhats <- lapply(resolutions_vec, 
                function(x) func.wave_pred(x_dat, y_dat, func.haar_basis, x)$pred) %>%
  bind_cols() %>%
  rename(J2 = names(.)[1],
         J4 = names(.)[2],
         J6 = names(.)[3]) %>%
  mutate(x = x_dat, y = y_dat)

##setwd('../Figs')

#pdf("J2_wavelet.#pdf")
ggplot(fhats, aes(x = x)) +
  geom_line(aes(y = y), col = "gray") +
  geom_line(aes(y = fhats$J2), col = "dodgerblue2", size = 1) +
  theme_university() +
  theme(legend.position = "right", 
        legend.direction = "vertical") +
  ylab("f(x)")
dev.off()

#pdf("J4_wavelet.#pdf")
ggplot(fhats, aes(x = x)) +
  geom_line(aes(y = y), col = "gray") +
  geom_line(aes(y = fhats$J4), col = "darkorange1", size = 1) +
  theme_university() +
  theme(legend.position = "right", 
        legend.direction = "vertical") +
  ylab("f(x)")
dev.off()

#pdf("J6_wavelet.#pdf")
ggplot(fhats, aes(x = x)) +
  geom_line(aes(y = y), col = "gray") +
  geom_line(aes(y = fhats$J6), col = "green4", size = 1) +
  theme_university() +
  theme(legend.position = "right", 
        legend.direction = "vertical") +
  ylab("f(x)")
dev.off()
#setwd('../Code')

mean((fhats$J2 - fhats$y)^2)
mean((fhats$J4 - fhats$y)^2)
mean((fhats$J6 - fhats$y)^2)

## 4.4 Using the bump function ----

v <- DJ.EX() # generate a signal
x_dat_qb <- (1:1024) / 1024 # x_data
v_noise <- DJ.EX(rsnr=sqrt(2), noisy=TRUE)$bumps # add noise

### Plot the Signals
qb_df <- data.frame(x = x_dat_qb, y = v$bumps,
                    y_w.noise =  v_noise)

init_signal <- ggplot(qb_df, aes(x = x, y = y)) +
  geom_line() +
  ggtitle("Signal without Noise") +
  theme_university()

noisy_signal <- ggplot(qb_df, aes(x = x, y = y_w.noise)) +
  geom_line() +
  theme_university() +
  ggtitle('Signal with Noise') +
  ylab("y")

#setwd('../Figs')
#pdf('wavelet_signals.#pdf')
grid.arrange(init_signal, noisy_signal, ncol = 2)
dev.off()

## 4.4.1 Wavelet transform and thresholding

Y <- qb_df$y_w.noise

#pdf('wavelet_DWT_haar.#pdf')
plot(wd(Y, filter.number = 1, family = 'DaubExPhase'))
dev.off()

#pdf('wavelet_DWT_thresh__haar.#pdf')
wd(Y, filter.number = 1, family = 'DaubExPhase') %>%
  threshold.wd(type = 'soft',
               policy = 'BayesThresh') %>%
  plot()
dev.off()


# sparsity(wavDWTMatrix(wavelet = 'haar', J = 10)) # 0.9892578
# sparsity(wavDWTMatrix(wavelet = 's8', J = 10)) # 0.9501648

thresh_hard <- wavShrink(Y, wavelet = 'haar', xform = 'dwt', 
                         n.level = 10, shrink.fun = 'hard')
thresh_soft <- wavShrink(Y, wavelet = 'haar', xform = 'dwt', 
                         n.level = 10, shrink.fun = 'soft')

df_thresh <- data.frame(Y = Y, x = x_dat_qb,
                        Hard = thresh_hard,
                        Soft = thresh_soft)

hard_thresh_plot <- ggplot(df_thresh, aes(x = x_dat_qb)) +
  geom_line(aes(y = Y)) +
  geom_line(aes(y = thresh_hard), col = 'dodgerblue3') +  
  theme_university() +
  xlab('x') + ylab('f(x)') +
  ggtitle('Hard Threshold')

soft_thresh_plot <- ggplot(df_thresh, aes(x = x_dat_qb)) +
  geom_line(aes(y = Y)) +
  geom_line(aes(y = thresh_soft), col = 'darkorange2') +
  theme_university() +
  xlab('x') + ylab('f(x)') +
  ggtitle('Soft Threshold')

#setwd('../Figs')
#pdf('wavelet_thresholds.#pdf')
grid.arrange(hard_thresh_plot, soft_thresh_plot, ncol = 2)
dev.off()

## 4.4.2 choose Hard or Soft
thresh_uni_hard <- wavShrink(Y, xform = 'dwt', n.level = 10, shrink.fun = 'hard',
                        thresh.fun = 'universal', wavelet = 'haar')
thresh_uni_soft <- wavShrink(Y, xform = 'dwt', n.level = 10, shrink.fun = 'soft',
                        thresh.fun = 'universal', wavelet = 'haar')
thresh_minimax_hard <- wavShrink(Y, xform = 'dwt', n.level = 10, shrink.fun = 'hard',
                            thresh.fun = 'minimax', wavelet = 'haar')
thresh_minimax_soft <- wavShrink(Y, xform = 'dwt', n.level = 10, shrink.fun = 'soft',
                                 thresh.fun = 'minimax', wavelet = 'haar')
thresh_adaptive_soft <- wavShrink(Y, xform = 'dwt', n.level = 10, shrink.fun = 'soft',
                             thresh.fun = 'adaptive', wavelet = 'haar')

df_shrink <- data.frame(UniversalHard = thresh_uni_hard,
                        UniversalSoft = thresh_uni_soft,
                        MiniMaxHard = thresh_minimax_hard,
                        MiniMaxSoft = thresh_minimax_soft,
                        AdaptiveSoft = thresh_adaptive_soft)

### By how much do we reduct the noisy signal?
MSEs <- apply(df_shrink, 2, function(x) mean((x - v$bumps)^2))
SNR <- apply(df_shrink, 2, function(x) sd(x)/3.5)

#setwd('../Figs')
#pdf('thrsh_universal_hard.#pdf')
ggplot(df_shrink, aes(x = x_dat_qb)) +
  geom_line(aes(y = Y)) +
  geom_line(aes(y = UniversalHard), col = 'dodgerblue3', size = 1) +
  theme_university() +
  xlab('x') + ylab('f(x)') 
dev.off()

#pdf('thrsh_minimax_hard.#pdf')
ggplot(df_shrink, aes(x = x_dat_qb)) +
  geom_line(aes(y = Y)) +
  geom_line(aes(y = MiniMaxHard), col = 'darkorange2', size = 1) +
  theme_university() +
  xlab('x') + ylab('f(x)')
dev.off()

#pdf('thrsh_adaptive_soft.#pdf')
ggplot(df_shrink, aes(x = x_dat_qb)) +
  geom_line(aes(y = Y)) +
  geom_line(aes(y = AdaptiveSoft), col = 'green4', size = 1) +
  theme_university() +
  xlab('x') + ylab('f(x)') 
dev.off()

#pdf('thrsh_universal_true.#pdf')
ggplot(df_shrink, aes(x = x_dat_qb)) +
  geom_line(aes(y = UniversalHard), col = 'dodgerblue3', size = 1) +
  geom_line(aes(y = v$bumps)) +
  theme_university() +
  xlab('x') + ylab('f(x)') 
dev.off()

#pdf('thrsh_minimax_true.#pdf')
ggplot(df_shrink, aes(x = x_dat_qb)) +
  geom_line(aes(y = MiniMaxHard), col = 'darkorange2', size = 1) +
  geom_line(aes(y = v$bumps)) +
  theme_university() +
  xlab('x') + ylab('f(x)')
dev.off()

#pdf('thrsh_adaptive_true.#pdf')
ggplot(df_shrink, aes(x = x_dat_qb)) +
  geom_line(aes(y = AdaptiveSoft), col = 'green4', size = 1) +
  geom_line(aes(y = v$bumps)) +
  theme_university() +
  xlab('x') + ylab('f(x)') 
dev.off()

# --- END --- #
