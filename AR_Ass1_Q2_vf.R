
# UCT ASSIGNMENT
# Author: Julian Albert
# Date:   18/03/2019

# Clean Environment
rm(list = ls()); dev.off()
# Assignment prelim script, useful packages and custom ggplot theme
source("../../../../PrelimScript/assignment_prelim_script.R")
## 0. Load all the packages for the assignment question
if (!require("pacman")) install.packages("pacman")
p_load(lattice, latticeExtra, viridis)

# Question 2 ----
# Thin-Plate Splines

## 2.1 Functions needed for tps data ----
true_function <- function(x, z, sx = 0.3, sz = 0.4)
{
  1.2*exp(-(x - 0.2)^2/sx^2-(z - 0.3)^2/sz^2) + 
    0.8*exp(-(x - 0.7)^2/sx^2-(z - 0.8)^2/sz^2)
}

## 2.1.1 thin-plate basis function
tpbf <- function(x, knot)
{
  euclid_dist <- sqrt(sum((x - knot)^2))
  tp <- ifelse(euclid_dist == 0, 0, 
               1/(8*pi) * euclid_dist^2 * log(euclid_dist))
  return(tp)
}

## 2.2 Function to get Thin-Plates for alternate Methods ----
func.tps <- function(init_n = 100, xrange = NULL, yrange = NULL, 
                     tpb_function = NULL, true_function = NULL, 
                     n.knots = 16, m = 2, d = 2, method, seed = 123)
{
  if(method == "Rank"){
    set.seed(seed) # reproducibility
    
    true_function <- match.fun(true_function)  
    tpbf <- match.fun(tpb_function)
    
    # data
    training_data <- matrix(runif(2*init_n), init_n, 2) # 100x2 matrix of random xvals
    x_knots <- training_data # every knot is an observation
    y_dat <- true_function(training_data[,1], training_data[,2]) # y_dat f(x, z)
    y_dat_w.noise <- y_dat + rnorm(init_n, 0, 0.1) # add noise
    
    # parameters
    n <- length(xrange) # number of points
    m <- m
    d <- d
    M <- choose((m+d-1), d) # number of linearly independent functions
    
    ## 2.2.1 Construct a Rank 16 approx matrix, wood2003thin Appendix ----
    
    ### Wood step 1 - Get the Enxn matrix
    E_mat_rank <- matrix(0, n, n)
    
    for (i in 1:n){
      for (j in 1:n){
        E_mat_rank[i, j] <- tpbf(x = training_data[i, ], knot = x_knots[j, ])
      }
    }
    
    ### Design matrix
    T_mat_rank <- matrix(0, n, M)
    T_mat_rank[, 1] <- 1
    T_mat_rank[, 2] <- training_data[, 1]
    T_mat_rank[, 3] <- training_data[, 2]
    
    ### Wood step 2 - Spectral Decomposition of E = U_D_Ut
    E_eigen.decomp <- eigen(E_mat_rank)
    U_vec <- E_eigen.decomp$vectors # U
    U_vec_truncated <- U_vec[, 1:n.knots] # rank = k so truncate
    D_mat <- diag(E_eigen.decomp$values) # D
    D_mat_truncated <-  diag(E_eigen.decomp$values[1:n.knots])
    Ut_vec_truncated <- t(U_vec_truncated) # U transpose
    
    ### Wood step 3 - QR decomp s.t QkRk <- t(U)kT
    QR.decomp <- qr(Ut_vec_truncated %*% T_mat_rank)
    Q_mat <- qr.Q(QR.decomp, complete = TRUE) # Q
    R_mat <- qr.R(QR.decomp, complete = TRUE) # R
    
    Zk <- Q_mat[, (n.knots-(n.knots-M)+1):n.knots] # Z basis of null
    
    ### Wood step 4 - X matrix
    X_rank <- cbind(U_vec_truncated %*% D_mat_truncated %*% Zk, T_mat_rank)
    ### where S is no penalty we have standard OLS for Beta
    beta_hat <- solve(crossprod(X_rank)) %*% t(X_rank) %*% y_dat_w.noise
    
    ### Wood step 5 - Predictions i.e. g(x)
    Uk.Zk <- U_vec_truncated %*% Zk
    
    func.predict <- function(data, # data matrix
                             betas, # beta vector
                             x_knots, # knot matrix, every obs = knots
                             M) # no. of linealry independent fx
    {
      n <- NROW(x_knots)
      k <- length(betas) - M
      a <- numeric(k)
      
      ## get the Enxn matrix
      for (i in 1:k){
        tmp.a <- numeric(n)
        for (j in 1:n){
          tmp.a[j] <- tpbf(x = data, knot = x_knots[j, ])*Uk.Zk[j, i]
        }
        a[i] <- sum(tmp.a)
      }
      
      b <- c(1, data)
      # rewrite the betas
      delta <- as.matrix(betas[1:k, ])
      alpha <- as.matrix(betas[(k+1):NROW(betas), ])
      # get predictions
      g_x <- t(a) %*% delta + t(b) %*% alpha
    }
    
    ### to get a surface need more points
    test_points <- expand.grid(x = xrange, y = yrange)
    
    z_rank <- apply(test_points, MARGIN = 1, func.predict, 
                    betas = beta_hat,
                    x_knots = x_knots,
                    M = M)
    
    fhat_dat_plot_rank <- matrix(z_rank, nrow = length(xrange))
    
    return(list(rank_pred = fhat_dat_plot_rank, 
                rank_xrange = xrange,
                rank_yrange = yrange))
  } else if(method == "Knot"){
    set.seed(seed)
    
    # Data
    x_dat_test <- expand.grid(x = xrange, y = yrange)
    n <- NROW(x_dat_test)
    y_dat_test <- true_function(x_dat_test[,1], x_dat_test[,2]) # y_dat f(x, z)
    y_dat_w.noise_test <- y_dat_test + rnorm(n, 0, 0.1) # add noise
    
    # Parameters
    m <- m
    d <- d
    M <- choose((m+d-1), d) ### number of linearly independent functions
    n.knots <- n.knots
    
    ## 2.2.2 Construct a 16 knot matrix, wood2003thin ----
    
    ### knots are in centre of the domain in 16 equal squares
    even_space <- 1/4
    knot_vector <- seq(even_space/2, (1-even_space/2), even_space)
    knot_matrix <- expand.grid(knot_vector, knot_vector)
    
    k <- NROW(knot_matrix)
    E_mat_knot <- matrix(0, n, k)
    
    ### get the Enxn matrix
    for (i in 1:n){
      for (j in 1:k){
        E_mat_knot[i, j] <- tpbf(x = x_dat_test[i, ], knot = knot_matrix[j, ])
      }
    }
    
    ### Same as Before
    T_mat_knot <- as.matrix(cbind(1, x_dat_test))
    X_knot <- as.matrix(cbind(E_mat_knot, T_mat_knot))
    ### normal ols
    beta_hat <- solve(t(X_knot) %*% X_knot) %*% t(X_knot) %*% y_dat_w.noise_test
    tmp_fhat <- X_knot %*% beta_hat
    fhat_dat_plot_rank <- matrix(tmp_fhat, nrow = length(xrange))
    
    return(list(knot_pred = fhat_dat_plot_rank, 
                knot_xrange = xrange,
                knot_yrange = yrange))
  } else(print("Method not valid!"))
}

## 2.3 Get the results ----

# parameters
set.seed(123)
init_n <- 100
xrange <- seq(0, 1, length = init_n)
yrange <- seq(0, 1, length = init_n)
m <- 2
d <- 2
n.knots <- 16

## 2.3.1 Surface Responses
true_y <- outer(xrange, yrange, true_function)
knot_method <- func.tps(xrange = xrange, yrange = yrange, method = "Knot")
rank_method <- func.tps(init_n = 100, xrange= xrange, 
                  yrange= yrange, tpb_function = tpbf, 
                  true_function = true_function, 
                  n.knots = 16, m = 2, d = 2, method = "Rank")

#setwd("../Figs")
#pdf("tps_true_surface.#pdf")
wireframe(true_y, row.values = xrange, col.values = yrange,
          drape = TRUE,
          scales = list(arrows = TRUE, cex = 0.5, 
                        z = list(arrows=F), distance =c(1, 1, 1)),
          col.regions = plasma(100),
          colorkey = FALSE,
          screen = list(z = 300, x = -75, y = 3),
          xlab = "x", ylab = "z", zlab = "y",
          main = "Thin-Plate Regression Spline Surface")
dev.off()

#pdf("tps_knot16_surface.#pdf")
wireframe(knot_method$knot_pred, row.values = xrange, col.values = yrange,
          drape = TRUE,
          scales = list(arrows = TRUE, cex = 0.5, 
                        z = list(arrows=F), distance =c(1, 1, 1)),
          col.regions = plasma(100),
          colorkey = FALSE,
          screen = list(z = 300, x = -75, y = 3),
          xlab = "x", ylab = "z", zlab = "y",
          main = "Thin-Plate Spline Surface")
dev.off()

#pdf("tps_rank16_surface.#pdf")
wireframe(rank_method$rank_pred, row.values = xrange, col.values = yrange,
          drape = TRUE,
          scales = list(arrows = TRUE, cex = 0.5, 
                        z = list(arrows=F), distance =c(1, 1, 1)),
          col.regions = plasma(100),
          colorkey = FALSE,
          screen = list(z = 300, x = -75, y = 3),
          xlab = "x", ylab = "z", zlab = "y",
          main = "Thin-Plate Regression Spline Surface")
dev.off()

#setwd('../Code')

# MSE
mean((knot_method$knot_pred - true_y)^2)
mean((rank_method$rank_pred - true_y)^2)

# --- END --- #





