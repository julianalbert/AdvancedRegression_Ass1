
# UCT ASSIGNMENT
# Author: Julian Albert
# Date:   11/03/2019

# Clean Environment
rm(list = ls()); dev.off()
# Assignment prelim script, useful packages and custom ggplot theme
source("../../../../PrelimScript/assignment_prelim_script.R")
## 0. Load all the packages for the assignment question
if (!require("pacman")) install.packages("pacman")
p_load(gamair, mgcv, mgcViz, earth, mda,
       splines, itsadug, vip, pdp, caret)

# Question 3 ---- 
# GAMS and MARS

## 3.1 Get the Data ----
set.seed(123)
data(chl)
chl_dat <- chl

## 3.2 GAMS ----

## 3.2.1 check the relationships between vars
#setwd('../Figs')
#pdf("VariablePairs.#pdf")
pairs(chl_dat, pch = ".")
dev.off()

## basic transformation
trans2_chl_dat <- chl_dat %>%
  mutate(chl.sw = sqrt(chl.sw),
         bath = sqrt(bath))

#pdf("VariablePairsTransformed.#pdf") # basic power
pairs(trans2_chl_dat, pch = ".")
dev.off()
#setwd('../Code')

## 3.2.2 Fitting a GAM
### split into training and testing
idx <- sample(0.60*NROW(trans2_chl_dat), replace = FALSE)
training_dat <- trans2_chl_dat[idx, ]
testing_dat <- trans2_chl_dat[-idx, ]

### Different GAMs
gam0_formula <- (chl ~ s(chl.sw, bs="cr", k=20) + s(bath, bs="cr", k=50) + 
                   s(jul.day, bs="cr", k=20))
gam0 <- gam(gam0_formula, data = training_dat)
gam0_pred <- predict.gam(gam0, testing_dat) 
gam0_mse <- sum((gam0_pred - testing_dat$chl)^2) #41578.47
summary(gam0); AIC(gam0); gam.check(gam0)
#R2 = 0.254, GCV = 2.031, AIC = 29450.91, 25.9%
#edf close for lon,lat and jul.day so increase k

gam1_formula <- (chl ~ s(chl.sw, bs="cr", k=20) + s(bath, bs="cr", k=50) + 
                   s(jul.day, bs="cr", k=20) + s(lon) + s(lat))
gam1 <- gam(gam1_formula, data = training_dat)
gam1_pred <- predict.gam(gam1, testing_dat) 
gam1_mse <- sum((gam1_pred - testing_dat$chl)^2) #43938.44
summary(gam1); AIC(gam1); gam.check(gam1)
#R2 = 0.293, GCV = 1.9297, AIC = 29025.97, 30%
#edf close for lon,lat and jul.day so increase k

gam2_formula <- (chl ~ s(chl.sw, bs="cr", k=20) + s(bath, bs="cr", k=50) + 
                   s(jul.day, bs="cr", k=50) + s(lon, lat, k=50))
gam2 <- gam(gam2_formula, data = training_dat)
gam2_pred <- predict.gam(gam2, testing_dat)
gam2_mse <- sum((gam2_pred - testing_dat$chl)^2) #46702.35
summary(gam2); AIC(gam2); gam.check(gam2)
#R2 = 0.348, GCV = 1.7912, AIC = 28405.97, 35.8%
#edf close for lon,lat and jul.day so reduce k

gam3_formula <- (chl ~ s(chl.sw, bs="cr", k=20) + s(bath, bs="cr", k=50) + 
                   s(jul.day, bs="cr", k=50) + ti(lon, lat, k = 20))
gam3 <- gam(gam3_formula, data = training_dat)
gam3_pred <- predict.gam(gam3, testing_dat)
gam3_mse <- sum((gam3_pred - testing_dat$chl)^2) #320436.6
summary(gam3); AIC(gam3); gam.check(gam3)
#R2 = 0.388, GCV = 1.7196, AIC = 28057.22, 41.1%
#edf close for lon,lat and jul.day so reduce k

gam4_formula <- (chl ~ s(chl.sw, bs="cr", k=20) + s(bath, bs="cr", k=20) + 
                   s(jul.day, bs="cr", k=100) + te(lon, lat, k = 20))
gam4 <- gam(gam4_formula, data = training_dat)
gam4_pred <- predict.gam(gam4, testing_dat)
gam4_mse <- sum((gam4_pred - testing_dat$chl)^2) #288011.3
summary(gam4); AIC(gam4); gam.check(gam4)
#R2 = 0.399, GCV = 1.6947, AIC = 27933.77, 42.4%
#edf close for lon,lat and jul.day so reduce k

gam5_formula <- (chl ~ s(chl.sw, bs="cr", k=20) + s(bath, bs="cr", k=20) + 
                   s(jul.day, bs="cr", k=100) + te(lon, lat, k = 20) +
                   te(chl.sw, bath))
gam5 <- gam(gam5_formula, data = training_dat)
gam5_pred <- predict.gam(gam5, testing_dat)
gam5_mse <- sum((gam5_pred - testing_dat$chl)^2) #325631.1
summary(gam5); AIC(gam5); gam.check(gam5)
#R2 = 0.404, GCV = 1.6833, AIC = 27877.57, 42.8%
#edf close for lon,lat and jul.day so reduce k

model_list <- list(model0 = gam0, model1 = gam1, 
                   model2 = gam2, model3 = gam3, 
                   model4 = gam4, model5 = gam5)

### Comparing the GCV & AIC scores for each model
model_gcvs <- lapply(model_list, function(x) x$gcv.ubre) %>%
  bind_rows()
model_aics <- lapply(model_list, AIC) %>%
  bind_rows()

### select optimal model
optimal_model <- model_list$model5
optimal_model.Viz <- getViz(optimal_model)

## 3.2.3 Check the residuals
#setwd('../Figs')
#pdf("gam_check.#pdf")
check(optimal_model.Viz,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))
dev.off()

# 3.2.4 Optimal GAM performance
plot_chlsw <- plot(sm(optimal_model.Viz, 1)) +
  l_fitLine(colour = "darkorange1") +
  l_ciLine(mul = 5, colour = "dodgerblue3", linetype = 2) + 
  theme_university()

plot_bath <- plot(sm(optimal_model.Viz, 2)) +
  l_fitLine(colour = "darkorange1") +
  l_ciLine(mul = 5, colour = "dodgerblue3", linetype = 2) + 
  theme_university()

plot_julday <- plot(sm(optimal_model.Viz, 3)) +
  l_fitLine(colour = "darkorange1") +
  l_ciLine(mul = 5, colour = "dodgerblue3", linetype = 2) + 
  theme_university()

#setwd('../Figs')
#pdf("gam_chlsw.#pdf")
plot_chlsw
dev.off()
#pdf("gam_bath.#pdf")
plot_bath
dev.off()
#pdf("gam_julday.#pdf")
plot_julday
dev.off()

# pvisgam(optimal_model.Viz,  n.grid = 50, theta = 310, phi = 0,
#         select = 4, plot.type = "persp", zlim = c(-30, 30)) # flat

#pdf("gam_chlswbath.#pdf")
pvisgam(optimal_model.Viz,  n.grid = 50, theta = 130, phi = 32,
        select = 5, plot.type = "persp")
dev.off()
#setwd('../Code')

## 3.3 Fitting a MARS ----

mars1 <- earth(chl ~ ., data = training_dat)
summary(mars1) %>% .$coefficients %>% head(10)

mars2 <- earth(chl ~ ., data = training_dat, degree = 2)
summary(mars2) %>% .$coefficients %>% head(10)

## 3.3.1 create a tuning grid for degree and nprune
hyper_grid <- expand.grid(degree = 1:2, 
  nprune = seq(2, 100, length.out = 10) %>% floor())

#### cross validated model
set.seed(123) # for reproducibiity
tuned_mars <- train(
  x = subset(training_dat, select = -chl),
  y = training_dat$chl,
  method = "earth",
  metric = "RMSE",
  trControl = trainControl(method = "cv", number = 10),
  tuneGrid = hyper_grid)

### Best model
tuned_mars$bestTune
mars_optimal <- tuned_mars$finalModel
summary(mars_optimal) %>% .$coefficients %>% head(10)

## 3.3.1 plot results

#setwd('../Figs')
#pdf('MARS_RMSEplot.#pdf')
ggplot(tuned_mars) + 
  theme_university() +
  theme(legend.position = "right", 
        legend.direction = "vertical") +
  xlab("No. of Terms")
dev.off()

### variabe importance plots
plot_vdp_GCV <- vip(tuned_mars, num_features = 40, 
                    bar = FALSE, value = "gcv") + 
  ggtitle("GCV") +
  theme_university()

plot_vdp_RSS <- vip(tuned_mars, num_features = 40, 
                    bar = FALSE, value = "rss") + 
  ggtitle("RSS") +
  theme_university()

#pdf('MARS_vimp_plot.#pdf')
grid.arrange(plot_vdp_GCV, plot_vdp_RSS, ncol = 2)
dev.off()

### partial dependence plots
pdp_chl.sw <- partial(tuned_mars, pred.var = "chl.sw", 
                      grid.resolution = 10) %>% autoplot() +
  theme_university()
pdp_jul.day <- partial(tuned_mars, pred.var = "jul.day",
                       grid.resolution = 10) %>% autoplot() +
  theme_university()
pdp_bath <- partial(tuned_mars, pred.var = "bath",
                    grid.resolution = 10) %>% autoplot() +
  theme_university()

#pdf('MARS_pdp_plot.#pdf')
grid.arrange(pdp_chl.sw, pdp_jul.day, pdp_bath, ncol = 3)
dev.off()

#setwd('../Code')

## 3.4 Compare OOS models ----
MARS_pred <- predict(mars_optimal, testing_dat) # GCV = 1.991194
MARS_mse <- mean((testing_dat$chl - MARS_pred)^2) #7.765728

gam5_pred <- predict.gam(gam5, testing_dat) #GCV = 1.683286
gam5_mse <- mean((testing_dat$chl - gam5_pred)^2) #58.82065




# --- END --- #
