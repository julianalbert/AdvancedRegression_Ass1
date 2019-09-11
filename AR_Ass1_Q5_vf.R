
# UCT ASSIGNMENT
# Author: Julian Albert
# Date:   18/03/2019

# Clean Environment
rm(list = ls()); dev.off()
# Assignment prelim script, useful packages and custom ggplot theme
source("../../../../PrelimScript/assignment_prelim_script.R")
## 0. Load all the packages for the assignment question
if (!require("pacman")) install.packages("pacman")
p_load(splines, ElemStatLearn)

# Question 5 ---- 

## 5.1 Data ----

set.seed(123) # reproducible
data(phoneme)
dat_phoneme <- phoneme

dat_aa <- dat_phoneme[phoneme$g == "aa", ]
dat_ao <- dat_phoneme[phoneme$g == "ao", ]

## 5.2 Replicate image in ESL fig 5.5 top panel ----
###log-periodograms for each of the two phonemes "aa" and "ao" measured at 256 frequencies.
total_sample_data <- rbind(dat_aa, dat_ao)
idx_sample <- sample((1:NROW(total_sample_data)), 1000)

dat_sample_train <- total_sample_data[idx_sample, ] %>%
  dplyr::select(-c(speaker))

dat_sample_train_aa <- dat_sample_train %>%
  dplyr::filter(g == "aa")
dat_sample_train_ao <- dat_sample_train %>%
  dplyr::filter(g == "ao")

#setwd("../Figs")
#pdf("phoneme_fig1.#pdf",width=10,height=5)
plot(as.numeric(dat_sample_train_aa[1, 1:256]), type = 'l', 
     col = 'green', ylim = c(0, 27), ylab = "Log-periodogram",
     xlab = "Frequency", main = "Phoneme Examples")
for(i in 2:7) lines(as.numeric(dat_sample_train_aa[i, 1:256]),
                     col = 'green')
for(i in 1:8) lines(as.numeric(dat_sample_train_ao[i, 1:256]),
                     col = "orange")
for(i in 7:15) lines(as.numeric(dat_sample_train_aa[i, 1:256]),
                     col = 'green')
for(i in 8:15) lines(as.numeric(dat_sample_train_ao[i, 1:256]),
                     col = "orange")
legend("topright", inset=.02, legend = c("aa", "ao"), 
       col=c("green", "orange"), lty = c(1,1), box.lty = 0)
dev.off()
#setwd("../Code")

## 5.3 Replicate image in ESL fig 5.5 bottom panel ----
### function bs(1:256) basis on frequencies!!!

dat_sample_train$g <- as.factor(as.character(dat_sample_train$g))

model1 <- glm(g ~ ., family = binomial(link = "logit"), data = dat_sample_train)
plot(coef(model1), type = 'l', ylim = c(-0.5, 0.5), 
     col = 'gray', ylab = "Logistic Regression Coefficients", xlab = "Frequency", 
     main = "Phoneme Classification: Raw and Restricted Logistic Regression")

q5_B_basis <- ns(1:(dim(dat_sample_train)[2]-1), df = 12, intercept = FALSE) # H_pxM

x_test1 <- (dat_sample_train)[,-NCOL(dat_sample_train)]
testfun <- function(x_vector) t(q5_B_basis) %*% as.numeric(x_vector)

### get x* for all rows 
tmp <- apply(x_test1, 1, testfun) %>%
  t() %>% as.data.frame() %>%
  mutate(g = dat_sample_train$g)

model2 <- glm(g ~ .,family = binomial(link = "logit"), data = tmp)

### beta is now B*theta
new_beta_hat <- q5_B_basis %*% coef(model2)[-1]

#setwd("../Figs")
#pdf("phoneme_fig2.#pdf",width=10,height=5)
plot(coef(model1), type = 'l', ylim = c(-0.5, 0.5), 
     col = 'gray', ylab = "Logistic Regression Coefficients", xlab = "Frequency", 
     main = "Phoneme Classification: Raw and Restricted Logistic Regression")
lines(new_beta_hat, col = "red", lwd = 1.5)
abline(h = 0)
dev.off()
#setwd("../Code")

# --- END --- #