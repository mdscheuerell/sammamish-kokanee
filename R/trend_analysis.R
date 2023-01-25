#### setup ####

library(here)
library(readr)
library(dplyr)
library(MARSS)

## read data
dat_counts <- read_csv(here("data", "samm_kokanee_counts.csv"))

## log-transform counts
log_counts <- dat_counts %>%
  select(count) %>%
  log() %>%
  unlist()


#### model fitting ####

control_list <- list(maxit = 5000, allow.degen = TRUE)
  
## 1) unbiased random walk

mod_list <- list(
  B = matrix(1),
  U = matrix(0),
  Q = matrix("q"),
  Z = matrix(1),
  A = matrix(0),
  R = matrix(1.5) # matrix("r")
)

mod_fit_unbiased <- MARSS(y = matrix(log_counts, nrow = 1),
                          model = mod_list, control = control_list)


## 2) biased random walk

mod_list <- list(
  B = matrix(1),
  U = matrix("u"),
  Q = matrix("q"),
  Z = matrix(1),
  A = matrix(0),
  R = matrix(1.5) # matrix("r")
)

mod_fit_biased <- MARSS(y = matrix(log_counts, nrow = 1),
                        model = mod_list, control = control_list)

## boostrapped CI's
mod_fit_biased_CI <- MARSSparamCIs(mod_fit_biased, method = "hessian", alpha = 0.05, nboot = 1000)


## 3) stationary AR(1)

mod_list <- list(
  B = matrix("b"),
  U = matrix(0),
  Q = matrix("q"),
  Z = matrix(1),
  A = matrix(0),
  R = matrix(1.5) # matrix("r")
)

mod_fit_ar1 <- MARSS(y = matrix(log_counts, nrow = 1),
                     model = mod_list, control = control_list)


#### Popn viability analysis ####

## set seed for random numbers

## number of MC simulations
n_sims <- 1000

## time horizon
n_years <- 50

## quasi-extinction threshold
qet <- log(10)

## initial state = log(count) in 2022
x0 <- log_counts[length(log_counts)]

## future population vector
# popn_vec <- c(x0, rep(NA, n_years))
popn_mat <- matrix(c(x0, rep(NA, n_years)),
                   nrow = n_sims, ncol = n_years + 1,
                   byrow = TRUE)

## simulations for unbiased
for(i in 1:n_sims) {
  for(t in 2:(n_years + 1)) {
    popn_mat[i, t] <- popn_mat[i, t - 1] + rnorm(1, 0, sqrt(mod_fit_unbiased$par$Q))
  }
}

pva_unbiased <- apply(popn_mat[, -1] < qet, 2, sum) / n_sims


## simulations for biased
for(i in 1:n_sims) {
  for(t in 2:(n_years + 1)) {
    popn_mat[i, t] <- popn_mat[i, t - 1] + mod_fit_biased$par$U + rnorm(1, 0, sqrt(mod_fit_biased$par$Q))
  }
}

pva_biased <- apply(popn_mat[, -1] < qet, 2, sum) / n_sims

