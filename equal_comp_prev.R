# Author: Khue-Tu Tran 

## Set working directory to project root if not already there
# setwd("GestALT")

## --- Load function definitions ----------------------------------------------
source("R/00_setup_equaprev.R")
source("R/01_simulate_data.R")
source("R/02_fit_models.R")
source("R/03_predict_auc.R")
source("R/04_sim_helpers.R")

## ------------------- Simulation with equal component prev --------------------
## Then run try_pred_new_ali() multiple times and compute average AUC 


#* ----------------------------- Coefficients sandbox --------------------------
# Outcome prevelance (admit:, ed: )
# Try parameters to see outcome prevelance before running simulations 
try_beta_admit <- 3 * my_beta_admit
try_beta_admit[1] <- -5
try_dat <- sim_data(N = 10000, 
                    beta_admit = try_beta_admit, 
                    beta_ed = 3 * my_beta_ed,
                    beta_ed2 = 3 * my_beta_ed2)
table(try_dat$Y_admit)

## Simulation
run_sim_print_res(sim_data, 
                  N = 1000,
                  beta_admit = try_beta_admit, 
                  beta_ed = 3 * my_beta_ed,
                  beta_ed2 = 3 * my_beta_ed2)

# covariates: age/sex (60F:40M and mean aged 48, range 18-65) 
try_beta_as_admit <- rbind(try_beta_admit, log(1.1) / 15, log(1))
try_beta_as_ed <- rbind(3 * my_beta_ed, log(0.9) / 15, log(0.8))
try_beta_as_ed2 <- rbind(3 * my_beta_ed2, log(0.9) / 15, log(0.8))
try_dat <- sim_data_agesex(N = 10000, 
                           beta_admit = try_beta_as_admit, 
                           beta_ed = try_beta_as_ed,
                           beta_ed2 = try_beta_as_ed2)
table(try_dat$Y_admit)

run_sim_print_res(sim_data_agesex, 
                  N = 1000,
                  beta_admit = try_beta_as_admit, 
                  beta_ed = try_beta_as_ed,
                  beta_ed2 = try_beta_as_ed2)


# Ordering of coeffs 

base_coefs <- seq(0, 1, by = 0.1) # magnitudes 0 to 1 (10 values)

## Scenario 1: complete opposite 
comp_order_admit1 <- 1:11          # A1C, ALB, BMI, CHOL, CRP, CREAT_C, HCST, TRIG, DIASTOLIC, SYSTOLIC
comp_order_ed1    <- c(1, 11:2)

or1_beta_admit <- matrix(c(3 * base_coefs[comp_order_admit1]), ncol = 1) 
or1_beta_admit[1] <- -5
# scaled down, reversed order
or1_beta_ed <- matrix(c(1.5 * base_coefs[comp_order_ed1]), ncol = 1) 
or1_beta_ed2 <- matrix(append(or1_beta_ed, value = 3 * 0.55, after = 1), ncol = 1) 

# without age and sex
or1_dat <- sim_data(N = 10000, 
                    beta_admit = or1_beta_admit, 
                    beta_ed = or1_beta_ed,
                    beta_ed2 = or1_beta_ed2)
table(or1_dat$Y_admit)

### auc_admit for new ali works better compared to old ali, but auc_ed and auc_ed2 are worse 
run_sim_print_res(sim_data, 
                  N = 1000,
                  beta_admit = or1_beta_admit, 
                  beta_ed = 3 * or1_beta_ed,
                  beta_ed2 = 3 * or1_beta_ed2)

# with age and sex
or1_beta_as_admit <- rbind(or1_beta_admit, log(1.1) / 15, log(1))
or1_beta_as_ed <- rbind(3 * or1_beta_ed, log(0.9) / 15, log(0.8))
or1_beta_as_ed2 <- rbind(3 * or1_beta_ed2, log(0.9) / 15, log(0.8))
or1_dat_as <- sim_data_agesex(N = 10000, 
                              beta_admit = or1_beta_as_admit, 
                              beta_ed = or1_beta_as_ed,
                              beta_ed2 = or1_beta_as_ed2)
table(or1_dat_as$Y_admit)

### auc_admit for new ali works better compared to old ali, but auc_ed and auc_ed2 are worse 
run_sim_print_res(sim_data_agesex, 
                  N = 1000,
                  beta_admit = or1_beta_as_admit, 
                  beta_ed = or1_beta_as_ed,
                  beta_ed2 = or1_beta_as_ed2)


## Scenario 2: partially overlapping ranking
comp_admit_2 <- 1:11
comp_ed_2    <- c(1, 2, 3, 6, 5, 4, 7, 8, 11, 10, 9) # swap ranks 4/6 and last two

or2_beta_admit <- matrix(3 * base_coefs[comp_admit_2], ncol = 1) 
or2_beta_admit[1] <- -5
# scaled down, reversed order
or2_beta_ed <- matrix(1.5 * base_coefs[comp_ed_2], ncol = 1) 
or2_beta_ed2 <- matrix(append(or2_beta_ed, value = 3 * 0.55, after = 8), ncol = 1) 

# without age and sex
or2_dat <- sim_data(N = 10000, 
                    beta_admit = or2_beta_admit, 
                    beta_ed = or2_beta_ed,
                    beta_ed2 = or2_beta_ed2)
table(or2_dat$Y_admit)

### new ali is better than old ali 
run_sim_print_res(sim_data, 
                  N = 1000,
                  beta_admit = or2_beta_admit, 
                  beta_ed = or2_beta_ed,
                  beta_ed2 = or2_beta_ed2)

# with age and sex
or2_beta_as_admit <- rbind(or2_beta_admit, log(1.1) / 15, log(1))
or2_beta_as_ed <- rbind(or2_beta_ed, log(0.9) / 15, log(0.8))
or2_beta_as_ed2 <- rbind(or2_beta_ed2, log(0.9) / 15, log(0.8))
or2_dat_as <- sim_data_agesex(N = 10000, 
                              beta_admit = or2_beta_as_admit, 
                              beta_ed = or2_beta_as_ed,
                              beta_ed2 = or2_beta_as_ed2)
table(or2_dat_as$Y_admit)

### new ali is better than old ali, new ali with sex and age accounted for also
### performs slightly better than when it's without age and sex 
run_sim_print_res(sim_data_agesex, 
                  N = 1000,
                  beta_admit = or2_beta_as_admit, 
                  beta_ed = or2_beta_as_ed,
                  beta_ed2 = or2_beta_as_ed2)



## Scenario 3: A1C dominant, different tails 
or3_beta_admit <- matrix(3 * c(0, 1.0, rep(0.4, 9)), ncol = 1)   # A1C dominant, rest small+equal
or3_beta_admit[1] <- -5
or3_beta_ed    <- matrix(3 * c(0, 1.0, rep(0.25, 9)), ncol = 1)  # A1C dominant, rest near-zero
or3_beta_ed2 <- matrix(3 * c(0, 1.0, rep(0.25, 9), 0.55), ncol = 1)
# without age and sex
or3_dat <- sim_data(N = 10000, 
                    beta_admit = or3_beta_admit, 
                    beta_ed = or3_beta_ed,
                    beta_ed2 = or3_beta_ed2)
table(or3_dat$Y_admit)

### new ali is better than old ali 
run_sim_print_res(sim_data, 
                  N = 1000,
                  beta_admit = or3_beta_admit, 
                  beta_ed = or3_beta_ed,
                  beta_ed2 = or3_beta_ed2)

# with age and sex
or3_beta_as_admit <- rbind(or3_beta_admit, log(1.1) / 15, log(1))
or3_beta_as_ed <- rbind(or3_beta_ed, log(0.9) / 15, log(0.8))
or3_beta_as_ed2 <- rbind(or3_beta_ed2, log(0.9) / 15, log(0.8))
or3_dat_as <- sim_data_agesex(N = 10000, 
                              beta_admit = or3_beta_as_admit, 
                              beta_ed = or3_beta_as_ed,
                              beta_ed2 = or3_beta_as_ed2)
table(or3_dat_as$Y_admit)

### new ali is better than old ali, new ali with sex and age accounted for also
### performs slightly better than when it's without age and sex 
run_sim_print_res(sim_data_agesex, 
                  N = 1000,
                  beta_admit = or3_beta_as_admit, 
                  beta_ed = or3_beta_as_ed,
                  beta_ed2 = or3_beta_as_ed2)

