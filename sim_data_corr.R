# Opt1: Multinormmal and pick correlated structure

mu_bio <- c(NUM_A1C = 6.286808, NUM_ALB = 4.277495, NUM_BMI = 30.729902,
            NUM_CHOL = 186.095427, NUM_CRP = 20, NUM_CREAT_C = 107,
            NUM_HCST = 10.511111, NUM_TRIG = 146.642879,
            NUM_BP_DIASTOLIC = 77.033838, NUM_BP_SYSTOLIC = 126.507262)

# original 
# var_bio <- c(NUM_A1C = 2.377849, NUM_ALB = 0.1010428, NUM_BMI = 58.64435,
#              NUM_CHOL = 1465.124, NUM_CRP = 2842.631, NUM_CREAT_C = 2835.045,
#              NUM_HCST = 7.641111, NUM_TRIG = 11536.85,
#              NUM_BP_DIASTOLIC = 70.23248, NUM_BP_SYSTOLIC = 162.4087)

var_bio <- c(NUM_A1C = 2.377849, NUM_ALB = 0.1010428, NUM_BMI = 58.64435,
             NUM_CHOL = 500, NUM_CRP = 16, NUM_CREAT_C = 144,
             NUM_HCST = 7.641111, NUM_TRIG = 500,
             NUM_BP_DIASTOLIC = 70.23248, NUM_BP_SYSTOLIC = 162.4087)
sd_bio <- sqrt(var_bio)

vars <- c("NUM_A1C", "NUM_ALB", "NUM_BMI", "NUM_CHOL", "NUM_CRP",
          "NUM_CREAT_C", "NUM_HCST", "NUM_TRIG",
          "NUM_BP_DIASTOLIC", "NUM_BP_SYSTOLIC")

R_bio <- matrix(c(
  # A1C          ALB           BMI           CHOL          CRP           CREAT_C  HCST          TRIG          BP_DIA        BP_SYS
  1.00000000,  -0.19846846,   0.20581931,  -0.05092428,  -0.50856668,   NA,      0.02343128,   0.32424369,   0.00556571,   0.24252539,  # A1C
  -0.19846846,   1.00000000,  -0.28579417,   0.08579578,  -0.04093183,   1,      -0.16064065,  -0.02887856,   0.02996051,  -0.07464056,  # ALB
  0.20581931,  -0.28579417,   1.00000000,  -0.00516256,   0.96173205,  -1,       0.09779848,   0.22172368,   0.22551612,   0.34733560,  # BMI
  -0.05092428,   0.08579578,  -0.00516256,   1.00000000,   0.36868348,   NA,      0.16742086,   0.32540027,   0.08608552,   0.09205011,  # CHOL
  -0.50856668,  -0.04093183,   0.96173205,   0.36868348,   1.00000000,   NA,      NA,          -0.02548241,   0.03669368,   0.11645110,  # CRP
  NA,           1.00000000,  -1.00000000,   NA,           NA,           1,       NA,          -1.00000000,  -1.00000000,  -1.00000000,  # CREAT_C
  0.02343128,  -0.16064065,   0.09779848,   0.16742086,   NA,           NA,      1.00000000,  -0.04092368,  -0.64630785,  -0.01319507,  # HCST
  0.32424369,  -0.02887856,   0.22172368,   0.32540027,  -0.02548241,  -1,      -0.04092368,   1.00000000,   0.12859438,   0.15821455,  # TRIG
  0.00556571,   0.02996051,   0.22551612,   0.08608552,   0.03669368,  -1,      -0.64630785,   0.12859438,   1.00000000,   0.64415706,  # BP_DIASTOLIC
  0.24252539,  -0.07464056,   0.34733560,   0.09205011,   0.11645110,  -1,      -0.01319507,   0.15821455,   0.64415706,   1.00000000   # BP_SYSTOLIC
), nrow = 10, ncol = 10, byrow = TRUE,
dimnames = list(vars, vars))

## sanity check: should be symmetric (allowing for NA)
### all.equal(R_bio, t(R_bio))

n_cor_comp <- matrix(
  c(468, 430, 468, 426,   5,   1,   4, 426,  468,  468,
    430, 810, 810, 703,   7,   2,   7, 703,  810,  810,
    468, 810, 998, 744,   7,   2,   9, 744,  998,  998,
    426, 703, 744, 744,   7,   0,   7, 744,  744,  744,
    5,   7,   7,   7,   7,   0,   0,   7,    7,    7,
    1,   2,   2,   0,   0,   2,   0,   0,    2,    2,
    4,   7,   9,   7,   0,   0,   9,   7,    9,    9,
    426, 703, 744, 744,   7,   0,   7, 744,  744,  744,
    468, 810, 998, 744,   7,   2,   9, 744, 1000, 1000,
    468, 810, 998, 744,   7,   2,   9, 744, 1000, 1000),
  nrow = 10, ncol = 10, byrow = TRUE
)

rownames(n_cor_comp) <- vars
colnames(n_cor_comp) <- vars


## Build correlation for components with insufficient data points
## Threshold for correlation estimate
min_n <- 10

## Cap for how large a shrunk correlation can be (only kicks in for LOW-N pairs)
shrink_cap <- 0.3

cor_comp_adj <- R_bio # copy 

## Find all off-diagonal cells estimated from too few complete pairs
low_n_idx <- which(n_cor_comp < min_n & row(n_cor_comp) != col(n_cor_comp), arr.ind = TRUE)

for (k in seq_len(nrow(low_n_idx))) {
  i <- low_n_idx[k, "row"]
  j <- low_n_idx[k, "col"]
  r <- R_bio[i, j]
  
  if (is.na(r)) {
    ## Zero overlapping obs -> no information -> assume independence
    cor_comp_adj[i, j] <- 0
  } else {
    ## Some obs, but too few to trust magnitude -> keep sign, cap size
    cor_comp_adj[i, j] <- sign(r) * min(abs(r), shrink_cap)
  }
}

## Diagonal should stay exactly 1
diag(cor_comp_adj) <- 1

## Sanity checks
## sum(is.na(cor_comp_adj))          # should be 0 now
## range(cor_comp_adj[row(cor_comp_adj) != col(cor_comp_adj)])  # no more exact ±1 from low-n pairs


## Build Sigma from sd's and the (subset) correlation matrix
D <- diag(sd_bio)
Sigma_bio <- D %*% cor_comp_adj %*% D
dimnames(Sigma_bio) <- list(names(sd_bio), names(sd_bio))

## Confirm/repair positive semi-definiteness
ev <- eigen(Sigma_bio, symmetric = TRUE, only.values = TRUE)$values
min(ev)  # if negative, nearPD below fixes it
if (any(ev < -1e-8)) {
  Sigma_bio <- as.matrix(nearPD(Sigma_bio, corr = FALSE)$mat)
}




# Approach 1: Threshold based on actual threshold in real-world cases 
## All components: unhealthy = value exceeds cutoff (direction = +1 for all)
## ------------------------------------------------------------------
thresholds <- c(NUM_A1C = 6.5, NUM_ALB = 3.5, NUM_BMI = 30,
                NUM_CHOL = 200, NUM_CRP = 10, NUM_CREAT_C = 110,
                NUM_HCST = 50, NUM_TRIG = 150,
                NUM_BP_DIASTOLIC = 90, NUM_BP_SYSTOLIC = 140)
thresholds <- thresholds[names(mu_bio)]  # enforce order

## >= vs > per the table (mostly >=, a couple strict >)
use_strict <- c(NUM_A1C = FALSE, NUM_ALB = FALSE, NUM_BMI = TRUE,
                NUM_CHOL = FALSE, NUM_CRP = FALSE, NUM_CREAT_C = FALSE,
                NUM_HCST = TRUE, NUM_TRIG = FALSE,
                NUM_BP_DIASTOLIC = TRUE, NUM_BP_SYSTOLIC = TRUE)
use_strict <- use_strict[names(mu_bio)]


## Simulate Function 
sim_data_corr <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
  X_bio <- MASS::mvrnorm(n = N, mu = mu_bio, Sigma = Sigma_bio)
  X_bio <- X_bio[, names(mu_bio), drop = FALSE]

  S_mat <- matrix(NA_real_, nrow = N, ncol = ncol(X_bio),
                  dimnames = list(NULL, factor_ALI_comp))
  for (j in seq_along(thresholds)) {
    bio_col <- names(thresholds)[j]
    cutoff  <- thresholds[bio_col]
    S_mat[, j] <- if (use_strict[bio_col]) {
      as.numeric(X_bio[, bio_col] > cutoff)
    } else {
      as.numeric(X_bio[, bio_col] >= cutoff)
    }
  }

  mu_admit <- cbind(1, S_mat) %*% beta_admit
  Y_admit  <- rbinom(N, 1, plogis(mu_admit))

  mu_ed <- cbind(1, S_mat) %*% beta_ed
  Y_ed  <- rbinom(N, 1, plogis(mu_ed))

  mu_ed2 <- cbind(1, S_mat, Y_admit) %*% beta_ed2
  Y_ed2  <- rbinom(N, 1, plogis(mu_ed2))

  data.frame(cbind(S_mat, Y_admit, Y_ed, Y_ed2))
}


# Approach 2: use quantile? 
sim_data_corr <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
  X_bio <- MASS::mvrnorm(n = N, mu = mu_bio, Sigma = Sigma_bio)
  X_bio <- X_bio[, names(mu_bio), drop = FALSE]
  
  S_mat <- matrix(NA_real_, nrow = N, ncol = ncol(X_bio),
                  dimnames = list(NULL, factor_ALI_comp))
  for (j in seq_along(pS)) {
    bio_col <- names(mu_bio)[j]
    thresh  <- quantile(X_bio[, bio_col], probs = 1 - pS[j])
    S_mat[, j] <- as.numeric(X_bio[, bio_col] > thresh)
  }
  
  mu_admit <- cbind(1, S_mat) %*% beta_admit
  Y_admit  <- rbinom(N, 1, plogis(mu_admit))
  
  mu_ed <- cbind(1, S_mat) %*% beta_ed
  Y_ed  <- rbinom(N, 1, plogis(mu_ed))
  
  mu_ed2 <- cbind(1, S_mat, Y_admit) %*% beta_ed2
  Y_ed2  <- rbinom(N, 1, plogis(mu_ed2))
  
  data.frame(cbind(S_mat, Y_admit, Y_ed, Y_ed2))
}


# -------------------------------- Main Script --------------------------------
## --- Load function definitions ----------------------------------------------
source("R/00_setup.R")
source("R/03_predict_auc.R")
source("R/04_sim_helpers.R")
## Check what prevalence mvrnorrm actually produces (won't exactly match pS anymore, 
## now driven by real cutoffs + mu/Sigma, not forced to hit pS); compare to original pS as a sanity check
check_dat <- sim_data_corr(N = 1000, beta_admit = my_beta_admit,
                           beta_ed = my_beta_ed, beta_ed2 = my_beta_ed2)
round(colMeans(check_dat[, factor_ALI_comp]), 3)
## For comparison, vector of "unhealthy" probabilities were component 
### See huge mismatch in a few components -> problem with multivariate normal? 
pS                                                                                                                                                                                                                 

table(check_dat$Y_admit)
table(check_dat$Y_ed)
table(check_dat$Y_ed2)

mean(check_dat$Y_admit)
mean(check_dat$Y_ed)
mean(check_dat$Y_ed2)

## Run simulation 
run_sim_print_res(sim_data_corr,
                  N = 1000,
                  beta_admit = 3 * my_beta_admit,
                  beta_ed = 3 * my_beta_ed,
                  beta_ed2 = 3 * my_beta_ed2)
