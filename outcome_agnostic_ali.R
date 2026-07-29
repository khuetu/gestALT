# Outcome Agnostic ALI 
library("pROC")

#* -------------------------------- Simulate data ------------------------------
## Fixed constants
N <- 1000 ### sample size
## Vector of "unhealthy" probabilities were component (taken from deidentified_ali_hospitalizations.csv)
pS <- c(0.125, 0.799, 0.458, 0.239, 0.002,
        0.000, 0.000, 0.251, 0.064, 0.146)
names(pS) <- c("A1C_F", "ALB_F", "BMI_F", "CHOL_F", "CRP_F",
               "CREAT_C_F", "HCST_F", "TRIG_F", "BP_DIASTOLIC_F", "BP_SYSTOLIC_F")
pS <- pmax(pS, rep(0.01, times = length(pS))) ### force 1%+ prob per component
sim_data <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
  ## Begin with ALI components (10 per person) from Bernoulli (pS)
  S <- rbinom(n = N * 10,
              size = 1,
              prob = rep(x = pS, times = N))
  S_mat <- matrix(data = S,
                  nrow = N,
                  ncol = 10,
                  byrow = TRUE)
  colnames(S_mat) <- factor_ALI_comp
  ## Simulate outcomes
  ### Hospital admission Y1 | S1,...,S10
  mu_admit <- cbind(1, S_mat) %*% beta_admit
  Y_admit <- rbinom(n = N,
                    size = 1,
                    prob = (1 + exp(-mu_admit)) ^ (- 1))
  ### Emergency department visit Y2 | S1,...,S10
  mu_ed <- cbind(1, S_mat) %*% beta_ed
  Y_ed <- rbinom(n = N,
                 size = 1,
                 prob = (1 + exp(-mu_ed)) ^ (- 1))
  ### Alt: Emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  mu_ed2 <- cbind(1, S_mat, Y_admit) %*% beta_ed2
  Y_ed2 <- rbinom(n = N,
                  size = 1,
                  prob = (1 + exp(-mu_ed2)) ^ (- 1))
  ## Combine all columns and return as dataframe
  dat <- data.frame(cbind(S_mat, Y_admit, Y_ed, Y_ed2))
  return(dat)
}
## Define ALI components
bin_ALI_comp <- c("A1C", "ALB", "BMI", "CHOL", "CRP", "CREAT_C", "HCST", "TRIG", "BP_DIASTOLIC", "BP_SYSTOLIC")
factor_ALI_comp <- paste0(bin_ALI_comp, "_F")

## Specify model coefficients (for now, assumed to all be different)
### For hospital admission Y1 | S1,...,S10
my_beta_admit <- seq(0, 1, by = 0.1)
my_beta_admit <- matrix(data = my_beta_admit,
                        ncol = 1)
### For emergency department visit Y2 | S1,...,S10
my_beta_ed <- seq(0, 0.5, by = 0.05)
my_beta_ed <- matrix(data = my_beta_ed,
                     ncol = 1)
### For alt: emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
my_beta_ed2 <- seq(0, 0.55, by = 0.05)
my_beta_ed2 <- matrix(data = my_beta_ed2,
                      ncol = 1)

#* -----------------------------------------------------------------------------
#* Write a function to check that models' data generation is as we expect (i.e., average estimated coefficients = true coefficients)
check_models <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
  ### Simulate data
  dat = sim_data(N = N, #### sample size
                 beta_admit = beta_admit, #### coefficients for hospital admission Y1 | S1,...,S10
                 beta_ed = beta_ed, #### coefficients for emergency department visit Y2 | S1,...,S10
                 beta_ed2 = beta_ed2) #### coefficients for alt: emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  
  ### Fit models
  #### Hospital admissions
  fit_admit = glm(formula = as.formula(paste0("Y_admit ~", paste0(factor_ALI_comp, collapse = "+"))),
                  data = dat,
                  family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed = glm(formula = as.formula(paste0("Y_ed ~", paste0(factor_ALI_comp, collapse = "+"))),
               data = dat,
               family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2 = glm(formula = as.formula(paste0("Y_ed2 ~", paste0(factor_ALI_comp, collapse = "+"), " + Y_admit")),
                data = dat,
                family = "binomial")
  
  ### Return coefficients
  return(list(
    mod_admit = fit_admit$coefficients,
    mod_ed = fit_ed$coefficients,
    mod_ed2 = fit_ed2$coefficients)
  )
}

#* -----------------------------------------------------------------------------
## Replicate data generation 1000 times, fit all three models, and compare average coefficients to truth
set.seed(918)
reps <- 1000
mod_admit <- mod_ed <- matrix(data = NA, nrow = reps, ncol = 11)
mod_ed2 <- matrix(data = NA, nrow = reps, ncol = 12)
for (r in 1:reps) {
  mods <- check_models(N = 1000,
                       beta_admit = my_beta_admit,
                       beta_ed = my_beta_ed,
                       beta_ed2 = my_beta_ed2)
  mod_admit[r, ] <- mods$mod_admit
  mod_ed[r, ] <- mods$mod_ed
  mod_ed2[r, ] <- mods$mod_ed2
}
colMeans(mod_admit, na.rm = TRUE)
colMeans(mod_ed, na.rm = TRUE)
colMeans(mod_ed2, na.rm = TRUE)

#* -----------------------------------------------------------------------------
## Try predicting Y1 and Y2 from correctly-specified models using full sample and return separate models' AUC
try_pred_full_sample <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
  ### Simulate data
  dat <- sim_data(N = N, #### sample size
                  beta_admit = beta_admit, #### coefficients for hospital admission Y1 | S1,...,S10
                  beta_ed = beta_ed, #### coefficients for emergency department visit Y2 | S1,...,S10
                  beta_ed2 = beta_ed2) #### coefficients for alt: emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  
  ### Fit models
  #### Hospital admissions
  fit_admit <- glm(formula = as.formula(paste0("Y_admit ~", paste0(factor_ALI_comp, collapse = "+"))),
                   data = dat,
                   family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed <- glm(formula = as.formula(paste0("Y_ed ~", paste0(factor_ALI_comp, collapse = "+"))),
                data = dat,
                family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2 <- glm(formula = as.formula(paste0("Y_ed2 ~", paste0(factor_ALI_comp, collapse = "+"), " + Y_admit")),
                 data = dat,
                 family = "binomial")
  
  ### Predict from models
  #### Hospital admissions
  pred_admit <- predict(fit_admit, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed <- predict(fit_ed, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2 = predict(fit_ed2, type = "response")
  
  ### Calculate AUC from models
  #### Hospital admissions
  auc_admit <- pROC::auc(pROC::roc(dat$Y_admit, pred_admit))
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  auc_ed <- pROC::auc(pROC::roc(dat$Y_ed, pred_ed))
  #### Emergency department visits (dependent on ANY_ADMIT)
  auc_ed2 = pROC::auc(pROC::roc(dat$Y_ed2, pred_ed2))
  
  ### Return coefficients
  return(
    data.frame(
      auc_admit = auc_admit,
      auc_ed = auc_ed,
      auc_ed2 = auc_ed2
    )
  )
}

try_pred_full_sample(
  N = 1000,
  beta_admit = my_beta_admit,
  beta_ed = my_beta_ed,
  beta_ed2 = my_beta_ed2
)

## Doubled all the coefficient vectors to make AUC bigger (0.6ish --> 0.7ish)
try_pred_full_sample(
  N = 1000,
  beta_admit = 2 * my_beta_admit,
  beta_ed = 2 * my_beta_ed,
  beta_ed2 = 2 * my_beta_ed2
)

## Tripled all the coefficient vectors to make AUC bigger (0.6ish --> 0.8ish)
try_pred_full_sample(
  N = 1000,
  beta_admit = 3 * my_beta_admit,
  beta_ed = 3 * my_beta_ed,
  beta_ed2 = 3 * my_beta_ed2
)

#* -----------------------------------------------------------------------------
## Try predicting Y1 and Y2 from correctly-specified models using full sample and return separate models' AUC
try_pred_new_ali <- function(sim_func, N = 1000, beta_admit, beta_ed, beta_ed2) {
  ### Simulate data
  dat <- sim_func(N = N, #### sample size
                  beta_admit = beta_admit, #### coefficients for hospital admission Y1 | S1,...,S10
                  beta_ed = beta_ed, #### coefficients for emergency department visit Y2 | S1,...,S10
                  beta_ed2 = beta_ed2) #### coefficients for alt: emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  
  ### Fit models using Y1, Y2 and separate ALI components
  #### Hospital admissions
  fit_admit <- glm(formula = as.formula(paste0("Y_admit ~", paste0(factor_ALI_comp, collapse = "+"))),
                   data = dat,
                   family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed <- glm(formula = as.formula(paste0("Y_ed ~", paste0(factor_ALI_comp, collapse = "+"))),
                data = dat,
                family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2 <- glm(formula = as.formula(paste0("Y_ed2 ~", paste0(factor_ALI_comp, collapse = "+"), " + Y_admit")),
                 data = dat,
                 family = "binomial")
  
  ### Predict from models using Y1, Y2 and separate ALI components
  #### Hospital admissions
  pred_admit <- predict(fit_admit, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed <- predict(fit_ed, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2 = predict(fit_ed2, type = "response")
  
  ### Calculate AUC from models
  #### Hospital admissions
  suppressMessages(
    auc_admit <- pROC::auc(pROC::roc(dat$Y_admit, pred_admit))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed <- pROC::auc(pROC::roc(dat$Y_ed, pred_ed))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2 <- pROC::auc(pROC::roc(dat$Y_ed2, pred_ed2))
  )
  
  ### Aggregate (faux-"ensemble") the predicted probabilities from models using Y1, Y2 and separate ALI components
  #### Sarah acknowledges that this is ad-hoc and could formally be derived using probability rules
  #### Aggregate by averaging over predicted probabilities
  dat$avg_pred_admit_ed <- (pred_admit + pred_ed) / 2
  dat$avg_pred_admit_ed2 <- (pred_admit + pred_ed2) / 2
  #### Aggregate by summing over predicted probabilities
  dat$sum_pred_admit_ed <- (pred_admit + pred_ed)
  dat$sum_pred_admit_ed2 <- (pred_admit + pred_ed2)
  
  ### Fit models using Y1, Y2 and new ALI instead of separate components
  #### Hospital admissions
  fit_admit_new <- glm(formula = Y_admit ~ avg_pred_admit_ed,
                       data = dat,
                       family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed_new <- glm(formula = Y_ed ~ avg_pred_admit_ed,
                    data = dat,
                    family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2_new <- glm(formula = Y_ed2 ~ avg_pred_admit_ed2 + Y_admit,
                     data = dat,
                     family = "binomial")
  
  ### Predict from models using Y1, Y2 and new ALI instead of separate components
  #### Hospital admissions
  pred_admit_new <- predict(fit_admit_new, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed_new <- predict(fit_ed_new, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2_new = predict(fit_ed2_new, type = "response")
  
  ### Calculate AUC from new ALI models
  #### Hospital admissions
  suppressMessages(
    auc_admit_new <- pROC::auc(pROC::roc(dat$Y_admit, pred_admit_new))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed_new <- pROC::auc(pROC::roc(dat$Y_ed, pred_ed_new))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2_new <- pROC::auc(pROC::roc(dat$Y_ed2, pred_ed2_new))
  )
  
  
  ### Fit models using Y1, Y2 and old ALI instead of separate components
  #### Calculate "old" ALI (equal weights for all components, all models)
  dat$old_ali <- rowSums(dat[, 1:10]) / 10
  
  #### Hospital admissions
  fit_admit_old <- glm(formula = Y_admit ~ old_ali,
                       data = dat,
                       family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed_old <- glm(formula = Y_ed ~ old_ali,
                    data = dat,
                    family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2_old <- glm(formula = Y_ed2 ~ old_ali + Y_admit,
                     data = dat,
                     family = "binomial")
  
  ### Predict from models using Y1, Y2 and old ALI instead of separate components
  #### Hospital admissions
  pred_admit_old <- predict(fit_admit_old, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed_old <- predict(fit_ed_old, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2_old = predict(fit_ed2_old, type = "response")
  
  ### Calculate AUC from models
  #### Hospital admissions
  suppressMessages(
    auc_admit_old <- pROC::auc(pROC::roc(dat$Y_admit, pred_admit_old))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed_old <- pROC::auc(pROC::roc(dat$Y_ed, pred_ed_old))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2_old <- pROC::auc(pROC::roc(dat$Y_ed2, pred_ed2_old))
  )
  
  ### Return coefficients
  return(
    list(
      sep_comp = data.frame(
        auc_admit = auc_admit,
        auc_ed = auc_ed,
        auc_ed2 = auc_ed2
      ),
      new_ali = data.frame(
        auc_admit = auc_admit_new,
        auc_ed = auc_ed_new,
        auc_ed2 = auc_ed2_new
      ),
      old_ali = data.frame(
        auc_admit = auc_admit_old,
        auc_ed = auc_ed_old,
        auc_ed2 = auc_ed2_old
      )
    )
  )
}

try_pred_new_ali(
  sim_data, 
  N = 1000,
  beta_admit = 3 * my_beta_admit,
  beta_ed = 3 * my_beta_ed,
  beta_ed2 = 3 * my_beta_ed2
)

#* -------------------------------- Simulation #* ------------------------------
## Simulation: run try_pred_new_ali() multiple times and compute average AUC
set.seed(22)
n_sims <- 1000

results <- lapply(1:n_sims, function(i) {
  suppressWarnings(suppressMessages(
    try_pred_new_ali(sim_data, 
                     N = 1000,
                     beta_admit = 3 * my_beta_admit,
                     beta_ed = 3 * my_beta_ed,
                     beta_ed2 = 3 * my_beta_ed2)
  ))
})

## Average AUC across reps for a given model type
avg_auc <- function(model_name) {
  mat <- do.call(rbind, lapply(results, function(r) as.numeric(r[[model_name]])))
  colnames(mat) <- names(results[[1]][[model_name]])
  colMeans(mat)
}

cat("Mean AUC - Separate components:\n"); print(avg_auc("sep_comp"))
cat("\nMean AUC - New ALI:\n");           print(avg_auc("new_ali"))
cat("\nMean AUC - Old ALI:\n");           print(avg_auc("old_ali"))


#* -----------------------------------------------------------------------------
## Try predicting Y1 and Y2 from correctly-specified models using full sample and return separate models' AUC
### Fit models using training data, evaluate using testing data (50/50 test/train split)
try_pred_new_ali_cv <- function(sim_func, N = 1000, train_prop = 0.5, beta_admit, beta_ed, beta_ed2) {
  ### Simulate data
  dat <- sim_func(N = N, #### sample size
                  beta_admit = beta_admit, #### coefficients for hospital admission Y1 | S1,...,S10
                  beta_ed = beta_ed, #### coefficients for emergency department visit Y2 | S1,...,S10
                  beta_ed2 = beta_ed2) #### coefficients for alt: emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  
  ### Split into 50/50 test/train <- New! As of 06/23/2026
  which_train <- sample(x = 1:N, size = N * train_prop, replace = FALSE) #### sample rows for training data
  train_dat <- dat[which_train, ] #### separate training data
  test_dat <- dat[-which_train, ] #### separate testing data
  
  ### Fit models using Y1, Y2 and separate ALI components
  #### Hospital admissions
  fit_admit <- glm(formula = as.formula(paste0("Y_admit ~", paste0(factor_ALI_comp, collapse = "+"))),
                   data = train_dat,
                   family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed <- glm(formula = as.formula(paste0("Y_ed ~", paste0(factor_ALI_comp, collapse = "+"))),
                data = train_dat,
                family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2 <- glm(formula = as.formula(paste0("Y_ed2 ~", paste0(factor_ALI_comp, collapse = "+"), " + Y_admit")),
                 data = train_dat,
                 family = "binomial")
  
  ### Predict from models using Y1, Y2 and separate ALI components
  #### Hospital admissions
  pred_admit <- predict(fit_admit, newdata = test_dat, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed <- predict(fit_ed, newdata = test_dat, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2 = predict(fit_ed2, newdata = test_dat, type = "response")
  
  ### Calculate AUC from models
  #### Hospital admissions
  suppressMessages(
    auc_admit <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed <- pROC::auc(pROC::roc(test_dat$Y_ed, pred_ed))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2 <- pROC::auc(pROC::roc(test_dat$Y_ed2, pred_ed2))
  )
  
  ### Predict from models using outcome-specific ALI (= pred_admit, pred_ed, pred_ed2)
  #### << Copy and paste from here to next marker and repeat for ed and ed 2 >>
  #### Fit models using Y1, Y2 and admission-specific ALI instead of separate components
  #### Hospital admissions
  fit_admit_admit <- glm(formula = Y_admit ~ pred_admit,
                         data = test_dat,
                         family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed_admit <- glm(formula = Y_ed ~ pred_admit,
                      data = test_dat,
                      family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2_admit <- glm(formula = Y_ed2 ~ pred_admit + Y_admit,
                       data = test_dat,
                       family = "binomial")
  
  #### Predict from models using Y1, Y2 and admission-specific ALI instead of separate components
  #### Hospital admissions
  pred_admit_admit <- predict(fit_admit_admit, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed_admit <- predict(fit_ed_admit, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2_admit = predict(fit_ed2_admit, type = "response")
  
  #### Calculate AUC from admission-specific ALI models
  #### Hospital admissions
  suppressMessages(
    auc_admit_admit <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_admit))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed_admit <- pROC::auc(pROC::roc(test_dat$Y_ed, pred_ed_admit))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2_admit <- pROC::auc(pROC::roc(test_dat$Y_ed2, pred_ed2_admit))
  )
  
  #### Fit models using Y1, Y2 and admission-specific ALI instead of separate components
  #### Emergency Department
  fit_admit_ed <- glm(formula = Y_admit ~ pred_ed,
                         data = test_dat,
                         family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed_ed <- glm(formula = Y_ed ~ pred_ed,
                      data = test_dat,
                      family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2_ed <- glm(formula = Y_ed2 ~ pred_ed + Y_admit,
                       data = test_dat,
                       family = "binomial")
  
  #### Predict from models using Y1, Y2 and admission-specific ALI instead of separate components
  #### Emeregency Dept
  pred_admit_ed <- predict(fit_admit_ed, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed_ed <- predict(fit_ed_ed, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2_ed = predict(fit_ed2_ed, type = "response")
  
  #### Calculate AUC from admission-specific ALI models
  #### Hospital admissions
  suppressMessages(
    auc_admit_ed <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_ed))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed_ed <- pROC::auc(pROC::roc(test_dat$Y_ed, pred_ed_ed))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2_ed <- pROC::auc(pROC::roc(test_dat$Y_ed2, pred_ed2_ed))
  )
  
  #### Fit models using Y1, Y2 and admission-specific ALI instead of separate components
  #### Emergency Department2
  #### NOTE: pred_ed2 is predicted from the mode: fit_ed2, where its formula is: 
  #### "Y_ed2 ~", paste0(factor_ALI_comp, collapse = "+"), " + Y_admit"
  fit_admit_ed2 <- glm(formula = Y_admit ~ pred_ed2,
                      data = test_dat,
                      family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed_ed2 <- glm(formula = Y_ed ~ pred_ed2,
                   data = test_dat,
                   family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2_ed2 <- glm(formula = Y_ed2 ~ pred_ed2 + Y_admit,
                    data = test_dat,
                    family = "binomial")
  
  #### Predict from models using Y1, Y2 and admission-specific ALI instead of separate components
  #### Hospital admissions
  pred_admit_ed2 <- predict(fit_admit_ed2, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed_ed2 <- predict(fit_ed_ed2, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2_ed2 <- predict(fit_ed2_ed2, type = "response")
  
  #### Calculate AUC from admission-specific ALI models
  #### Hospital admissions
  suppressMessages(
    auc_admit_ed2 <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_ed2))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed_ed2 <- pROC::auc(pROC::roc(test_dat$Y_ed, pred_ed_ed2))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2_ed2 <- pROC::auc(pROC::roc(test_dat$Y_ed2, pred_ed2_ed2))
  )
  
  #------------------------------------------------------------------
  ### Aggregate (faux-"ensemble") the predicted probabilities from models using Y1, Y2 and separate ALI components
  #### Sarah acknowledges that this is ad-hoc and could formally be derived using probability rules
  #### Aggregate by averaging over predicted probabilities
  test_dat$avg_pred_admit_ed <- (pred_admit + pred_ed) / 2
  test_dat$avg_pred_admit_ed2 <- (pred_admit + pred_ed2) / 2
  #### Aggregate by summing over predicted probabilities
  test_dat$sum_pred_admit_ed <- (pred_admit + pred_ed)
  test_dat$sum_pred_admit_ed2 <- (pred_admit + pred_ed2)
  
  ### Fit models using Y1, Y2 and new ALI instead of separate components
  #### Hospital admissions
  fit_admit_new <- glm(formula = Y_admit ~ avg_pred_admit_ed,
                       data = test_dat,
                       family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed_new <- glm(formula = Y_ed ~ avg_pred_admit_ed,
                    data = test_dat,
                    family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2_new <- glm(formula = Y_ed2 ~ avg_pred_admit_ed2 + Y_admit,
                     data = test_dat,
                     family = "binomial")
  
  ### Predict from models using Y1, Y2 and new ALI instead of separate components
  #### Hospital admissions
  pred_admit_new <- predict(fit_admit_new, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed_new <- predict(fit_ed_new, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2_new = predict(fit_ed2_new, type = "response")
  
  ### Calculate AUC from new ALI models
  #### Hospital admissions
  suppressMessages(
    auc_admit_new <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_new))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed_new <- pROC::auc(pROC::roc(test_dat$Y_ed, pred_ed_new))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2_new <- pROC::auc(pROC::roc(test_dat$Y_ed2, pred_ed2_new))
  )
  
  ### Fit models using Y1, Y2 and old ALI instead of separate components
  #### Calculate "old" ALI (equal weights for all components, all models)
  test_dat$old_ali <- rowSums(test_dat[, 1:10]) / 10
  
  #### Hospital admissions
  fit_admit_old <- glm(formula = Y_admit ~ old_ali,
                       data = test_dat,
                       family = "binomial")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  fit_ed_old <- glm(formula = Y_ed ~ old_ali,
                    data = test_dat,
                    family = "binomial")
  #### Emergency department visits (dependent on ANY_ADMIT)
  fit_ed2_old <- glm(formula = Y_ed2 ~ old_ali + Y_admit,
                     data = test_dat,
                     family = "binomial")
  
  ### Predict from models using Y1, Y2 and old ALI instead of separate components
  #### Hospital admissions
  pred_admit_old <- predict(fit_admit_old, type = "response")
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  pred_ed_old <- predict(fit_ed_old, type = "response")
  #### Emergency department visits (dependent on ANY_ADMIT)
  pred_ed2_old = predict(fit_ed2_old, type = "response")
  
  ### Calculate AUC from models
  #### Hospital admissions
  suppressMessages(
    auc_admit_old <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_old))
  )
  #### Emergency department visits (conditionally independent of ANY_ADMIT)
  suppressMessages(
    auc_ed_old <- pROC::auc(pROC::roc(test_dat$Y_ed, pred_ed_old))
  )
  #### Emergency department visits (dependent on ANY_ADMIT)
  suppressMessages(
    auc_ed2_old <- pROC::auc(pROC::roc(test_dat$Y_ed2, pred_ed2_old))
  )
  
  ### Return coefficients
  return(
    list(
      sep_comp = data.frame(
        auc_admit = auc_admit,
        auc_ed = auc_ed,
        auc_ed2 = auc_ed2
      ),
      new_ali = data.frame(
        auc_admit = auc_admit_new,
        auc_ed = auc_ed_new,
        auc_ed2 = auc_ed2_new
      ),
      old_ali = data.frame(
        auc_admit = auc_admit_old,
        auc_ed = auc_ed_old,
        auc_ed2 = auc_ed2_old
      ),
      admit_ali = data.frame(
        auc_admit = auc_admit_admit,
        auc_ed = auc_ed_admit,
        auc_ed2 = auc_ed2_admit
      ),
      ed_ali = data.frame(
        auc_admit = auc_admit_ed,
        auc_ed = auc_ed_ed,
        auc_ed2 = auc_ed2_ed
      ),
      ed2_ali = data.frame(
        auc_admit = auc_admit_ed2,
        auc_ed = auc_ed_ed2,
        auc_ed2 = auc_ed2_ed2
      )
    )
  )
}

#* --------------------------------- Simulation -------------------------------- 
## Simulation: run try_pred_new_ali() multiple times and compute average AUC
set.seed(22)
n_sims <- 1000

## Helper function: Average AUC across reps for a given model type
avg_auc <- function(model_name, results) {
  mat <- do.call(rbind, lapply(results, function(r) as.numeric(r[[model_name]])))
  colnames(mat) <- names(results[[1]][[model_name]])
  colMeans(mat)
}

run_sim_print_res <- function(sim_func, N, beta_admit, beta_ed, beta_ed2) {
  results <- lapply(1:n_sims, function(i) {
    suppressWarnings(suppressMessages(
      try_pred_new_ali_cv(sim_func, 
                          N = 1000,
                          beta_admit = beta_admit,
                          beta_ed = beta_ed,
                          beta_ed2 = beta_ed2)
    ))
  })
  
  cat("Mean AUC - Separate components:\n"); print(avg_auc("sep_comp", results))
  cat("\nMean AUC - New ALI:\n");           print(avg_auc("new_ali", results))
  cat("\nMean AUC - Old ALI:\n");           print(avg_auc("old_ali", results))
  cat("\nMean AUC - Admission-specific: Admit:\n"); print(avg_auc("admit_ali", results))
  cat("\nMean AUC - Admission-specific: Ed:\n"); print(avg_auc("ed_ali", results))
  cat("\nMean AUC - Admission-specific: Ed2:\n"); print(avg_auc("ed2_ali", results)) 
}

###Note: Average auc_admit under ed2_ali is surprisingly high!? 
run_sim_print_res(sim_data, 
                  N = 1000,
                  beta_admit = 3 * my_beta_admit,
                  beta_ed = 3 * my_beta_ed,
                  beta_ed2 = 3 * my_beta_ed2)

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
# Create new sim data function to include age and sex 
sim_data_agesex <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
  ## Begin with ALI components (10 per person) from Bernoulli (pS)
  S <- rbinom(n = N * 10,
              size = 1,
              prob = rep(x = pS, times = N))
  S_mat <- matrix(data = S,
                  nrow = N,
                  ncol = 10,
                  byrow = TRUE)
  colnames(S_mat) <- factor_ALI_comp
  ## Simulate outcomes
  ### Age: mean aged 48, range 18-65
  age <- rnorm(N, mean = 48, sd = 15)
  age <- pmin(pmax(age, 18), 65)
  
  ### Sex: binary (0 = male, 1 = female)
  sex <- rbinom(N, size = 1, prob = 0.6)
  
  ### Hospital admission Y1 | S1,...,S10
  mu_admit <- cbind(1, S_mat, age, sex) %*% beta_admit
  Y_admit <- rbinom(n = N,
                    size = 1,
                    prob = (1 + exp(-mu_admit)) ^ (- 1))
  ### Emergency department visit Y2 | S1,...,S10
  mu_ed <- cbind(1, S_mat, age, sex) %*% beta_ed
  Y_ed <- rbinom(n = N,
                 size = 1,
                 prob = (1 + exp(-mu_ed)) ^ (- 1))
  ### Alt: Emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  mu_ed2 <- cbind(1, S_mat, Y_admit, age, sex) %*% beta_ed2
  Y_ed2 <- rbinom(n = N,
                  size = 1,
                  prob = (1 + exp(-mu_ed2)) ^ (- 1))
  ## Combine all columns and return as dataframe
  dat <- data.frame(cbind(S_mat, age = age, sex = sex, Y_admit, Y_ed, Y_ed2))
  return(dat)
}

# beta coefficients with age and sex 
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
                  beta_ed = or1_beta_ed,
                  beta_ed2 = or1_beta_ed2)

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
comp_ed_2    <- c(1, 2, 3, 6, 5, 4, 7, 8, 11, 10, 9) # swap ranks 4/6 and last three

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


#*------------------------------- Possible Next Steps? -------------------------
#* pS equal for all components 
#* Real biomarkers are correlated (A1C and BMI, blood pressure components, etc.). Adding a correlation structure
#* Opt1: Multinormmal and pick correlated structure -> discretize instead of binary (bernoulli components)
#* mvrnorm {MASS}: mu = complete case mean, cov matrix, na.remove = TRUE, correlation component are missing 
#* Opt2: Conditional bernoulli on one another (component wise)
#* Aggregate probabilities: simplex? but minimize what? or constrained MLE? 
#* Run simulations where some components truly don't matter (zero true coefficient) and verify the new ALI down-weights them appropriately
#* Multinomial 
#* Down the road: missingness: pattern submodel 

                               