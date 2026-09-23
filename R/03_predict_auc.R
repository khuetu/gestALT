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

## Try predicting Y1 and Y2 from correctly-specified models using full sample and return aggregated AUC
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

## Try predicting Y1 and Y2 from correctly-specified models using training sample and return aggregated AUC
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

try_pred_new_ali_3folds <- function(sim_func, N = 1000, 
                                train_prop = 0.5, 
                                valid_prop = 0.25, # remaining 0.25 goes to test
                                beta_admit, beta_ed, beta_ed2) {
  ### Simulate data
  dat <- sim_func(N = N,
                  beta_admit = beta_admit,
                  beta_ed = beta_ed,
                  beta_ed2 = beta_ed2)
  
  ### Three-way split: train / valid / test
  n <- nrow(dat)
  which_train <- sample(x = 1:n, size = round(n * train_prop), replace = FALSE)
  remaining   <- setdiff(1:n, which_train)
  which_valid <- sample(x = remaining, size = round(n * valid_prop), replace = FALSE)
  which_test  <- setdiff(remaining, which_valid)
  
  train_dat <- dat[which_train, ]
  valid_dat <- dat[which_valid, ]
  test_dat  <- dat[which_test, ]
  
  ### STAGE 1: Fit first-stage models (separate ALI components) on TRAIN data  --------------------------------
  fit_admit <- glm(formula = as.formula(paste0("Y_admit ~", paste0(factor_ALI_comp, collapse = "+"))),
                   data = train_dat,
                   family = "binomial")
  fit_ed <- glm(formula = as.formula(paste0("Y_ed ~", paste0(factor_ALI_comp, collapse = "+"))),
                data = train_dat,
                family = "binomial")
  fit_ed2 <- glm(formula = as.formula(paste0("Y_ed2 ~", paste0(factor_ALI_comp, collapse = "+"), " + Y_admit")),
                 data = train_dat,
                 family = "binomial")
  
  #### Predict on new "VALIDATE dataset" ----------------------------------------------------------------------
  valid_dat$pred_admit <- predict(fit_admit, newdata = valid_dat, type = "response")
  valid_dat$pred_ed    <- predict(fit_ed,    newdata = valid_dat, type = "response")
  valid_dat$pred_ed2   <- predict(fit_ed2,   newdata = valid_dat, type = "response")
  
  #### First-stage predictions on TEST data (for final evaluation) --------------------------------------------
  test_dat$pred_admit <- predict(fit_admit, newdata = test_dat, type = "response")
  test_dat$pred_ed    <- predict(fit_ed,    newdata = test_dat, type = "response")
  test_dat$pred_ed2   <- predict(fit_ed2,   newdata = test_dat, type = "response")
  
  #### AUC for separate components (first-stage), evaluated on TEST
  suppressMessages(auc_admit <- pROC::auc(pROC::roc(test_dat$Y_admit, test_dat$pred_admit)))
  suppressMessages(auc_ed    <- pROC::auc(pROC::roc(test_dat$Y_ed,    test_dat$pred_ed)))
  suppressMessages(auc_ed2   <- pROC::auc(pROC::roc(test_dat$Y_ed2,   test_dat$pred_ed2)))
  
  ### STAGE 2: Fit second-stage models on VALID, evaluate on TEST --------------------------------------------
  #### ---- admit-specific ALI ----
  fit_admit_admit <- glm(Y_admit ~ pred_admit, data = valid_dat, family = "binomial")
  fit_ed_admit    <- glm(Y_ed    ~ pred_admit, data = valid_dat, family = "binomial")
  fit_ed2_admit   <- glm(Y_ed2   ~ pred_admit + Y_admit, data = valid_dat, family = "binomial")
  
  pred_admit_admit <- predict(fit_admit_admit, newdata = test_dat, type = "response")
  pred_ed_admit    <- predict(fit_ed_admit,    newdata = test_dat, type = "response")
  pred_ed2_admit   <- predict(fit_ed2_admit,   newdata = test_dat, type = "response")
  
  suppressMessages(auc_admit_admit <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_admit)))
  suppressMessages(auc_ed_admit    <- pROC::auc(pROC::roc(test_dat$Y_ed,    pred_ed_admit)))
  suppressMessages(auc_ed2_admit   <- pROC::auc(pROC::roc(test_dat$Y_ed2,   pred_ed2_admit)))
  
  #### ---- ed-specific ALI ----
  fit_admit_ed <- glm(Y_admit ~ pred_ed, data = valid_dat, family = "binomial")
  fit_ed_ed    <- glm(Y_ed    ~ pred_ed, data = valid_dat, family = "binomial")
  fit_ed2_ed   <- glm(Y_ed2   ~ pred_ed + Y_admit, data = valid_dat, family = "binomial")
  
  pred_admit_ed <- predict(fit_admit_ed, newdata = test_dat, type = "response")
  pred_ed_ed    <- predict(fit_ed_ed,    newdata = test_dat, type = "response")
  pred_ed2_ed   <- predict(fit_ed2_ed,   newdata = test_dat, type = "response")
  
  suppressMessages(auc_admit_ed <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_ed)))
  suppressMessages(auc_ed_ed    <- pROC::auc(pROC::roc(test_dat$Y_ed,    pred_ed_ed)))
  suppressMessages(auc_ed2_ed   <- pROC::auc(pROC::roc(test_dat$Y_ed2,   pred_ed2_ed)))
  
  #### ---- ed2-specific ALI ----
  fit_admit_ed2 <- glm(Y_admit ~ pred_ed2, data = valid_dat, family = "binomial")
  fit_ed_ed2    <- glm(Y_ed    ~ pred_ed2, data = valid_dat, family = "binomial")
  fit_ed2_ed2   <- glm(Y_ed2   ~ pred_ed2 + Y_admit, data = valid_dat, family = "binomial")
  
  pred_admit_ed2 <- predict(fit_admit_ed2, newdata = test_dat, type = "response")
  pred_ed_ed2    <- predict(fit_ed_ed2,    newdata = test_dat, type = "response")
  pred_ed2_ed2   <- predict(fit_ed2_ed2,   newdata = test_dat, type = "response")
  
  suppressMessages(auc_admit_ed2 <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_ed2)))
  suppressMessages(auc_ed_ed2    <- pROC::auc(pROC::roc(test_dat$Y_ed,    pred_ed_ed2)))
  suppressMessages(auc_ed2_ed2   <- pROC::auc(pROC::roc(test_dat$Y_ed2,   pred_ed2_ed2)))
  
  #### "New" ALI: use the pred fitted on VALID from the outcome-specìic section above, eval on TEST ----------
  valid_dat$avg_pred_admit_ed  <- (valid_dat$pred_admit + valid_dat$pred_ed)  / 2
  valid_dat$avg_pred_admit_ed2 <- (valid_dat$pred_admit + valid_dat$pred_ed2) / 2
  
  test_dat$avg_pred_admit_ed  <- (test_dat$pred_admit + test_dat$pred_ed)  / 2
  test_dat$avg_pred_admit_ed2 <- (test_dat$pred_admit + test_dat$pred_ed2) / 2
  
  fit_admit_new <- glm(Y_admit ~ avg_pred_admit_ed, data = valid_dat, family = "binomial")
  fit_ed_new    <- glm(Y_ed    ~ avg_pred_admit_ed, data = valid_dat, family = "binomial")
  fit_ed2_new   <- glm(Y_ed2   ~ avg_pred_admit_ed2 + Y_admit, data = valid_dat, family = "binomial")
  
  pred_admit_new <- predict(fit_admit_new, newdata = test_dat, type = "response")
  pred_ed_new    <- predict(fit_ed_new,    newdata = test_dat, type = "response")
  pred_ed2_new   <- predict(fit_ed2_new,   newdata = test_dat, type = "response")
  
  suppressMessages(auc_admit_new <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_new)))
  suppressMessages(auc_ed_new    <- pROC::auc(pROC::roc(test_dat$Y_ed,    pred_ed_new)))
  suppressMessages(auc_ed2_new   <- pROC::auc(pROC::roc(test_dat$Y_ed2,   pred_ed2_new)))
  
  # equal ALI: fit on VALID, eval on TEST -----------------------------------------------------------------
  valid_dat$old_ali <- rowSums(valid_dat[, 1:10]) / 10
  test_dat$old_ali  <- rowSums(test_dat[, 1:10]) / 10 

  fit_admit_old <- glm(Y_admit ~ old_ali, data = valid_dat, family = "binomial")
  fit_ed_old    <- glm(Y_ed    ~ old_ali, data = valid_dat, family = "binomial")
  fit_ed2_old   <- glm(Y_ed2   ~ old_ali + Y_admit, data = valid_dat, family = "binomial")
  
  pred_admit_old <- predict(fit_admit_old, newdata = test_dat, type = "response")
  pred_ed_old    <- predict(fit_ed_old,    newdata = test_dat, type = "response")
  pred_ed2_old   <- predict(fit_ed2_old,   newdata = test_dat, type = "response")
  
  suppressMessages(auc_admit_old <- pROC::auc(pROC::roc(test_dat$Y_admit, pred_admit_old)))
  suppressMessages(auc_ed_old    <- pROC::auc(pROC::roc(test_dat$Y_ed,    pred_ed_old)))
  suppressMessages(auc_ed2_old   <- pROC::auc(pROC::roc(test_dat$Y_ed2,   pred_ed2_old)))
  
  return(
    list(
      sep_comp  = data.frame(auc_admit = auc_admit,       auc_ed = auc_ed,       auc_ed2 = auc_ed2),
      new_ali   = data.frame(auc_admit = auc_admit_new,   auc_ed = auc_ed_new,   auc_ed2 = auc_ed2_new),
      old_ali   = data.frame(auc_admit = auc_admit_old,   auc_ed = auc_ed_old,   auc_ed2 = auc_ed2_old),
      admit_ali = data.frame(auc_admit = auc_admit_admit, auc_ed = auc_ed_admit, auc_ed2 = auc_ed2_admit),
      ed_ali    = data.frame(auc_admit = auc_admit_ed,    auc_ed = auc_ed_ed,    auc_ed2 = auc_ed2_ed),
      ed2_ali   = data.frame(auc_admit = auc_admit_ed2,   auc_ed = auc_ed_ed2,   auc_ed2 = auc_ed2_ed2)
    )
  )
}