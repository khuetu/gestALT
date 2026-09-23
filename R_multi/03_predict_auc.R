## PREDICT AUC FOR DIFFERFENT ALIs USING MULTIPLE OUTCOMES (K) -----------------
try_pred_new_ali_3folds_multi <- function(sim_data_multi, N = 1000, K = 5, 
                                          train_prop = 0.5, 
                                          valid_prop = 0.25, # remaining 0.25 goes to test
                                          my_betas) {
  
  ### Simulate data
  dat <- sim_data_multi(N, K, my_betas) 
  
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
  K_outcomes <- outcomes_reference[1:K, ]
  
  clean_name <- gsub(" ", "", tolower(K_outcomes$outcomes))
  dynamic_names_auc <- paste0("auc_", clean_name)
  
  auc_outcomes <- data.frame(matrix(NA, nrow = K + 3, ncol = K))
  colnames(auc_outcomes) <- dynamic_names_auc
  rownames(auc_outcomes) <- c("sep_comp", "new_ali", "old_ali", dynamic_names_auc)
  
  #### Predict on new "VALIDATE dataset" ----------------------------------------------------------------------
  for (outcome in K_outcomes$outcomes) {
    index <- match(outcome, K_outcomes$outcomes)
    clean_outcome <- gsub(" ", "", tolower(outcome))
    
    ##### 
    fit_temp <- glm(formula = as.formula(paste0(paste0("Y_", clean_outcome, "~"), paste0(factor_ALI_comp, collapse = "+"))),
                    data = train_dat,
                    family = "binomial")
    ##### 
    valid_dat[[paste0("pred_", clean_outcome)]] <- predict(fit_temp, newdata = valid_dat, type = "response")
    test_dat[[paste0("pred_", clean_outcome)]] <- predict(fit_temp, newdata = test_dat, type = "response")
    
    #####
    suppressMessages(auc_outcomes["sep_comp", index] <- 
                       pROC::auc(pROC::roc(test_dat[[paste0("Y_", clean_outcome)]], 
                                           test_dat[[paste0("pred_", clean_outcome)]])))
  }
  
  ### STAGE 2: Fit second-stage models on VALID, evaluate on TEST --------------------------------------------
  for (outcome_specific in K_outcomes$outcomes) {
    clean_outcome_specific <- gsub(" ", "", tolower(outcome_specific))
    i <- match(outcome_specific, K_outcomes$outcomes)
    
    for (outcome in K_outcomes$outcomes) {
      clean_outcome <- gsub(" ", "", tolower(outcome))
      j <- match(outcome, K_outcomes$outcomes)
      
      fit_temp <- glm(as.formula(paste0(paste0("Y_", clean_outcome), paste0("~ pred_", clean_outcome_specific))), 
                      data = valid_dat, family = "binomial")
      
      ##### We dont need to store this within valid or test data because we only need this once to calculate the auc below
      pred_temp <- predict(fit_temp, newdata = test_dat, type = "response")
      
      suppressMessages(auc_outcomes[i + 3, j] <- 
                         pROC::auc(pROC::roc(test_dat[[paste0("Y_", clean_outcome)]], 
                                             pred_temp)))
    }
  }
  
  #### "New" and "Old" ALI: use the pred fitted on VALID from the outcome-specìic section above, eval on TEST ----------
  pred_cols <- (11 + K):(10 + 2 * K) 
  ##### 10 first elements are ALI comps, K next elements are Y outcomes, then the pred
  ##### print(head(valid_dat)) 
  
  valid_dat$new_ali <- rowSums(valid_dat[, pred_cols]) / K
  test_dat$new_ali  <- rowSums(test_dat[, pred_cols]) / K
  
  valid_dat$old_ali <- rowSums(valid_dat[, 1:10]) / 10
  test_dat$old_ali <- rowSums(test_dat[, 1:10]) / 10
  
  for (outcome in K_outcomes$outcomes) {
    clean_outcome <- gsub(" ", "", tolower(outcome))
    index <- match(outcome, K_outcomes$outcomes)
    
    ##### "New" 
    fit_temp <- glm(as.formula(paste0(paste0("Y_", clean_outcome), "~ new_ali")), data = valid_dat, family = "binomial")
    pred_temp <- predict(fit_temp, newdata = test_dat, type = "response")
    suppressMessages(auc_outcomes["new_ali", index] <- pROC::auc(pROC::roc(test_dat[[paste0("Y_", clean_outcome)]], pred_temp)))
    
    ##### "Old" 
    fit_temp <- glm(as.formula(paste0(paste0("Y_", clean_outcome), "~ old_ali")), data = valid_dat, family = "binomial")
    pred_temp <- predict(fit_temp, newdata = test_dat, type = "response")
    suppressMessages(auc_outcomes["old_ali", index] <- pROC::auc(pROC::roc(test_dat[[paste0("Y_", clean_outcome)]], pred_temp)))
  }
  
  return(auc_outcomes)
}