## Set working directory to project root if not already there
# setwd("GestALT")

## --- Load function definitions -----------------------------------------------
source("R/00_setup.R")

## --- Else if too inconvenient, run this setup --------------------------------
## Fixed constants
N <- 1000 ### sample size
## Vector of "unhealthy" probabilities were component (taken from deidentified_ali_hospitalizations.csv)
pS <- c(0.125, 0.799, 0.458, 0.239, 0.002,
        0.000, 0.000, 0.251, 0.064, 0.146)
names(pS) <- c("A1C_F", "ALB_F", "BMI_F", "CHOL_F", "CRP_F",
               "CREAT_C_F", "HCST_F", "TRIG_F", "BP_DIASTOLIC_F", "BP_SYSTOLIC_F")
pS <- pmax(pS, rep(0.01, times = length(pS))) ### force 1%+ prob per component

## Define ALI components
bin_ALI_comp <- c("A1C", "ALB", "BMI", "CHOL", "CRP", "CREAT_C", "HCST", "TRIG", "BP_DIASTOLIC", "BP_SYSTOLIC")
factor_ALI_comp <- paste0(bin_ALI_comp, "_F")

## Define outcomes for each K = number of outcomes to aggregate from -----------
outcomes_reference <- data.frame(
  outcomes = c(
    "CVD", "Stroke", "Heart Failure", "Type II Diabetes", "Obesity", 
    "Schizophrenia", "Depression", "Anxiety", "Dementia", "Frailty", 
    "Chronic Kidney Disease", "Asthma", "Osteoarithritis", "Migraines", "Stomach Cancer", 
    "Liver Cancer", "Kidney Cancer", "Escphageal Cancer", "Lung Cancer", "Breast Cancer"
  ), 
  prev = c(1:20)
)


## Function to simulate multiple outcomes data ---------------------------------
sim_data_multi <- function(N, K, betas_outcomes) {
  ## Check dimension of betas_outcomes matches with what we need 
  stopifnot(ncol(betas_outcomes) == K, nrow(betas_outcomes) == 11)
  
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
  K_outcomes <- outcomes_reference |> slice(1:K)
  
  clean_name <- gsub(" ", "", tolower(K_outcomes$outcomes))
  dynamic_names_mu <- paste0("mu_", clean_name)
  dynamic_names_Y <- paste0("Y_", clean_name)
  
  mu_outcomes <- lapply(dynamic_names_mu, function(x) 0)
  names(mu_outcomes) <- dynamic_names_mu
  Y_outcomes <- lapply(dynamic_names_Y, function(x) 0)
  names(Y_outcomes) <- dynamic_names_Y
  
  for (outcome in K_outcomes$outcomes) {
    index <- match(outcome, K_outcomes$outcomes)
    clean_outcome <- gsub(" ", "", tolower(outcome))
    beta_outcome <- betas_outcomes[, index]
    
    mu_outcomes[[paste0("mu_", clean_outcome)]] <- cbind(1, S_mat) %*% beta_outcome
    Y_outcomes[[paste0("Y_", clean_outcome)]] <- rbinom(n = N,
                                               size = 1,
                                               prob = (1 + exp(-mu_outcomes[[paste0("mu_", clean_outcome)]])) ^ (- 1))
  }
  
  Y_outcomes <- as.data.frame(Y_outcomes)
  dat <- data.frame(cbind(S_mat, Y_outcomes))
  return(dat)
}


## Test simulation function above ----------------------------------------------
set.seed(123)
my_betas <- matrix(rnorm(11 * 5), nrow = 11, ncol = 5)
colnames(my_betas) <- paste0("beta", 1:5)
my_betas
temp <- sim_data_multi(1000, 5, my_betas)
head(temp)



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
  K_outcomes <- outcomes_reference |> slice(1:K)
  
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


## Test run --------------------------------------------------------------------
aucs <- try_pred_new_ali_3folds_multi(sim_data_multi, N = 1000, K = 5, 
                              train_prop = 0.5, 
                              valid_prop = 0.25,
                              my_betas)
aucs

## Try radar plots
install.packages("devtools")
devtools::install_github("ricardo-bion/ggradar")
library(ggradar)

install.packages("fmsb")
library(fmsb)

radar_data <- rbind(
  rep(1, ncol(aucs)),   # max
  rep(0, ncol(aucs)),   # min
  aucs
)

# Select which rows to actually plot (e.g., the first 3 methods)
plot_data <- radar_data[c("sep_comp", "new_ali", "old_ali") %in% rownames(radar_data) |
                          rownames(radar_data) %in% c(NA), ]

# Simpler: subset explicitly
plot_data <- rbind(radar_data[1:2, ], radar_data[c("sep_comp", "new_ali", "old_ali"), ])

colors_border <- c("#E41A1C", "#377EB8", "#4DAF4A")
colors_in <- adjustcolor(colors_border, alpha.f = 0.2)

radarchart(plot_data,
           axistype = 1,
           pcol = colors_border,
           plwd = 2,
           plty = 1,
           cglcol = "grey70",
           cglty = 1,
           axislabcol = "grey70",
           calcex = 0.8,
           caxislabels = seq(0, 1, 0.25),  # adjust to match your min/max
           cglwd = 0.8,
           vlcex = 0.9)

legend("topright", legend = rownames(plot_data)[3:5],
       col = colors_border, lty = 1, lwd = 2, bty = "n")





