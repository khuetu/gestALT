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
  K_outcomes <- outcomes_reference[1:K, ]
  
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

## equal prev ------------------------------------------------------------------
sim_data_multi_eqprev <- function(N, K, betas_outcomes) {
  ## Check dimension of betas_outcomes matches with what we need 
  stopifnot(ncol(betas_outcomes) == K, nrow(betas_outcomes) == 11)
  
  ## Begin with ALI components (10 per person) from Bernoulli (pS)
  S <- rbinom(n = N * 10,
              size = 1,
              prob = rep(x = pS_eq, times = N))
  
  S_mat <- matrix(data = S,
                  nrow = N,
                  ncol = 10,
                  byrow = TRUE)
  colnames(S_mat) <- factor_ALI_comp
  ## Simulate outcomes
  K_outcomes <- outcomes_reference[1:K, ]
  
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

# ------------- Similar simulation but for correlated outcomes -----------------
sim_data_multi_corr <- function(N, K, betas_outcomes) {
  ## Check dimension of betas_outcomes matches with what we need 
  stopifnot(ncol(betas_outcomes) == K, nrow(betas_outcomes) == 11)
  
  ## --- Correlated ALI components -------------------------------------------
  ## Draw correlated continuous biomarkers, then threshold each at its
  ## real-world cutoff to get the binary "unhealthy" indicators (S_mat),
  ## replacing the independent rbinom(pS) draws.
  X_bio <- MASS::mvrnorm(n = N, mu = mu_bio, Sigma = Sigma_bio)
  X_bio <- X_bio[, names(mu_bio), drop = FALSE]
  
  S_mat <- matrix(NA_real_, nrow = N, ncol = length(thresholds),
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
  
  ## --- Simulate outcomes (unchanged generalized logic) ----------------------
  K_outcomes <- outcomes_reference[1:K, ]
  
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
                                                        prob = plogis(mu_outcomes[[paste0("mu_", clean_outcome)]]))
  }
  
  Y_outcomes <- as.data.frame(Y_outcomes)
  dat <- data.frame(cbind(S_mat, Y_outcomes))
  return(dat)
}