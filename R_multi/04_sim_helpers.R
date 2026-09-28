## ---------------------------- AUC averaging across reps ---------------------
average_auc_table <- function(results) {
  # results: list of data.frames/matrices, all same dim and dimnames
  mats <- lapply(results, as.matrix)
  arr <- simplify2array(mats)              # rows x cols x n_reps
  avg_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  avg_df <- as.data.frame(avg_mat)
  rownames(avg_df) <- rownames(mats[[1]])
  colnames(avg_df) <- colnames(mats[[1]])
  avg_df
}

run_sim_print_res <- function(sim_func, N, K, betas) {
  n_sims <- 1000
  results <- lapply(1:n_sims, function(i) {
    out <- tryCatch(
      suppressWarnings(suppressMessages(
        try_pred_new_ali_3folds_multi(sim_func,
                                      N = N,
                                      K = K,
                                      my_betas = betas)
      )),
      error = function(e) { message("Rep failed: ", conditionMessage(e)); NULL }
    )
    out
  })
  
  n_failed <- sum(sapply(results, is.null))
  if (n_failed > 0) {
    cat(sprintf("Note: %d of %d reps skipped due to degenerate test-fold outcomes\n\n", n_failed, n_sims))
  }
  results <- Filter(Negate(is.null), results)
  
  final_table <- average_auc_table(results)
  print(final_table)
  invisible(final_table)
}


## ---------------------------- Prevalence check helper ------------------------
## Generalized version of table(try_dat$Y_admit) for K outcomes
check_prev <- function(dat, K_outcomes) {
  y_cols <- paste0("Y_", gsub(" ", "", tolower(K_outcomes$outcomes)))
  sapply(dat[y_cols], function(x) c(n0 = sum(x == 0), n1 = sum(x == 1), prev = mean(x)))
}

## ---------------------------- Prevalence scenario helper ------------------------
check_scenario <- function(my_betas, label, K = 5) {
  dat <- sim_data_multi(N = 200000, K = K, betas_outcomes = my_betas)
  
  K_outcomes <- outcomes_reference |> slice(1:K)
  K_target_prev <- target_prev[1:K] 
  
  simulated <- check_prev(dat, K_outcomes)["prev", ]

  cat(sprintf("\n--- %s ---\n", label))
  print(rbind(target = K_target_prev, simulated = simulated))
}

## ---------------------------- Intercept calibration ---------------------------
## Calc marginal prevalence given an intercept + a component-coefficient assignment
## We know the P(Y = 1 | S = s) is equal plogis(β₀ + sᵀβ), where s is s vector 10 components 
## We want to find P(Y=1), using the total law of prob, P(Y=1) = \sum_s P(Y = 1, S = s) over all possible combination of s
calc_prevalence <- function(beta0, betas_comp, pS_vec) {
  # grid 2^10 = 1024, S assume to be independent 
  grid  <- as.matrix(expand.grid(rep(list(c(0, 1)), length(pS_vec))))
  # P(S = s) = \prod_{k=1}^{K} p_k^{s_k} * (1 - p_k)^{1 - s_k} (Bernoulli)
  probs <- apply(grid, 1, function(s) prod(ifelse(s == 1, pS_vec, 1 - pS_vec)))
  # mu = β₀ + sᵀβ
  mu    <- beta0 + grid %*% betas_comp
  # P(Y=1) = \sum_s P(Y = 1, S = s) = \sum_s P(Y = 1 | S = s) * P(S = s) 
  sum(probs * plogis(mu))
}

## Solve for the intercept that enable this ordering to achieve target prevalence 
solve_intercept <- function(target_prev, betas_comp, pS_vec = pS, interval = c(-25, 25)) {
  uniroot(function(beta0) calc_prevalence(beta0, betas_comp, pS_vec) - target_prev,
          interval = interval)$root
}

## ---------------------------- Beta-matrix builders -----------------------------
## Build an 11 x K beta matrix from a list of K component-orderings (each a
## permutation of 1:10), with beta0 is so that every ordering outputs the same target marginal prevalence.
build_my_betas <- function(order_list, base_coefs, scale, target_prev, pS_vec = pS) {
  stopifnot(length(order_list) == length(scale), length(scale) == length(target_prev))
  
  # scale 
  comp_mat <- sapply(seq_along(order_list), function(k) scale[k] * base_coefs[order_list[[k]]])
  # solve for beta0 given prev for each outcome 
  beta0 <- sapply(seq_along(order_list), function(k)
    solve_intercept(target_prev[k], comp_mat[, k], pS_vec))
  
  my_betas <- rbind(beta0, comp_mat)
  colnames(my_betas) <- paste0("beta", seq_along(order_list))
  rownames(my_betas) <- NULL
  my_betas
}

## Build an 11 x K beta matrix for the "A1C dominant" scenario
build_a1c_dominant <- function(scale, rest_mag, target_prev, pS_vec = pS) {
  stopifnot(length(rest_mag) == length(scale), length(scale) == length(target_prev))
  
  # scale 
  comp_mat <- sapply(seq_along(scale), function(k) scale[k] * c(1.0, rep(rest_mag[k], 9)))
  # solve for beta0 given prev for each outcome 
  b0 <- sapply(seq_along(scale), function(k)
    solve_intercept(target_prev[k], comp_mat[, k], pS_vec))
  
  my_betas <- rbind(b0, comp_mat)
  colnames(my_betas) <- paste0("beta", seq_along(scale))
  rownames(my_betas) <- NULL
  my_betas
}