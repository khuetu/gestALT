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

## ---------------------------- Intercept calibration ---------------------------
## Exact marginal prevalence given an intercept + a specific component-coefficient
## assignment, using the true independent-Bernoulli distribution of S (2^10 = 1024
## patterns -- exact)
exact_prevalence <- function(beta0, betas_comp, pS_vec) {
  grid  <- as.matrix(expand.grid(rep(list(c(0, 1)), length(pS_vec))))
  probs <- apply(grid, 1, function(s) prod(ifelse(s == 1, pS_vec, 1 - pS_vec)))
  mu    <- beta0 + grid %*% betas_comp
  sum(probs * plogis(mu))
}

## Solve for the intercept that makes THIS ordering hit the target prevalence exactly
calibrate_intercept <- function(target_prev, betas_comp, pS_vec = pS, interval = c(-25, 25)) {
  uniroot(function(b0) exact_prevalence(b0, betas_comp, pS_vec) - target_prev,
          interval = interval)$root
}

## ---------------------------- Beta-matrix builders (calibrated) ---------------
## Build an 11 x K beta matrix from a list of K component-orderings (each a
## permutation of 1:10), with beta0 calibrated per-column so every ordering
## reproduces the SAME target marginal prevalence.
build_my_betas <- function(order_list, base_coefs, scale, target_prev, pS_vec = pS) {
  stopifnot(length(order_list) == length(scale), length(scale) == length(target_prev))
  comp_mat <- sapply(seq_along(order_list), function(k) scale[k] * base_coefs[order_list[[k]]])
  b0 <- sapply(seq_along(order_list), function(k)
    calibrate_intercept(target_prev[k], comp_mat[, k], pS_vec))
  my_betas <- rbind(b0, comp_mat)
  colnames(my_betas) <- paste0("beta", seq_along(order_list))
  rownames(my_betas) <- NULL
  my_betas
}

## Build an 11 x K beta matrix for the "A1C dominant" scenario, calibrated the same way
build_a1c_dominant <- function(scale, rest_mag, target_prev, pS_vec = pS) {
  stopifnot(length(rest_mag) == length(scale), length(scale) == length(target_prev))
  comp_mat <- sapply(seq_along(scale), function(k) scale[k] * c(1.0, rep(rest_mag[k], 9)))
  b0 <- sapply(seq_along(scale), function(k)
    calibrate_intercept(target_prev[k], comp_mat[, k], pS_vec))
  my_betas <- rbind(b0, comp_mat)
  colnames(my_betas) <- paste0("beta", seq_along(scale))
  rownames(my_betas) <- NULL
  my_betas
}