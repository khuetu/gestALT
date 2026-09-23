## Average AUC across reps for a given model type
avg_auc <- function(model_name, results) {
  mat <- do.call(rbind, lapply(results, function(r) as.numeric(r[[model_name]])))
  colnames(mat) <- names(results[[1]][[model_name]])
  colMeans(mat)
}

run_sim_print_res <- function(sim_func, N, beta_admit, beta_ed, beta_ed2) {
  n_sims <- 1000
  results <- lapply(1:n_sims, function(i) {
    out <- tryCatch(
      suppressWarnings(suppressMessages(
        try_pred_new_ali_3folds(sim_func, 
                                N = N,
                                beta_admit = beta_admit,
                                beta_ed = beta_ed,
                                beta_ed2 = beta_ed2)
      )),
      error = function(e) NULL
    )
    out
  })
  
  n_failed <- sum(sapply(results, is.null))
  if (n_failed > 0) {
    cat(sprintf("Note: %d of %d reps skipped due to degenerate test-fold outcomes\n\n", n_failed, n_sims))
  }
  results <- Filter(Negate(is.null), results)
  
  auc_outcomes <- data.frame(matrix(NA, nrow = 6, ncol = 3))
  colnames(auc_outcomes) <- c("auc_admit", "auc_ed", "auc_ed2")
  row_names <- c("sep_comp", "new_ali", "old_ali", "admit_ali", "ed_ali", "ed2_ali")
  rownames(auc_outcomes) <- row_names
  
  for (name in row_names) {
    auc_outcomes[name, ] <- avg_auc(name, results) 
  }
  
  cat("Mean AUC - Separate components:\n"); print(avg_auc("sep_comp", results))
  cat("\nMean AUC - New ALI:\n");           print(avg_auc("new_ali", results))
  cat("\nMean AUC - Old ALI:\n");           print(avg_auc("old_ali", results))
  cat("\nMean AUC - Admission-specific: Admit:\n"); print(avg_auc("admit_ali", results))
  cat("\nMean AUC - Admission-specific: Ed:\n"); print(avg_auc("ed_ali", results))
  cat("\nMean AUC - Admission-specific: Ed2:\n"); print(avg_auc("ed2_ali", results)) 
  
  return(auc_outcomes)
}