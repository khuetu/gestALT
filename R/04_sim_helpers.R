## Average AUC across reps for a given model type
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
