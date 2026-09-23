## Set working directory to project root if not already there
# setwd("GestALT")
library(tidyverse)

## --- Load function definitions ----------------------------------------------
source("R_multi/00_setup.R")
source("R_multi/01_simulate_data.R")
source("R_multi/03_predict_auc.R")
source("R_multi/04_sim_helpers.R")



#* ----------------------------- Coefficients sandbox (multi-outcome, K = 5) ---
set.seed(123)

K <- 5
K_outcomes <- outcomes_reference |> slice(1:K)
target_prev <- K_outcomes$prev   # the prevalence every scenario below must match

## Magnitudes 0.1 to 1.0 -- one per ALI component
base_coefs <- seq(0.1, 1, by = 0.1)

## Sign/scale per outcome -- matches your original "Test" my_betas pattern
outcome_scale <- c(0.5, 0.05, -0.5, -0.5, 0.5)

#* --- Default ------------------------------------------------------------------
my_betas <- matrix(c(seq(0, 0.5, length.out = 11), 
                     seq(0, 0.5, length.out = 11), 
                     seq(0, 0.5, length.out = 11), 
                     seq(0, -0.5, length.out = 11), 
                     seq(0, 0.5, length.out = 11)
), nrow = 11, ncol = 5)
beta0s <- qlogis(outcomes_reference$prev + 0.02)
my_betas[1, ] <- beta0s[1:5]

#* --- Ordering -----------------------------------------------------------------
order_ref <- 1:10   # A1C -> SYSTOLIC, increasing magnitude
order_opp <- 10:1   # SYSTOLIC -> A1C, decreasing magnitude

#* --- Scenario 1: only outcome 1 gets the opposite ordering --------------------
order_list_opp <- list(order_opp, order_ref, order_ref, order_ref, order_ref)
my_betas_opp   <- build_my_betas(order_list_opp, base_coefs, outcome_scale, target_prev)

#* --- Scenario 2: only outcome 1 gets a partial swap ---------------------------
order_part      <- c(1, 2, 3, 6, 5, 4, 7, 8, 10, 9)   # swap ranks 4/6 and last two
order_list_part <- list(order_part, order_ref, order_ref, order_ref, order_ref)
my_betas_part   <- build_my_betas(order_list_part, base_coefs, outcome_scale, target_prev)

#* --- Scenario 3: only outcome 1 is A1C dominant -------------------------------
## Outcome 1 column from build_a1c_dominant, outcomes 2-5 from build_my_betas
a1c_col  <- build_a1c_dominant(outcome_scale[1], rest_mag = 0.4, target_prev[1])
ref_cols <- build_my_betas(list(order_ref, order_ref, order_ref, order_ref),
                           base_coefs, outcome_scale[-1], target_prev[-1])

my_betas_a1cdom <- cbind(a1c_col, ref_cols)
colnames(my_betas_a1cdom) <- paste0("beta", 1:5)

#* ------------------------ Sanity check: prevalence match ----------------------
check_scenario(my_betas_opp,    "Outcome 1: Opposite Ordering")
check_scenario(my_betas_part,   "Outcome 1: Partial Ordering")
check_scenario(my_betas_a1cdom, "Outcome 1: A1C Dominant")

#* -------------------------------- Orderings -----------------------------------
auc_default <- run_sim_print_res(sim_data_multi, N = 10000, K = K, betas = my_betas)
auc_opp <- run_sim_print_res(sim_data_multi, N = 10000, K = K, betas = my_betas_opp)
auc_part <- run_sim_print_res(sim_data_multi, N = 10000, K = K, betas = my_betas_part)
auc_a1cdom <- run_sim_print_res(sim_data_multi, N = 10000, K = K, betas = my_betas_a1cdom)

#* ---------------------- Correlation vs No correlation ------------------------
auc_corr <- run_sim_print_res(sim_data_multi_corr, N = 10000, K = K, betas = my_betas) # same my_betas
#* ---------------------- Equal Prevalence vs Empirical ------------------------
auc_eq_prev <- run_sim_print_res(sim_data_multi_eqprev, N = 10000, K = K, betas = my_betas) # same my_betas

#* ------------------------------- Store Data ------------------------------------
library(writexl)
library(openxlsx)

aucs_list <- list(
  "Default Ordering" = auc_default,
  "Opposite Ordering" = auc_opp,
  "Partial Ordering" = auc_part,
  "A1C Dominant Ordering" = auc_a1cdom, 
  "Correlated Components (Default)" = auc_corr, 
  "Equal Component Prev (Default)" = auc_eq_prev
)

write.xlsx(aucs_list, file = "aucs_K5.xlsx", rowNames = TRUE, colWidths = "auto")

#* ------------------------------ Radar/Spider Plot -------------------------------
library(fmsb)

aucs <- aucs_list[["A1C Dominant Ordering"]]

radar_data <- rbind(rep(1, ncol(aucs)), rep(0.5, ncol(aucs)),
                    aucs["sep_comp", ], aucs["new_ali", ], aucs["old_ali", ])
rownames(radar_data)[1:2] <- c("max", "min")

colors_border <- c("#E41A1C", "#377EB8", "#4DAF4A")
colors_in <- adjustcolor(colors_border, alpha.f = 0.2)

radarchart(radar_data,
           axistype = 1, pcol = colors_border, plwd = 2, plty = 1,
           cglcol = "grey70", cglty = 1, axislabcol = "grey70", calcex = 0.8,
           caxislabels = round(seq(0.5, 1, length = 5), 2), cglwd = 0.8, vlcex = 0.9)

legend("topright", legend = c("sep_comp", "new_ali", "old_ali"),
       col = colors_border, lty = 1, lwd = 2, bty = "n")




