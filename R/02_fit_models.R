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
