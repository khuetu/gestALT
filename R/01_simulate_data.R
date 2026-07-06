sim_data <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
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
  ### Hospital admission Y1 | S1,...,S10
  mu_admit <- cbind(1, S_mat) %*% beta_admit
  Y_admit <- rbinom(n = N,
                    size = 1,
                    prob = (1 + exp(-mu_admit)) ^ (- 1))
  ### Emergency department visit Y2 | S1,...,S10
  mu_ed <- cbind(1, S_mat) %*% beta_ed
  Y_ed <- rbinom(n = N,
                 size = 1,
                 prob = (1 + exp(-mu_ed)) ^ (- 1))
  ### Alt: Emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  mu_ed2 <- cbind(1, S_mat, Y_admit) %*% beta_ed2
  Y_ed2 <- rbinom(n = N,
                  size = 1,
                  prob = (1 + exp(-mu_ed2)) ^ (- 1))
  ## Combine all columns and return as dataframe
  dat <- data.frame(cbind(S_mat, Y_admit, Y_ed, Y_ed2))
  return(dat)
}

# Create new sim data function to include age and sex
sim_data_agesex <- function(N = 1000, beta_admit, beta_ed, beta_ed2) {
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
  ### Age: mean aged 48, range 18-65
  age <- rnorm(N, mean = 48, sd = 15)
  age <- pmin(pmax(age, 18), 65)
  
  ### Sex: binary (0 = male, 1 = female)
  sex <- rbinom(N, size = 1, prob = 0.6)
  
  ### Hospital admission Y1 | S1,...,S10
  mu_admit <- cbind(1, S_mat, age, sex) %*% beta_admit
  Y_admit <- rbinom(n = N,
                    size = 1,
                    prob = (1 + exp(-mu_admit)) ^ (- 1))
  ### Emergency department visit Y2 | S1,...,S10
  mu_ed <- cbind(1, S_mat, age, sex) %*% beta_ed
  Y_ed <- rbinom(n = N,
                 size = 1,
                 prob = (1 + exp(-mu_ed)) ^ (- 1))
  ### Alt: Emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
  mu_ed2 <- cbind(1, S_mat, Y_admit, age, sex) %*% beta_ed2
  Y_ed2 <- rbinom(n = N,
                  size = 1,
                  prob = (1 + exp(-mu_ed2)) ^ (- 1))
  ## Combine all columns and return as dataframe
  dat <- data.frame(cbind(S_mat, age = age, sex = sex, Y_admit, Y_ed, Y_ed2))
  return(dat)
}
