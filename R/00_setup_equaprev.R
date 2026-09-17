library("pROC")

#* -------------------------------- Simulate data ------------------------------
## Fixed constants
N <- 1000 ### sample size
## Vector of "unhealthy" probabilities were component: take everything to be equal as a base case 
pS <- rep(0.208, 10)
names(pS) <- c("A1C_F", "ALB_F", "BMI_F", "CHOL_F", "CRP_F",
               "CREAT_C_F", "HCST_F", "TRIG_F", "BP_DIASTOLIC_F", "BP_SYSTOLIC_F")
pS <- pmax(pS, rep(0.01, times = length(pS))) ### force 1%+ prob per component

## Define ALI components
bin_ALI_comp <- c("A1C", "ALB", "BMI", "CHOL", "CRP", "CREAT_C", "HCST", "TRIG", "BP_DIASTOLIC", "BP_SYSTOLIC")
factor_ALI_comp <- paste0(bin_ALI_comp, "_F")

## Specify model coefficients (for now, assumed to all be different)
### For hospital admission Y1 | S1,...,S10
my_beta_admit <- seq(0, 1, by = 0.1)
my_beta_admit <- matrix(data = my_beta_admit,
                        ncol = 1)
### For emergency department visit Y2 | S1,...,S10
my_beta_ed <- seq(0, 0.5, by = 0.05)
my_beta_ed <- matrix(data = my_beta_ed,
                     ncol = 1)
### For alt: emergency department visit (dependent on hospital admission) Y2 | S1,...,S10,Y1
my_beta_ed2 <- seq(0, 0.55, by = 0.05)
my_beta_ed2 <- matrix(data = my_beta_ed2,
                      ncol = 1)
