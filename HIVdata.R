# get outcome info
# get correlation info 

hivdat <- read.csv("~/Downloads/hiv_6624_final.csv")
hivdat$lVload <- log(hivdat$VLOAD, base=10)
hivdat$hard_drugs <- as.factor(hivdat$hard_drugs)
#Create variables for differences in outcomes
#only consider observations with all 4 outcomes
hivdat_0 <- hivdat[hivdat$years==0 & complete.cases(hivdat$AGG_MENT, hivdat$lVload, hivdat$AGG_PHYS, hivdat$LEU3N),]
hivdat_2 <- hivdat[hivdat$years==2 & complete.cases(hivdat$AGG_MENT, hivdat$lVload, hivdat$AGG_PHYS, hivdat$LEU3N),]

hivdat_clean <- merge(
  hivdat_0,
  hivdat_2,
  by = "newid",
  suffixes = c("_0", "_2")
)




hivdat_clean$diff_LEU3N    <- hivdat_clean$LEU3N_2    - hivdat_clean$LEU3N_0
hivdat_clean$diff_lvload   <- hivdat_clean$lVload_2   - hivdat_clean$lVload_0
hivdat_clean$diff_AGG_MENT <- hivdat_clean$AGG_MENT_2 - hivdat_clean$AGG_MENT_0
hivdat_clean$diff_AGG_PHYS <- hivdat_clean$AGG_PHYS_2 - hivdat_clean$AGG_PHYS_0

# adherence at year 2
hivdat_clean$ADH <- if_else(hivdat_clean$ADH_2>2,1,0)
hivdat_clean$ADH_cat <- if_else(hivdat_clean$ADH_2>2,"95% Adherence or Higher","Less than 95% Adherence")


#education:  5-7 or 1-4 for college educ or higher vs. no college
hivdat_clean$educ_cat <- if_else(hivdat_clean$EDUCBAS_0>4, "College Degree or Higher", "No College Degree")

hivdat_clean$race_cat <- if_else(hivdat_clean$RACE_0==1, "White, Non-Hispanic", "Other")

# SMOKE Smoking status 1= Never smoked
# 2= Former smoker
# 3= Current smoker
# Blank= Missing
hivdat_clean$smoke_cat <- if_else(hivdat_clean$SMOKE_0==1, "Never Smoked", 
                                  if_else(hivdat_clean$SMOKE_0==2, "Former Smoker", 
                                          if_else(!is.na(hivdat_clean$SMOKE_0), "Current Smoker", NA)))

#remove BMI outliers
hivdat_clean$BMI <- if_else((hivdat_clean$BMI_0 > 0 & hivdat_clean$BMI_0 < 250), hivdat_clean$BMI_0, NA)

# BMI, age_0, hard_drugs_0

hivdat_sub <- hivdat_clean |>select(BMI, age_0, DIAB_0, hard_drugs_0, diff_lvload)
hivdat_sub$hard_drugs_0 <- as.numeric(hivdat_sub$hard_drugs_0)
#mm <- model.matrix(~ . - 1, data = hivdat_clean)   # dummy-codes factors/characters
C <- cov(hivdat_sub, use = "pairwise.complete.obs")

vars5 <- c("BMI", "age_0", "DIAB_0", "hard_drugs_01", "diff_lvload")

# Correlations with a different pattern from the original:
# older age -> slightly lower BMI, less hard drug use; hard drug use -> larger viral load change
R5 <- matrix(1, 5, 5, dimnames = list(vars5, vars5))
R5[upper.tri(R5)] <- c(
  -0.10,                 # BMI ~ age
  0.25,  0.20,          # DIAB ~ BMI, age
  -0.05, -0.15,  0.10,   # hard_drugs ~ BMI, age, DIAB
  0.10, -0.05,  0.15,  0.20   # diff_lvload ~ BMI, age, DIAB, hard_drugs
)
R5[lower.tri(R5)] <- t(R5)[lower.tri(R5)]

# SDs similar in scale to your original
sds <- c(BMI = 4.34, age_0 = 8.74, DIAB_0 = 3.87,
         hard_drugs_01 = 0.268, diff_lvload = 1.23)

C_trial <- diag(sds) %*% R5 %*% diag(sds)
dimnames(C_trial) <- list(vars5, vars5)
mod <-lm(diff_lvload ~ BMI + age_0 + DIAB_0+hard_drugs_0, data = hivdat_clean)
summary(mod)

# so in the treatment group, your outcome model is 
# this is your outcome model 
#-2.8 * trt + -.1 * X1 (BMI) + .9 * age_0 + .5 * X3 (diab) + .4 * X4 (hard_drugs) 

#pi_trial <- plogis(selection_lp(x1, x2, scenario))
#baseline = 0.1 + 1.4 * x1 - 1.0 * x2,

# here is your trial participation model
#pi_trial <- plogis(.1 + .05 * X1  +.01 * X2+ .5 * X3 + .2 * X4)



