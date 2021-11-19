####################################################################################

### Summary table demographics and plots 

####################################################################################

EWAS <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_EWAS_722_demos.csv")
GWAS <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_GWAS_769_demos.csv")


# Create a demographics table with an EWAS and GWAS column 
# Get summary stats here for the table 

# n 

dim(EWAS)

dim(GWAS)


# Mean s100b and sd 

mean(EWAS$s100b, na.rm = T)
mean(GWAS$s100b, na.rm = T)

sd(EWAS$s100b, na.rm = T)
sd(GWAS$s100b, na.rm = T)


# Sex


table(EWAS$sex.x)

340 / 781 * 100


> table(GWAS$sex)

#   1   2
# 402 367
# 1 = M 2 = F

367 / 769 * 100 

# [1] 47.72432

# Age 

GWAS$age <- GWAS$age_W2 / 365.25

mean(EWAS$age, na.rm = T)
gmean <- mean(GWAS$age, na.rm = T)

sd(EWAS$age, na.rm = T)
gsd <- sd(GWAS$age, na.rm = T)


# BMI

emean <- mean(EWAS$bmi_w2, na.rm = T)
gmean <- mean(GWAS$bmi_w2, na.rm = T)

esd <- sd(EWAS$bmi_w2, na.rm = T)
gsd <- sd(GWAS$bmi_w2, na.rm = T)



# EpiSmokEr

mean(EWAS$smokingScore, na.rm = T)
sd(EWAS$smokingScore, na.rm = T)


# > mean(EWAS$smokingScore, na.rm = T)
# [1] 0.8426407
# > sd(EWAS$smokingScore, na.rm = T)
# [1] 5.148012
# >



# Smoking


table(GWAS$smokcurr_w2)

#   0   1   2
# 376 333  60


# > 376 / 769 * 100
# [1] 48.89467
# > 333 / 769 * 100
# [1] 43.30299
# > 60 / 769 * 100
# [1] 7.802341








