####################################################################################

### Phenotype file preps

####################################################################################

library(tidyverse)

# s100b protein measures 
s100b <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/W2_s100b_6_outliers_above_4_sd_from_mean_removed.csv")
names(s100b)[1] <- "ID"

# Target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
target <- target %>% filter(WAVE == "2") # Filter to wave 2
target <- target %>% filter(cohort == "LBC36") # Filter to LBC36 - 801 people
target$lbc36no <- as.character(target$ID) # Make sure no factor issues when joining by IDs 

# Create the dataset with both DNAm and s100b at W2
overlap <- which(s100b$lbc36no %in% target$lbc36no) # 790 
join <- merge(target, s100b, by = "ID") # 790

# Load IDs for people that will need to be removed from LBC and remove from dataset
ID_remove <- c("LBC360558", "LBC360760", "LBC360536", "LBC361272",
"LBC360213", "LBC360262", "LBC361264", "LBC361030", "LBC360412",
"LBC361076", "LBC360721")
overlap <- which(join$lbc36no %in% ID_remove)
length(overlap) # 9
data <- join[-overlap,]
dim(data) # 781

# Read in epismoker scores and merge with phenotypes file 
smok <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/lbc_epismoker.rds")
data3 <- merge(data, smok, by.x = "Basename", by.y = "SampleName")

# Add BMI
library(haven)
phenos <- read_sav("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/LBC1936_Blood_DNAm_And_Brain_DNAm_RM_22APR2020.sav")
test <- as.data.frame(phenos)
bmi <- test[c(1,24)]
names(bmi)[1] <- "ID"
data4 <- merge(data3, bmi, by = "ID")

# Read in principal components for LBC1936 and merge with phenotypes file
pcs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/LBC1936_PCA_from_Gail_rob_sent.csv")
names(pcs)[1] <- "ID"
data5 <- merge(data4, pcs, by = "ID") # 722 people 

## Save phenotypes file 
write.csv(data5, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_EWAS_722_demos.csv",row.names = F)

## Log Transform 
phenotypes <- data5
for(i in colnames(phenotypes)[22]){ 
  phenotypes[,i]<- log(phenotypes[,i])
}

## Regress Proteins onto Covariates 
phenotypes_residualised <- phenotypes
for(i in colnames(phenotypes_residualised)[22]){ 
  phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex.x) + bmi_w2 +
                                     C1 + C2 + C3 + C4, 
                                    na.action = na.exclude, data = phenotypes_residualised)$residuals
}

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(phenotypes_residualised)[22]){ 
  phenotypes_residualised[,i]<- orderNorm(phenotypes_residualised[,i])$x.t
}

## Scale 
phenotypes_residualised$s100b <- scale(phenotypes_residualised$s100b)

## Save phenotypes file 
write.csv(phenotypes_residualised, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_EWAS_722.csv",row.names = F)

######################################################################################

### Genetic preps 

library(tidyverse)

# s100b protein measures 
s100b <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/W2_s100b_6_outliers_above_4_sd_from_mean_removed.csv")

# Read in principal components for LBC1936 and merge with phenotypes file
pcs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/LBC1936_PCA_from_Gail_rob_sent.csv")
data2 <- merge(s100b, pcs, by = "lbc36no") # 769 people 
names(data2)[1] <- "ID"

# Add BMI and smoking status 
library(haven)
phenos <- read_sav("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/LBC1936_Blood_DNAm_And_Brain_DNAm_RM_22APR2020.sav")
test <- as.data.frame(phenos)
bmi <- test[c(1,24,63)]
names(bmi)[1] <- "ID"
data4 <- left_join(data2, bmi, by = "ID")

data4$age <- data4$age_W2 / 365.25

## Save phenotypes file 
write.csv(data4, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_GWAS_769_demos.csv",row.names = F)

## Log Transform 
phenotypes <- data4
for(i in colnames(phenotypes)[4]){ 
  phenotypes[,i]<- log(phenotypes[,i])
}

## Regress Proteins onto Covariates 
phenotypes_residualised <- phenotypes
for(i in colnames(phenotypes_residualised)[4]){ 
  phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + bmi_w2 +
                                     C1 + C2 + C3 + C4, 
                                    na.action = na.exclude, data = phenotypes_residualised)$residuals
}

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(phenotypes_residualised)[4]){ 
  phenotypes_residualised[,i]<- orderNorm(phenotypes_residualised[,i])$x.t
}

## Scale 
phenotypes_residualised$s100b <- scale(phenotypes_residualised$s100b)

## Save phenotypes file 
write.csv(phenotypes_residualised, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_GWAS_769.csv",row.names = F)

         
## Make phenotype file using same phenotype file as reference as above for the reml analyses as a .phen file 
## 3 columns: IID, FID, Phenotype
phenos <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_GWAS_769.csv")

pheno2 <- data.frame(FID = phenos$ID,
                    IID = phenos$ID,
                    pheno = phenos$s100b)
                    
#save out phenotype file.  
write.table(pheno2, file = "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/GREML_080921/REML/s100b_GWAS_769.phen", col.names = F, quote = F, row.names=F, sep=' ')    

######################################################################################





















