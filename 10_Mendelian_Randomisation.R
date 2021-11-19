
#########################################################################################

### MR WITH S100B

#########################################################################################

library(devtools)
library(TwoSampleMR)
library(data.table)

# # Format s100b file for MR base upload 
# # Format s100b exposure data as per robs example input 
# ex <- read.delim("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/fastGWA_Dec/MR_Inputs/O15232.txt") # rob example input 
# file <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits.csv") # s100b hits 
# file <- file[c(1,3,2,4,5,9,6,10,11,13,3,14)] # set order to match example input rob 
# names(file) <- names(ex) # set names to match example input rob 
# up <- file[c(12,8,9,4,5,7)]
# names(up) <- c("SNP", "beta", "se", "effect_allele", "other_allele", "eaf")
# write.table(up, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/MR_coloc_090921/MR/top_hits_mr_base.txt", quote = F, row.names = F)


## S100B to AD direction 

# Format s100b exposure data as per robs example input 
ex <- read.delim("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/fastGWA_Dec/MR_Inputs/O15232.txt") # rob example input 
file <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits.csv") # s100b hits 
file <- file[c(1,3,2,4,5,9,6,10,11,13,3,14)] # set order to match example input rob 
names(file) <- names(ex) # set names to match example input rob 
write.table(file, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits.txt", quote = F, row.names = F)

# file <- "/Cluster_Filespace/Marioni_Group/Rob/Somalogic/fastGWA_Dec/MR_Inputs/O15232.txt" - this one reads in fine 
file <-"/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits.txt"

exposure_dat <- read_exposure_data( 
filename = file, 
sep = " ", 
snp_col = 'rsid',
beta_col = 'BETA',
se_col = 'SE',
effect_allele_col = 'A1',
phenotype_col = '',
units_col = '', 
other_allele_col = 'A2',
eaf_col = 'AF1', 
samplesize_col = 'N', 
ncase_col = '',
ncontrol_col = '',
gene_col = '',
pval_col = 'P' 
)


# Format outcome data as per robs example input 
ad = fread("/Cluster_Filespace/Marioni_Group/Rob/AD_sumstats_Jansenetal.txt")
ad =as.data.frame(ad)
names(ad) <- c("Marker", "CHR", "BP", "A1", "A2", "SNP", "Z", "P", "N", "Neff", "dir", "MAF", "BETA", "SE")
write.table(ad, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/AD_janssen_sumstats.txt", quote = F, row.names = F)

# ex2 <- fread("/Cluster_Filespace/Marioni_Group/Rob/Neurology Proteins/MR/4_UK_Biobank_IGAP_17May2018.txt") - this one reads in fine 
ad <- fread("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/AD_janssen_sumstats.txt")

file2 <- "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/AD_janssen_sumstats.txt"

outcome_dat <- read_outcome_data(
filename = file2, 
snp_col = 'SNP',
beta_col = 'BETA',
se_col = 'SE',
effect_allele_col = 'A1',
phenotype_col = '',
units_col = '',
other_allele_col = 'A2',
eaf_col = '',
samplesize_col = '',
ncase_col = '',
ncontrol_col = '',
gene_col = '',
pval_col = 'P'
)


# Clump the exposure 
exposure_dat <- clump_data(exposure_dat)

# Check SNPs match 
length(which(exposure_dat$SNP %in% outcome_dat$SNP)) # yes, our SNP is present 

# Harmonise alleles if needed 
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

# Check prior to MR
length(dat$mr_keep) == length(which(dat$mr_keep %in% "FALSE"))

## Run MR 
mr_results <- mr(dat) 

generate_odds_ratios(mr_results)

# write.csv(mr_results, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/MR_coloc_090921/MR/MR_s100b_to_AD_100921.csv", row.names = F)
write.csv(mr_results, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/MR_coloc_090921/MR/MR_s100b_to_AD_211021_ORs.csv", row.names = F)

#################################################################################

### NOW DO THE REVERSE 

library(devtools)
library(TwoSampleMR)
library(data.table)

## AD to S100B direction

# Format AD exposure data as per robs example input 
ad = fread("/Cluster_Filespace/Marioni_Group/Rob/AD_sumstats_Jansenetal.txt")
ad =as.data.frame(ad)
names(ad) <- c("Marker", "CHR", "BP", "A1", "A2", "SNP", "Z", "P", "N", "Neff", "dir", "MAF", "BETA", "SE")
ad <- ad[which(ad$P < 5e-8),] # 2357
ad$EAF <- 1 - ad$MAF
write.table(ad, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/AD_janssen_sumstats_reversed.txt", quote = F, row.names = F)

file <- "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/AD_janssen_sumstats_reversed.txt"

exposure_dat <- read_exposure_data(
filename = file, 
snp_col = 'SNP',
beta_col = 'BETA',
se_col = 'SE',
effect_allele_col = 'A1',
phenotype_col = '',
units_col = '',
other_allele_col = 'A2',
eaf_col = 'EAF',
samplesize_col = '',
ncase_col = '',
ncontrol_col = '',
gene_col = '',
pval_col = 'P',
)


# Format s100b outcome data as per robs example input 
ex <- read.delim("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/fastGWA_Dec/MR_Inputs/O15232.txt") # rob example input 
file <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/s100b_GWAS_output.txt") # s100b hits (all)
file <- file[c(1,3,2,4,5,9,6,10,11,13,3,14)] # set order to match example input rob 
names(file) <- names(ex) # set names to match example input rob 
file[file==""] <- NA
table(is.na(file))
file <- na.omit(file)
write.table(file, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits_reversed.txt", quote = F, row.names = F)

# file <- "/Cluster_Filespace/Marioni_Group/Rob/Somalogic/fastGWA_Dec/MR_Inputs/O15232.txt" - this one reads in fine 
file2 <-"/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits_reversed.txt"

outcome_dat <- read_outcome_data( 
filename = file2, 
sep = " ", 
snp_col = 'rsid',
beta_col = 'BETA',
se_col = 'SE',
effect_allele_col = 'A1',
phenotype_col = '',
units_col = '', 
other_allele_col = 'A2',
eaf_col = 'AF1', 
samplesize_col = 'N', 
ncase_col = '',
ncontrol_col = '',
gene_col = '',
pval_col = 'P' 
)



# Clump the exposure 
exposure_dat <- clump_data(exposure_dat) # 2327 of 2357 variants removed - leaving 30 post-clumping 

# Check SNPs match 
length(which(exposure_dat$SNP %in% outcome_dat$SNP)) # yes, 29/30 present 

# Harmonise alleles if needed 
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)


# Harmonising exposure (5O2T0B) and outcome (XbXEkR)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
# rs10933431, rs111278892, rs148601586, rs346757, rs41289512, rs6931277


# Check prior to MR
# length(dat$mr_keep) == length(which(dat$mr_keep %in% "FALSE"))

# Remove non harmonised SNPs 
dat2 <- dat[which(dat$mr_keep %in% TRUE),]

## Run MR 
mr_results <- mr(dat2) 
write.csv(mr_results, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/MR_coloc_090921/MR/MR_AD_to_s100b_100921.csv", row.names = F)

# Check in mr eggar for pleitorpohy 
mr_pleiotropy_test(dat2)





# #################################################################################

# ### ROB ORIGINAL EXAMPLE CODE FOR REFERENCE 

# ## Set up loop

# library(devtools)

# library(TwoSampleMR)

# snps <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/GWAS_Hits_Final_Dec_with_betas.csv")

# ## Gets SOMAMers/proteins with significant pQTLs 
# list = unique(snps$Somamer)

# output <- list()

# i <- "O15232"

# for(i in list){ 

# ## Protein as exposure 

# print(i)

# ## Open exposure - Protein
# file = paste0("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/fastGWA_Dec/MR_Inputs/",i, ".txt")


# exposure_dat <- read_exposure_data( 
# filename = file, 
# sep = "\t", 
# snp_col = 'rsid',
# beta_col = 'BETA',
# se_col = 'SE',
# effect_allele_col = 'A1',
# phenotype_col = '',
# units_col = '', 
# other_allele_col = 'A2',
# eaf_col = 'AF1', 
# samplesize_col = 'N', 
# ncase_col = '',
# ncontrol_col = '',
# gene_col = '',
# pval_col = 'P' 
# )


# ## Clump Exposures 

# exposure_dat <- clump_data(exposure_dat)

# ## Open Outcome  ### Example of Alzheimer's Disease trait from before - should update to newer Alzheimer's GWAS if you're running this, like Jansen et al. (mine got deleted in last server update - it was in /home/robert)####


# outcome_dat <- read_outcome_data(
# filename = "/Cluster_Filespace/Marioni_Group/Rob/Neurology Proteins/MR/4_UK_Biobank_IGAP_17May2018.txt", 
# snp_col = 'SNP',
# beta_col = 'BETA',
# se_col = 'SE',
# effect_allele_col = 'A1',
# phenotype_col = '',
# units_col = '',
# other_allele_col = 'A2',
# eaf_col = '',
# samplesize_col = '',
# ncase_col = '',
# ncontrol_col = '',
# gene_col = '',
# pval_col = 'P'
# )


# ## Check to see if SNPs match 
# if(length(which(exposure_dat$SNP %in% outcome_dat$SNP)) == 0) { 
# NULL 
# } else { 


# ## Harmonisation of Alleles 

# dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#   if(length(dat$mr_keep) == length(which(dat$mr_keep %in% "FALSE"))){ 
# NULL

# } else { 
# ## Run MR 

# mr_results <- mr(dat) 
# mr_results$Somamer <- i

# output[[i]] <- mr_results

# }
# } 

# }


# output <- do.call("rbind", output)


# write.csv(output, "/Cluster_Filespace/Marioni_Group/Rob/Somalogic/MR/Proteins_on_AD.csv", row.names = F)

##########################################################################################