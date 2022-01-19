
#########################################################################################

### MR WITH S100B

#########################################################################################

library(devtools)
library(TwoSampleMR)
library(data.table)

## S100B to AD direction 

# Format s100b exposure data as per robs example input 
ex <- read.delim("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/fastGWA_Dec/MR_Inputs/O15232.txt") # rob example input 
file <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits.csv") # s100b hits 
file <- file[c(1,3,2,4,5,9,6,10,11,13,3,14)] # set order to match example input
names(file) <- names(ex) # set names to match example input
write.table(file, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits.txt", quote = F, row.names = F)

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
eaf_col = 'MAF',
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
write.table(ad, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/AD_janssen_sumstats_rev.txt", quote = F, row.names = F)

file <- "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/AD_janssen_sumstats_rev.txt"

exposure_dat <- read_exposure_data(
filename = file, 
snp_col = 'SNP',
beta_col = 'BETA',
se_col = 'SE',
effect_allele_col = 'A1',
phenotype_col = '',
units_col = '',
other_allele_col = 'A2',
eaf_col = 'MAF',
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
write.table(file, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits_rev.txt", quote = F, row.names = F)

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

# Remove non harmonised SNPs 
dat2 <- dat[which(dat$mr_keep %in% TRUE),]

## Run MR 
mr_results <- mr(dat2) 
write.csv(mr_results, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/MR_coloc_090921/MR/MR_AD_to_s100b_rev2.csv", row.names = F)


