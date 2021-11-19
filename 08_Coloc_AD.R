library(coloc)
library(data.table)

## Read in AD sumstats - Jansen et al 2019 

ad = fread("/Cluster_Filespace/Marioni_Group/Rob/AD_sumstats_Jansenetal.txt")
ad =as.data.frame(ad)
names(ad) <- c("Marker", "CHR", "BP", "A1", "A2", "SNP", "Z", "P", "N", "Neff", "dir", "MAF", "BETA", "SE")

## Read in s100 beta data

s100b = fread("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/s100b_GWAS_output.txt")
s100b = as.data.frame(s100b)


## Extract top SNP - rs8128872. We want all SNPs +200kb and -200kb from this region
chr <- s100b[which(s100b$rsid %in% "rs8128872"), "Chr"]
pos <- s100b[which(s100b$rsid %in% "rs8128872"), "Position"]
## subset to chr of interest
s100b_coloc <- s100b[s100b$Chr %in% chr, ]
## subset to positions of interest 
s100b_coloc <- s100b_coloc[which(s100b_coloc$Position > (pos-2e5)), ]
s100b_coloc <- s100b_coloc[which(s100b_coloc$Position < (pos+2e5)), ]
## check range
range(s100b_coloc$Position)
max(s100b_coloc$Position) - min(s100b_coloc$Position) ## not quite 400k but could be due to small no. of SNPs

## Extract the same region of SNPs in AD sum stats 

ad_coloc <- ad[ad$CHR %in% chr, ]
## subset to positions of interest 
ad_coloc <- ad_coloc[which(ad_coloc$BP > (pos-2e5)), ]
ad_coloc <- ad_coloc[which(ad_coloc$BP < (pos+2e5)), ]
## check range
range(ad_coloc$BP)
max(ad_coloc$BP) - min(ad_coloc$BP) ## not quite 400k but could be due to small no. of SNPs


## Clean the files 

## make sure MAF actually corresponds to MAF i.e. should be < 0.5 
ad_coloc$MAF <- ifelse(ad_coloc$MAF > 0.5, 1-ad_coloc$MAF, ad_coloc$MAF)
s100b_coloc$MAF <- ifelse(s100b_coloc$MAF > 0.5, 1-s100b_coloc$MAF, s100b_coloc$MAF)

## make sure MAF is numeric 
ad_coloc$MAF <- as.numeric(ad_coloc$MAF)
s100b_coloc$MAF <- as.numeric(s100b_coloc$MAF)

## trim rare variants in AD GWAS 
ad_coloc = ad_coloc[which(ad_coloc$MAF > 0.01),]

## Coloc doesn't handle NAs so trim these 

## MAF missigness 
if(length(which(is.na(s100b_coloc$MAF))) >= 1) { 
  s100b_coloc <- s100b_coloc[-which(is.na(s100b_coloc$MAF)),]
}

if(length(which(is.na(ad_coloc$MAF))) >= 1) { 
  ad_coloc <- ad_coloc[-which(is.na(ad_coloc$MAF)),]
} 


## Duplicated SNPs
if(length(which(duplicated(ad_coloc$SNP))) >= 1) { 
  ad_coloc <- ad_coloc[-which(duplicated(ad_coloc$SNP)),]
} 

if(length(which(duplicated(s100b_coloc$rsid))) >= 1) { 
  s100b_coloc <- s100b_coloc[-which(duplicated(s100b_coloc$rsid)),]
} 


## Missing SNPs
if(length(which(is.na(ad_coloc$SNP))) >= 1) { 
  ad_coloc <- ad_coloc[-which(is.na(ad_coloc$SNP)),]
} 

if(length(which(is.na(s100b_coloc$rsid))) >= 1) { 
  s100b_coloc <- s100b_coloc[-which(is.na(s100b_coloc$rsid)),]
} 


# > dim(ad_coloc)
# [1] 1346   14
# > dim(s100b_coloc)
# [1] 1010   14


## Coloc step 

dataset1 = list(snp = as.character(ad_coloc$SNP), N = as.numeric(ad_coloc$N), beta= as.numeric(ad_coloc$BETA), MAF = as.numeric(ad_coloc$MAF), varbeta = as.numeric(ad_coloc$SE)^2,  pvalues = as.numeric(ad_coloc$P), type = "quant", sdY =1 )
dataset2 = list(snp = as.character(s100b_coloc$rsid), N = as.numeric(s100b_coloc$N), pvalues = as.numeric(s100b_coloc$p), beta = as.numeric(s100b_coloc$beta), varbeta = as.numeric(s100b_coloc$se)^2, MAF = as.numeric(s100b_coloc$MAF), type = "quant", sdY= 1)
coloc = coloc.abf(dataset1, dataset2)


# Dataset 1 = Ad, dataset 2 = s100b
# > coloc = coloc.abf(dataset1, dataset2)
# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf
#  1.25e-10  1.37e-12  9.88e-01  1.08e-02  1.51e-03
# [1] "PP abf for shared variant: 0.151%"

# H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only),
# H3 (two distinct causal variants) and H4 (one common causal variant)

data <- coloc$summary 
data <- as.data.frame(data)

write.csv(data, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/MR_coloc_090921/AD.csv")

