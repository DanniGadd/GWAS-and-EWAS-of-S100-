####################################################################################

### OSCA EWAS S100b

####################################################################################

setwd("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/")

## Read in methylation data  
meth = load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")
dim(dat)

## Read in target file - with all LBC information
target = read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

## Subset to LBC36 - excluding LBC21 individuals 
target = target[target$cohort %in% "LBC36",]

## Subset to the wave with S100B 
target = target[target$WAVE == 2,]
dim(dat)

## Read in target file - with all LBC information
target = read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

## Subset to LBC36 - excluding LBC21 individuals 
target = target[target$cohort %in% "LBC36",]

## Subset to the wave with S100B 
target = target[target$WAVE == 2,]

# Read in unprepped pheno file to get IDs for subsetting
phenos <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_EWAS_722.csv")

## Subset methylation data to just those in s100b file 
dat <- dat[,which(colnames(dat) %in% phenos$Basename)]

## Make sure IDs match up in both files 
ids = phenos$Basename 
dat <- dat[,match(ids, colnames(dat))]

# Match ID order between phenotype and meth files 
dat2 = dat[,match(as.character(phenos$Basename), colnames(dat))]
table(colnames(dat2) == phenos$Basename)

# Transpose methylation data
meth = t(dat2)

# Resdualise meth matrix 
for(i in 1:(ncol(meth))){
  print(i)
  meth[,i] <- resid(lm(meth[,i] ~ phenos$age + as.factor(phenos$sex.x) + as.factor(phenos$set) + as.factor(phenos$plate) + phenos$smokingScore + phenos$bmi_w2 + phenos$eosin + phenos$baso + phenos$lymph + phenos$neut + phenos$mono, na.action = na.exclude))
} 

meth <-as.data.frame(meth)

osca_dat <- data.frame(IID=rownames(meth), FID=rownames(meth))
osca_dat <- cbind(osca_dat, meth)

write.table(osca_dat, file="/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921.txt", row.names=F, sep=' ')

##########################################################################################################
            
## Make phenotype file using same phenotype file as reference as above 
## 3 columns: IID, FID, Phenotype
phenos <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_EWAS_722.csv")

pheno2 <- data.frame(FID = phenos$Basename,
                    IID = phenos$Basename,
                    pheno = phenos$s100b)
                    
#save out phenotype file.  
write.table(pheno2, file = "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/s100b.phen", col.names = F, row.names=F, sep=' ')    

########################################################################################################

# osca_dat <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921.txt")

## Must have to manually make an .opi annotation file - annotates the CpG data 
## Manual Annotation 
anno = readRDS("/Cluster_Filespace/Marioni_Group/Daniel/450k_annotation.rds")
opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
## get rid of 'chr' before chromosomes 
opi$chr <- gsub("chr", "", opi$chr)
## change X chromosome to 23 
opi$chr <- gsub("X", "23", opi$chr)
## change Y chromosome to 24 
opi$chr <- gsub("Y", "24", opi$chr)
## Change all to numeric 
opi$chr <- as.numeric(opi$chr)
## change gene names to factor variables
opi$UCSC_RefGene_Name <- as.factor(opi$UCSC_RefGene_Name)
## change the strand to variables 
opi$strand <- as.factor(opi$strand)
## In case there are blank gene names, set to NA 
opi[which(opi$UCSC_RefGene_Name==""), "UCSC_RefGene_Name"] <- NA
## write out result 
write.table(opi, file="/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921.opi", 
            col.names=F, 
            row.names=F, 
            quote=F, sep='\t')
                    
##########################################################################################################

# Make Binary Methylation File
osca_Linux --efile /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921.txt --methylation-m --make-bod --out /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921

head /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921.opi # annotations are present 

# Make ORM from methylation data 
# osca_Linux --befile /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921 --make-orm --out /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921

##########################################################################################

# BATCH EWAS 

cd /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/

osca_Linux --moment --befile /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/PREPS/methylation_722_060921 --pheno /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/s100b.phen --out /Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/s100b_MOMENT --methylation-beta

# Save results file 
result <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/s100b_MOMENT.mlma")
result <- as.data.frame(result)
result <- result[order(result$p),]
write.csv(result, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/result_s100b.csv", row.names = F)

##########################################################################################

# Do QC plots and check lambdas

# PROBEID, BETA, SE, and P_VAL - include as inputs 
setwd("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/")
data <- result
data <- data[c(2,6,7,8)]
names(data) <- c("PROBEID", "BETA", "SE", "P_VAL")

map <- result
map <- map[c(2,1,3)]
names(map) <- c("TARGETID", "CHR", "MAPINFO")

library(QCEWAS)

EWAS_QC(data = data,
map = map,
outputname = "QCEWAS")


# Define lambda function - can also be used from daniel 
# Needs A table with a column called "P.Value" or "P.value"

summstats <- result
summstats <- summstats[8]
names(summstats) <- "P.Value"

getlambda = function(x){
pvalue <- x[,grep("P.Value|P.value", colnames(x))]
chisq <- qchisq(1-pvalue,1)
# For z-scores as association, just square them
    # chisq <- data$z^2
        #For chi-squared values, keep as is
        #chisq <- data$chisq
    lambda = median(chisq)/qchisq(0.5,1)
    return(signif(lambda,3))
}
getlambda(summstats) # 0.941


source("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/qqunif.plot.R") 
names(summstats)[1] <- "pval"

pdf("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/OSCA_EWAS_060921/OSCA/QQ_EWAS.pdf")
qqunif.plot(summstats$pval)
dev.off()

##########################################################################################################

