## Run in command line 

## Loop through all 22 chromosomes 
## Swap in path to .ped file from preparation code 
## Saves a result for each chromosome - will have to concatenate later 

# First, make sure we can replicate the results Rob passed on when he left the lab - run the same code with same files but in my cluster area 
screen -S gwas 

cd /Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/pheno/
for chr in {1..22} 

do  

plink2 --allow-no-sex --double-id --mac 25 --glm allow-no-covars  --out outfile_DG_150721${chr} --pheno /Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/pheno/s100b_GWAS_prep.ped --vcf /GWAS_Source/GWAS_Data/LBC_36/HRC_v1.1/LBC1936_HRCv1.1_VCF/chr${chr}.dose.vcf.gz

done


## Open R, concatenate results 
#install.packages("data.table")
library(data.table) ## install.packages("data.table") if you don't have

# setwd("/Cluster_Filespace/Marioni_Group/Robert_Mc/s100b_GWAS/")
setwd("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/pheno/")

## Take out all results files - all 22 chromosomes 
loop = list.files(".", ".linear")


files = lapply(loop, read.table)

combined.df = do.call("rbind", files)

## Next we need to filter out any poorly imputed SNPs - Rsq < 0.6, these come from info files

setwd("/Cluster_Filespace/Marioni_Group/Rob/LBC_GWAS/HRC_INFO")

## Read the files in for all chromosomes and combine 
loop1 = list.files(".", ".")
files = lapply(loop1, fread)

files = do.call("rbind", files)
## Take out the SNPs which fail to show an Rsq > 0.6, a sign of decent imputation 
snps_to_remove = files[-which(files$Rsq > 0.6),]

gwas = combined.df[-which(combined.df$V3 %in% snps_to_remove$SNP),]

## Now we need to put in the rsid format of the SNPs, these are commonly used e.g. rs12345...
library(R.utils)
anno = as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Cohort_Summary_Stats/EasyQC/Allele_Freq_and_Mapping_Info/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map.gz"))
anno$MarkerName = paste0(anno$chr, ":", anno$pos)

## Make our gwas file fit the same format 
gwas$MarkerName <-sub('^([^:]+:[^:]+).*', '\\1', gwas$V3)

## match up the MarkerName column in each and put in rsid name 

gwas$rsid = anno[match(gwas$MarkerName, anno$MarkerName), "rsid"]

## Format column names 
gwas$MarkerName <- NULL 
gwas$V13 <- NULL 
names(gwas) <- c("Chr", "Position", "Marker", "REF", "ALT", "A1", "Test", "N", "beta", "se", "t_stat", "p", "rsid")

## Add in allele frequencies 

freq = as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/LBC_GWAS/LBC_allele_frequencies.txt"))

# Subset freq so that it has the GWAS markers included now that MAC has been filtered 
freq <- freq[which(freq$ID %in% gwas$Marker),]

## check that order of SNPs is the same 
table(freq2$ID == gwas$Marker)
## check that order of ALT allele is the same 
table(freq2$ALT == gwas$ALT)

## if both are TRUE, proceed to this step, otherwise align the orders of SNPs and ALT alleles if wrong
## Add in frequencies of Alternative Allele
gwas$Alt_Allele_Freq <- freq$ALT_FREQS

## A1 is the 'effect allele' - we want to get the frequency of this allele, however what we have
## at the moment is the 'alt allele' frequency 

## from this step 
table(gwas$ALT == gwas$A1)
## we can see that not Alt alleles are A1 alleles

## for those that ALT == A1, we already know their frequencies - Alt_Allele_Freq 
alt = gwas[which(gwas$ALT == gwas$A1), ]
alt$A1_Freq <- alt$Alt_Allele_Freq

## identify the other allele - A2, because all 'ALT' column == A1, we know A2 must be in 'REF' 
alt$A2 <- alt$REF 

## for those that don't match up, we can just take 1-Alt_Allele_Freq as this will be the allele frequency for 
## the other ('ref' allele)

ref = gwas[-which(gwas$ALT == gwas$A1), ]
ref$A1_Freq <- 1-ref$Alt_Allele_Freq

## identify the other allele - A2, because all 'REF' column == A1, we know A2 must be in 'ALT' 
ref$A2 <- ref$ALT 

## combine together 
gwas_freq = rbind(alt,ref)

## Restore original order of rsids 
ids = gwas$Marker 
gwas_freq = gwas_freq[match(ids, gwas_freq$Marker), ]
table(gwas$Marker == gwas_freq$Marker)

gwas_freq$Alt_Allele_Freq <- NULL

## Get minor allele frequency  
gwas_freq$MAF <- ifelse(gwas_freq$A1_Freq > 0.5, 1-gwas_freq$A1_Freq, gwas_freq$A1_Freq)

## Reorder columns 
gwas <- gwas_freq[,c(1,2,3,4,5,6,15,14,16,7,8,9,10,11,12,13)]
gwas$REF <- NULL 
gwas$ALT <- NULL 

## Write out gwas 
setwd("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/")
fwrite(gwas, "s100b_GWAS_output.txt", row.names = F)

write.csv(gwas, "s100b_output.csv", row.names = F)

## write out COJO table 

cojo <- gwas[,c("Marker", "A1", "A2", "A1_Freq", "beta", "se", "p", "N")]
names(cojo) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
write.table(cojo, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/COJO_table.txt", row.names = F, quote = F)

#see most significant SNP
gwas[which.min(gwas$p),]

#see where hits are 
gwas1 <- gwas[gwas$p < 5e-8,]
table(gwas1$Chr)

# Write out the top hits below the threshold for a suppl table 
write.csv(gwas1, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/top_hits.csv", row.names = F)


# ### PREP FUMA TABLE 
# # SNP, pval, A1, A2, Effect, SE and N
# tmp <- gwas 
# tmp ,- tmp[c()]
# names(tmp) = c("SNP", "P.value", "A1", "A2", "Effect", "StdErr", "TotalSampleSize")
# write.table(tmp, file=paste0(i, "_for_Fuma_13022021.txt"), sep='\t', quote=F, row.names=F)


# # Read in GWAS prior to MAC 25 filtering for quick comparison

# og <- read.delim("/Cluster_Filespace/Marioni_Group/Robert_Mc/s100b_GWAS/s100b_GWAS_output.txt", sep = ",")

# #see most significant SNP
# og[which.min(og$p),]

# #see where hits are 
# og1 <- og[og$p < 5e-8,]
# table(og1$Chr)



#################################################################

# Get lambdas and QQ

# Define lambda function - can also be used from daniel 
# Needs A table with a column called "P.Value" or "P.value"

summstats <- gwas
summstats <- summstats[13]
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
getlambda(summstats) # 0.991


source("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/qqunif.plot.R") 
names(summstats)[1] <- "pval"

pdf("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/QQ_GWAS.pdf")
qqunif.plot(summstats$pval)
dev.off()

####################################################

### COJO 

## Recode chromosome 21 to binary format - reference dataset 
plink19 --vcf /GWAS_Source/GWAS_Data/LBC_36/HRC_v1.1/LBC1936_HRCv1.1_VCF/chr21.dose.vcf.gz --make-bed --out /Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/COJO/test

## COJO step - change file paths to s100b directory 
gcta64  --bfile /Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/COJO/  --chr 21 --cojo-file /Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/COJO_table.txt --cojo-slct --cojo-actual-geno  --out /Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/COJO/cojo_file

cojo1 <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/COJO/cojo_file.ldr.cojo")

cojo2 <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/COJO/cojo_file.jma.cojo")


######################################################

### COLOCALISATIONS

# MR with one hit is ok just can't interpret much
# it would be very worth running coloc I think in this case, because it takes a 400kb region around your top hit and assesses shared causal signals between AD and s100b
# i'd do both, it helps to tease apart the results particularly because we only have one instrument for MR we can't test its assumptions which makes it hard to interpret. Coloc on the other hand we are ok to run and interpret
# Coloc is insanely simple to run and I could put together that code no problem
# It gives probabilities for 5 outcomes. it seems to be throwing out a high probability that in the locus 400 kb around the lead s100b variant, there is only a causal variant for s100b levels only and not for AD.





























