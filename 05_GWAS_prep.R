## In the .ped file, there are 3 columns:
## Columns 1 and 2: IID and FID - these are both the IDs for LBC participants i.e. same entry in both columns 
## normally FID is family id and IID is individual id, which is present in Generation Scotland - family study
## but this isn't present in LBC 

## Column 3: scaled phenotype i.e. s100b 

# getwd()
# setwd("M:/")

d <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/BayesR_T2/s100b_bayes_phenotypes_GWAS_769.csv")

out <- d[,c("IID", "FID", "s100b")]
head(out)

out$IID <- paste(out$IID, out$IID, sep= "_")
out$FID <- paste(out$FID, out$FID, sep="_")
head(out)

## Save .ped file, in the new directory you just created   
# getwd()
write.table(out, file="/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/pheno/s100b_GWAS_prep.ped", quote=F, col.names=F, row.names=F)




