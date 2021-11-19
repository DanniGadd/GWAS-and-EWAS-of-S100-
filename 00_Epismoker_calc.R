
####################################################################################

### Calculate the epismoker measure for LBC1936

####################################################################################

source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")
install.packages("IlluminaHumanMethylation450kmanifest")

library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)

suppressPackageStartupMessages({
library(EpiSmokEr)  
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)
})

Read in target file with matched IDs for methylation data 
load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")
result2 <- epismoker(dataset=dat, method = "SSc")

Save out data scores for smoking
write.rds("/Cluster_Filespace/Marioni_Group/Danni/lbc_epismoker.rds") 

