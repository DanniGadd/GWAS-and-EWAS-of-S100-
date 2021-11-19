## eQTLGen Analyses - s100b

## load data.table library - helps read in large fastGWA files 
library(data.table)

## Read in eQTL Gen data 
eqtls = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Neurology Proteins/eqtls/eQTLgen_data.rds")

eqtls$BETA <- eqtls$Zscore/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))
eqtls$SE <- 1/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))

## Read in s100b data
s100b = fread("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/s100b_GWAS_output.txt")
s100b = as.data.frame(s100b)

snps1 <- s100b[which(s100b$rsid %in% "rs8128872"), ]

# ## install/load coloc library for colocalisation analyses 
# if(!require("remotes"))
#    install.packages("remotes") # if necessary
# library(remotes)
# install_github("chr1swallace/coloc")

library(coloc)

## Extract SNP 
tmp = snps1
tmp$Gene_of_Protein <- "S100B"
## Extract the SNPs coordinates 
position = as.numeric(tmp$Position)
chr = as.numeric(tmp$Chr)

## Subset eqtl file to the chr of the queried SNP 
tmp_eqtl <- eqtls[which(eqtls$SNPChr %in% chr), ]
## Subset eqtl file to region +/- 200 kb from queried SNP 
tmp_eqtl <- tmp_eqtl[which(tmp_eqtl$SNPPos >= (position - 2e5)),]
tmp_eqtl <- tmp_eqtl[which(tmp_eqtl$SNPPos <= (position + 2e5)),]

tmp_eqtl = tmp_eqtl[tmp_eqtl$GeneSymbol %in% tmp$Gene_of_Protein,]  


## Subset gwas file to chr + position 
fastGWA <- s100b[which(s100b$Chr %in% chr), ]
## Subset eqtl file to region +/- 200 kb from queried SNP 
fastGWA <- fastGWA[which(fastGWA$Position >= (position - 2e5)),]
fastGWA <- fastGWA[which(fastGWA$Position <= (position + 2e5)),]
## Get MAF - coloc needs this 
fastGWA$MAF <- as.numeric(fastGWA$MAF)
fastGWA$MAF <- ifelse(fastGWA$MAF > 0.5, 1-fastGWA$MAF, fastGWA$MAF)



## Clean up datasets 
## Missing MAF 
if(length(which(is.na(tmp_eqtl$MAF))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(is.na(tmp_eqtl$MAF)),]
} 

if(length(which(is.na(fastGWA$MAF))) >= 1) { 
  fastGWA <- fastGWA[-which(is.na(fastGWA$MAF)),]
} 

## Duplicated SNPs
if(length(which(duplicated(tmp_eqtl$SNP))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(duplicated(tmp_eqtl$SNP)),]
} 

if(length(which(duplicated(fastGWA$rsid))) >= 1) { 
  fastGWA <- fastGWA[-which(duplicated(fastGWA$rsid)),]
} 

## Missing SNPs
if(length(which(is.na(tmp_eqtl$SNP))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(is.na(tmp_eqtl$SNP)),]
} 

if(length(which(is.na(fastGWA$rsid))) >= 1) { 
  fastGWA <- fastGWA[-which(is.na(fastGWA$rsid)),]
} 

dataset1 = list(snp = as.character(tmp_eqtl$SNP), N = as.numeric(tmp_eqtl$NrSamples), beta = as.numeric(tmp_eqtl$BETA), varbeta = as.numeric(tmp_eqtl$SE)^2, pvalues = as.numeric(tmp_eqtl$Pvalue), MAF = as.numeric(tmp_eqtl$MAF), type = "quant", sdY =1 )
dataset2 = list(snp = as.character(fastGWA$rsid), N = as.numeric(fastGWA$N), pvalues = as.numeric(fastGWA$p), beta = as.numeric(fastGWA$beta), varbeta = as.numeric(fastGWA$se)^2, MAF = as.numeric(fastGWA$MAF), type = "quant", sdY= 1)
coloc = coloc.abf(dataset2, dataset1)

data <- coloc$summary 
data <- as.data.frame(data)

# Dataset 1 = eqtl, dataset 2 = s100b
# > coloc = coloc.abf(dataset2, dataset1)
# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf
#  0.00e+00  0.00e+00  1.27e-10  1.00e+00  3.38e-06
# [1] "PP abf for shared variant: 0.000338%"
# >

# H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only),
# H3 (two distinct causal variants) and H4 (one common causal variant)

# Save out 
write.csv(data, "/Cluster_Filespace/Marioni_Group/Danni/s100b_project/MR_coloc_090921/eQTL.csv")





#############

### POSITIONAL INFO LOOKUP

## Read in s100b data
s100b = fread("/Cluster_Filespace/Marioni_Group/Danni/s100b_project/PLINK_GWAS_060921/s100b_GWAS_output.txt")
s100b = as.data.frame(s100b)

snps1 <- s100b[which(s100b$rsid %in% "rs8128872"), ]

# Load in s100b gene info TSS 

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_circus/")

# Libraries
library(tidyverse)
library(readxl)

# Load in annotational information for the STRADL proteins
library(readxl)
anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/heritability_050321/outputs_combined/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,18,4,13)]
# comb <- left_join(osca, anno, by = "SeqId")

# Load in the significant EWAS results with proteins listed to get unique proteins
# cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_240521_results/pQTL_kept.csv")
# cpgs1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_240521_results/pQTL_WBC_covs_EWAS_thr_bonfer.csv")
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080821_full/results/no_eGFR_FULL_EWAS_thr_bonfer.csv")
unique <- unique(cpgs$SeqId)

length(unique)
# [1] 153

unique <- cpgs[c(1,10,11)] %>% unique() # xtracted into dataset with gene names 
tab <- unique
names(tab)[2] <- "gene" 

# get gene info set up from external source (Use Ensembl to Extract Relevant Features)
library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl, page = "structure")
attributes[grep("transcript", attributes$description, ignore.case = TRUE), ]

# get list of proteins gene names 
list <- tab$gene

list <- "S100B"

## Get Transcription Start Sites for soma proteins 
tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values = c(list),
             mart = ensembl)



## set up code to create biomart dataframe with gene, start and end columns
biomart_df <- tss

names(biomart_df)[8] <- "gene"
names(biomart_df)[3] <- "start"
names(biomart_df)[4] <- "end"
names(biomart_df)[2] <- "chromosome_ensembl"

# Again, lets do a check to see if each protein has its own TSS extracted that is unique 

subset <- biomart_df[,c("transcription_start_site", "gene")]
subset <- unique(subset$gene)

## separate data into genes on positive strand and genes on negative strand

biomart_df_pos <- biomart_df[biomart_df$strand == 1,]
biomart_df_neg <- biomart_df[biomart_df$strand == -1,]

# sanity check to see if overlap is present/correct number of genes are present

a = unique(biomart_df_pos$gene)
b = unique(biomart_df_neg$gene)
which(a %in% b)
length(a)
length(b)

## create output data frames for genes on positive strand and those on negative strand to be r-bound later
   
   
    ##out_df1 == genes on positive strand
   
out_df1 <- matrix(nrow=length(unique(biomart_df_pos$gene)), ncol=2)
colnames(out_df1) <- c("start", "end")
rownames(out_df1) <- unique(biomart_df_pos$gene)


    ##out_df2 == genes on negative strand
   
out_df2 <- matrix(nrow=length(unique(biomart_df_neg$gene)), ncol=2)
colnames(out_df2) <- c("start", "end")
rownames(out_df2) <- unique(biomart_df_neg$gene)

## loop through gene names to fill in min 5'end (start) and max 3'end (end) as on positive strand you want furthest away 5'(tss) and biggest distance to 3' end
## fill out out_df1 based on the above

for(gene in unique(biomart_df_pos$gene)) {
tmp <- biomart_df_pos[which(biomart_df_pos$gene==gene), ]
out_df1[gene,"start"] <-  min(tmp$start)
out_df1[gene,"end"] <-  max(tmp$end)
}

## loop through gene names to fill in min 3'end (end) and max 5'end (start) as on negative strand you want furthest away 3'(tss (relative to positive strand)) and biggest distance to 5' end
## fill out out_df2 based on the above

for(gene in unique(biomart_df_neg$gene)) {
tmp <- biomart_df_neg[which(biomart_df_neg$gene==gene), ]
out_df2[gene,"start"] <-  max(tmp$end)
out_df2[gene,"end"] <-  min(tmp$start)
}

## row bind the positive and negative strand dataframes to give one out_df
out_df <- rbind(out_df1, out_df2)


## Tidy the files before saving for new row with protein names (out_df), then one for length of transcripts in biomart_df for future reference

    ## set as data frames
    out_df <- as.data.frame(out_df)
    biomart_df <- as.data.frame(biomart_df)

out_df$Gene <- row.names(out_df)
biomart_df$transcript_length <- abs(biomart_df$start - biomart_df$end)

# Check unique start sites here 
unique(out_df$start) # 149
plot(out_df$start)

## Get Chromosome Number for Protein 
biomart_df$chromosome_ensembl2 <-  gsub("CHR_HSCHR", "", biomart_df$chromosome_ensembl)
biomart_df$chromosome_ensembl2 <- gsub("_.*", "", biomart_df$chromosome_ensembl2)
biomart_df <- biomart_df[,c(2,8,10)]
names(biomart_df) <- c("Chromosome_Biomarker", "Gene", "Chromosome_edit")
biomart_df = biomart_df[-which(duplicated(biomart_df$Gene)),]

# Join biomart_df file to the out_df file 
out_df <- merge(out_df, biomart_df, by = "Gene")

# Now join in the location of the SNP
out_df$SNP <- 48025097

# Now work out the TSS difference to the rsid 
out_df$diff <- out_df$SNP - out_df$Biomarker_TSS_Start