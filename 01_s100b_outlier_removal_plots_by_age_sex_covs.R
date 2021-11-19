####################################################################################

### s100b processing 

####################################################################################

# Read in data 
library(tidyverse)
library(haven)
s1 <- read_sav("Y:/Danni/s100b_project/LBC1936_S100b_WithGeneticEpigeneticCognitiveNeuroImagingPhenotypes_RM_DG_02NOV2020.sav")

s2 <- s1 %>% as.data.frame()

# Create a summary table with key info about the dataset 
names <- colnames(s2)
output <- matrix(nrow = length(names), ncol = 2)
output <- as.data.frame(output)
names(output) <- c("Variable", "Count")
output$Variable <- names

for(i in 1:length(names)){
	variable <- as.character(names[i])
	df <- s2[variable]
	complete <- na.omit(df)
	count <- nrow(complete)
	output[i,2] <- count
}

# write out initial variable counts 
write.csv(output, "Y:/Danni/s100b_project/summary_of_counts.csv", row.names = F)


# Plot s100b levels to identify outliers across waves 
data <- s2
data <- data[c("lbc36no", "sex", "agedays_w2", "bld_S100_w2", "bld_S100_w3", "bld_S100B_w4")]

W2 <- data[c("lbc36no", "bld_S100_w2")]
W2$wave <- 2
names(W2) <- c("lbc36no", "s100b", "wave")
W2$s100b <- as.numeric(W2$s100b)

W3 <- data[c("lbc36no", "bld_S100_w3")]
W3$wave <- 3
names(W3) <- c("lbc36no", "s100b", "wave")
W3$s100b <- as.numeric(W3$s100b)

W4 <- data[c("lbc36no", "bld_S100B_w4")]
W4$wave <- 4
names(W4) <- c("lbc36no", "s100b", "wave")
W4$s100b <- as.numeric(W4$s100b)

join <- rbind(W2, W3)
join <- rbind(join, W4)
join$wave <- as.character(join$wave)

pdf("Y:/Danni/s100b_project/s100b_raw_data_plotted_for_each_wave.pdf")
ggplot(join, aes(x=wave, y=s100b)) +
geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
geom_boxplot(width=0.1) + theme_minimal()
dev.off()

# Point to consider:

# Extreme outlier points greater than 4 sd above mean excluded (w2 = 6, W3 = 4)
# 7 points at W3 below sensitivity assay threshold - < 0.02 ug/L
# 776 and 619 - for waves 2 and 3, ages 73 and 76 years 

data <- s2
W2 <- data[c("lbc36no", "sex", "agedays_w2", "bld_S100_w2")]
W2$wave <- 2
names(W2) <- c("lbc36no", "sex", "age_W2", "s100b", "wave")
W2$s100b <- as.numeric(W2$s100b)
W2 <- na.omit(W2) # 834 from W2 without exclusion of anyone 

write.csv(W2, "Y:/Danni/s100b_project/W2_s100b_n_834_no_exclusions_full_set.csv", row.names = F)


# Remove outliers from wave 2
mean <- mean(W2$s100b)
sd <- sd(W2$s100b)
lower <- mean - (4 * sd)
upper <- mean + (4 * sd)

outliers <- which(W2$s100b > upper)
out <- W2[outliers,]
write.csv(out,  "Y:/Danni/s100b_project/W2_s100b_6_outliers_above_4_sd_from_mean.csv", row.names = F) # saved 6 points considered to be outliers (matches simons paper)

lower_outliers <- which(W2$s100b < lower) # no outliers on the lower end 

non_outliers <- which(W2$s100b < upper)
int <- W2[non_outliers,]
write.csv(int,  "Y:/Danni/s100b_project/W2_s100b_6_outliers_above_4_sd_from_mean_removed.csv", row.names = F) # 824 remaining data points with outliers removed 


# Replot with outliers removed 
pdf("Y:/Danni/s100b_project/s100b_W2_outliers_removed.pdf")
ggplot(int, aes(x=wave, y=s100b)) +
geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
geom_boxplot(width=0.1) + theme_minimal()
dev.off()

#################################################################################

# Examine s100b by covariates in W2
# Age at W2 and sex 

library(ggpubr)

setwd("Y:/Danni/s100b_project/Submission_prep/Plots/Age_sex/")

int$age <- int$age_W2 / 325.25

pdf("s100b_W2_outliers_removed_plotted_by_sex_V2.pdf")
ggplot(int, aes(x=as.character(sex), y=s100b)) +
  geom_violin(trim=FALSE, fill='cyan4', color="darkcyan")+
  geom_boxplot(width=0.1) + theme_minimal() + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + ylab("S100β Wave 2") + xlab("Sex") + 
  scale_x_discrete(breaks=c("1", "2"),
        labels=c("Male", "Female"))
dev.off()

fem <- int[which(int$sex == 1),]
mean1 <- mean(fem$s100b)

mal <- int[which(int$sex == 2),]
mean2 <- mean(mal$s100b)

mod <- lm(scale(int$s100b) ~ scale(int$sex))

library(ggpubr)
pdf("s100b_W2_outliers_removed_plotted_by_age_at_W2_V2.pdf")
ggplot(int, aes(x=age, y=s100b)) +
geom_point(color = "darkgoldenrod2") + geom_smooth(method = 'lm', color = "darkgoldenrod3") + theme_minimal() + 
   stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + ylab("S100β Wave 2") + xlab("Age")
dev.off()



############################################################################

### TEST ASSOCIATIONS WITH SMOKING AND BMI AND PLOT 

setwd("Y:/Danni/s100b_project/Submission_prep/Plots/Age_sex/")

EWAS <- read.csv("s100b_bayes_phenotypes_EWAS_781.csv")

df <- EWAS

library(ggpmisc)
library(ggpubr)  

pdf("s100b_W2_outliers_removed_plotted_by_smokingScore_V2.pdf")
model1=lm(s100b~smokingScore, df) # no association
ggplot(df, aes(x=smokingScore, y=s100b)) + geom_point(color = "dodgerblue2") + geom_smooth(method=lm, color = "dodgerblue2", se=TRUE, level=0.95) +
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7) + theme_minimal() + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + ylab("S100β Wave 2") + xlab("EpiSmokEr")
dev.off()

pdf("s100b_W2_outliers_removed_plotted_by_BMI_V2.pdf")
model1=lm(s100b~bmi_w2, df) # an association
ggplot(df, aes(x=bmi_w2, y=s100b)) + geom_point(color = "firebrick2") + geom_smooth(method=lm, color = "firebrick2", se=TRUE, level=0.95) +
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7) + theme_minimal() + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 14), 
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + ylab("S100β Wave 2") + xlab("BMI")
dev.off()









































