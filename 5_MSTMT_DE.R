rm(list = ls())
library(MSstatsTMT)

setwd('/mnt/8tb/DARPA/Proteomics_R2')

input.pd_1 <- read.csv('./IRS/MSTMT_Brain.csv')

quant.median_1<-list(ProteinLevelData=input.pd_1,FeatureLevelData=input.pd_1)

test.pairwise <- groupComparisonTMT(quant.median_1, moderated = TRUE)
write.csv(test.pairwise$ComparisonResult,'./DE/Brain.csv', row.names = FALSE)

input.pd_1 <- read.csv('./IRS/MSTMT_Kidney.csv')
quant.median_1<-list(ProteinLevelData=input.pd_1,FeatureLevelData=input.pd_1)
test.pairwise <- groupComparisonTMT(quant.median_1, moderated = TRUE)
write.csv(test.pairwise$ComparisonResult,'./DE/Kidney.csv', row.names = FALSE)

