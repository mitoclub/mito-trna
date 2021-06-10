################################################
rm(list=ls(all=TRUE))
library(dplyr)
library(plyr)
library(tidyr)
codonusage = read.table("../../Body/2Derived/CodonUsageinMammalsbyTRNA.txt",sep="\t", header = TRUE)
remove = c('Start','Stop')
codonusage = codonusage[,!colnames(codonusage) %in% remove]
tbbsandgibbs = read.table("../../Body/2Derived/TRNASpeciesandGibbs.txt",sep="\t", header = TRUE)
remove = c('id','parameter','anticodon','sequence','SecondaryStructure','GenerationLength_d','LongLived')
tbbsandgibbs = tbbsandgibbs[,!colnames(tbbsandgibbs) %in% remove]
