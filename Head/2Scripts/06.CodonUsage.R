###################################
rm(list=ls(all=TRUE))
library(dplyr)
library(plyr)
library(tidyr)
data = read.table("../../Body/2Derived/AllGenesCodonUsageNoOverlap.txt", header = TRUE, fill = TRUE)
Mammals = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", sep = '\t', header = TRUE)
Mammals$species = gsub(' ','_',Mammals$Scientific_name)
Mammals = Mammals[13]
data = merge(data,Mammals, by = 'species')
colnames(data)
result = ddply(data,"species",numcolwise(sum))
remove = c('GeneEnd','CodonsToDeleteInTheBeginning','CodonsToDeleteAtTheEnd','NeutralA','NeutralT','NeutralG','NeutralC')
result = result[,!colnames(result) %in% remove]
result$Start = result$ATG
##remove = c('ATG') DO NOT REMOVE/ATG = MET
result = result[,!colnames(result) %in% remove]
result$Stop = result$AGA + result$AGG + result$TAA + result$TAG
remove = c('AGA','AGG','TAA','TAG')
result = result[,!colnames(result) %in% remove]
result$Ala = result$GCT + result$GCC + result$GCA + result$GCG 
remove = c('GCT','GCC','GCA','GCG')
result = result[,!colnames(result) %in% remove]
result$Arg = result$CGT + result$CGC + result$CGA + result$CGG
remove = c('CGT','CGC','CGA','CGG')
result = result[,!colnames(result) %in% remove]
result$Asn = result$AAT + result$AAC
remove = c('AAT','AAC')
result = result[,!colnames(result) %in% remove]
result$Asp = result$GAT + result$GAC
remove = c('GAT','GAC')
result = result[,!colnames(result) %in% remove]
result$Cys = result$TGT + result$TGC
remove = c('TGT','TGC')
result = result[,!colnames(result) %in% remove]
result$Gln = result$CAA + result$CAG 
remove = c('CAA','CAG')
result = result[,!colnames(result) %in% remove]
result$Glu = result$GAA + result$GAG
remove = c('GAA','GAG')
result = result[,!colnames(result) %in% remove]
result$Gly = result$GGT + result$GGC + result$GGA  + result$GGG
remove = c('GGT','GGC','GGA','GGG')
result = result[,!colnames(result) %in% remove]
result$His = result$CAT + result$CAC
remove = c('CAT','CAC')
result = result[,!colnames(result) %in% remove]
result$Ile = result$ATT + result$ATC
remove = c('ATT','ATC')
result = result[,!colnames(result) %in% remove]
result$Leu = result$CTT + result$CTC + result$CTA + result$CTG + result$TTA + result$TTG
remove = c('CTT','CTC','CTA','CTG','TTA','TTG')
result = result[,!colnames(result) %in% remove]
result$Lys = result$AAA + result$AAG
remove = c('AAA','AAG')
result = result[,!colnames(result) %in% remove]
result$Met = result$ATG + result$ATA
remove = c('ATG','ATA')
result = result[,!colnames(result) %in% remove]
result$Phe = result$TTT + result$TTC
remove = c('TTT','TTC')
result = result[,!colnames(result) %in% remove]
result$Pro = result$CCT + result$CCC + result$CCA + result$CCG
remove = c('CCT','CCC','CCA','CCG')
result = result[,!colnames(result) %in% remove]
result$Ser = result$TCT + result$TCC + result$TCA + result$TCG + result$AGT + result$AGC
remove = c('TCT','TCC','TCA','TCG','AGT','AGC')
result = result[,!colnames(result) %in% remove]
result$Thr = result$ACT + result$ACC + result$ACA + result$ACG
remove = c('ACT','ACC','ACA','ACG')
result = result[,!colnames(result) %in% remove]
result$Trp = result$TGG + result$TGA
remove = c('TGG','TGA')
result = result[,!colnames(result) %in% remove]
result$Tyr = result$TAT + result$TAC
remove = c('TAT','TAC')
result = result[,!colnames(result) %in% remove]
result$Val = result$GTT + result$GTC + result$GTA + result$GTG
remove = c('GTT','GTC','GTA','GTG')
result = result[,!colnames(result) %in% remove]
write.table(result, "../../Body/2Derived/CodonUsageinMammalsbyTRNA.txt",sep="\t", quote = FALSE, row.names=FALSE, col.names=TRUE)
