################################################
rm(list=ls(all=TRUE))
library(dplyr)
library(plyr)
library(tidyr)

codonusage = read.table("../../Body/2Derived/CodonUsageinMammalsbyTRNA.txt",sep="\t", header = TRUE)
remove = c('Start','Stop')
codonusage = codonusage[,!colnames(codonusage) %in% remove]
colnames(codonusage)

tbbsandgibbs = read.table("../../Body/2Derived/TRNASpeciesandGibbs.txt",sep="\t", header = TRUE)
remove = c('id','parameter','anticodon','sequence','SecondaryStructure','GenerationLength_d','LongLived')
tbbsandgibbs = tbbsandgibbs[,!colnames(tbbsandgibbs) %in% remove]
table(tbbsandgibbs$trna)

##### derive fraction of aminoacids (from total) specific for each "tRNA x species":
Final=data.frame()
for (i in 1:nrow(codonusage))
{ # i = 1
  OneSpecies = codonusage[i,]
  Species = OneSpecies$species
  OneSpecies = OneSpecies[,-1]
  OneSpecies = data.frame(t(OneSpecies))
  names(OneSpecies)=c('AaUsageNumber')
  OneSpecies$AaUsageFraction =  OneSpecies$AaUsageNumber/sum(OneSpecies$AaUsageNumber)
  OneSpecies$trna = row.names(OneSpecies)
  OneSpecies$species = Species
  Final=rbind(Final,OneSpecies)
}
# nrow(Final)/20 = number of species

nrow(Final)
#### merge CodonUsage with GibbsEnergy

nrow(tbbsandgibbs)
CodonUsageGibbs = merge(tbbsandgibbs,Final,by=c('species','trna'))
nrow(CodonUsageGibbs)

## species without codon usage in our Final dataset
setdiff(unique(tbbsandgibbs$species),unique(Final$species))

## species without GibbsEnergey (There are hundreds!!!!)
setdiff(unique(Final$species),unique(tbbsandgibbs$species))

## how many specise in our current merged dataset
length(unique(CodonUsageGibbs$species)) # 206
nrow(CodonUsageGibbs)/length(unique(CodonUsageGibbs$species))  # 17.9 < 20! some species do not have all 20 tRNAs (marsupials?)

## how many species have all tRNAs:
freq = data.frame(table(CodonUsageGibbs$species))
nrow(freq[freq$Freq == 18,]) # 186 < 206
freq[freq$Freq < 18,]$Var1 # ~ 20 potentil marsupiales + monotremata


###### ANALYSES:
#####
AveragePerTrna = aggregate(list(-CodonUsageGibbs$GibbsEnergy,CodonUsageGibbs$AaUsageFraction,CodonUsageGibbs$TimeBeingSingleStrangedForAll,CodonUsageGibbs$GeneCodedOnLightChain),by=list(CodonUsageGibbs$trna),FUN = median)
names(AveragePerTrna) = c('trna','MedianStability','MedianAaUsageFraction','TimeBeingSingleStrangedForAll','GeneCodedOnLightChain')
plot(AveragePerTrna$MedianAaUsageFraction,AveragePerTrna$MedianStability)
# text()

##### FOR ALL TRNAS TOGETHER:
summary(lm(-CodonUsageGibbs$GibbsEnergy ~ log2(CodonUsageGibbs$TimeBeingSingleStrangedForAll) + CodonUsageGibbs$GeneCodedOnLightChain + CodonUsageGibbs$AaUsageFraction))


##### FOR EACH TRNA SEPARATELY (REMINDER: some tRNAs show negative, other - positive correlations):
summary(lm(-GibbsEnergy ~ AaUsageFraction, data = CodonUsageGibbs[CodonUsageGibbs$trna == 'Ala',])) # pos...
summary(lm(-GibbsEnergy ~ AaUsageFraction, data = CodonUsageGibbs[CodonUsageGibbs$trna == 'Ala',]))
summary(lm(-GibbsEnergy ~ AaUsageFraction, data = CodonUsageGibbs[CodonUsageGibbs$trna == 'Ala',]))
summary(lm(-GibbsEnergy ~ AaUsageFraction, data = CodonUsageGibbs[CodonUsageGibbs$trna == 'Ala',]))
summary(lm(-GibbsEnergy ~ AaUsageFraction, data = CodonUsageGibbs[CodonUsageGibbs$trna == 'Ala',]))





