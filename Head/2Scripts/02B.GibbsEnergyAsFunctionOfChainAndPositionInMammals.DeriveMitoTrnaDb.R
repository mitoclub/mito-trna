###################################

rm(list=ls(all=TRUE))

data = read.table("../../Body/2Derived/02A.GibbsEnergyAsFunctionOfChainAndPositionInMammals.DeriveMitoTrnaDb.txt", header = TRUE)
names(data)

### 1 keep only mammals and derive Longlived (!!!!!! not ideal!!!!)
Mammals = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", sep = '\t', header = TRUE)
Mammals$species = gsub(' ','_',Mammals$Scientific_name)
Mammals = Mammals[names(Mammals) %in% c('species','GenerationLength_d')]
ThresholdGenLength = quantile(Mammals$GenerationLength_d,0.25) # this work better than 0.5 or 0.75 => don't know why
Mammals$LongLived = 1
for (i in 1:nrow(Mammals))  { if (Mammals$GenerationLength_d[i] <= ThresholdGenLength) {Mammals$LongLived[i] = 0} }
table(Mammals$LongLived)
nrow(data)
data = merge(data,Mammals, by = 'species')
nrow(data)/22 # ~ 432 species
data = data[order(data$species),]
table(data$trna) # Phe, Thr, Met and super rare!!!

#### 2 check if each mammalian species has 22 tRNA: !!!!!!!!!!!
SpeciesFreq = data.frame(table(data$species))
SpeciesFreq = SpeciesFreq[order(-SpeciesFreq$Freq),]
summary(SpeciesFreq$Freq)
hist(SpeciesFreq$Freq, breaks = 50)
Ideal = SpeciesFreq[SpeciesFreq$Freq == 22,]
nrow(Ideal) # 0 !!!!!! why we don't have species with standard set of tRNAs???

### 3 read annotation and merge with data
Annot = read.table("../../Body/2Derived/MammalianTrnaManualAnnotation.txt", sep = '\t', header = TRUE)
Annot$GeneCodedOnLightChain = as.numeric(as.factor(Annot$Chain))-1 # heavy = 0; light = 1
Annot = Annot[colnames(Annot) %in% c('trna','GeneCodedOnLightChain','TimeBeingSingleStrangedForAll')]
nrow(data)
data = merge(data,Annot, by = 'trna')
nrow(data)

### 4 average Gibbs Energy for each tRNA
agg = aggregate(data$GibbsEnergy, by = list(data$trna,data$GeneCodedOnLightChain,data$TimeBeingSingleStrangedForAll), FUN = median)
names(agg) = c('tRna','GeneCodedOnLightChain','TimeBeingSingleStrangedForAll','MedianGibbs')
agg = agg[order(agg$MedianGibbs),]

### 5 analysis with 22 values
summary(lm(-MedianGibbs ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll, data = agg))
# if GeneCodedOnLightChain == 1, it means that tRNA after transcription is identical to heavy chain => G and T rich (due to to common transitions: A>G, C>T)
# althougth non-significant, TimeBeingSingleStrangedForAll has positive coefficient, which is expected
summary(lm(-MedianGibbs ~ scale(GeneCodedOnLightChain) + scale(TimeBeingSingleStrangedForAll), data = agg)) # no interaction

### 6 analysis with whole dataset
summary(lm(-GibbsEnergy ~ scale(GeneCodedOnLightChain) + scale(TimeBeingSingleStrangedForAll), data = data)) # strong effect of all

# to do 1: introduce species factor into the model. glm 
# to do 2: introduce species as species-specific generation length 
# to do 3: PICs

summary(lm(-GibbsEnergy ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll + log2(GenerationLength_d), data = data)) # PIC
summary(lm(-GibbsEnergy ~ (GeneCodedOnLightChain + TimeBeingSingleStrangedForAll)*log2(GenerationLength_d), data = data))
summary(lm(-GibbsEnergy ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll + LongLived, data = data))  #### !!!!!!!!!!
summary(lm(-GibbsEnergy ~ scale(GeneCodedOnLightChain) + scale(TimeBeingSingleStrangedForAll) + scale(LongLived), data = data))  #### !!!!!!!!!!

summary(lm(-GibbsEnergy ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll*LongLived, data = data))

summary(lm(-GibbsEnergy ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll*log2(GenerationLength_d), data = data))

summary(lm(-scale(GibbsEnergy) ~ scale(GeneCodedOnLightChain) + scale(TimeBeingSingleStrangedForAll)*scale(log2(GenerationLength_d)), data = data))
summary(lm(-scale(GibbsEnergy) ~ 0 + scale(GeneCodedOnLightChain) + scale(TimeBeingSingleStrangedForAll)*scale(log2(GenerationLength_d)), data = data))

### 9 species by species
#VecOfSpecies = unique(data$species); length(VecOfSpecies)



### 10 drawing
#boxplot(data[data$trna == 'Pro',]$gibbs, data[data$trna == 'Tyr',]$gibbs, notch = TRUE, outline = FALSE, names = c('Pro','Tyr'))




