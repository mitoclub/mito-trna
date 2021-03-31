###################################

rm(list=ls(all=TRUE))

data = read.table("../../Body/2Derived/harvest.txt", header = TRUE)
names(data)

### 1 keep only 37
table(data$temp)
nrow(data) # 243476
data = data[data$temp == 37,]
nrow(data) # 81157

### 2 keep only normal tRNA
length(table(data$trna)) # 23 becaus there is 'local' ???
data[data$trna == 'local',]
nrow(data)
data = data[data$trna != 'local',]
nrow(data)

### 3 keep only tRNA with estimated Gibbs Energy (not NA)
str(data)
nrow(data)
data = data[!is.na(data$gibbs),]
nrow(data)
summary(data$gibbs)

### 4 keep only mammals (!!!!!! not ideal!!!!)
Mammals = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", sep = '\t', header = TRUE)
Mammals$species = gsub(' ','_',Mammals$Scientific_name)
Mammals = Mammals[names(Mammals) %in% c('species','GenerationLength_d')]
nrow(data)
data = merge(data,Mammals, by = 'species')
nrow(data)/22 # ~ 432 species

### 5 read annotation and merge with data
Annot = read.table("../../Body/2Derived/MammalianTrnaManualAnnotation.txt", sep = '\t', header = TRUE)
Annot$GeneCodedOnLightChain = as.numeric(as.factor(Annot$Chain))-1 # heavy = 0; light = 1
Annot = Annot[colnames(Annot) %in% c('trna','GeneCodedOnLightChain','TimeBeingSingleStrangedForAll')]
nrow(data)
data = merge(data,Annot, by = 'trna')
nrow(data)

### 6 average Gibbs Energy for each tRNA
agg = aggregate(data$gibbs, by = list(data$trna,data$GeneCodedOnLightChain,data$TimeBeingSingleStrangedForAll), FUN = median)
names(agg) = c('trna','GeneCodedOnLightChain','TimeBeingSingleStrangedForAll','MedianGibbs')
agg = agg[order(agg$MedianGibbs),]

### 7 analysis with 22 values
summary(lm(-MedianGibbs ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll, data = agg))
# if GeneCodedOnLightChain == 1, it means that tRNA after transcription is identical to heavy chain => G and T rich (due to to common transitions: A>G, C>T)
# althougth non-significant, TimeBeingSingleStrangedForAll has positive coefficient, which is expected

summary(lm(-MedianGibbs ~ GeneCodedOnLightChain*TimeBeingSingleStrangedForAll, data = agg)) # no interaction

### 8 analysis with whole dataset
summary(lm(-gibbs ~ GeneCodedOnLightChain+TimeBeingSingleStrangedForAll, data = data)) # strong effect of all

# to do 1: introduce species factor into the model. glm 
# to do 2: introduce species as species-specific generation length 
# to do 3: PICs

summary(lm(-gibbs ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll + GenerationLength_d, data = data)) # PIC
summary(lm(-gibbs ~ (GeneCodedOnLightChain + TimeBeingSingleStrangedForAll)*GenerationLength_d, data = data))

#  Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)                                       1.299e+01  1.037e-01 125.361  < 2e-16 ***
#  GeneCodedOnLightChain                             1.229e+00  9.487e-02  12.959  < 2e-16 ***
#  TimeBeingSingleStrangedForAll                     5.241e-05  7.589e-06   6.906 5.31e-12 ***
#  GenerationLength_d                                9.767e-05  2.985e-05   3.272  0.00107 ** 
#  GeneCodedOnLightChain:GenerationLength_d         -4.071e-06  2.687e-05  -0.152  0.87957    
#  TimeBeingSingleStrangedForAll:GenerationLength_d -1.052e-08  2.176e-09  -4.836 1.35e-06 ***

summary(lm(-gibbs ~ GeneCodedOnLightChain + TimeBeingSingleStrangedForAll*GenerationLength_d, data = data))
#  Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)                                       1.300e+01  9.755e-02 133.270  < 2e-16 ***
#  GeneCodedOnLightChain                             1.218e+00  5.658e-02  21.524  < 2e-16 ***
#  TimeBeingSingleStrangedForAll                     5.240e-05  7.589e-06   6.905 5.33e-12 ***
#  GenerationLength_d                                9.577e-05  2.711e-05   3.533 0.000412 ***
#  TimeBeingSingleStrangedForAll:GenerationLength_d -1.052e-08  2.176e-09  -4.834 1.36e-06 ***

# -gibbs ~ +GeneCodedOnLightChain + +TimeBeingSingleStrangedForAll  +GenerationLength_d -TimeBeingSingleStrangedForAll*GenerationLength_d
# -gibbs ~ +TimeBeingSingleStrangedForAll + GenerationLength_d -TimeBeingSingleStrangedForAll*GenerationLength_d
# -gibbs ~ +TimeBeingSingleStrangedForAll # if short lived => GT = 0
# -gibbs ~ +BIG*TimeBeingSingleStrangedForAll + BIG+GenerationLength_d -SMALL*TimeBeingSingleStrangedForAll*GenerationLength_d # if long lived

summary(lm(-scale(gibbs) ~ scale(GeneCodedOnLightChain) + scale(TimeBeingSingleStrangedForAll)*scale(GenerationLength_d), data = data))
summary(lm(-scale(gibbs) ~ 0 + scale(GeneCodedOnLightChain) + scale(TimeBeingSingleStrangedForAll)*scale(GenerationLength_d), data = data))

### 9 continue till all 22 (color by chain)
#boxplot(data[data$trna == 'Pro',]$gibbs, data[data$trna == 'Tyr',]$gibbs, notch = TRUE, outline = FALSE, names = c('Pro','Tyr'))




