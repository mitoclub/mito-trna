###################################

rm(list=ls(all=TRUE))

### harvest Gibbs Energies from individual files /../../Body/1Raw/5MitoTrnaDb/GibbsEnergyFromVictor/
AllFiles = list.files(path = "../../Body/1Raw/5MitoTrnaDb/GibbsEnergyFromVictor")
length(AllFiles)
Final=data.frame()
for (i in 1:length(AllFiles))
{ # i = 1
  # AllFiles[i]
  infile = paste("../../Body/1Raw/5MitoTrnaDb/GibbsEnergyFromVictor/",AllFiles[i],sep='')
  data = read.table(infile, header = FALSE,sep = '\t')
  id = unlist(strsplit(data$V1[1],"\\|"))[1]; id = gsub(">",'',id);
  species = unlist(strsplit(data$V1[1],"\\|"))[2]
  position = unlist(strsplit(data$V1[1],"\\|"))[3]
  trna = unlist(strsplit(data$V1[1],"\\|"))[4]
  anticodon = unlist(strsplit(data$V1[1],"\\|"))[5]
  sequence = data$V1[2]
  SecondaryStructure = unlist(strsplit(data$V1[3]," "))[1]
  GibbsEnergy = unlist(strsplit(data$V1[3]," \\("))[2]
  GibbsEnergy = gsub("\\(",'',GibbsEnergy); GibbsEnergy = gsub("\\)",'',GibbsEnergy)
  GibbsEnergy = as.numeric(GibbsEnergy)
  
  OneLine = c(id,species,position,trna,anticodon,sequence,SecondaryStructure,GibbsEnergy)
  Final = rbind(Final,OneLine)
}
names(Final)=c("id","species","parameter","trna","anticodon","sequence","SecondaryStructure","GibbsEnergy")
nrow(Final[is.na(Final$GibbsEnergy),]) # 617
table(Final$trna)
# Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile Leu1 Leu2  Lys  Met  Phe  Pro Ser1 Ser2  Thr  Trp  Tyr  Val 
# 282  282  282  282  282  282  282  282  282  282  282  282  259  282  282  282  282  282  282  282  282  282 
length(table(Final$trna))
write.table(Final, "../../Body/2Derived/02A.GibbsEnergyAsFunctionOfChainAndPositionInMammals.DeriveMitoTrnaDb.txt", quote = FALSE)



