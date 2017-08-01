source("X:/Labs/Chapple/0-Cole Wunderlich/R Files/Code/Methods/scaleData3.r")
pred = read.csv("X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/snyderRun4/OPLSDA Analysis/MasterPredictors_range_Standardization_12C_Only_OSCPLS.csv",check.names = F,stringsAsFactors = F)
fullData = read.csv("X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/snyderRun4/4cl1_2_3_WT_atomt_cadC-D_ccr1_fah1-2_med5ab_ref3-3_ref8-1_med5ab_tt4-2_valid_Clusters_tTestFilter_into.csv",check.names = F,stringsAsFactors = F)
phenoData = read.csv("X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/snyderRun4/phenoData.csv",row.names=1,stringsAsFactors = F,check.names = F)
pred = pred[!grepl("Global",pred$Predictive),]
# highpred = pred[grepl("High",pred$Predictive),]


dataCols = rownames(phenoData)[which(phenoData[,2]=="12C"&!grepl("[Ww][Tt]",phenoData[,1]))]
dataCols = which(colnames(fullData)%in%dataCols)

scaleddata = scaleData(fullData[,dataCols],byCol=F,"range")

#must match phenoData[1,]
pred$`Predictive`[grepl("4cl",pred$`Predictive`)] = "4cl1_2_3"
pred$`Predictive`[grepl("atomt",pred$`Predictive`)] = "atomt"
pred$`Predictive`[grepl("ref8-1",pred$`Predictive`)] = "ref8-1_med5ab"
pred$`Predictive`[grepl("cad",pred$`Predictive`)] = "cadC-D"
pred$`Predictive`[grepl("fah1-2",pred$`Predictive`)] = "fah1-2"
pred$`Predictive`[grepl("ccr1",pred$`Predictive`)] = "ccr1"
pred$`Predictive`[grepl("ref3-3",pred$`Predictive`)] = "ref3-3"
pred$`Predictive`[grepl("tt4-2",pred$`Predictive`)] = "tt4-2"
pred$`Predictive`[grepl("WT",pred$`Predictive`)] = "WT"
pred$`Predictive`[grepl("(?<!ref8-1)(.med5a-?5?b)",pred$`Predictive`,perl=T)] = "med5ab"

metabs = unique(pred$name)
genos = unique(pred$Predictive)

#create matrix of predictors and genotypes
predmatx =t(vapply(metabs,function (x,genos,pred){

  out = vector("double",length(genos))
  preddata = pred[which(pred$name==x),]
  for(geno in genos)
  {
    hit=geno%in%preddata$Predictive
    if(hit)
    {
      out[which(genos==geno)] = mean(as.vector(scaleddata[fullData$groupname==x,colnames(scaleddata)%in%rownames(phenoData)[phenoData[,1]==geno]],mode="double"))
    }else {
      out[which(genos==geno)] = 0
    }
  }
  # rownames(out)=x
  # names(out)=genos
  return(out)

},vector("double",length(genos)),genos,pred,USE.NAMES = T))

colnames(predmatx)=genos

write.csv(predmatx,"X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/snyderRun4/OPLSDA Analysis/clusteringDataAll_stdizedData.csv")

