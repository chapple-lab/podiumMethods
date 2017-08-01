pred = read.csv("X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/snyderRun4/OPLSDA Analysis/MasterPredictors_range_Standardization_12C_Only_OSCPLS.csv",check.names = F,stringsAsFactors = F)
#clean out global
pred = pred[!grepl("Global",pred$Predictive),]
highpred = pred[grepl("High",pred$Predictive),]


#replace names
highpred$`Predictive`[grepl("4cl",highpred$`Predictive`)] = "4cl1,2,3"
highpred$`Predictive`[grepl("atomt",highpred$`Predictive`)] = "atomt"
highpred$`Predictive`[grepl("ref8-1",highpred$`Predictive`)] = "ref8-1 med5a-5b"
highpred$`Predictive`[grepl("cad",highpred$`Predictive`)] = "cadC-D"
highpred$`Predictive`[grepl("fah1-2",highpred$`Predictive`)] = "fah1-2"
highpred$`Predictive`[grepl("ccr1",highpred$`Predictive`)] = "ccr1"
highpred$`Predictive`[grepl("ref3-3",highpred$`Predictive`)] = "ref3-3"
highpred$`Predictive`[grepl("tt4-2",highpred$`Predictive`)] = "tt4-2"
highpred$`Predictive`[grepl("WT",highpred$`Predictive`)] = "WT"
highpred$`Predictive`[grepl("(?<!ref8-1)(.med5a-?5?b)",highpred$`Predictive`,perl=T)] = "med5a-5b"
# highpred$`Predictive`[grepl("(?<!ref8-1)[, ]med5a-5b",highpred$`Predictive`,perl=T)] = "med5a-5b"
# highpred$`Predictive`[grepl("(?<!ref8-1)*med5ab",highpred$`Predictive`,perl=T)] = "med5a-5b"


#binary predictors
# predin = highpred[,c("name","Comp1","Predictive")]
# predin$Comp1 = as.numeric(predin$Comp1<0)
#up, down, not there 1,-1,0

metabs = unique(highpred$name)
genos = unique(highpred$Predictive)

finaldata =t(vapply(metabs,function (x,genos,highpred){

  out = vector("numeric",length(genos))
  preddata = highpred[which(highpred$name==x),]
  for(geno in genos)
  {
    hit=geno%in%preddata$Predictive
    if(hit)
    {
      if(preddata$Comp1[preddata$Predictive==geno]<0)
        out[which(genos==geno)] = 1
      if(preddata$Comp1[preddata$Predictive==geno]>0)
        out[which(genos==geno)] = -1
    }else {
      out[which(genos==geno)] = 0
    }
  }
  # rownames(out)=x
  names(out)=genos
  return(out)

},vector("numeric",length(genos)),genos,highpred,USE.NAMES = T))
colnames(finaldata)=genos

#check to make sure its working
finaldata[1:10,1:9]
pred[which(pred$name==pred$name[1]),]
pred[which(pred$name==pred$name[10]),]

write.csv(finaldata,"X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/snyderRun4/OPLSDA Analysis/clusteringDataHighOnly.csv")

