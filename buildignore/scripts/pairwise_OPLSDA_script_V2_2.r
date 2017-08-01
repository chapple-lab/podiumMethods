#OPLS-DA test script, uses devium functions that are built off of pls package
# source("http://pastebin.com/raw.php?i=UyDBTA57")
##Note: The following is setup to deal with 12 samples total (2 genos with 3 reps each)
#Requires that WT be in the phenoData for the WT sample (the one everything is pairwise compared to)
#Uses global CAMERA data for all individual pairwise comparisons when full data output is turned on.
#Always uses pairwise CAMREA data in the masterFile output when full data output is turned on.
source('X:/Labs/Chapple/0-Cole Wunderlich/R files/Devium OSC-PLS/Devium PLS and OPLS.r')
source("X:/Labs/Chapple/0-Cole Wunderlich/R Files/Code/Methods/scaleData3.r")
source('X:/Labs/Chapple/0-Cole Wunderlich/R files/Devium OSC-PLS/Modified Plot Function.R') #Overides Default Plot Function

plots=F #currently not supported for pairwise comparison
C12=T
camera12C=T #was camera run in 12C mode? (note: NOT C12Compat mode but actual 12C mode)
reshape=T
normType="range"
masterFile=T
global=T
flag=F
fullDataOut=F

# resultsPath="X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/GWAS_Comparison_20sec_PooledVariance_NoPrefilter"
# resultsPath = "X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/GWAScomparison_20sec_highResp"

if(!file.exists(file.path(resultsPath,"OPLSDA Analysis")))
  dir.create(file.path(resultsPath,"OPLSDA Analysis"))

phenoData = read.csv(file=file.path(resultsPath,"phenoData.csv"),row.names=1,stringsAsFactors = F,check.names=F)
if(camera12C)
{
  fullData = read.csv(file=file.path(resultsPath,"CAMERA Analysis","12C_CAMERA_FinalResults_MassAnnosAndHits.csv"),check.names=F,stringsAsFactors=F)
  pairwiseAnnos = read.csv(file=file.path(resultsPath,"CAMERA Analysis","12C_CAMERA_AllAnnotations_AllGenotypesPairedWithWT.csv"),check.names=F,stringsAsFactors = F)
}else
{
  fullData = read.csv(file=file.path(resultsPath,"CAMERA Analysis","CAMERA_FinalResults_MassAnnosAndHits.csv"),check.names=F,stringsAsFactors=F)
  pairwiseAnnos = read.csv(file=file.path(resultsPath,"CAMERA Analysis","CAMERA_AllAnnotations_AllGenotypesPairedWithWT.csv"),check.names=F,stringsAsFactors = F)
}

resultsPath2 = file.path(resultsPath,"OPLSDA Analysis")

##ID desired phenotypes
wtDataCols = grepl("[Ww][Tt]",phenoData[,1]) #change for desired base phenotype
if(C12)
{
  phenos=phenoData[which(phenoData[,2]=="12C"),1]
  wtDataCols = which(colnames(fullData)%in%rownames(phenoData)[phenoData[,2]=="12C"&wtDataCols])
} else
{
  phenos=phenoData[,1]
  wtDataCols = which(colnames(fullData)%in%rownames(phenoData)[which(wtDataCols)])
}

LoopPhenos = unique(phenos)
LoopPhenos = LoopPhenos[!grepl("[Ww][Tt]",LoopPhenos)] #remove WT from pheno list

if(global)
  LoopPhenos = c(LoopPhenos,"Global")

if(masterFile)
{
  masterScores = NULL
  masterLoadings = NULL
  masterPredictors = NULL
}

for(pheno in LoopPhenos)
{
  cat("\nProcessing Phenotype: ",pheno,"\n")
  if(pheno=="Global")
  {
    pheno=LoopPhenos[!LoopPhenos=="Global"]
    flag=T
  }


  if(C12)
    dataCols = rownames(phenoData)[which(phenoData[,2]=="12C"&phenoData[,1]%in%pheno)] else
      dataCols = rownames(phenoData)[phenoData[,1]%in%pheno]

  dataCols = which(colnames(fullData)%in%dataCols)
  dataCols = c(dataCols,wtDataCols) #order of arguments is important!!
  if(reshape)
  {
    testData = fullData[ ,dataCols]
    testData = do.call(rbind,testData[,1:ncol(testData)])
    colnames(testData) = fullData$groupname
    #generate numeric codes
    count=1
    numCodes=NULL
    for(x in pheno)
    {
      if(C12)
      numSamps = length(which(phenoData[,2]=="12C"&phenoData[,1]==x)) else
        numSamps=length(which(phenoData[,1]%in%x)) #number of occurences of pheno

      numCodes = c(numCodes,rep(count,numSamps))
      count=count+1
    }
    numCodes = c(numCodes,rep(count,length(wtDataCols))) # WT order important??
    testData = cbind(numCodes,testData) #add numeric codes
    testData = data.frame(testData)
  }


  pls.data = testData[,-1]
  pls.y = testData[,1]

  groupLabels = do.call(rbind,list(c(phenos[phenos%in%pheno],phenos[grepl("[Ww][Tt]",phenos)])))

  scaledData = scaleData(pls.data,byCol=T,normType)

  # opls.results<-make.OSC.PLS.model(pls.y,pls.data,
  #                                  comp=2,
  #                                  OSC.comp=1,
  #                                  validation = "LOO",
  #                                  cv.scale = TRUE,
  #                                  train.test.index=NULL,
  #                                  progress=FALSE)

  results = make.OSC.PLS.model(pls.y,scaledData,
                               comp=2,
                               OSC.comp=1,
                               validation = "none",
                               cv.scale = F,
                               method = "oscorespls",
                               train.test.index=NULL,
                               progress=F)

  if(plots)
  {
    plot.OSC.results(results,plot="summary",groups=groupLabels,lgroups=fullData$iso,alpha=0.5)
    plot.OSC.results(results,plot="multiloadings")
    plot.OSC.results(results,plot="scores",groups=pls.y,nOSCcomp=1)
  }

  if(C12)
    pheno = paste(pheno,".12C",sep="")

  if(flag)
    pheno="Global"

  data = results$loadings[[2]]
  data = cbind(as.numeric(data[,1]),as.numeric(data[,2]))
  data= as.data.frame(data)
  data = cbind(fullData$groupname,fullData$iso,fullData$clusterNumber,data)
  colnames(data) = c("name","iso","clusterNumber","Comp1","Comp2")

  #get predictive vars. High predictive is +/- 1/2 stdev,  Med predictive is +/- 4/5 stdev
  Predictability = vector("character",length=nrow(fullData))
  stdev = sd(data$Comp1)
  minLoading = min(data$Comp1)
  maxLoading = max(data$Comp1)
  highPredVars = which(data$Comp1<=minLoading+stdev/2|data$Comp1>=maxLoading-stdev/2)
  medPredVars = which(data$Comp1<=minLoading+stdev/1.25|data$Comp1>=maxLoading-stdev/1.25)
  medPredVars = setdiff(medPredVars,highPredVars)

  Predictability[highPredVars] = "High"
  Predictability[medPredVars] = "Medium"
  Predictability[-c(medPredVars,highPredVars)] = "None"
  predictors = union(highPredVars,medPredVars)
  if(fullDataOut)
  {
    predData = cbind(fullData[predictors,1:7],Predictive=Predictability[predictors],Component_1=data$Comp1[predictors],fullData[predictors,-(1:7)])
    data = cbind(fullData[,1:7],Predictive=Predictability,data$Comp1,fullData[ ,-(1:7)])
  } else
  {
    predData = cbind(data[predictors,],Database_Matches=fullData[predictors,"Database_Matches"],Predictive=Predictability[predictors])
    data=cbind(data,Predictive=Predictability)
  }


  if(C12)
    write.csv(predData,file=file.path(resultsPath2,paste(pheno,"vs_WT",normType,"Standardization_12C_Only_OSCPLS_Predictors.csv",sep="_")),row.names=F) else
      write.csv(predData,file=file.path(resultsPath2,paste(pheno,"vs_WT",normType,"Standardization_OSCPLS_Predictors.csv",sep="_")),row.names=F)


  if(C12)
    write.csv(data,file=file.path(resultsPath2,paste(pheno,"vs_WT",normType,"Standardization_12C_Only_OSCPLS_Loadings.csv",sep="_")),row.names=F) else
      write.csv(data,file=file.path(resultsPath2,paste(pheno,"vs_WT",normType,"Standardization_OSCPLS_Loadings.csv",sep="_")),row.names=F)


  if(masterFile)
  {
    if(fullDataOut) #load pairwise data for master files
    {
      pwaIdx= which(pairwiseAnnos$Type==pheno)
      data= cbind(data[,1:2],pairwiseAnnos[pwaIdx,1:5],data[,-(1:6)])
      predData= cbind(predData[,1:2],pairwiseAnnos[pwaIdx[predictors],1:5],predData[,-(1:2)])
    } else	
    {
      data = cbind(Type=rep(pheno,nrow(data)),data)
      predData = cbind(predData[,-ncol(predData)],Predictive=paste(predData$Predictive,pheno,sep="_"))
    }
    masterLoadings = rbind(masterLoadings,data)
    masterPredictors = rbind(masterPredictors,predData)
  }

  data = results$scores[[2]]
  data = cbind("Sample Name"=colnames(fullData)[dataCols],data)
  if(C12)
    write.csv(data,file=file.path(resultsPath2,paste(pheno,"vs_WT",normType,"Standardization_12C_Only_OSCPLS_Scores.csv",sep="_")),row.names=F) else
      write.csv(data,file=file.path(resultsPath2,paste(pheno,"vs_WT",normType,"Standardization_OSCPLS_Scores.csv",sep="_")),row.names=F)
  if(masterFile)
  {
    data = cbind(Type=rep(pheno,nrow(data)),data)
    masterScores = rbind(masterScores,data)
  }

  if(plots)
    plotLoadingClusters(Group2results,clusters=fullData$clusterNumber,isos=fullData$iso,alpha=1/3,path="X:/Labs/Chapple/0-Cole Wunderlich/R files/Devium OSC-PLS/Loadings By Cluster Iso Coded",ftype=".png")

}

if(masterFile)
{
  if(C12)
  {
    write.csv(masterScores,file=file.path(resultsPath2,paste("MasterScores",normType,"Standardization_12C_Only_OSCPLS.csv",sep="_")),row.names=F)
    write.csv(masterLoadings,file=file.path(resultsPath2,paste("MasterLoadings",normType,"Standardization_12C_Only_OSCPLS.csv",sep="_")),row.names=F)
    write.csv(masterPredictors,file=file.path(resultsPath2,paste("MasterPredictors",normType,"Standardization_12C_Only_OSCPLS.csv",sep="_")),row.names=F)
  }else
  {
    write.csv(masterScores,file=file.path(resultsPath2,paste("MasterScores",normType,"Standardization_OSCPLS.csv",sep="_")),row.names=F)
    write.csv(masterLoadings,file=file.path(resultsPath2,paste("MasterLoadings",normType,"Standardization_OSCPLS.csv",sep="_")),row.names=F)
    write.csv(masterPredictors,file=file.path(resultsPath2,paste("MasterPredictors",normType,"Standardization_OSCPLS.csv",sep="_")),row.names=F)
  }
}
cat("\n-----OPLSDA Analysis Complete-----\n")


###12C data only
# testData = testData[c(-4,-5,-6,-10,-11,-12),]
# groupLabels = rbind('ref3-3','ref3-3','ref3-3','WT','WT','WT')
# write.csv(data,file=file.path(resultsPath,"12C_Only_OSCPLSDA_Loadings.csv"),row.names=F)
# write.csv(data,file=file.path(resultsPath,"12C_Only_OSCPLSDA_Scores.csv"),row.names=F)




# final.opls.results<-get.OSC.model(obj=results,OSC.comp=1)
#
#
# loadings = results$loadings[[2]]
# mx= cbind(as.double(loadings[,1]),as.double(loadings[,2]))
# df = as.data.frame(mx)
# colnames(df)=c("Comp1","Comp2")
# rownames(df) = rownames(loadings)
# plot = ggplot(df,aes(Comp1,Comp2))
# plot + geom_point(alpha=1/3)
#
# scores = results$scores[[2]]
# scoresDf = as.data.frame(cbind(as.double(scores[,1]), as.double(scores[,2])))
# colnames(scoresDf)=c("Comp1","Comp2")


