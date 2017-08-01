##Pairwise CAMERA Script
#to be run in environment after standard analysis
#switched order of calling libraries

require(CAMERA)
#source('X:/Labs/Chapple/0-Cole Wunderlich/R files/CAMERA/CAMERA_Prep.R')
#source("X:/Labs/Chapple/0-Cole Wunderlich/R files/Code/Methods/LogFunctions.r")
#source('X:/Labs/Chapple/0-Cole Wunderlich/R files/CAMERA/getPeakAnnotations.R')
#source("X:/Labs/Chapple/0-Cole Wunderlich/R files/Database Search/findMassMatches.R")

raw=T
global=T
flag=F
C12=T #only uses 12C fed samples
C12Compat=F #for compatibility with OPLSDA ran in C12 mode, only needed when all samples are used for CAMERA but only C12 samples in OPLSDA
masterData = NULL
phenoData = dataSet2@phenoData #may substitute xs@xcmsSet@phenoData if using a previously run file with dataSet2 removed
phenoData = cbind(as.character(levels(phenoData[,1])[phenoData[,1]]),as.character(levels(phenoData[,2])[phenoData[,2]]))
rownames(phenoData) = rownames(dataSet2@phenoData)

if(raw)
{
  cameraData = CAMERA_Prep(dataSet2,resultsPath,isos="12C",pheno=phenoTag,type=type,tTest=T,dataVal=val)
  cameraData = convert(cameraData)
}else
{
  cameraData = xs@xcmsSet
}


if(!file.exists(file.path(resultsPath,"CAMERA Analysis")))
  dir.create(file.path(resultsPath,"CAMERA Analysis"))


labelData = read.csv(file=file.path(resultsPath,paste(phenoTag,type,"Clusters_tTestFilter",paste(val,".csv",sep=""),sep="_")),check.names=F,stringsAsFactors=F)
labelData12C = labelData[which(labelData$iso=="12C"),]

dataBase = data(Database_Final_V4)
#dataBase = read.csv("X:/Labs/Chapple/0-Cole Wunderlich/Database/Database_Final_V4.csv",stringsAsFactors = F,check.names = F)

phenos=phenoData[,1]
wtindx = which(grepl("[Ww][Tt]",phenoData[,1])) #change for desired base phenotype
phenos = phenos[-wtindx] #remove WT from pheno list
phenos = unique(phenos)
if(C12)
{
  wtindx = wtindx[which(phenoData[wtindx,ncol(phenoData)]=="12C")] #get only WT 12C samples
}


if(global)
  phenos = c(phenos,"Global")


for(pheno in phenos)
{
  cat("\nProcessing Phenotype: ",pheno,"\n")
  if(pheno=="Global")
  {
    pheno=phenos[!phenos=="Global"]
    flag=T
  }

    phenosindx = which(phenoData[,1]%in%pheno)
    if(C12)
        phenosindx=which(phenoData[phenosindx,ncol(phenoData)]=="12C")
    phenosindx = c(phenosindx,wtindx) #order of arguments is important!!

    xs = xsAnnotate(cameraData,nSlaves=7,polarity="negative",sample=phenosindx)
    xs=groupFWHM(xs,intval="into")#perfwhm=1
    xs = findIsotopes(xs,intval="into",minfrac=0.15)
    xs = groupCorr(xs,calcIso=T,calcCaS = T, cor_eic_th = 0.6,cor_exp_th = 0.6) #cor_exp_th = 0.6,cor_eic_th = 0.6
    xs = findAdducts(xs,polarity="negative")
    cleanParallel(xs)


    annos=getPeakAnnotations(xs)

    if(C12Compat|C12) #for compatibility
        pheno = paste(pheno,".12C",sep="")

    if(flag)
      pheno="Global"

    if(pheno=="Global")#Run database search
    {
      endData2 = cbind(labelData12C[,1:2],annos,labelData12C[,-(1:2)])

      #Create final file with mass hit information, flatten multiple hits into one if needed
      hits = findMassMatches(endData2,dataBase,ppm=30,rawmz="mzmed")

      hitList = sapply(endData2$mzmed,function(x){
        indx = which(hits$Peak_Mz==x) ##NOTE: will need to change 'Peak_Mz' depending on type of findMassMatches search
        if(length(indx)>1)
          matches = paste(hits[indx,"Name"],collapse="_or_")
        else if(length(indx)==0)
          matches = "None"
        else
          matches = hits[indx,"Name"]

        return(matches)
      })

      endData2 = cbind(endData2[,1:6],Database_Matches=hitList,endData2[,-(1:6)])
      if(C12)
        write.csv(endData2,file.path(resultsPath,"CAMERA Analysis",paste("12C_CAMERA_FinalResults_MassAnnosAndHits.csv",sep="_")),row.names = F)
      else write.csv(endData2,file.path(resultsPath,"CAMERA Analysis",paste("CAMERA_FinalResults_MassAnnosAndHits.csv",sep="_")),row.names = F)

    }

    annos = cbind(Type=rep(pheno,nrow(annos)),annos,Groupname=labelData12C$groupname,stringsAsFactors=F)
    masterData = rbind(masterData,annos)

}

if(C12)
  write.csv(masterData,file=file.path(resultsPath,"CAMERA Analysis",paste("12C_CAMERA_AllAnnotations_AllGenotypesPairedWithWT.csv",sep="_")),row.names=F) else
  write.csv(masterData,file=file.path(resultsPath,"CAMERA Analysis",paste("CAMERA_AllAnnotations_AllGenotypesPairedWithWT.csv",sep="_")),row.names=F)


cat("\n-----CAMERA Analysis Complete-----\n")
