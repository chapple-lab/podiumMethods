#!/depot/chapple/apps/R-3.0.2/bin/Rscript --vanilla
# Cole Wunderlich
# 10/13/2015
# MS type = Multiple
# EIC ppm = 50
# Alpha = 0.05
## Manual Inspections Notes:
## Large BPC are around 15 seconds full width, small BPC are around X sec
## Background is around X
## Noise is around X throughout
##scans taken about every .7sec on our machine

require(coleXcmsMethods,quietly=T)
#LOAD COLOR LIST!!!! (for EIC generation)
data(LightColorList)

##Enable Logging? Append or overwrite existing log?
Log=T
Append=F

#Run CAMERA and OPLSDA afterwards?
scripts=T

##Minimize retained data? If = 1 will just remove setList, if = 2 will remove dataSet and setList
##SaveData, 0 = do not save, 1=save final dataSet2, 2=save dataSet and dataSet2
datMin=1
saveData =2
readRawData=T

##Set Run Mode
##Compare 2 phenotypes, run multiple (more than 2) phenotypes together, or run each individually?
##Note: To run each individually, set both comparison and multi to false.  Also, only one can be set to ture.
comparison=F
multi =T
if(comparison&multi)
{
  stop("Error: Only one run mode at most can be set to true.")
}


##List of phenotypes to test. If "All" all phenotypes will be tested. Must be a list of strings matching the names
##of the folders in the highest level of the start directory.  (ie. must match phenotypes contained in the first row
##of the xcmsSet @phenoData slot). eg. phenos=c("ref3-3","fah1-2","cadC/D")
##NOTE: This parameter is ignored if running in single/batch mode. (will run all of them automatically)
selectedPhenos = c("All")

##Use Statistical or Heuristic Based Pairing Algorithm (P)? Use Statistical or Heuristic Based Validation Algorithm (V)?
tTestPairing = T
tTestFilter = T

##What value should be used for peak intensity? (into=integrated peak intensity, maxo=maximum peak intensity)
#Note: diffval is the peak intensity value used for the diffreport() function.  This parmeter allows for one
#value to be used for pairing and validation purposes and another to be used for the differential experssion
#analysis.  THE FINAL REPORT FILE WILL CONTAIN RAW VALUES ACCORDING TO DIFFVAL
##diffval currently turned off
val="into"
diffval=val

##Final statistical report?
Freport = T

##Remove groups from initial column elution period?
##If so, what is the retention time threshold (minimum time after which to keep)?
rem=T
rtmin= 110

#define paths to source data and results folder
dataPath = "C:/0-TempData/Stem Feeding T1/Start"
resultsPath = "X:/Labs/Chapple/0-Cole Wunderlich/Results/0-Stem Feeding/Test2"
if(!file.exists(resultsPath))
  dir.create(resultsPath)
if(!comparison & !multi)
{
  baseResultsPath = resultsPath
}

#Log console output
#NOTE: Cannot Seem to capture warnings in either RunLog or ErrorLog
if(Log)
  startLog(resultsPath,Append)

cat("\nData Path: ",dataPath)
cat("\nResults Path: ",resultsPath,"\n")
if(comparison)
{
  cat("\nRun Mode: Comparison\n")
}else if(multi)
{
  cat("\nRun Mode: Multiple\n")
} else
{
  cat("\nRun Mode: Single/Batch\n")
}

#-----------
###Load Data
#Params for centWave ONLY
ppm = 40
peakwidth = c(4,20)
snthresh = 5
prefilter = c(3,150)
integrate = 1
nslaves = 16

#print centWave params
cat("\n\ncentWave Parameters\n\nROI max M/z width in ppm:",ppm,"\nPeak width tolerance (min,max) in sec:",peakwidth,"\nS/N Threshold:",snthresh,"\nROI Prefilter (minNumPeaks,IntensityThresh):",prefilter)
cat("\nIntegration Technique (1=Mexican Hat, 2=Raw Data):",integrate,"\n","nSlaves:",nslaves,"\n")

if(readRawData) {dataSet = xcmsSet2(files=dataPath,profmethod="binlin",profparam=list(profstep=0.01, step=0.01),nSlaves=nslaves,method="centWave",ppm=ppm,peakwidth=peakwidth,snthresh=snthresh,prefilter=prefilter,integrate=integrate,mzdiff=-0.01,fitgauss=T,verbose.columns=T)
cat("\n----Peak Picking Complete----\n")
  print(Sys.time())
  cat("\n")
  print(proc.time())
  cat("\n")

}else 
{
 load(file.path(resultsPath,"dataSet.RData"))
 cat("\nUsing Preloaded Data\n")
}
 

#----------
##Separate out each phenotype (based on highest level phenotype)
#phenoTag is an identifier for the given run based on the run mode and selected phenotypes
phenotypes = dataSet@phenoData[ ,1]
phenotypes = as.character(unique(phenotypes))
write.csv(dataSet@phenoData,file=file.path(resultsPath,"phenoData.csv"))

if("12C"%in%phenotypes | "13C"%in%phenotypes) #if only one phenotype, generate phenotype name using common elements in each sample's file path
{
  phenotypes = genPheno(dataSet)
  dataSet@phenoData = cbind(phenotype=rep(phenotypes,nrow(dataSet@phenoData)),dataSet@phenoData) #add phenotype information to dataSet
  setList = list(dataSet)
  names(setList) = phenotypes
  cat("\n\nDetected Phenotype:",phenotypes,"\n")
  runList = phenotypes
  #ToDo: Code writeAllPeaks for just one phenotype
}else if(comparison|multi) #If analyzing more than one phenotype at a time
{
  if(selectedPhenos[1]=="All"|selectedPhenos[1]=="all")#run all detected phenotypes, check to make sure only "all" is being specified.
  {
    if(length(selectedPhenos)>1)
      stop("Error: Both \"All\" and other phenotypes specified. Please use either \"All\" or specify individual phenotypes for analysis. (Not both)")

    cat("\nAnalyzing all phenotypes:",as.character(phenotypes))
    if(comparison&(length(phenotypes)>2|length(phenotypes)==1))
      stop("Error: Incorrect number of phenotypes selected.  Comparison mode requires excatly two phenotypes.")

    phenoTag = paste(phenotypes,collapse="_")
  }else{
    #Select desired phenotypes
    #first check to ensure all selected phenotypes have been detected
    if(length(which(selectedPhenos%in%phenotypes))!=length(selectedPhenos))
      stop(paste("Error: Not all specified phenotypes where detected\n\nDetected:",paste(phenotypes,collapse = " "),"\nSpecified:",paste(selectedPhenos,collapse = " "),
                 "\nMissing:",selectedPhenos[which(!selectedPhenos%in%phenotypes)]))

    if(comparison&(length(phenotypes)>2|length(phenotypes)==1))
      stop("Error: Incorrect number of phenotypes selected.  Comparison mode requires excatly two phenotypes.")

    phenoTag = paste(selectedPhenos,collapse="_")
    if(comparison)
      phenoTag = paste("Comparison",phenoTag,sep="_")

    phenotypes = phenotypes[which(phenotypes%in%selectedPhenos)]
  }
  setList = list(dataSet)
  names(setList) = phenoTag
  if(!file.exists(file.path(resultsPath,paste(phenoTag,"AllPeaks.csv",sep="_"))))
    writeAllPeaks(dataSet,resultsPath,name=paste(phenoTag,"AllPeaks.csv",sep="_"))
  runList = phenoTag

  } else {
    setList = split(dataSet,dataSet@phenoData[ ,1])
    cat("\n\nDetected Phenotypes:",phenotypes,"\n")
    runList = phenotypes
  }


  if(datMin==2)
    {rm(dataSet)
    gc()}


#----------
###Process Each Phenotype or process group of phenotypes depending on settings
##Will use all of loaded data except for peak pairing and validation output where only the selected phenotypes are used
##tag is a unique identifier for each run based on the run mode and selected phenotypes
for(phenoTag in runList)
{
  cat("\n\nProcessing:",phenoTag,"\n")
  if(!(comparison|multi)) #if single/batch mode
  {
    resultsPath=file.path(resultsPath,phenoTag)
    if(!file.exists(resultsPath))
      dir.create(resultsPath)
    phenotypes = phenoTag

  }

  ###Perform Rt correction and Grouping across ALL samples
  ##Rt Correction
  response=100
  #align using WT sample???  currently center=NULL which causes sample with most peaks to be chosen
  cat("\n\nAligning Retention Times\nObiwarp Response:",response,"\n")
  dataSet2 = retcor.obiwarp(setList[[phenoTag]],response=response,plottype="deviation",distFunc="cor_opt",gapInit=0.3,gapExtend=2.4) #may want to change back to 10 resp
  dev.copy(pdf,file=file.path(resultsPath,paste(phenoTag,"Resp",response,"ObiwarpPlot.pdf",sep="_")),width=17.50,height=9.71)
  dev.off()
  dev.off()

  if(!file.exists(file.path(resultsPath,paste(phenoTag,"Aligned_AllPeaks.csv",sep="_"))))
    writeAllPeaks(dataSet2,resultsPath,name=paste(phenoTag,"Aligned_AllPeaks.csv",sep="_"))


  if(!file.exists(file.path(resultsPath,paste(phenoTag,"Aligned_TIC.pdf",sep="_"))))
  {
    cat("\n\nPrinting TIC to PDF\n")
    printTIC(dataSet2,pdfname=file.path(resultsPath,paste(phenoTag,"Aligned_TIC.pdf",sep="_")),rt="corrected")
  }

  ##Grouping
  #  group(res10RtcorData,minfrac=0.25,minsamp=2,bw=3,max=50,sleep=0.01)
  #  group(res10RtcorData,minfrac=0.25,minsamp=2,bw=5,max=50,sleep=0.01)
  cat("\n\nGrouping Detected Peaks\n")
  minfrac=0
  minsamp=3
  mzwid=0.025
  bw=2
  max=50

  cat("\nMinimum Fraction of Samples Required for Peak Retention: ",minfrac,"\n")
  cat("Minimum Number of Samples Required for Peak Retention: ",minsamp,"\n")
  cat("M/z step: ",mzwid,"\nBandwidth of Guassian Kernel (stdev): ",bw,"\nMax groups per M/z step: ",max,"\n\n")

  dataSet2 = group(dataSet2,minsamp=minsamp,minfrac=minfrac,mzwid=mzwid,bw=bw,max=max) #minfrac=0.25, mzwid!!
  show(dataSet2)
  cat("\n")

  write.csv(cbind(groupnames(dataSet2),dataSet2@groups,groupval(dataSet2,method="medret",value=val)),file=file.path(resultsPath,paste(phenoTag,"InitialGroupings_AllGroups.csv",sep="_")),row.names=F)

  if(rem)
  {
    cat("\n\nRemoving Groups with Rt less than ",rtmin,"sec.\n")
    dataSet2=removeGroups(dataSet2,rtrng=c(0,rtmin))
    write.csv(cbind(groupnames(dataSet2),dataSet2@groups,groupval(dataSet2,method="medret",value=val)),file=file.path(resultsPath,paste(phenoTag,"_InitialGroupings_Under_",rtmin,"sec_Groups_Removed.csv",sep="")),row.names=F)
  }

  ###find C13 labelled peaks
  #parameters
  mzppm<-10 ## was 15
  mzabs<-0
  rterror<-2
  C12Fchng = 1.25
  C13Fchng = 1.40
  preA = 0.05
  postA = 0.05
  postfilter=F


  if(tTestPairing)
  {

    cat("\n\nFilling Peaks\n")
    dataSet2 = fillPeaks(dataSet2,"chrom",nSlaves=nslaves,expand.mz=1,expand.rt=1,min.width.mz=0.005,min.width.rt=1) #note: Minwidth for rt is in seconds, scans taken about every .7sec on our machine, exp is factor based (ie 3 = 3X expansion)->1 will do nothing
    cat("\n\nPairing Groups Using t-Test Based Pairing Algorithm\n")
    cat("\nPairing Parameters\nmzppm:",mzppm,"\nmzabs:",mzabs,"\nrterror:",rterror,"\nData Type:",val,"\n")
    cat("\nPrefilter alpha:",preA,"\n")
    cat("Postfilter:",postfilter)
    if(postfilter){
      cat("\nPostfilter alpha:",postA)
    }else cat("\n\n")
    C13peaks_GroupPairing_tTest(dataSet2,mzppm,mzabs,rterror,resultsPath,phenoTag,phenotypes,preA,postA,value=val)
  }else
  {
    cat("\n\nPairing Groups Using Heuristic Pairing Algorithm (V2.3)\n")
    cat("\nPairing Parameters\nmzppm:",mzppm,"\nmzabs:",mzabs,"\nrterror:",rterror,"\nBase Peak 12C Fchng Tolerance:",C12Fchng,"\n")
    cat("Base Peak 13C Fchng Tolerance",C13Fchng,"\nData Type:",val,"\n\n")
    C13peaks_GroupPairing_Heuristic(dataSet2,mzppm,mzabs,rterror,resultsPath,phenoTag,C12Fchng,C13Fchng,value=val)
    cat("\n\nFilling Peaks\n")
    dataSet2 = fillPeaks(dataSet2,"chrom",nSlaves=nslaves,expand.mz=1,expand.rt=1,min.width.mz=0.005,min.width.rt=1) #note: Minwidth for rt is in seconds, scans taken about every .7sec on our machine, exp is factor based (ie 3 = 3X expansion)->1 will do nothing
  }

  #Fill Peaks output
  if(!file.exists(file.path(resultsPath,paste(phenoTag,val,"FilledGroups.csv",sep="_"))))
    write.csv(cbind(groupnames(dataSet2),dataSet2@groups,groupval(dataSet2,method="medret",value=val)),file=file.path(resultsPath,paste(phenoTag,val,"FilledGroups.csv",sep="_")),row.names=F)

  if(!file.exists(file.path(resultsPath,paste(phenoTag,"Filled_Aligned_AllPeaks.csv",sep="_"))))
    writeAllPeaks(dataSet2,resultsPath,name=paste(phenoTag,"Filled_Aligned_AllPeaks.csv",sep="_"))

  ###Validate/Filter Results
  ##Initialize Validation Parameters, including thresholds for profiling suspect groups (applied only to clustered groups)
  Fchng = 1.30
  type = "valid"
  ppm = 15

  #Perform statistical filtration of paired groups
  if(postfilter&tTestFilter)
  {
    cat("\n\nFiltering groups using Welch's unpaired t-test\nData Type:",val,"\n")
    #13C_base > 12C_base?
    tresults = tFilter_V4(dataSet2,resultsPath,phenotypes=phenotypes,phenoTag=phenoTag,postA,peakType1="12C",peakType2="12C","greater",value=val,type=type,"13C","12C",returnType=2,mode="filter")
    clusters = tresults$clusters
    idx = tresults$idx
    #12C_6 > 12C_base?
    tresults = tFilter_V4(dataSet2,resultsPath,phenotypes=phenotypes,phenoTag=phenoTag,postA,peakType1="13C_6",peakType2="12C","greater",value=val,type=type,"12C","12C",returnType=2,mode="filter")
    clusters = intersect(clusters,tresults$clusters) #keep only clusters returned as valid/ok by both filters
    idx = union(idx,tresults$idx) #keep all groups returned as valid/ok by any filter
  } else if(!tTestFilter)
  {
    ##Heuristic Validation
    #determine invalid C13_6 peaks using a heuristic filter (12C_6*Fchng >= 13C_6)
    cat("\n\nValidating 13C_6 groups using heuristic filtration\n")
    cat("\nValidation Fold Change: ",Fchng,"\n")
    filterout = hfilter(dataSet2,resultsPath,"13C_6","13C_6","greater|=",value=val,returnType=type,Fchng,"13C","12C")
    idx = filterout[[1]]
    clusters = filterout[[2]]
  }else
  {
    pairedGroups_matrix = read.csv(file.path(resultsPath,paste(phenoTag,"Clusters_MatrixForm.csv",sep="_")),check.names=F,stringsAsFactors=F)
    clusters = unique(pairedGroups_matrix$clusterNumber)
    idx = 1:nrow(pairedGroups_matrix)
    rm(pairedGroups_matrix)
    gc()
  }

  #output EICs, MS, and Tables for clusters and groups according to the Validation Parameters selected above
  group_validationOutput(dataSet2,phenoTag,resultsPath,type,tTestFilter,value=val,clusters,Indx=idx,Rtindent=2,ppm)



  #Generate Statistical Report
  if(type=="valid"&&comparison&&Freport)
  {
    report = diffreport2(dataSet2,resultsPath,type,tTestFilter,phenoTag,filebase=file.path(resultsPath,paste(phenoTag,diffval,"StatResults",sep="_")),eicmax=100,dataVal=val,value=diffval)
  }

  cat("\n----Finished Processing Phenotype:",phenoTag,"----\n\n")
  print(proc.time())
  cat("\n")
  if(!(comparison|multi))
    resultsPath=baseResultsPath
}

cat("\n----Analysis Complete----\n")
print(Sys.time())
proc.time()
cat("\n")

if(datMin>=1)
  {rm(setList)
  gc()}

if(saveData>=2&datMin<=1)
{
  save(dataSet,file=file.path(resultsPath,"dataSet.RData"))
  save(dataSet2,file=file.path(resultsPath,"dataSet2.RData"))
  #save environment vars
  env = !grepl("data[Ss]et|[sS]et[Ll]ist",ls())
  save(list=ls()[env],file=file.path(resultsPath,"env.RData"))
  cat("\n----Data Save Complete----\n")
  print(Sys.time())
  cat("\n")
  print(proc.time())
  cat("\n")
}else if(saveData>=1)
{
  save(dataSet2,file=file.path(resultsPath,"dataSet2.RData"))
  #save environment vars
  env = !grepl("data[Ss]et|[sS]et[Ll]ist",ls())
  save(list=ls()[env],file=file.path(resultsPath,"env.RData"))
  cat("\n----Data Save Complete----\n")
  print(Sys.time())
  cat("\n")
  print(proc.time())
  cat("\n")
}
##Need to update non-cluster scripts before enabeling
#cat("\n----Running Postrun Analysis----\n\n")
#cat("----Running CAMERA----\n\n")
#source("/depot/chapple/data/coleXcms/scripts/CAMERA_Cluster_Script.r")
#cat("\n\n----Running OPLSDA----\n\n")
#source("/depot/chapple/data/coleXcms/scripts/OPLSDA_Cluster_Script.r")
#cat("\n----Postrun Analysis Complete----\n\n")


if(Log)
  stopLog()
