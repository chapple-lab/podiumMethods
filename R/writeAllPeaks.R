writeAllPeaks <-
function(xcmsSet2=NULL,resultsPath=NULL,name=NULL)
{
  peaks=xcmsSet2@peaks
  phenoData=xcmsSet2@phenoData

  ##Generate sample type (phenotype and variable name) and sample name columns
  #Initialize variables by performing first operation
  idx=which(peaks[,"sample"]==1)
  sampleType = cbind(rep(paste(phenoData[1,1],phenoData[1,2]),length(idx)))
  colnames(sampleType)="SampleType"

  sampleName = cbind(rep(rownames(phenoData)[1],length(idx)))
  colnames(sampleType)="SampleName"

  #if only one sample
  if(length(unique(peaks[,"sample"]))==1)
  {
    peaks=cbind(peaks,sampleType,sampleName)
    write.csv(peaks,file=file.path(resultsPath,name),row.names=F)
    return(paste(name,"File Write Sucessful"))
  }


  samps = unique(peaks[,"sample"])[-1] #remove first sample since already processed above
  for(x in samps)
  {
    idx=which(peaks[,"sample"]==x)
    sampleType = rbind(sampleType,cbind(rep(paste(phenoData[x,1],phenoData[x,2]),length(idx))))
    sampleName = rbind(sampleName,cbind(rep(rownames(phenoData)[x],length(idx))))
  }

  peaks=cbind(peaks,sampleType,sampleName)
  write.csv(peaks,file=file.path(resultsPath,name),row.names=F)
  return(paste(name,"File Write Sucessful"))

}
