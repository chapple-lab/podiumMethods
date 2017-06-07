CAMERA_Prep = function(xcmsSet2,resultsPath=NULL,isos="all",pheno,type,tTest,dataVal)
{

  #read in valid Labelled Group Data
  if(type=="suspect")
  {
    stop("Type must be \"valid\", make sure validation of type \"valid\" has already been run")
  }else if(type=="valid")
  {
    if(tTest==T)
    {
      labelData = read.csv(file=file.path(resultsPath,paste(pheno,type,"Clusters_tTestFilter",paste(dataVal,".csv",sep=""),sep="_")),check.names=F,stringsAsFactors=F)

    } else {
      labelData = read.csv(file=file.path(resultsPath,paste(pheno,type,"Clusters_HeuristicFilter",paste(dataVal,".csv",sep=""),sep="_")),check.names=F,stringsAsFactors=F)
    }
  }else
  {
    warning("invalid type argument")
  }
  #Remove all groups and related information from xcmsSet2 that are not labeled Groups

  if(isos!="all")
  {
    labeledGpIdx = which(labelData$iso%in%isos)
    labeledGpIdx = match(labelData[labeledGpIdx,"groupname"],groupnames(xcmsSet2))
  }
  else
  {
    labeledGpIdx = match(labelData[ ,"groupname"], groupnames(xcmsSet2))
  }


  xcmsSet2@groups = xcmsSet2@groups[labeledGpIdx, ]
  xcmsSet2@groupidx = xcmsSet2@groupidx[labeledGpIdx]
  return(xcmsSet2)
}


