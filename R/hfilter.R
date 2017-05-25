hfilter <-
function(xcmsSet2=NULL,resultsPath,iso1,iso2,type=c("greater|=","less|="),value="maxo",returnType=c("valid","suspect"),Fchng=1.15,sampType1=c("12C","13C"),sampType2=c("12C","13C"),gpvals=NULL)
{
  #get C12 and C13 sample names
  sampleNames = rownames(xcmsSet2@phenoData)
  lastCol = ncol(xcmsSet2@phenoData)
  C12sampleNames = sampleNames[grepl("12C",xcmsSet2@phenoData[ ,lastCol])] #may need to change  "12C" to user defined identifier
  C13sampleNames = sampleNames[grepl("13C",xcmsSet2@phenoData[ ,lastCol])] #may need to change  "13C" to user defined identifier
  pairedGroups_matrix = read.csv(file.path(resultsPath,paste(pheno,"Clusters_MatrixForm.csv",sep="_")),check.names=F,stringsAsFactors=F)

  #check input values, generate group value table if necessary
    if(is.null(gpvals)&&!is.null(xcmsSet2))
  {
    noDupe = which(duplicated(pairedGroups_matrix[,"groupname"]))
    if(length(noDupe)==0)
    {
      noDupe=pairedGroups_matrix
    }else
    {
      noDupe = pairedGroups_matrix[-noDupe,]
    }
    idx = match(noDupe[,"groupname"],groupnames(xcmsSet2))
    gpvals = cbind(noDupe,groupval(xcmsSet2,method="medret",value=value)[idx,]) #make sure works
#     gpvalsOut <<-gpvals
  }else if(is.null(gpvals)&&is.null(xcmsSet2))
  {
    stop("Error: No group values or xcmsSet2 input. Cannot generate group value table.")
  }

  idx1 = which(gpvals[,"iso"]==iso1)
  idx2 = which(gpvals[,"iso"]==iso2)

  #get data for both groups from the designated sample types
  if(length(C12sampleNames)==1)
  {
    if(sampType1=="12C")
      samp1 = gpvals[idx1,C12sampleNames]
    else if(sampType1=="13C")
      samp1 = gpvals[idx1,C13sampleNames]

    if(sampType2=="12C")
      samp2 = gpvals[idx2,C12sampleNames]
    else if(sampType2=="13C")
      samp2 = gpvals[idx2,C13sampleNames]
  }
  else
  {
    if(sampType1=="12C")
      samp1 = rowMeans(gpvals[idx1,C12sampleNames],na.rm=T)
    else if(sampType1=="13C")
      samp1 = rowMeans(gpvals[idx1,C13sampleNames],na.rm=T)

    if(sampType2=="12C")
      samp2 = rowMeans(gpvals[idx2,C12sampleNames],na.rm=T)
    else if(sampType2=="13C")
      samp2 = rowMeans(gpvals[idx2,C13sampleNames],na.rm=T)
  }


  #perform test based on type and fold change
  if(type=="greater|=")
  {
    result= (samp1 >= samp2*Fchng)
  }else if (type=="less|=")
  {
    result = (samp1*Fchng <= samp2)
  }
  resultOut <<- result

  #Order of definition is such that groups that are in mulitple clusters (due to duplication) will have all of their clusters returned
  if(returnType=="valid")
  {
    groups = unique(gpvals[idx1,"groupname"][result])
    idx = which(pairedGroups_matrix[,"groupname"]%in%groups)
    clusters = unique(pairedGroups_matrix[idx,"clusterNumber"])
  }else if(returnType=="suspect")
  {
    groups=unique(gpvals[idx1,"groupname"][-which(result)])
    idx = which(pairedGroups_matrix[,"groupname"]%in%groups)
    clusters = unique(pairedGroups_matrix[idx,"clusterNumber"])
  }else
  {
    stop("Invalid return type.  Must be either 'valid' or 'suspect'\n")
  }

  return(list(idx,clusters,groups))
}
