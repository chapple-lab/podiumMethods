singleHfilter <-
function(xcmsSet2=NULL,group1,group2,type=c("greater|=","less|="),value="maxo",Fchng=1.15,sampType1=c("12C","13C"),sampType2=c("12C","13C"),gpvals=NULL)
{
  #get C12 and C13 sample names
  sampleNames = rownames(xcmsSet2@phenoData)
  lastCol = ncol(xcmsSet2@phenoData)
  C12sampleNames = sampleNames[grepl("12C",xcmsSet2@phenoData[ ,lastCol])] #may need to change  "12C" to user defined identifier
  C13sampleNames = sampleNames[grepl("13C",xcmsSet2@phenoData[ ,lastCol])] #may need to change  "13C" to user defined identifier
  #check input values, generate group value table if necessary
  if(is.null(gpvals)&&!is.null(xcmsSet2))
  {
    gpvals = groupval(xcmsSet2,method="medret",value=value)
  }else if(is.null(gpvals)|is.null(xcmsSet2))
  {
    stop("Error: No xcmsSet2 or groupvals input. Cannot generate needed information.")
  }

  #get data for both groups from the designated sample types
  if(sampType1=="12C")
    samp1 = rowMeans(gpvals[match(group1,groupnames(xcmsSet2)),C12sampleNames],na.rm=T)
  else if(sampType1=="13C")
    samp1 = rowMeans(gpvals[match(group1,groupnames(xcmsSet2)),C13sampleNames],na.rm=T)

  if(sampType2=="12C")
    samp2 = rowMeans(gpvals[match(group2,groupnames(xcmsSet2)),C12sampleNames],na.rm=T)
  else if(sampType2=="13C")
    samp2 = rowMeans(gpvals[match(group2,groupnames(xcmsSet2)),C13sampleNames],na.rm=T)

  #perform test based on type and fold change
  if(type=="greater|=")
  {
    result= samp1 >= samp2*Fchng
  }else if (type=="less|=")
  {
    result = samp1*Fchng <= samp2
  }
  return(result)
}
