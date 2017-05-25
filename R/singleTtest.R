singleTtest <-
function(xcmsSet2=NULL,phenotypes=NULL,group1,group2,alternative=c("two.sided","greater","less"),value="maxo",iso1=c("12C","13C"),iso2=c("12C","13C"),gpvals=NULL,returnType=1,...)
{

  #if phenotypes not input, run based on highest level phenotype in @phenoData
  if(is.null(phenotypes))
  {
    phenotypes = xcmsSet2@phenoData[ ,1]
    phenotypes = unique(phenotypes)
  }

  #determine files corresponding to phenotype, then separate samples based on isotope
  sampleNames = rownames(xcmsSet2@phenoData)
  lastCol = ncol(xcmsSet2@phenoData)
  C12sampleidx = grepl("12C",xcmsSet2@phenoData[ ,lastCol]) #may need to change  "12C" to user defined identifier
  C13sampleidx = grepl("13C",xcmsSet2@phenoData[ ,lastCol]) #may need to change  "13C" to user defined identifier


  if(is.null(gpvals)&&!is.null(xcmsSet2))
  {
    gpvals = groupval(xcmsSet2,method="medret",value=value)
  }else if(is.null(gpvals)|is.null(xcmsSet2))
  {
    stop("Error: No xcmsSet or groupvals input. Cannot generate needed information.")
  }


  tTests = mat.or.vec(length(phenotypes),1)
  j=1
  for(phenotype in phenotypes)
  {
    phenoIdx = grepl(phenotype,xcmsSet2@phenoData[,1])
    C12sampleNames = sampleNames[phenoIdx&C12sampleidx]
    C13sampleNames = sampleNames[phenoIdx&C13sampleidx]

    #get data for both groups from the designated sample types
    if(iso1=="12C")
      samp1 = gpvals[match(group1,groupnames(xcmsSet2)),C12sampleNames]
    else if(iso1=="13C")
      samp1 = gpvals[match(group1,groupnames(xcmsSet2)),C13sampleNames]

    if(iso2=="12C")
      samp2 = gpvals[match(group2,groupnames(xcmsSet2)),C12sampleNames]
    else if(iso2=="13C")
      samp2 = gpvals[match(group2,groupnames(xcmsSet2)),C13sampleNames]


    #perform t-Test  mu=0,paired=F,var.equal=F
    tTests[j] = t.test(samp1,samp2,alternative=alternative,var.equal=T,...)$p.value
    j=j+1
  }

  if(returnType==1)
  {
    names(tTests)=phenotypes
    return(tTests)
  }
  else if(returnType==2)
  {
    return(min(tTests))
  }
  else if(returnType==3)
  {
    return(max(tTests))
  }else stop("Error: Return Type must be either 1 (all p-values), 2 (smallest p-value), or 3 (largest p-value)")

}
