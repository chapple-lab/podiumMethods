tFilter_V4 <-
function(xcmsSet2=NULL,resultsPath,phenotypes=NULL,phenoTag,mode=c("filter","validate"),alpha=0.05,peakType1,peakType2,alternative=c("two.sided","greater","less"),value="maxo",type=c("valid","suspect"),
                      sampleType1=c("12C","13C"),sampleType2=c("12C","13C"),fout=T,cout=T,gpvals=NULL,removeDupes=F,returnType=1,overide=F,...)
{
  #check to see if xcmsSet has been filled
  check = is.na(xcmsSet2@filled[1])
  if(check&!overide)
  {
    stop("\nError: xcmsSet not filled. Run fillPeaks() on xcmsSet first.\nTo run without filling, use overide=TRUE.\n")
  }else if(overide&check)
  {
    cat("\nWarning: xcmsSet appears to not be filled\n")
  }


  #get C12 and C13 sample names
  sampleNames = rownames(xcmsSet2@phenoData)
  lastCol = ncol(xcmsSet2@phenoData)
  C12sampleidx = grepl("12C",xcmsSet2@phenoData[ ,lastCol]) #may need to change  "12C" to user defined identifier
  C13sampleidx = grepl("13C",xcmsSet2@phenoData[ ,lastCol]) #may need to change  "13C" to user defined identifier
  pairedGroups_matrix = read.csv(file.path(resultsPath,paste(phenoTag,"Clusters_MatrixForm.csv",sep="_")),check.names=F,stringsAsFactors=F)


  #define file extensions
  filePattern = "\\.[Cc][Dd][Ff]$|\\.[Nn][Cc]$|\\.([Mm][Zz])?[Xx][Mm][Ll]$|\\.[Mm][Zz][Dd][Aa][Tt][Aa]$|\\.[Mm][Zz][Mm][Ll]$"
  datePattern = "(19|20)\\d\\d([- /._])(0[1-9]|1[012])\\2(0[1-9]|[12][0-9]|3[01])"

  #check input values, generate group value table if necessary
  if(is.null(gpvals)&&!is.null(xcmsSet2))
  {
    if(removeDupes)
    {
      noDupe = which(duplicated(pairedGroups_matrix[,"groupname"]))

      if(length(noDupe)==0) ## no duplicated groupnames found
      {
        noDupe = pairedGroups_matrix
        idx = match(noDupe[,"groupname"],groupnames(xcmsSet2))

      }else{

        noDupe = pairedGroups_matrix[-noDupe,] #removes all but first occurance of duplicated groups
        idx = match(noDupe[,"groupname"],groupnames(xcmsSet2))
      }
      gpvals = cbind(noDupe,groupval(xcmsSet2,method="medret",value=value)[idx,])
    } else
    {
      idx = match(pairedGroups_matrix[,"groupname"],groupnames(xcmsSet2))
      gpvals = cbind(pairedGroups_matrix,groupval(xcmsSet2,method="medret",value=value)[idx,])
    }

  }else if(is.null(gpvals)&&is.null(xcmsSet2))
  {
    stop("Error: No group values or xcmsSet2 input. Cannot generate group value table.")
  }



  #if phenotypes not input, run based on highest level phenotype in @phenoData
  if(is.null(phenotypes))
  {
    phenotypes = xcmsSet2@phenoData[ ,1]
    phenotypes = unique(phenotypes)
  }

  idx1 = which(gpvals[,"iso"]==peakType1)
  idx2 = which(gpvals[,"iso"]==peakType2)
  ##handle cases where there are uneven numver of iso occurances
  clusts = intersect(gpvals$clusterNumber[idx1],gpvals$clusterNumber[idx2])
  idx1 = idx1[which(gpvals$clusterNumber[idx1]%in%clusts)]
  idx2 = idx2[which(gpvals$clusterNumber[idx2]%in%clusts)]

  gps = gpvals[idx1,"groupname"]

  tTests = mat.or.vec(length(idx1),length(phenotypes))
  validGroups = vector(mode = "list", length=length(phenotypes))
  peakType1data = vector(mode = "list", length=length(phenotypes))
  peakType2data = vector(mode = "list", length=length(phenotypes))
  j=1
  for(pheno in phenotypes)
  {
    #determine files corresponding to phenotype, then separate samples based on isotope
    phenoIdx = grepl(pheno,xcmsSet2@phenoData[,1])
    C12sampleNames = sampleNames[phenoIdx&C12sampleidx]
    C13sampleNames = sampleNames[phenoIdx&C13sampleidx]


    #set variables according to input peak type and iso
    if(sampleType1=="12C")
    {
      sampT1=C12sampleNames
    }else if(sampleType1=="13C")
    {
      sampT1=C13sampleNames
    }else stop("\nInvalid Iso:",sampleType1,"entered. Iso must be either '12C' or '13C'.\n")

    if(sampleType2=="12C")
    {
      sampT2=C12sampleNames
    }else if(sampleType2=="13C")
    {
      sampT2=C13sampleNames
    }else stop("\nInvalid Iso:",sampleType2,"entered. Iso must be either '12C' or '13C'.\n")


    #get raw data according to peak type
    peakList1 = lapply(idx1,function(x) gpvals[x,sampT1])
    peakList2 = lapply(idx2, function(x) gpvals[x,sampT2])


    #perform unpaired, unequal variance, one tailed welches's t-test for each pair of peaks
    #replaced args with ... to allow for full flexibility "alternative=H0,mu=0,paired=FALSE,var.equal=FALSE"
    tTest = mapply(t.test,peakList1,peakList2, MoreArgs=list(alternative=alternative,var.equal=T,...),SIMPLIFY = F)
    tTest = sapply(tTest,function(x) x$p.value)


    #store/append results for this phenotype
    tTests[,j] = tTest
    validGroups[[j]] = as.numeric(tTest<=alpha)
    peakType1data[[j]] = do.call(rbind,peakList1)
    peakType2data[[j]] = do.call(rbind,peakList2)#cbind(peakType2data,matx)

    j=j+1
  }


  rm(peakList1)
  rm(peakList2)

  #Format for Output
  peakType1data = do.call(cbind,peakType1data)
  peakType2data = do.call(cbind,peakType2data)
  validGroups = do.call(cbind,validGroups)
  colnames(tTests)= paste("pval",phenotypes,sep="_")
  colnames(validGroups)= paste(phenotypes,"_valid(a=",alpha,")",sep="")
  colnames(peakType1data)=sub(filePattern,"",x=colnames(peakType1data))
  colnames(peakType2data)=sub(filePattern,"",x=colnames(peakType2data))
  colnames(peakType1data)=sub(datePattern,"",x=colnames(peakType1data))
  colnames(peakType2data)=sub(datePattern,"",x=colnames(peakType2data))


  peakSet1 = rep(paste(sampleType1,peakType1,sep="-"),length(gps))
  peakSet2 = rep(paste(sampleType2,peakType2,sep="-"),length(gps))
  # stuff<<-list(gps,tTests,validGroups,peakSet1,peakType1data,peakSet2,peakType2data)
  Tresults = cbind(groupname=gps,tTests,validGroups,peakType1=peakSet1,peakType1data,peakType2=peakSet2,peakType2data)

  #flatten valid groups and convert (both) to logical
  if(returnType==1)
  {
    flatValidGroups = rowSums(validGroups)>0
    validGroups = validGroups==1
  }else if(returnType==2)
  {
    flatValidGroups = rowSums(validGroups)==length(phenotypes)
    validGroups = validGroups==1
  }else
  {
    stop("Error: Return type must be either 1 (valid in at least one phenotype) or 2 (valid in all phenotypes).")
  }

  if(fout)
    write.csv(Tresults,file=file.path(resultsPath,paste(phenoTag,"tTestFilterValidation",sampleType1,peakType1,alternative,sampleType2,paste(peakType2,".csv",sep=""),sep="_")),row.names=F)
  if(cout)
  {
    cat("\n\n",paste(sampleType1,"_",peakType1," ",alternative," than ",sampleType2,"_",peakType2,"?",sep=""),"\n\nNumber of Valid Groups:\t",length(which(flatValidGroups)),
        "\nTotal Group Comparisons:",length(flatValidGroups),"\nPercent Valid:\t",round(length(which(flatValidGroups))/length(flatValidGroups)*100,digits=2),"%",sep="")
    if(returnType==1)
      cat("\n\n(Note: Valid indicates the group was valid in at least one phenotype)\n")
    else if(returnType==2)
    {
      cat("\n\n(Note: Valid indicates the group was valid in all phenotypes)\n")
    }else
      stop("Error: Return type must be either 1 (valid in at least one phenotype) or 2 (valid in all phenotypes).")

  }

  #Order of definition is such that groups that are in mulitple clusters (due to duplication) will have all of their clusters returned/removed
  #unique() removed when determining groups because duplicate groups are already removed prior to this step
  if(type=="valid")
  {
    if(mode=="validate")
    {
      groups = gps[flatValidGroups]
      idx = which(pairedGroups_matrix[,"groupname"]%in%groups&pairedGroups_matrix$iso==peakType1) #check
      clusters = unique(pairedGroups_matrix[idx,"clusterNumber"]) #may want to change to more stringent form: gpvals[idx1[flatValidGroups],"clusterNumber"]
    }else if(mode=="filter")
    {
      groups = gps[!flatValidGroups]
      idx = which(pairedGroups_matrix[,"groupname"]%in%groups&pairedGroups_matrix$iso==peakType1)
      clusters = unique(pairedGroups_matrix[idx,"clusterNumber"])
    }else
    {
      stop("Error: Mode must either be \'filter\' or \'validate\'.")
    }

  }else if(type=="suspect")
  {
    if(mode=="validate")
    {
      groups=gps[!flatValidGroups]
      idx = which(pairedGroups_matrix[,"groupname"]%in%groups&pairedGroups_matrix$iso==peakType1)
      clusters = unique(pairedGroups_matrix[idx,"clusterNumber"])
    }else if(mode=="filter")
    {
      groups = gps[flatValidGroups]
      idx = which(pairedGroups_matrix[,"groupname"]%in%groups&pairedGroups_matrix$iso==peakType1)
      clusters = unique(pairedGroups_matrix[idx,"clusterNumber"])
    }

  }else
  {
    stop("Invalid type.  Must be either 'valid' or 'suspect'\n")
  }

  return(list(idx=idx,clusters=clusters,groups=groups))
}
