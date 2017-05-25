tFilterPreProcess <-
function(xcmsSet2=NULL,plus6,resultsPath,phenoTag,phenotypes=NULL,alpha=0.05,peakType1=c("M","Mplus"),peakType2=c("M","Mplus"),iso1=c("12C","13C"),iso2=c("13C","12C"),
                             alternative=c("two.sided","greater","less"),value="maxo",cout=T,fout=T,returnType=1,vect1=NULL, vect2=NULL,...)
{

  if(is.null(xcmsSet2)|is.null(plus6))
  {
    stop("Error: No xcmsSet2 or 'plus6' group matrix input. Cannot generate needed information.")
  }


  #if vectors present, extract information
  if(length(vect1)>0)
  {
    peakType1 = as.character(vect1[1])
    iso1 = as.character(vect1[2])
  }

  if(length(vect2)>0)
  {
    peakType2 = as.character(vect2[1])
    iso2 = as.character(vect2[2])
  }

  #if phenotypes not input, run based on highest level phenotype in @phenoData
  if(is.null(phenotypes))
  {
    phenotypes = xcmsSet2@phenoData[ ,1]
    phenotypes = unique(phenotypes)
  }

  #define file extensions and date format regular expressions
  filePattern = "\\.[Cc][Dd][Ff]$|\\.[Nn][Cc]$|\\.([Mm][Zz])?[Xx][Mm][Ll]$|\\.[Mm][Zz][Dd][Aa][Tt][Aa]$|\\.[Mm][Zz][Mm][Ll]$"
  datePattern = "(19|20)\\d\\d([- /._])(0[1-9]|1[012])\\2(0[1-9]|[12][0-9]|3[01])"

  #get C12 and C13 sample names, read in group data, separate data by specified iso
  sampleNames = rownames(xcmsSet2@phenoData)
  lastCol = ncol(xcmsSet2@phenoData)
  C12sampleidx = grepl("12C",xcmsSet2@phenoData[ ,lastCol]) #may need to change  "12C" to user defined identifier
  C13sampleidx = grepl("13C",xcmsSet2@phenoData[ ,lastCol]) #may need to change  "13C" to user defined identifier
  groupNames = plus6[,peakType1]
  gpvals = groupval(xcmsSet2,method="medret",value=value)
  idx1 = match(plus6[,peakType1],groupnames(xcmsSet2)) #note: this indexing strategy only works becaue every M must have a corresponding Mplus
  idx2 = match(plus6[,peakType2],groupnames(xcmsSet2))

  #initialize vars
  #note: use =NULL instead of mat.or.vec(0,1) if you want TRUE/FALSE in the out put instead of 1/0
  tTests = mat.or.vec(length(idx1),length(phenotypes))
  validGroups = vector(mode = "list", length=length(phenotypes))
  peakType1data = vector(mode = "list", length=length(phenotypes))
  peakType2data = vector(mode = "list", length=length(phenotypes))
  j=1

  ##Process each phenotype
  for(phenotype in phenotypes)
  {

    #determine files corresponding to phenotype, then separate samples based on isotope
    phenoIdx = grepl(phenotype,xcmsSet2@phenoData[,1])
    C12sampleNames = sampleNames[phenoIdx&C12sampleidx]
    C13sampleNames = sampleNames[phenoIdx&C13sampleidx]


    #set variables according to input peak type and iso
    if(iso1=="12C")
    {
      sampT1=C12sampleNames
    }else if(iso1=="13C")
    {
      sampT1=C13sampleNames
    }else stop("\nInvalid Iso:",iso1,"entered. Iso must be either '12C' or '13C'.\n")

    if(iso2=="12C")
    {
      sampT2=C12sampleNames
    }else if(iso2=="13C")
    {
      sampT2=C13sampleNames
    }else stop("\nInvalid Iso:",iso2,"entered. Iso must be either '12C' or '13C'.\n")


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

  peakSet1 = rep(paste(iso1,peakType1,sep="_"),length(groupNames))
  peakSet2 = rep(paste(iso2,peakType2,sep="_"),length(groupNames))
  Tresults = cbind(groupname=groupNames,tTests,validGroups,peakType1=peakSet1,peakType1data,peakType2=peakSet2,peakType2data)

  #flatten valid groups and convert (both) to logical
  flatValidGroups = rowSums(validGroups)>0
  flatAllValidGroups = rowSums(validGroups)==length(phenotypes)
  validGroups = validGroups==1

  if(fout)
    write.csv(Tresults,file=file.path(resultsPath,paste(phenoTag,"tTestPrefilter",iso1,peakType1,alternative,iso2,paste(peakType2,".csv",sep=""),sep="_")),row.names=F)
  if(cout)
  {
    cat("\n\n",paste(iso1,"_",peakType1," ",alternative," than ",iso2,"_",peakType2,"?",sep=""),"\n\nNumber of Valid Groups:\t",length(which(flatValidGroups)),
        "\nTotal Group Comparisons:",length(flatValidGroups),"\nPercent Valid:\t",round(length(which(flatValidGroups))/length(flatValidGroups)*100,digits=2),"%",sep="")
    if(returnType<5)
      cat("\n\n(Note: Valid indicates the group was valid in at least one phenotype)\n")
    else
      cat("\n\n(Note: Valid indicates the group was valid in all phenotypes)\n")
  }

  #Return type 1 returns the indicies of valid groups for each phenotype in the plus6 data frame
  #Return type 2 returns a flattened index of valid groups (ie. a group's index will be returned if it's valid in at least one phenotype)
  #Return type 3 returns the names of the valid groups in the plus6 data frame
  #Retrun type 4 returns the data frame containing information and statistics regarding all of the input groups
  colnames(validGroups) = phenotypes
  if(returnType==1)
    return(apply(validGroups,2,which))
  else if(returnType==2)
    return(which(flatValidGroups))
  else if(returnType==3)
    return(apply(validGroups,2,function(x) groupNames[x]))
  else if (returnType==4)
      return(Tresults)
  else if (returnType==5)
	  return(which(flatAllValidGroups))
}
