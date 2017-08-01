
##DEPRECATED::calculates range using ppm of 12C peak, not 13C peak, may want to change
##changed to calculate based on the 13C peak
groupmzrange<-function(groups,mzppm=10,mzabs=0.005,nC13=6) {
  C13mzdif<-1.0033548   ##mass of 13C
  mztor<-(groups$mzmed+C13mzdif*nC13)*mzppm*1e-6+mzabs  ##calculate M/z tolerance
  uppermz<-groups$mzmed+C13mzdif*nC13+mztor
  lowermz<-groups$mzmed+C13mzdif*nC13-mztor
  allpeaks<-cbind(groups,lowermz=lowermz,uppermz=uppermz) ##add calculated upper and lower ranges as columns to allpeaks dataframe
  #   write.csv(allpeaks,paste("allpeaksC",nC13,".csv",sep=""),row.names=F)
}

findC13cmpd<-function(groupswithmzrange,rterror=1)
{
  hitlist<-lapply(sort(unique(groupswithmzrange$rtmed)),function(x) ## make a sorted vector of peaks with unique retention times and apply function to it
  {
    peakgp<-groupswithmzrange[which(abs(groupswithmzrange$rtmed-x)<rterror),] # return a filtered matrix of peaks that fall within Rterror of peaks X's retention time -> a peak group
    hit<-sapply (c(1:nrow(peakgp)),function(y) ## apply following function to all peaks in the peak group
    {
      hit<-which(peakgp$mzmed<=peakgp$uppermz[y] & peakgp$mzmed>=peakgp$lowermz[y]) # look in peak group and return all peaks that fall within peak y's M/z range
      if (length(hit)>0) #if at least one match was found, create and return a matrix of that peak and its matching peaks (using peak name and RT) ##Modified to now include all peaks M/z vals
      {
        cbind(M=as.character(peakgp[y,"groupname"]),Mplus=as.character(peakgp[hit,"groupname"]),mzmed_M=peakgp[y,"mzmed"],mzmed_Mplus=peakgp[hit,"mzmed"],rtmed_M=peakgp[y,"rtmed"],rtmed_Mplus=peakgp[hit,"rtmed"],C12avg_M=peakgp[y,"C12avg"],C12avg_Mplus=peakgp[hit,"C12avg"],C13avg_M=peakgp[y,"C13avg"],C13avg_Mplus=peakgp[hit,"C13avg"]) # Return each hit as its: Tag, M+ Tag, Rt, M+ Rt##Modified to include Mz vals and C12/C13avg vals
      }
    })
  })

  hitlist<-unlist(hitlist,recursive=F) ##Flatten hitlist
  ## Determine Null elements in hitlist and set them to NULL
  isnull<-unlist(lapply (hitlist,function(x)
  {
    is.null(x)
  }))
  hitlist[isnull]<-NULL
  #   print(hitlist)
  output<-unique(do.call(rbind,hitlist)) # Make a matrix by setting each element of the flattened hitlist as a row???
  return(output)
}

findC13cmpd_stat<-function(groupswithmzrange,rterror=1)
{
  hitlist<-lapply(sort(unique(groupswithmzrange$rtmed)),function(x) ## make a sorted vector of peaks with unique retention times and apply function to it
  {
    peakgp<-groupswithmzrange[which(abs(groupswithmzrange$rtmed-x)<rterror),] # return a filtered matrix of peaks that fall within Rterror of peaks X's retention time -> a peak group
    hit<-sapply (c(1:nrow(peakgp)),function(y) ## apply following function to all peaks in the peak group
    {
      hit<-which(peakgp$mzmed<=peakgp$uppermz[y] & peakgp$mzmed>=peakgp$lowermz[y]) # look in peak group and return all peaks that fall within peak y's M/z range
      if (length(hit)>0) #if at least one match was found, create and return a matrix of that peak and its matching peaks (using peak name and RT) ##Modified to now include all peaks M/z vals
      {
        cbind(M=as.character(peakgp[y,"groupname"]),Mplus=as.character(peakgp[hit,"groupname"]),mzmed_M=peakgp[y,"mzmed"],mzmed_Mplus=peakgp[hit,"mzmed"],rtmed_M=peakgp[y,"rtmed"],rtmed_Mplus=peakgp[hit,"rtmed"]) # Return each hit as its: Tag, M+ Tag, Rt, M+ Rt##Modified to include Mz vals
      }
    })
  })

  hitlist<-unlist(hitlist,recursive=F) ##Flatten hitlist
  ## Determine Null elements in hitlist and set them to NULL
  isnull<-unlist(lapply (hitlist,function(x)
  {
    is.null(x)
  }))
  hitlist[isnull]<-NULL
  #   print(hitlist)
  output<-unique(do.call(rbind,hitlist)) # Make a matrix by setting each element of the flattened hitlist as a row???
  return(output)
}

#----------

C13peaks_GroupPairing_tTest_Catabolism<-function(xcmsSet2,mzppm=15,mzabs=0.005,rterror=1,resultspath=NULL,pheno=NULL,preA=0.05,postA=0.05,value="maxo",n13C=6)
{
  groups = as.data.frame(xcmsSet2@groups) #read xcmspeaks
  groups = cbind(groupnames(xcmsSet2),groups)
  colnames(groups) = c("groupname",colnames(groups)[-1]) #clean up column names
  sampcls = levels(sampclass(xcmsSet2))
  filledSet = fillPeaks(xcmsSet2,"chrom",nSlaves=8,expand.mz=1,expand.rt=1,min.width.mz=0.005,min.width.rt=1)
  gpvals = groupval(filledSet,method="medret",value=value)
  #get names of 12C and 13C samples
  sampleNames = rownames(filledSet@phenoData)
  lastCol = ncol(filledSet@phenoData)
  C12sampleNames = sampleNames[grepl("12C",filledSet@phenoData[ ,lastCol])] #may need to change  "12C" to user defined identifier
  C13sampleNames = sampleNames[grepl("13C",filledSet@phenoData[ ,lastCol])] #may need to change  "13C" to user defined identifier
  #If unique folders were not used for the 12C and 13C samples, use the file name
  if(length(C12sampleNames)==0)
  {
    C12sampleNames = sampleNames[grepl("12C",sampnames(filledSet))] #may need to change  "12C" to user defined identifier
    C13sampleNames = sampleNames[grepl("13C",sampnames(filledSet))] #may need to change  "13C" to user defined identifier
  }

  #   rownames(groups)<-groupnames(xcmsSet2) ##set group Id tags as row names
  groupsplusNmzrange<-groupmzrange(groups,mzppm,mzabs,nC13=6) #add uppermz and lowermz of M+6 for each group
  plusN<-findC13cmpd_stat(groupsplusNmzrange,rterror) #find M+6 C13 peaks
#   plusNout <<-plusN

  
  ###Filter out spurrious M peaks using an unpaired t-test with unequal variances.
  ###NOTE: This will remove any valid compounds whose +6 peak coelutes with another compound of identical mass.
  ###      Also, may want to perform before the parameter optimization step

  #12C_base > 12C_6??
  t1 = c(iso="M",sampType="12C")
  t2 = c(iso="Mplus",sampType="12C")
  idx = tFilterPreProcess(filledSet,plusN,resultsPath,pheno,preA,H0="greater",value=value,cout=T,fout=T,returnType=1,vect1=t1,vect2=t2)

  #13C_base > 12C_base??
  t1 = c(iso="M",sampType="13C")
  t2 = c(iso="M",sampType="12C")
  kpidx1 = tFilterPreProcess(filledSet,plusN,resultsPath,pheno,preA,H0="greater",value=value,cout=T,fout=T,returnType=1,vect1=t1,vect2=t2)

  #13C_6 > 12C_6 ??
  t1 = c(iso="Mplus",sampType="13C")
  t2 = c(iso="Mplus",sampType="12C")
  kpidx2 = tFilterPreProcess(filledSet,plusN,resultsPath,pheno,preA,H0="greater",value=value,cout=T,fout=T,returnType=1,vect1=t1,vect2=t2)

  #Determine which groups are +6, +12 pairs and prevent them from being removed
  kpidx = intersect(kpidx1,kpidx2)
  kpidx = union(idx,kpidx)

  #remove invalid pairings
  plusN=plusN[kpidx,]
  rm(kpidx,idx)


#   t1 = c(iso="M",sampType="12C")
#   t2 = c(iso="Mplus",sampType="12C")
  C12peak<-plusN[plusN[,"M"] %in% setdiff(plusN[,"M"],plusN[,"Mplus"]),] #find C12 (M) peaks  ##returns the elements of "M" that are not found in "Mplus"
  ##Find +12 and up peaks, remove peaks from M list that are actually Mplus peaks
  count = 0
  holder = mat.or.vec(0,1)
  grouped<-lapply(c(1:nrow(C12peak)),function(x)
  {
    Mplus<-C12peak[x,"Mplus"]

    idx<-vector()
    gp<-vector()
    mult=F

    ## Filter out peaks that were found as "M" but are actually an M+6 peak
    while (length(which(plusN[,"M"]==Mplus))>0) ##For each peak in "M" that matches the "Mplus" of the current peak
    {
      if(length(which(plusN[,"M"]==Mplus))>=2) #if there are multiple dupes, print diagnositc info
      {
        mult = T
        cat("M",C12peak[x,"M"],"\t","M+",Mplus,"\n")
        cat("Dupe Index: ",which(plusN[,"M"]==Mplus),"\n")
        print(plusN[which(plusN[,"M"]==Mplus)])
        cat("idx_a: ",idx,"\n")
      }
      idx<-c(idx,which(plusN[,"M"]==Mplus)) ##Concatenate to idx
      if(mult)#if there are multiple dupes, print diagnositc info
        cat("idx_b: ",idx,"\n")
      Mplus<-plusN[plusN[,"M"]==Mplus,"Mplus"] ##Change Mplus to the "Mplus" peak of one of the M+6 peaks just found in the M column
      if(mult) #if there are multiple dupes, print diagnositc info
        cat("NewM+: ",Mplus,"\n\n")
      mult=F
    }
    l<-length(idx) #L
    if (l>0)
    {
      rtmin<-min(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(plusN[idx,"rtmed_M"]),as.numeric(plusN[idx[l],"rtmed_Mplus"])) ##idx[L], calculate min rt for group
      rtmax<-max(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(plusN[idx,"rtmed_M"]),as.numeric(plusN[idx,"rtmed_Mplus"])) ##idx[L], calculate max rt for group
      gp<-c(C12peak[x,"M"],unique(plusN[idx,"M"]),plusN[idx[l],"Mplus"],C12peak[x,"mzmed_M"],unique(plusN[idx,"mzmed_M"]),plusN[idx[l],"mzmed_Mplus"],rtmin,rtmax) #group peaks for return, Group is base peak, all peaks in Idx, and the +6 Peak of the last peak in Idx
      nc13gps = length(unique(plusN[idx,"M"]))
      names(gp)<-c("12C",paste("13C_",seq(n13C,(nc13gps+1)*n13C,n13C),sep=""),"12C_mz",paste("13C_",seq(n13C,(nc13gps+1)*n13C,n13C),"_mz",sep=""),"rtmedmin","rtmedmax") #add identifiers to each peak in the group and label columns
    } else
    {
      rtmin<-min(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      rtmax<-max(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      gp<-c(C12peak[x,"M"],C12peak[x,"Mplus"],C12peak[x,"mzmed_M"],C12peak[x,"mzmed_Mplus"],rtmin,rtmax)
      names(gp)<-c("12C",paste("13C_",n13C,sep=""),"12C_mz",paste("13C_",n13C,"_mz",sep=""),"rtmedmin","rtmedmax")
    }
    #If 13C_base peak is significantly greater than the 12C_base peak, discard the peak group
    if(singleTtest(filledSet,C12peak[x,"M"],C12peak[x,"M"],H0="greater",value=value,sampType1="13C",sampType2="12C",gpvals)<=postA)
    {
#       count=count+1
#       holder[count]=C12peak[x,"M"]
#       cat("\nC12:",C12peak[x,"C12avg_M"],"\tC13:",C12peak[x,"C13avg_M"],"\tC12Name:",C12peak[x,"M"],"Count:",count,"\n")
    }
    else
      return(gp)
  })
#   rm(count)
#   holder <<- holder
#   rm(holder)
  colnames<-unique(unlist(c(sapply(grouped,names)))) ## shorter version unique(unlist(sapply(grouped,names))), determine all of the unique col names
  out<-as.data.frame(do.call(rbind,lapply(grouped,"[",colnames)),stringsAsFactors = FALSE) #turn grouped into a dataframe with a column for each label in colnames
  names(out)<-colnames #fix column labels
  out<-out[rowSums(is.na(out))!=ncol(out),] #remove all "NA" rows
  out<-out[,c(setdiff(names(out),c("rtmedmin","rtmedmax")),"rtmedmin","rtmedmax")] #reorder columns
  clusterNumber = 1:length(out[ ,1])
  out=cbind(clusterNumber,out)
  write.csv(out,file=file.path(resultsPath,paste(pheno,"Clusters_GroupNamesOnly.csv",sep="_")),row.names=F)
  cat("\nBrief Output\n")
  print(head(out))
  #   out <<-out
  #matrix format
  ##apply to each peak group (row)
  #   groupsplusNmzrange <<- groupsplusNmzrange
  mat<-lapply(c(1:nrow(out)),function(x)
  {
    pkname<-out[x,-c(ncol(out)-1,ncol(out))] ##leave off last two cols (rt min and max) -> peak names
    #cat("\nPeakname:",pkname,"\n")
    idx_a<-which(!is.na(pkname)) ##return indicies of non-NULL col names in group (peak slots that have a match, ie. 12c, 13C_6, ect.)
    idx = grep("clusterNumber|_mz",names(pkname)[idx_a],invert=T) ##return the indicies for the name slots (as opposed to _mz slots) of non-Null Peaks
    #print(head(unlist(pkname[idx])))
    groups_idx = which(groupsplusNmzrange[ ,"groupname"] %in% unlist(pkname[idx]))
    cbind(clusterNumber=rep(pkname[1,1],length(idx)),iso=names(pkname)[idx], groupsplusNmzrange[groups_idx,c("groupname","mzmed","rtmed","mzmin","mzmax","rtmin","rtmax","npeaks",sampcls,"lowermz","uppermz")])
  })
  mat<-do.call(rbind,mat) #bind previously generated cols together
  cat("\nBrief Output Matrix Form\n")
  print(head(mat))
  write.csv(mat,file=file.path(resultsPath,paste(pheno,"Clusters_MatrixForm.csv",sep="_")),row.names=F)
}





C13peaks_GroupPairing_Heuristic_Catabolism<-function(xcmsSet2,mzppm=15,mzabs=0.005,rterror=1,resultspath=NULL,pheno=NULL,C12Fchng=1.25,C13Fchng=1.40,value="maxo",n13C=6)
{
  groups = as.data.frame(xcmsSet2@groups) #read xcmspeaks
  groups = cbind(groupnames(xcmsSet2),groups)
  colnames(groups) = c("groupname",colnames(groups)[-1]) #clean up column names
  sampcls = levels(sampclass(xcmsSet2))
  filledSet = fillPeaks(xcmsSet2,"chrom",nSlaves=8,expand.mz=1,expand.rt=1,min.width.mz=0.005,min.width.rt=1)
  gpvals = groupval(filledSet,method="medret",value=value)
  #get names of 12C and 13C samples
  sampleNames = rownames(xcmsSet2@phenoData)
  lastCol = ncol(xcmsSet2@phenoData)
  C12sampleNames = sampleNames[grepl("12C",xcmsSet2@phenoData[ ,lastCol])] #may need to change  "12C" to user defined identifier
  C13sampleNames = sampleNames[grepl("13C",xcmsSet2@phenoData[ ,lastCol])] #may need to change  "13C" to user defined identifier
  #If unique folders were not used for the 12C and 13C samples, use the file name
  if(length(C12sampleNames)==0)
  {
    C12sampleNames = sampleNames[grepl("12C",sampnames(xcmsSet2))] #may need to change  "12C" to user defined identifier
    C13sampleNames = sampleNames[grepl("13C",sampnames(xcmsSet2))] #may need to change  "13C" to user defined identifier
  }
  #need to consider how to handle this senario in the case of multiple genotpes
  C12avg = rowMeans(gpvals[,C12sampleNames,drop=F],na.rm=T)
  C13avg = rowMeans(gpvals[,C13sampleNames,drop=F],na.rm=T)
  groups = cbind(groups,C12avg,C13avg)

  #   rownames(groups)<-groupnames(xcmsSet2) ##set group Id tags as row names
  groupsplusNmzrange<-groupmzrange(groups,mzppm,mzabs,n13C) #add uppermz and lowermz of M+6 for each group
  plusN<-findC13cmpd(groupsplusNmzrange,rterror) #find M+6 C13 peaks
#   plusNout <<-plusN

  
  #filter out spurrious M peaks using a maximum fold change of 1.  May want to also consider using statistics.
  #NOTE: This will remove any valid compounds whose +6 peak coelutes with another compound of the same mass.
  #      Also, may want to perform before the parameter optimization step
  rmidx = which(as.numeric(plusN[,"C12avg_M"])<as.numeric(plusN[,"C12avg_Mplus"]))
  #Added fold change of 1.25 for threshold
  kpidx = which(as.numeric(plusN[,"C12avg_M"])*C12Fchng<as.numeric(plusN[,"C13avg_M"]) & as.numeric(plusN[,"C12avg_Mplus"])*C12Fchng<as.numeric(plusN[,"C13avg_Mplus"]))  #Guard against removal of valid 13C_12 peaks, if 13C signal is greater than 12C singal by > 12C*fchng at both masses, keep peak
  rmidx = setdiff(rmidx,kpidx)
  plusN=plusN[-rmidx,]
  # plusNOut <<-plusN

  C12peak<-plusN[plusN[,"M"] %in% setdiff(plusN[,"M"],plusN[,"Mplus"]),] #find C12 (M) peaks  ##returns the elements of "M" that are not found in "Mplus"
  ##Find +12 and up peaks, remove peaks from M list that are actually Mplus peaks
  count = 0
  holder = mat.or.vec(0,1)
  grouped<-lapply(c(1:nrow(C12peak)),function(x)
  {
    Mplus<-C12peak[x,"Mplus"]

    idx<-vector()
    gp<-vector()
    mult=F

    ## Filter out peaks that were found as "M" but are actually an M+6 peak
    while (length(which(plusN[,"M"]==Mplus))>0) ##For each peak in "M" that matches the "Mplus" of the current peak
    {
      if(length(which(plusN[,"M"]==Mplus))>=2) #if there are multiple dupes, print diagnositc info
      {
        mult = T
        cat("M",C12peak[x,"M"],"\t","M+",Mplus,"\n")
        cat("Dupe Index: ",which(plusN[,"M"]==Mplus),"\n")
        print(plusN[which(plusN[,"M"]==Mplus)])
        cat("idx_a: ",idx,"\n")
      }
      idx<-c(idx,which(plusN[,"M"]==Mplus)) ##Concatenate to idx
      if(mult)#if there are multiple dupes, print diagnositc info
        cat("idx_b: ",idx,"\n")
      Mplus<-plusN[plusN[,"M"]==Mplus,"Mplus"] ##Change Mplus to the "Mplus" peak of one of the M+6 peaks just found in the M column
      if(mult) #if there are multiple dupes, print diagnositc info
        cat("NewM+: ",Mplus,"\n\n")
      mult=F
    }
    l<-length(idx) #L
    if (l>0)
    {
      rtmin<-min(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(plusN[idx,"rtmed_M"]),as.numeric(plusN[idx[l],"rtmed_Mplus"])) ##idx[L], calculate min rt for group
      rtmax<-max(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(plusN[idx,"rtmed_M"]),as.numeric(plusN[idx,"rtmed_Mplus"])) ##idx[L], calculate max rt for group
      gp<-c(C12peak[x,"M"],unique(plusN[idx,"M"]),plusN[idx[l],"Mplus"],C12peak[x,"mzmed_M"],unique(plusN[idx,"mzmed_M"]),plusN[idx[l],"mzmed_Mplus"],rtmin,rtmax) #group peaks for return, Group is base peak, all peaks in Idx, and the +6 Peak of the last peak in Idx
      nc13gps = length(unique(plusN[idx,"M"]))
      names(gp)<-c("12C",paste("13C_",seq(n13C,(nc13gps+1)*n13C,n13C),sep=""),"12C_mz",paste("13C_",seq(n13C,(nc13gps+1)*n13C,n13C),"_mz",sep=""),"rtmedmin","rtmedmax") #add identifiers to each peak in the group and label columns
    } else
    {
      rtmin<-min(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      rtmax<-max(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      gp<-c(C12peak[x,"M"],C12peak[x,"Mplus"],C12peak[x,"mzmed_M"],C12peak[x,"mzmed_Mplus"],rtmin,rtmax)
      names(gp)<-c("12C",paste("13C_",n13C,sep=""),"12C_mz",paste("13C_",n13C,"_mz",sep=""),"rtmedmin","rtmedmax")
    }
    if(as.numeric(C12peak[x,"C13avg_M"])>as.numeric(C12peak[x,"C12avg_M"])*C13Fchng) #If 12C_C12 peak is less than the 13C_C12 (with a 40% error tolerance as default), discard the peak group
    {
#       count=count+1
#       holder[count]=C12peak[x,"M"]
#       count <<-count
#       cat("\nC12:",C12peak[x,"C12avg_M"],"\tC13:",C12peak[x,"C13avg_M"],"\tC12Name:",C12peak[x,"M"],"Count:",count,"\n")
    }
    else
      return(gp)
  })
#   rm(count)
#   holder <<- holder
#   rm(holder)
  colnames<-unique(unlist(c(sapply(grouped,names)))) ## shorter version unique(unlist(sapply(grouped,names))), determine all of the unique col names
  out<-as.data.frame(do.call(rbind,lapply(grouped,"[",colnames)),stringsAsFactors = FALSE) #turn grouped into a dataframe with a column for each label in colnames
  names(out)<-colnames #fix column labels
  out<-out[rowSums(is.na(out))!=ncol(out),] #remove all "NA" rows
  out<-out[,c(setdiff(names(out),c("rtmedmin","rtmedmax")),"rtmedmin","rtmedmax")] #reorder columns
  clusterNumber = 1:length(out[ ,1])
  out=cbind(clusterNumber,out)
  write.csv(out,file=file.path(resultsPath,paste(pheno,"Clusters_GroupNamesOnly.csv",sep="_")),row.names=F)
  cat("\nBrief Output\n")
  print(head(out))
  #   out <<-out
  #matrix format
  ##apply to each peak group (row)
  #   groupsplusNmzrange <<- groupsplusNmzrange
  mat<-lapply(c(1:nrow(out)),function(x)
  {
    pkname<-out[x,-c(ncol(out)-1,ncol(out))] ##leave off last two cols (rt min and max) -> peak names
    #cat("\nPeakname:",pkname,"\n")
    idx_a<-which(!is.na(pkname)) ##return indicies of non-NULL col names in group (peak slots that have a match, ie. 12c, 13C_6, ect.)
    idx = grep("clusterNumber|_mz",names(pkname)[idx_a],invert=T) ##return the indicies for the name slots (as opposed to _mz slots) of non-Null Peaks
    #print(head(unlist(pkname[idx])))
    groups_idx = which(groupsplusNmzrange[ ,"groupname"] %in% unlist(pkname[idx]))
    cbind(clusterNumber=rep(pkname[1,1],length(idx)),iso=names(pkname)[idx], groupsplusNmzrange[groups_idx,c("groupname","mzmed","rtmed","mzmin","mzmax","rtmin","rtmax","npeaks",sampcls,"lowermz","uppermz")])
  })
  mat<-do.call(rbind,mat) #bind previously generated cols together
  cat("\nBrief Output Matrix Form\n")
  print(head(mat))
  write.csv(mat,file=file.path(resultsPath,paste(pheno,"Clusters_MatrixForm.csv",sep="_")),row.names=F)
}



