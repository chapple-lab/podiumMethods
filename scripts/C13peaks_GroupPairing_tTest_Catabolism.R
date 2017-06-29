labeledPeaks_GroupPairing_tTest <-
function(xcmsSet2,mLabel=1.0033548,nLabel=6,mzppm=15,mzabs=0.005,rterror=1,resultspath=NULL,phenoTag=NULL,phenotypes=NULL,preA=0.05,postA=0.05,value="maxo",overide=F)
{
  check = is.na(xcmsSet2@filled[1])
  if(check&!overide)
  {
    stop("\nError: xcmsSet not filled. Run fillPeaks() on xcmsSet first.\nTo run without filling, use overide=TRUE.\n")
  }else if(overide&check)
  {
    cat("\nWarning: xcmsSet appears to not be filled\n")
  }

  groups = as.data.frame(xcmsSet2@groups) #read xcmspeaks
  groups = cbind(groupnames(xcmsSet2),groups)
  colnames(groups) = c("groupname",colnames(groups)[-1]) #clean up column names
  sampcls = levels(sampclass(xcmsSet2))
  # filledSet = fillPeaks(xcmsSet2,"chrom",nSlaves=8,expand.mz=1,expand.rt=1,min.width.mz=0.005,min.width.rt=1)
  gpvals = groupval(xcmsSet2,method="medret",value=value)
  #   rownames(groups)<-groupnames(xcmsSet2) ##set group Id tags as row names
  groupsplusNmzrange<-groupmzrange(groups,mzppm,mzabs,nLabel=nLabel,mLabel=mLabel) #add uppermz and lowermz of M+6 for each group
  plusN<-findLablcmpd_stat(groupsplusNmzrange,rterror) #find M+6 C13 peaks


  ##--parameter optimization
  groupsplus2Nmzrange<-groupmzrange(groups,mzppm,mzabs,nC13=12) #add uppermz and lowermz of M+12 for each peak
  plus2N<-findLablcmpd_stat(groupsplus2Nmzrange,rterror)
  #optimize mzppm and mzabs so serialpath is identical to parallelpath
  conflictpair6<-sum(duplicated(plusN[,"M"]))+sum(duplicated(plusN[,"Mplus"])) ##determine number of peaks that occur multiple times in the same list (M list or M+6), does NOT look at multiple occurrences across lists  (ie. One in both M and M+6 does NOT count)
  conflictpair12<-sum(duplicated(plus2N[,"M"]))+sum(duplicated(plus2N[,"Mplus"])) ##repeat above for +12 path
  cat("conflictpair6: ",conflictpair6,"\n","conflictpair12: ",conflictpair12, "\n")
  cat("If any conflictpair is non-zero, reduce rterror or mzerror. \n")
  serialpath<-plusN[which(plusN[,"Mplus"] %in% intersect(plusN[,"M"],plusN[,"Mplus"])),"M"] #M -> M+6 -> M+12 ##Return M peaks for Mplus peaks that are found in both M and Mplus cols (should represent M->M6->m12). Will fail if one M peak mapps to two distinct Mplus peaks
  parallelpath<-intersect(plusN[,"M"],plus2N[,"M"]) #M -> M+6; M -> M+12 ##Returns all M peaks that are in both plusN and plus2N M list
  cat("Results from two paths are identical?", setequal(serialpath,parallelpath),"\n") #If FALSE, increase the rterror or mzerror.
  ##---end parameter optimization

  ###Filter out spurrious M peaks using an unpaired t-test with (OLD:unequal) equal variances.

  #13C_6 > 12C_6 ??
  t1 = c(peakType="Mplus",iso="13C")
  t2 = c(peakType="Mplus",iso="12C")
  kpidx = tFilterPreProcess(xcmsSet2,plusN,resultsPath,phenoTag,phenotypes,preA,alternative="greater",value=value,cout=T,fout=T,returnType=2,vect1=t1,vect2=t2)

  #remove invalid pairings
  plusN=plusN[kpidx,]
  rm(kpidx)


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
      nLablgps = length(unique(plusN[idx,"M"]))
      names(gp)<-c("12C",paste("13C_",seq(nLabel,(nLablgps+1)*nLabel,nLabel),sep=""),"12C_mz",paste("13C_",seq(nLabel,(nLablgps+1)*nLabel,nLabel),"_mz",sep=""),"rtmedmin","rtmedmax") #add identifiers to each peak in the group and label columns
    } else
    {
      rtmin<-min(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      rtmax<-max(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      gp<-c(C12peak[x,"M"],C12peak[x,"Mplus"],C12peak[x,"mzmed_M"],C12peak[x,"mzmed_Mplus"],rtmin,rtmax)
      names(gp)<-c("12C",paste("13C_",nLabel,sep=""),"12C_mz",paste("13C_",nLabel,"_mz",sep=""),"rtmedmin","rtmedmax")
    }
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
  write.csv(out,file=file.path(resultsPath,paste(phenoTag,"nLabel",nLabel,"Clusters_GroupNamesOnly.csv",sep="_")),row.names=F)
  cat("\nBrief Output\n")
  print(head(out))
  #   out <<-out
  #matrix format
  ##apply to each peak group (row)
  #   groupsplus6mzrange <<- groupsplus6mzrange
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
  write.csv(mat,file=file.path(resultsPath,paste(phenoTag,"nLabel",nLabel,"Clusters_MatrixForm.csv",sep="_")),row.names=F)
}
