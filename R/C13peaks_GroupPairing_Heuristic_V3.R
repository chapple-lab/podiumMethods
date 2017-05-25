C13peaks_GroupPairing_Heuristic_V3 <-
function(xcmsSet2,mzppm=15,mzabs=0.005,rterror=1,resultspath=NULL,pheno=NULL,C12Fchng=1.25,C13Fchng=1.40,value="maxo")
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
  groupsplus6mzrange<-groupmzrange(groups,mzppm,mzabs,nC13=6) #add uppermz and lowermz of M+6 for each group
  plus6<-findC13cmpd(groupsplus6mzrange,rterror) #find M+6 C13 peaks
  #   plus6out <<-plus6

  ##--parameter optimization
  groupsplus12mzrange<-groupmzrange(groups,mzppm,mzabs,nC13=12) #add uppermz and lowermz of M+12 for each peak
  plus12<-findC13cmpd(groupsplus12mzrange,rterror)
  #optimize mzppm and mzabs so serialpath is identical to parallelpath
  conflictpair6<-sum(duplicated(plus6[,"M"]))+sum(duplicated(plus6[,"Mplus"])) ##determine number of peaks that occur multiple times in the same list (M list or M+6), does NOT look at multiple occurrences across lists  (ie. One in both M and M+6 does NOT count)
  conflictpair12<-sum(duplicated(plus12[,"M"]))+sum(duplicated(plus12[,"Mplus"])) ##repeat above for +12 path
  cat("conflictpair6: ",conflictpair6,"\n","conflictpair12: ",conflictpair12, "\n")
  cat("If any conflictpair is non-zero, reduce rterror or mzerror. \n")
  serialpath<-plus6[which(plus6[,"Mplus"] %in% intersect(plus6[,"M"],plus6[,"Mplus"])),"M"] #M -> M+6 -> M+12 ##Returns the M+6 peaks that are also in the M peak list, think there might be a mistake here
  parallelpath<-intersect(plus6[,"M"],plus12[,"M"]) #M -> M+6; M -> M+12 ##Returns all M peaks that are in both plus6 and plus12 M list
  cat("Results from two paths are identical?", setequal(serialpath,parallelpath),"\n") #If FALSE, increase the rterror or mzerror.
  ##---end parameter optimization

  #filter out spurrious M peaks using a maximum fold change of 1.  May want to also consider using statistics.
  #NOTE: This will remove any valid compounds whose +6 peak coelutes with another compound of the same mass.
  #      Also, may want to perform before the parameter optimization step
  rmidx = which(as.numeric(plus6[,"C13avg_Mplus"])*C13Fchng<as.numeric(plus6[,"C12avg_Mplus"]))
  #Added fold change of 1.25 for threshold
  # kpidx = which(as.numeric(plus6[,"C12avg_M"])*C12Fchng<as.numeric(plus6[,"C13avg_M"]) & as.numeric(plus6[,"C12avg_Mplus"])*C12Fchng<as.numeric(plus6[,"C13avg_Mplus"]))  #Guard against removal of valid 13C_12 peaks, if 13C signal is greater than 12C singal by > 12C*fchng at both masses, keep peak
  # rmidx = setdiff(rmidx,kpidx)
  plus6=plus6[-rmidx,]
  plus6Out<<-plus6

  C12peak<-plus6[plus6[,"M"] %in% setdiff(plus6[,"M"],plus6[,"Mplus"]),] #find C12 (M) peaks  ##returns the elements of "M" that are not found in "Mplus"
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
    while (length(which(plus6[,"M"]==Mplus))>0) ##For each peak in "M" that matches the "Mplus" of the current peak
    {
      if(length(which(plus6[,"M"]==Mplus))>=2) #if there are multiple dupes, print diagnositc info
      {
        mult = T
        cat("M",C12peak[x,"M"],"\t","M+",Mplus,"\n")
        cat("Dupe Index: ",which(plus6[,"M"]==Mplus),"\n")
        print(plus6[which(plus6[,"M"]==Mplus)])
        cat("idx_a: ",idx,"\n")
      }
      idx<-c(idx,which(plus6[,"M"]==Mplus)) ##Concatenate to idx
      if(mult)#if there are multiple dupes, print diagnositc info
        cat("idx_b: ",idx,"\n")
      Mplus<-plus6[plus6[,"M"]==Mplus,"Mplus"] ##Change Mplus to the "Mplus" peak of one of the M+6 peaks just found in the M column
      if(mult) #if there are multiple dupes, print diagnositc info
        cat("NewM+: ",Mplus,"\n\n")
      mult=F
    }
    l<-length(idx) #L
    if (l>0)
    {
      rtmin<-min(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(plus6[idx,"rtmed_M"]),as.numeric(plus6[idx[l],"rtmed_Mplus"])) ##idx[L], calculate min rt for group
      rtmax<-max(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(plus6[idx,"rtmed_M"]),as.numeric(plus6[idx,"rtmed_Mplus"])) ##idx[L], calculate max rt for group
      gp<-c(C12peak[x,"M"],unique(plus6[idx,"M"]),plus6[idx[l],"Mplus"],C12peak[x,"mzmed_M"],unique(plus6[idx,"mzmed_M"]),plus6[idx[l],"mzmed_Mplus"],rtmin,rtmax) #group peaks for return, Group is base peak, all peaks in Idx, and the +6 Peak of the last peak in Idx
      nc13gps = length(unique(plus6[idx,"M"]))
      names(gp)<-c("12C",paste("13C_",seq(6,(nc13gps+1)*6,6),sep=""),"12C_mz",paste("13C_",seq(6,(nc13gps+1)*6,6),"_mz",sep=""),"rtmedmin","rtmedmax") #add identifiers to each peak in the group and label columns
    } else
    {
      rtmin<-min(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      rtmax<-max(as.numeric(C12peak[x,"rtmed_M"]),as.numeric(C12peak[x,"rtmed_Mplus"]))
      gp<-c(C12peak[x,"M"],C12peak[x,"Mplus"],C12peak[x,"mzmed_M"],C12peak[x,"mzmed_Mplus"],rtmin,rtmax)
      names(gp)<-c("12C","13C_6","12C_mz","13C_6_mz","rtmedmin","rtmedmax")
    }
    # if(as.numeric(C12peak[x,"C13avg_M"])>as.numeric(C12peak[x,"C12avg_M"])*C13Fchng) #If 12C_C12 peak is less than the 13C_C12 (with a 40% error tolerance as default), discard the peak group
    # {
      #       count=count+1
      #       holder[count]=C12peak[x,"M"]
      #       count <<-count
      #       cat("\nC12:",C12peak[x,"C12avg_M"],"\tC13:",C12peak[x,"C13avg_M"],"\tC12Name:",C12peak[x,"M"],"Count:",count,"\n")
    # }
    # else
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
  #   groupsplus6mzrange <<- groupsplus6mzrange
  mat<-lapply(c(1:nrow(out)),function(x)
  {
    pkname<-out[x,-c(ncol(out)-1,ncol(out))] ##leave off last two cols (rt min and max) -> peak names
    #cat("\nPeakname:",pkname,"\n")
    idx_a<-which(!is.na(pkname)) ##return indicies of non-NULL col names in group (peak slots that have a match, ie. 12c, 13C_6, ect.)
    idx = grep("clusterNumber|_mz",names(pkname)[idx_a],invert=T) ##return the indicies for the name slots (as opposed to _mz slots) of non-Null Peaks
    #print(head(unlist(pkname[idx])))
    groups_idx = which(groupsplus6mzrange[ ,"groupname"] %in% unlist(pkname[idx]))
    cbind(clusterNumber=rep(pkname[1,1],length(idx)),iso=names(pkname)[idx], groupsplus6mzrange[groups_idx,c("groupname","mzmed","rtmed","mzmin","mzmax","rtmin","rtmax","npeaks",sampcls,"lowermz","uppermz")])
  })
  mat<-do.call(rbind,mat) #bind previously generated cols together
  cat("\nBrief Output Matrix Form\n")
  print(head(mat))
  write.csv(mat,file=file.path(resultsPath,paste(pheno,"Clusters_MatrixForm.csv",sep="_")),row.names=F)
}
