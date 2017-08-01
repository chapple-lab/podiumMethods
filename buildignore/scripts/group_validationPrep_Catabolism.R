group_composite_Catabolism = function(xcmsSet2=NULL,pheno=NULL,resultsPath=NULL,n13C=6)
{
  if(type != "suspect" && type != "valid")
    warning("Type does not match accepted inputs \"suspect\" or \"valid\".\n")

  pairedGroups = read.csv(file.path(resultsPath,paste(pheno,"Clusters_GroupNamesOnly.csv",sep="_")),check.names=F,stringsAsFactors=F)
  pairedGroups_matrix = read.csv(file.path(resultsPath,paste(pheno,"Clusters_MatrixForm.csv",sep="_")),check.names=F,stringsAsFactors=F)

  #prepair for statistical analysis by calculating and storing metrics in a data frame
  #get groups
  cat("\n\nCreating Composite Statistical Profile\n")
  groupNames = groupnames(xcmsSet2)
  C12Groups = which(groupNames %in% pairedGroups[ ,"12C"])
  C13Groups = which(groupNames %in% pairedGroups[ ,paste("13C_",n13C,sep="")])
  #get peaks in each group
  C12groupList = xcmsSet2@groupidx[C12Groups]
  C13groupList = xcmsSet2@groupidx[C13Groups]

  ##create data frame of all peaks from each group
  peakData = mat.or.vec(0,1)
  groupidx=1
  groupHolder = list()
  for(x in C12groupList)
  {
    peakData = rbind(peakData,xcmsSet2@peaks[x, ])
    groupName = groupNames[C12Groups[groupidx]]
    groupHolder[length(groupHolder)+1] = list(rep(groupName,length(x)))
    groupidx = groupidx+1
  }

  groupTags = unlist(groupHolder)
  iso = rep("12C",length(groupTags))
  peakData = data.frame(groupTags,iso,peakData,stringsAsFactors=F)

  peakData2 = mat.or.vec(0,1)
  groupidx=1
  groupHolder = list()
  for(x in C13groupList)
  {
    peakData2 = rbind(peakData2,xcmsSet2@peaks[x, ])
    groupHolder[length(groupHolder)+1] = list(rep(groupNames[C13Groups[groupidx]],length(x)))
    groupidx = groupidx+1
  }

  groupTags = unlist(groupHolder)
  iso = rep(paste("13C_",n13C,sep=""),length(groupTags))
  peakData2 = data.frame(groupTags,iso,peakData2,stringsAsFactors=F)

  peakData = merge(peakData,peakData2,all=T)
  write.csv(peakData,file=file.path(resultsPath,paste(pheno,"Clusters_AllPeaksMatrix.csv",sep="_")), row.names=F)

  ##Compute mean of each stat across each group's peak set to create a composite statistical profile for each group
  allGroups = union(C12Groups,C13Groups)
  allGroups = groupNames[allGroups]

  meanHolder = list()
  isoHolder = list()
  groupComposite = mat.or.vec(0,1)
  for(x in allGroups)
  {
    group = peakData[which(peakData$groupTags==x), ]
    isoHolder = append(isoHolder,group[1,"iso"])
    for(z in (3:(length(group)-1)))
    {
      meanHolder[length(meanHolder)+1] = mean(group[ ,z])
    }
    groupComposite = rbind(groupComposite,unlist(meanHolder))
    meanHolder = list()
  }
  iso = unlist(isoHolder)
  groupComposite = data.frame(allGroups,iso,groupComposite,stringsAsFactors=F)
  #adjust col names
  names = colnames(peakData[c(-1,-2,-length(peakData))])
  names = unlist(lapply(names,"paste","mean",sep="_"))
  names = c("group","iso",names)
  colnames(groupComposite) = names
  #variable clean up
  rm(meanHolder)
  rm(isoHolder)
  rm(groupHolder)
  write.csv(groupComposite,file=file.path(resultsPath,paste(pheno,"GroupComposite_unfilled.csv",sep="_")),row.names=F)

}



group_validationPrep_Catabolism = function(xcmsSet2=NULL,pheno=NULL,resultsPath=NULL,intThresh=100, type=c("valid","suspect"),snThresh=5,widthThresh=9,ppmThresh=35,n13C=6)
{
    if(type != "suspect" && type != "valid")
		warning("Type does not match accepted inputs \"suspect\" or \"valid\".\n")

    pairedGroups = read.csv(file.path(resultsPath,paste(pheno,"Clusters_GroupNamesOnly.csv",sep="_")),check.names=F,stringsAsFactors=F)
    pairedGroups_matrix = read.csv(file.path(resultsPath,paste(pheno,"Clusters_MatrixForm.csv",sep="_")),check.names=F,stringsAsFactors=F)


    #prepair for statistical analysis by calculating and storing metrics in a data frame
    #get groups
    cat("\n\nCreating Composite Statistical Profile\n")
    groupNames = groupnames(xcmsSet2)
    C12Groups = which(groupNames %in% pairedGroups[ ,"12C"])
    C13Groups = which(groupNames %in% pairedGroups[ ,paste("13C_",n13C,sep="")])
    #get peaks in each group
    C12groupList = xcmsSet2@groupidx[C12Groups]
    C13groupList = xcmsSet2@groupidx[C13Groups]

    ##create data frame of all peaks from each group
    peakData = mat.or.vec(0,1)
    groupidx=1
    groupHolder = list()
    for(x in C12groupList)
    {
      peakData = rbind(peakData,xcmsSet2@peaks[x, ])
      groupName = groupNames[C12Groups[groupidx]]
      groupHolder[length(groupHolder)+1] = list(rep(groupName,length(x)))
      groupidx = groupidx+1
    }

    groupTags = unlist(groupHolder)
    iso = rep("12C",length(groupTags))
    peakData = data.frame(groupTags,iso,peakData,stringsAsFactors=F)

    peakData2 = mat.or.vec(0,1)
    groupidx=1
    groupHolder = list()
    for(x in C13groupList)
    {
      peakData2 = rbind(peakData2,xcmsSet2@peaks[x, ])
      groupHolder[length(groupHolder)+1] = list(rep(groupNames[C13Groups[groupidx]],length(x)))
      groupidx = groupidx+1
    }

    groupTags = unlist(groupHolder)
    iso = rep(paste("13C_",n13C,sep=""),length(groupTags))
    peakData2 = data.frame(groupTags,iso,peakData2,stringsAsFactors=F)

    peakData = merge(peakData,peakData2,all=T)
    write.csv(peakData,file=file.path(resultsPath,paste(pheno,"Clusters_AllPeaksMatrix.csv",sep="_")), row.names=F)

    ##Compute mean of each stat across each group's peak set to create a composite statistical profile for each group
    allGroups = union(C12Groups,C13Groups)
    allGroups = groupNames[allGroups]

    meanHolder = list()
    isoHolder = list()
    groupComposite = mat.or.vec(0,1)
    for(x in allGroups)
    {
      group = peakData[which(peakData$groupTags==x), ]
      isoHolder = append(isoHolder,group[1,"iso"])
      for(z in (3:(length(group)-1)))
      {
        meanHolder[length(meanHolder)+1] = mean(group[ ,z])
      }
      groupComposite = rbind(groupComposite,unlist(meanHolder))
      meanHolder = list()
    }
    iso = unlist(isoHolder)
    groupComposite = data.frame(allGroups,iso,groupComposite,stringsAsFactors=F)
    #adjust col names
    names = colnames(peakData[c(-1,-2,-length(peakData))])
    names = unlist(lapply(names,"paste","mean",sep="_"))
    names = c("group","iso",names)
    colnames(groupComposite) = names
    #variable clean up
    rm(meanHolder)
    rm(isoHolder)
    rm(groupHolder)
    write.csv(groupComposite,file=file.path(resultsPath,paste(pheno,"GroupComposite_unfilled.csv",sep="_")),row.names=F)
#     groupcomp_out <<-groupComposite


    #Get indicies of groups which fail suspect criteria, merge so that validation is performed when either the 12C or 13C peak is suspect
    #may want to revert to old from (drop parens around group of or statments)
    indx_12C = which((groupComposite[ ,"iso"]=="12C")&((groupComposite [ ,"maxo_mean"]<intThresh) | (groupComposite[ ,"sn_mean"]<snThresh) | (groupComposite[,"dppm_mean"]>ppmThresh) | ((groupComposite[,"rtmax_mean"]-groupComposite[,"rtmin_mean"])<widthThresh)))
    indx_13C = which((groupComposite[,"iso"]==paste("13C_",n13C,sep=""))&((groupComposite[,"maxo_mean"]<intThresh) | (groupComposite[,"sn_mean"]<snThresh) | (groupComposite[,"dppm_mean"]>ppmThresh) | ((groupComposite[,"rtmax_mean"]-groupComposite[,"rtmin_mean"])<widthThresh)))
    indx = union(indx_12C,indx_13C)
    indx = indx[order(indx)] #sort index

    #get clusterNumber for each suspect group
    C12Clusters = pairedGroups[which(pairedGroups[ ,"12C"] %in% groupComposite[indx,"group"]),"clusterNumber"]
    C13Clusters = pairedGroups[which(pairedGroups[ ,paste("13C_",n13C,sep="")] %in% groupComposite[indx,"group"]),"clusterNumber"]
    clusters = unique(C12Clusters,C13Clusters)
    if(type=="valid")
    {
      clusters = pairedGroups[!(pairedGroups[ ,"clusterNumber"]%in%clusters),"clusterNumber"] #get all non-suspect groups for printing valid group clusters (fixed 6/13/2014)
	}
	return (list(clusters,indx))
}
