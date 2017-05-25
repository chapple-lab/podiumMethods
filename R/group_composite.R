group_composite <-
function(xcmsSet2=NULL,pheno=NULL,resultsPath=NULL)
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
  C13Groups = which(groupNames %in% pairedGroups[ ,"13C_6"])
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
  iso = rep("13C_6",length(groupTags))
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
