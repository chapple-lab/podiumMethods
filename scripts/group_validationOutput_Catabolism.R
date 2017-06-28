group_validationOutput_Catabolism = function(xcmsSet2=NULL,pheno=NULL,resultsPath=NULL,type=c("valid","suspect"),tTest=F,value=NULL,clusters=NULL,Indx=NULL,Rtindent=2,ppm=NULL,n13C=6)
#   tresults=NULL
{
	#output EIC, MS, and Table
    #Rtindent=2

	if(length(clusters)==0)
	{
		stop("\nNo ",type," clusters were found.\n\n")
	}


  filter = "HeuristicFilter"
  if(tTest)
    filter="tTestFilter"


  if(!is.null("ppm"))
  {
    path=file.path(resultsPath,paste(type,filter,pheno,value,"ppm",ppm,sep="_"))
  }
  else
  {
    path=file.path(resultsPath,paste(type,filter,pheno,value,sep="_"))
  }
  if(file.exists(path))
    unlink(file.path(path,"*")) #delete contents of existing folder in prep for writing new data
  else
    dir.create(path)

# 	pairedGroups = read.csv(file.path(resultsPath,paste(pheno,"Clusters_GroupNamesOnly.csv",sep="_")),check.names=F,stringsAsFactors=F)
  pairedGroups_matrix = read.csv(file.path(resultsPath,paste(pheno,"Clusters_MatrixForm.csv",sep="_")),check.names=F,stringsAsFactors=F)
	groupComposite = read.csv(file=file.path(resultsPath,paste(pheno,"GroupComposite_unfilled.csv",sep="_")),check.names=F,stringsAsFactors=F)


	#print EIC
	cat("\n\nGenerating EICs for",type,"Clusters in",pheno,"\n")
    group_printEICs(xcmsSet2,pairedGroups_matrix,groupComposite,clusters,path,rtrng=24,ppm=ppm,rtindent=Rtindent,output=F,n13C) #Default ppm=50,rtrng used to be 5

	#print MS
	cat("\n\nGenerating MS for",type,"Clusters in",pheno,"\n")
    group_printMS(xcmsSet2,pairedGroups_matrix,groupComposite,clusters,path,type="multiple",rtindent=Rtindent,n13C)

	#print tables according to analysis parameters
	cat("\n\nGenerating Tables for",type,"Clusters in",pheno,"\n")
	matx = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]%in%clusters),  ] #get group info for all clusters contained in 'clusters'
	groupNames = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]%in%clusters),"groupname"]
	gpIdx = match(groupNames,groupnames(xcmsSet2))
	matx = cbind(matx,groupval(xcmsSet2,method="medret",value=value)[gpIdx, ])

# 	if(type=="suspect")
# 	{
    gmatx = pairedGroups_matrix[Indx,  ]
    if(tTest==T)
		{
			write.csv(matx,file=file.path(resultsPath,paste(pheno,type,"Clusters_tTestFilter",paste(value,".csv",sep=""),sep="_")),row.names=F)
			write.csv(gmatx,file=file.path(resultsPath,paste(pheno,type,"Groups_tTestFilter",paste(value,".csv",sep=""),sep="_")),row.names=F)
		} else {
			write.csv(matx,file=file.path(resultsPath,paste(pheno,type,"Clusters_HeuristicFilter",paste(value,".csv",sep=""),sep="_")),row.names=F)
			write.csv(gmatx,file=file.path(resultsPath,paste(pheno,type,"Groups_HeuristicFilter",paste(value,".csv",sep=""),sep="_")),row.names=F)
		}
# 	}else if(type=="valid")
# 	{
# 		if(tTest==T)
# 		{
# 			write.csv(matx,file=file.path(resultsPath,paste(pheno,"ValidClusters_tTestFilter.csv",sep="_")),row.names=F)
# 			gmatx = pairedGroups_matrix[which(pairedGroups_matrix[ ,"groupname"]%in%tresults),  ]
# 			write.csv(gmatx,file=file.path(resultsPath,paste(pheno,"ValidGroups_tTestFilter.csv",sep="_")),row.names=F)
# 		} else {
# 			write.csv(matx,file=file.path(resultsPath,paste(pheno,"ValidClusters_HeuristicFilter.csv",sep="_")),row.names=F)
# 			gmatx = pairedGroups_matrix[Indx,  ]
# 			write.csv(gmatx,file=file.path(resultsPath,paste(pheno,"ValidGroups_HeuristicFilter.csv",sep="_")),row.names=F)
# 		}
# 	}
}
