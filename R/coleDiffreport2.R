coleDiffreport2 <-
function(object,resultsPath=NULL,type,tTest,pheno,
											 class1=NULL, class2=NULL, filebase = character(),
											 eicmax = 0, eicwidth = 200, sortpval = TRUE, classeic =NULL,
                                             value = c("into","maxo","intb"), dataVal = c("into","maxo","intb"),
											 metlin = FALSE, h = 480, w = 640, mzdec=2,ppm=50, ...)

{

	##Prep for calling diffreport
	classes = sampclass(object)
	classes = levels(classes)#[as.vector(unclass(classes))]

	#if not specified, set parameters to default values
	if(length(class1)==0)
		class1=classes[1]

	if(length(class2)==0)
		class2=classes[2]

	if(length(classeic)==0)
		classeic = c(class1,class2)

	#read in valid Labelled Group Data
	if(type=="suspect")
	{
		stop("Type must be \"valid\", make sure validation of type \"valid\" has already been run")
	}else if(type=="valid")
	{
		if(tTest==T)
		{
			labelData = read.csv(file=file.path(resultsPath,paste(pheno,type,"Clusters_tTestFilter",paste(dataVal,".csv",sep=""),sep="_")),check.names=F,stringsAsFactors=F)
			fname = paste(pheno,value,class1,"vs",class2,"tTestFilter_FinalReport.csv",sep="_")

		} else {
			labelData = read.csv(file=file.path(resultsPath,paste(pheno,type,"Clusters_HeuristicFilter",paste(dataVal,".csv",sep=""),sep="_")),check.names=F,stringsAsFactors=F)
			fname = paste(pheno,value,class1,"vs",class2,"HeuristicFilter_FinalReport.csv",sep="_")
		}
	}else
	{
		warning("invalid type argument")
	}
	#Remove all groups and related information from xcmsSet that are not labelled Groups
	labeledGpIdx = match(labelData[ ,"groupname"], groupnames(object))
	object@groups = object@groups[labeledGpIdx, ]
	object@groupidx = object@groupidx[labeledGpIdx]
	#if(modPhenoData&&(ncol(object@phenoData)>1))
	#{
		#Remove all sample class info except highest level phenotype
	#	object@phenoData[ ,-1] = NULL
	#}


	cat("\n\nGenerating Report for:",class1,"and",class2,".\n")
	report = coleDiffreport(object,groupnames = labelData[,"groupname"],class1=class1,class2=class2,filebase = filebase, eicmax = eicmax, eicwidth = eicwidth,
                                             sortpval = sortpval, classeic = classeic,value = value, metlin = FALSE,
                                             h = 480, w = 640, mzdec=2,ppm=50,...)

	##Add labeling data to diffreport
	reportIdx = match(report[ ,"name"],labelData[ ,"groupname"])
	finalResults = cbind(labelData[reportIdx,1:2],report)

	write.csv(finalResults,file=file.path(resultsPath,fname),row.names=F)
	return (finalResults)

}
