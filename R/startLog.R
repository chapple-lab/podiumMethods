startLog <-
function(resultsPath,Append=F)
{
	cat("Logging Console Output (but not warnings).\n")
	logFile = file.path(resultsPath,"RunLog.txt")
	if(file.exists(logFile)&Append)
	{ 
		sink(file=logFile,append=T,type="output",split=T)
		cat("\n\n-----------New Run-----------\n\n")
		print(Sys.time())
		cat("\n")
	}else{
	sink(file=file.path(resultsPath,"RunLog.txt"),append=F,type="output",split=T)
	cat("\n\n")
	print(Sys.time())
	cat("\n")
	}
	
}
