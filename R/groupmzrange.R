groupmzrange <-
function(groups,mzppm=10,mzabs=0.005,nC13=6) {
  C13mzdif<-1.0033548   ##mass of 13C
  mztor<-(groups$mzmed+C13mzdif*nC13)*mzppm*1e-6+mzabs  ##calculate M/z tolerance
  uppermz<-groups$mzmed+C13mzdif*nC13+mztor
  lowermz<-groups$mzmed+C13mzdif*nC13-mztor
  allpeaks<-cbind(groups,lowermz=lowermz,uppermz=uppermz) ##add calculated upper and lower ranges as columns to allpeaks dataframe
  #   write.csv(allpeaks,paste("allpeaksC",nC13,".csv",sep=""),row.names=F)
}
