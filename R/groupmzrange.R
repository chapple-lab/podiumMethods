##DEPRECATED::calculates range using ppm of 12C peak, not 13C peak, may want to change
##changed to calculate based on the 13C peak
##For ring-labeled phe use mLabel = 1.0033548 and nLabel = 6
groupmzrange <-
function(groups,mzppm=10,mzabs=0,nLabel=6,mLabel=1.0033548) {
  ##mass of 13C C13mzdif<-1.0033548
  ##old mzabs default was 0.005
  mztor<-(groups$mzmed+mLabel*nLabel)*mzppm*1e-6+mzabs  ##calculate M/z tolerance
  uppermz<-groups$mzmed+mLabel*nLabel+mztor
  lowermz<-groups$mzmed+mLabel*nLabel-mztor
  allpeaks<-cbind(groups,lowermz=lowermz,uppermz=uppermz) ##add calculated upper and lower ranges as columns to allpeaks dataframe
  #   write.csv(allpeaks,paste("allpeaksC",nLabel,".csv",sep=""),row.names=F)
}
