findC13cmpd <-
function(groupswithmzrange,rterror=1)
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
