group_printMS <-
function(xcmsSet2=NULL,pairedGroups_matrix=NULL,clusters=NULL,path,type=c("single","average","multiple"),mzexp=10,perstdev=1.25,rtindent=1.25)
{
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
  C12sampleidx = which(sampnames(xcmsSet2) %in% C12sampleNames)
	C13sampleidx = which(sampnames(xcmsSet2) %in% C13sampleNames)

    #get m/z and rt values for each group
    ##switched to using my groupComposite instead of the groups() values stored in paired groups matrix as they do not contain
    ##the info I need
    #gpvals = groups(xcmsSet2)
    # gpvals = groupComposite #read.csv(file=file.path(resultsPath,paste(pheno,"groupComposite_unfilled.csv",sep="_")),check.names=F,stringsAsFactors=F)
    rtminMatx = groupval(xcmsSet2,value="rtmin")
    rtmaxMatx = groupval(xcmsSet2,value="rtmax")

    for(j in clusters)
    {
      #determine 12C group
      theseGroups = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]==j), "groupname"]
      basegp = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]==j & pairedGroups_matrix[ ,"iso"]=="12C"), ]
      baseidx = which(groupnames(xcmsSet2)==basegp$groupname)
      # basegp = gpvals[which(gpvals[ ,"group"]==basegp$groupname), ]
      C13gp = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]==j & pairedGroups_matrix[ ,"iso"]=="13C_6"), ]
      # C13gp = gpvals[which(gpvals[ ,"group"]==C13gp$groupname), ]

      #determine mz range
      mzmin = basegp["mzmed"]-mzexp
      mzmax = basegp["mzmed"]+12+mzexp


      #determine rt bounds as rtmin_med and rtmax_med across one representative peak from each group
      #determine MS slice bounds as +/- perstdev*stdev arround median rt
      #note: for a normal distribution, one standard deviation can be approximated as 1/6 of the range of x values
      rtmin_med = median(rtminMatx[baseidx,])
      rtmax_med = median(rtmaxMatx[baseidx,])
      stdev = (rtmax_med-rtmin_med)/6
      msRtmin = basegp$rtmed - stdev*perstdev
      msRtmax = basegp$rtmed + stdev*perstdev

      #check to ensure that rt window is > 1 sec
      if(msRtmax-msRtmin<1)
      {
        msRtmax = msRtmax+0.5
        msRtmin = msRtmin+0.5
      }

#       #determine rt range adjusting for rtindent if possible, ensure minimum rt window of 1sec
#       if((basegp$rtmax-basegp$rtmin)<=(rtindent*2+1)| ((basegp$rtmax-basegp$rtmin) <= 2 ))#check to make sure rtindent will not make rtmax smaller than rtmin
#       {
#         if((basegp$rtmax-basegp$rtmin) < 1 ) #ensure 1 sec rt window
#         {
#           rtmin = basegp$rtmin-0.5
#           rtmax = basegp$rtmax+0.5
#         }
#         else
#         {
#           rtmin = basegp$rtmin
#           rtmax = basegp$rtmax
#         }
#
#       }
#       else
#       {
#         rtmin = basegp$rtmin+rtindent
#         rtmax = basegp$rtmax-rtindent
#       }

      # cat("\n\nCluster: ",j,"\nMz range: ",as.vector(c(mzmin, mzmax),mode="numeric"),"\nRt range:",as.vector(c(rtmin, rtmax),mode="numeric"),"\n")

      #get spectra
      if(type=="average") #print MS as an averaged composite of all Samples,  probably want to reprogram to use a spline for interpolation instead of averaging.
      {
        C12ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C12sampleidx,rtrange=as.vector(c(msRtmin, msRtmax),mode="numeric"),rt="corrected")
        C12ms[which(is.na(C12ms[,2])),2]=0 #change missing data to zeros
        C13ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C13sampleidx,rtrange=as.vector(c(msRtmin, msRtmax),mode="numeric"),rt="corrected")
        C13ms[which(is.na(C13ms[,2])),2]=0 #change missing data to zeros
      }
      else if(type=="single") ##Probably need to make sure Samples line up, print MS using info from only one sample
      {
        C12ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C12sampleidx[1],rtrange=as.vector(c(msRtmin, msRtmax),mode="numeric"),rt="corrected")
        C12ms[which(is.na(C12ms[,2])),2]=0 #change missing data to zeros
        C13ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C13sampleidx[1],rtrange=as.vector(c(msRtmin, msRtmax),mode="numeric"),rt="corrected")
        C13ms[which(is.na(C13ms[,2])),2]=0 #change missing data to zeros
      }
      else if(type=="multiple")
      {
        C12ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C12sampleidx,rtrange=as.vector(c(msRtmin, msRtmax),mode="numeric"),rt="corrected")
        #change missing data to zeros
        for(x in 1:length(C12ms))
        {
          C12ms[[x]][which(is.na(C12ms[[x]][,2])),2]=0
        }

        C13ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C13sampleidx,rtrange=as.vector(c(msRtmin, msRtmax),mode="numeric"),rt="corrected")
        #change missing data to zeros
        for(x in 1:length(C13ms))
        {
          C13ms[[x]][which(is.na(C13ms[[x]][,2])),2]=0
        }

      }


      #print to file
      name = paste("Cluster",j,"Groups",do.call(paste,as.list(theseGroups)),sep="_") #generate descriptive name for group cluster
      png(file.path(path,paste(name,"MS.png",sep="_")),width=1200,height=800,res=170) #generate file name, path is global variable set with setPath()
      #calc nums for scaling x and y axis, scales y axis based on intensities found in the middle 2 quartiles of the MS spectra
      if(type=="single"|type=="average")
      {
        C12rng = length(C12ms[ ,1])/4
        C13rng = length(C13ms[ ,1])/4
        xlim = range(C12ms[,1],C13ms[,1],na.rm=T)
        ylim = c(0,max(C12ms[c(C12rng:(length(C12ms[,1])-C12rng)),2],C13ms[c(C13rng:(length(C13ms[,1])-C13rng)),2])+20)
      }

      if(type=="multiple")
      {

        C12rng = 0
        C13rng = 0
        xlim = 0
        ylim = 0

        for(x in 1:min(length(C12ms),length(C13ms)))
        {
          C12rng1 = length(C12ms[[x]][ ,1])/4
          C13rng1 = length(C13ms[[x]][ ,1])/4
          xlim1 = range(C12ms[[x]][,1],C13ms[[x]][,1],na.rm=T)
          ylim1 = c(0,max(C12ms[[x]][c(C12rng1:(length(C12ms[[x]][,1])-C12rng1)),2],C13ms[[x]][c(C13rng1:(length(C13ms[[x]][,1])-C13rng1)),2])+20)
          if(x==1)
          {
            xlim = xlim1
            ylim = ylim1
          }
          else
          {
            xlim=range(xlim,xlim1,na.rm=T)
            ylim=c(0,max(ylim1,ylim))
          }
        }
      }
      # cat("Cluster: ",j,"\txlim: ",xlim,"\n")
      plot(0,0,xlim=xlim,ylim=ylim,main=paste(name,"Rt:",round(msRtmin,digits=2),"-",round(msRtmax,digits=2),"Validation Plot",sep=" "),xlab="M/z",ylab="Intensity")
      #plot vertical lines for C12 and C13 compound base peaks
      C12_Mz = as.vector(basegp$mzmed,mode="numeric")
      C13_Mz = as.vector(C13gp$mzmed,mode="numeric")
      #plot MS
      if(type=="single"|type=="average")
      {
        lines(C13ms,col="red") #order changed, used to be black first
        lines(C12ms,col="black")
      }
      else if(type=="multiple")
      {
        for(x in 1:min(length(C12ms),length(C13ms)))
        {
          lines(C13ms[[x]],col="red")
          lines(C12ms[[x]],col="black") #order changed, used to be red first
        }
      }

#       print(C12_Mz)
#       print(C13_Mz)
      abline(v=C12_Mz,col="black",lty=2,lwd=1)
      abline(v=C13_Mz,col="red",lty=2,lwd=1) #peaks[which((peaks[,"group"]==j)&(peaks[,"iso"]=="13C_6")),"mz"]
      text(xlim[2]-.9,ylim[2]-.015*ylim[2],paste("Shift:",round((C13_Mz-C12_Mz),digits=4)),cex=.75)
      legend("topright",c("C12","C13"),col = c("black","red"),lty=1,horiz=T,xpd=T,seg.len=.5,cex=.75,x.intersp=.45,y.intersp=.75,inset=c(-0.03,-.06))
      dev.off()

    }

}
