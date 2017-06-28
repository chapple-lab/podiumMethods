source("x:/Labs/Chapple/0-Cole Wunderlich/R files/Code/Methods/group_validationPrep_Catabolism.r")
source("x:/Labs/Chapple/0-Cole Wunderlich/R files/Code/Methods/group_validationOutput_Catabolism.r")
source("x:/Labs/Chapple/0-Cole Wunderlich/R files/Code/Methods/removeGroups.r")
source("x:/Labs/Chapple/0-Cole Wunderlich/R files/Code/Methods/hFilter.r")
source("x:/Labs/Chapple/0-Cole Wunderlich/R files/Code/Methods/tFilter_V3.r")
# DEPRICATED source("x:/Labs/Chapple/0-Cole Wunderlich/R files/Code/Methods/tFilter.r")

group_printEICs = function(xcmsSet2=NULL,pairedGroups_matrix=NULL,groupComposite=NULL,clusters=NULL,path,rtrng=100,ppm=50,rtindent=NULL,output=F,n13C=6) #groupComposite=Null
{

  ##generate indicies for groups with suspect peaks, then get EIC's for all of the peaks in each group
  clusters=sort(clusters)
  groups = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]%in%clusters),"groupname"] #names of all groups in suspect clusters
  EIC = getEIC(xcmsSet2,rtrange=rtrng,groupidx=groups,rt="corrected",mzExpMeth="ppm",ppm=ppm)
  if(output)
  {
    cat("\n\nEICs and Groups Output to Global Enviroment\n")
    EIC_out <<- EIC
    group_out <<- groups
  }

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

  peaks = peaks(xcmsSet2)
  pkindx = groupval(xcmsSet2)


  ##write EIC png's for each suspect cluster
  x=1
  for(j in clusters)
  {
    theseGroups = which(pairedGroups_matrix[ ,"clusterNumber"]==j)
    numGroups = length(theseGroups)
    names = pairedGroups_matrix[theseGroups,"groupname"]
#     groupPeaks_idx = which(groupComposite[ ,"group"]%in%names)
    group_idx = which(groupnames(xcmsSet2)%in%names)
    #cat("Peak_indxs: ",thesePeaks,"\nNumPeaks: ",numPeaks,"\n\n")
    name = paste("Cluster",j,"Groups",do.call(paste,as.list(pairedGroups_matrix[theseGroups,"groupname"])),sep="_") #generate descriptive name for group cluster
    png(file.path(path,paste(name,"EIC.png",sep="_")),width=1200,height=800,res=170) #generate file name, path is global variable set with setPath()

    #setup plot window
    #par(oma=c(0,0,0,2.7))
    xlim= range(EIC@rtrange[x:(x+numGroups-1),])
    #         cat("\n\nXlim: ",xlim,"\nEIC Range")#\n12C: ",range(C12eic@eic$xcmsRaw[[x]][,1]),"\n13C: ",range(C13eic@eic$xcmsRaw[[x]][,1]))
    #         for(r in (x:(x+numPeaks)))
    #         {
    #           cat("\n12C:",range(C12eic@eic$xcmsRaw[[r]][,1]),"\n13C: ",range(C13eic@eic$xcmsRaw[[r]][,1]))
    #         }
    #expand x axis
    shift=2
    xlim[1]=xlim[1]-shift
    xlim[2]=xlim[2]+shift
    #set Y axis
#     ylim=range(groupComposite[which(groupComposite[ ,"group"]%in%names),"maxo_mean"])#,C12eic@eic$xcmsRaw[[x]][,2]) #may need to adjust
#
#     highGroup = max(groupComposite[groupPeaks_idx,"maxo_mean"]) #find max intensity in group cluster
#     highGroup = groupComposite[which(groupComposite[groupPeaks_idx,"maxo_mean"]==highGroup), "group"] #get name of group with max intensity
#     highGroup_indx = which(pairedGroups_matrix[theseGroups,"groupName"]==)
    intenRange = list()
    for(sname in C12sampleNames)
    {
      intenRange=append(intenRange,range(EIC@eic[[sname]][[x]][ ,2]))#C12_C12
      intenRange=append(intenRange,range(EIC@eic[[sname]][[x+1]][ ,2]))#C12_6
    }
    for(sname in C13sampleNames)
    {
      intenRange=append(intenRange,range(EIC@eic[[sname]][[x]][ ,2]))#C12_C12
      intenRange=append(intenRange,range(EIC@eic[[sname]][[x+1]][ ,2]))#C12_6
    }
    ylim = range(unlist(intenRange))
#     range_out <<-intenRange
#     print(intenRange)
	  ylim[1]=0
    ylim[2]=ylim[2]+25 #may need to adjust
# 	  print(ylim)
    plot(0,0,xlim=xlim,ylim=ylim,main=paste(name,"Validation Plot",sep=" "),xlab="Rt (sec)",ylab="Intensity")

    #plot rt lines for generation of MS if rtindent is specified, otherwise, plot a reference line
    #need to decide where to take MS slice for group validation, this will affect where bounds are set
    #FIX CODE BELOW
    if(F)
    {


    if(!is.null(rtindent))
    {
      #based on C12 peak
      if((peaks[thesePeaks[1],"rtmax"]-peaks[thesePeaks[1],"rtmin"]) <= (rtindent*2+1)) #check to make sure rtindent will not make rtmax smaller than rtmin
      {
        rtmin = peaks[thesePeaks[1],"rtmin"]
        rtmax = peaks[thesePeaks[1],"rtmax"]
      }
      else
      {
        rtmin = peaks[thesePeaks[1],"rtmin"]+rtindent
        rtmax = peaks[thesePeaks[1],"rtmax"]-rtindent
      }

      abline(v=rtmin,col="gray",lty=2)
      abline(v=rtmax,col="gray",lty=2)
    }
    else
    {
      abline(v=peaks[thesePeaks[1],"rt"],col="gray",lty=2)
    }

    }
    #END FIX SECTION

    #NEED TO CHECK TO MAKE SURE X IS WORKING (referencing the correct EICs for each cluster)
    #plot 12C_C12, 12C_6, and 12C_more
    for(sname in C12sampleNames)
    {
      peakrngs = peaks[pkindx[group_idx,sname], c("rtmin","rtmax")]#get Rt range over which actual peaks occur
#       cat("\n\ncluster ",j,"\tSampName: ",sname,"\n")
#       print(peakrngs)
      pts = EIC@eic[[sname]][[x]]
      lines(pts,col="gray")#C12_C12 entire EIC

      int_idx = (pts[ ,"rt"] >= peakrngs[1,1] & pts[ ,"rt"] <= peakrngs[1,2])
      points(pts[int_idx, ],type="l", col="black")#C12_C12 peak

      pts = EIC@eic[[sname]][[x+1]]
      lines(pts,col=lightColorList$lightRed)#C12_6 entire EIC
      int_idx = (pts[ ,"rt"] >= peakrngs[2,1] & pts[ ,"rt"] <= peakrngs[2,2])
      points(pts[int_idx, ],type="l", col="red")#C12_6 peak

      if(numGroups >= 3)#if extra peaks (beyond +6)
      {
        cntr = 3
        for(y in ((x+2):(x+numGroups-1)))#handle rest of peaks
        {
          pts = EIC@eic[[sname]][[y]]
          lines(pts,col=lightColorList$lightPurp)
          int_idx = (pts[ ,"rt"] >= peakrngs[cntr,1] & pts[ ,"rt"] <= peakrngs[cntr,2])
          points(pts[int_idx, ],type="l", col="purple")
          cntr = cntr+1
        }
      }
    }

    #plot 13C_C12 and 13C_6
    for(sname in C13sampleNames)
    {

      peakrngs = peaks[pkindx[group_idx, sname], c("rtmin","rtmax")]#get Rt range over which actual peaks occur
#       cat("\n\ncluster ",j,"\tSampName: ",sname,"\n")
#       print(peakrngs)
      pts = EIC@eic[[sname]][[x]]
      lines(pts,col=lightColorList$lightGrn)#C13_C12 entire EIC

      int_idx = (pts[ ,"rt"] >= peakrngs[1,1] & pts[ ,"rt"] <= peakrngs[1,2])
      points(pts[int_idx, ],type="l", col="green")#C13_C12 peak

      pts = EIC@eic[[sname]][[x+1]]
      lines(pts,col=lightColorList$lightBlue)#C13_6 entire EIC

      int_idx = (pts[ ,"rt"] >= peakrngs[1,1] & pts[ ,"rt"] <= peakrngs[1,2])
      points(pts[int_idx, ],type="l", col="blue")#C13_6 peak

      if(numGroups >=3)#if extra peaks (beyond +6)
      {
        cntr = 3
        for(y in ((x+2):(x+numGroups-1)))#handle rest of peaks
        {
          pts = EIC@eic[[sname]][[y]]
          lines(pts,col=lightColorList$lightBrwn)
          int_idx = (pts[ ,"rt"] >= peakrngs[cntr,1] & pts[ ,"rt"] <= peakrngs[cntr,2])
          points(pts[int_idx, ],type="l", col="brown")
          cntr = cntr+1
        }
      }
    }


    legend("topright",c("C12","C13",paste("C12_+",n13C,sep=""),paste("C13_+",n13C,sep=""),"C12_+More","C13_+More"),col = c("black","green","red","blue","purple","brown"),lty=1,ncol=3,xpd=T,seg.len=.5,cex=.75,x.intersp=.45,y.intersp=.75,inset=c(-0.03,-.06)) #lty=1,ncol=3,xpd=T,seg.len=.5,cex=.75,xjust=0.5,text.width=1.6,x.intersp=.75,y.intersp=.75,inset=c(0,-.02)
    dev.off()
    x=x+numGroups
	#cat("Cluster:",j,"\n")
    #cat("New x: ",x,"\n")
  }
  graphics.off()#just in case
}


group_printMS = function(xcmsSet2=NULL,pairedGroups_matrix=NULL,groupComposite=NULL,clusters=NULL,path,type=c("single","average","multiple"),mzexp=10,rtindent=1.25)
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
    gpvals = groupComposite #read.csv(file=file.path(resultsPath,paste(pheno,"groupComposite_unfilled.csv",sep="_")),check.names=F,stringsAsFactors=F)


    for(j in clusters)
    {
      #determine 12C group
      theseGroups = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]==j), "groupname"]
      basegp = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]==j & pairedGroups_matrix[ ,"iso"]=="12C"), ]
      basegp = gpvals[which(gpvals[ ,"group"]==basegp$groupname), ]
      C13gp = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]==j & pairedGroups_matrix[ ,"iso"]==paste("13C_",n13C,sep="")), ]
      C13gp = gpvals[which(gpvals[ ,"group"]==C13gp$groupname), ]


      #determine mz range
      mzmin = basegp["mz_mean"]-mzexp
      mzmax = basegp["mz_mean"]+12+mzexp


      #determine rt range adjusting for rtindent if possible
      if(((basegp$rtmax_mean-basegp$rtmin_mean)<=(rtindent*2+1)) | ((basegp$rtmax_mean-basegp$rtmin_mean) <= 2 ))#check to make sure rtindent will not make rtmax smaller than rtmin
      {
        rtmin = basegp$rtmin_mean
        rtmax = basegp$rtmax_mean
      }
      else
      {
        rtmin = basegp$rtmin_mean+rtindent
        rtmax = basegp$rtmax_mean-rtindent
      }

#       cat("\n\nCluster: ",j,"\nMz range: ",as.vector(c(mzmin, mzmax),mode="numeric"),"\n")

      #get spectra
      if(type=="average") #print MS as an averaged composite of all Samples,  probably want to reprogram to use a spline for interpolation instead of averaging.
      {
        C12ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C12sampleidx,rtrange=as.vector(c(rtmin, rtmax),mode="numeric"),rt="corrected")
        C12ms[which(is.na(C12ms[,2])),2]=0 #change missing data to zeros
        C13ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C13sampleidx,rtrange=as.vector(c(rtmin, rtmax),mode="numeric"),rt="corrected")
        C13ms[which(is.na(C13ms[,2])),2]=0 #change missing data to zeros
      }
      else if(type=="single") ##Probably need to make sure Samples line up, print MS using info from only one sample
      {
        C12ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C12sampleidx[1],rtrange=as.vector(c(rtmin, rtmax),mode="numeric"),rt="corrected")
        C12ms[which(is.na(C12ms[,2])),2]=0 #change missing data to zeros
        C13ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C13sampleidx[1],rtrange=as.vector(c(rtmin, rtmax),mode="numeric"),rt="corrected")
        C13ms[which(is.na(C13ms[,2])),2]=0 #change missing data to zeros
      }
      else if(type=="multiple")
      {
        C12ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C12sampleidx,rtrange=as.vector(c(rtmin, rtmax),mode="numeric"),rt="corrected")
        #change missing data to zeros
        for(x in 1:length(C12ms))
        {
          C12ms[[x]][which(is.na(C12ms[[x]][,2])),2]=0
        }

        C13ms= getSpec(xcmsSet2,type=type,mzrange=as.vector(c(mzmin, mzmax),mode="numeric"),rawidx=C13sampleidx,rtrange=as.vector(c(rtmin, rtmax),mode="numeric"),rt="corrected")
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

      plot(0,0,xlim=xlim,ylim=ylim,main=paste(name,"Rt:",round(rtmin,digits=2),"-",round(rtmax,digits=2),"Validation Plot",sep=" "),xlab="M/z",ylab="Intensity")
      #plot vertical lines for C12 and C13 compound base peaks
      C12_Mz = as.vector(basegp$mz_mean,mode="numeric")
      C13_Mz = as.vector(C13gp$mz_mean,mode="numeric")
      #plot MS
      if(type=="single"|type=="average")
      {
        lines(C12ms,col="black")
        lines(C13ms,col="red")
      }
      else if(type=="multiple")
      {
        for(x in 1:min(length(C12ms),length(C13ms)))
        {
          lines(C12ms[[x]],col="black")
          lines(C13ms[[x]],col="red")
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
