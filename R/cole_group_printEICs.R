cole_group_printEICs <-
function(xcmsSet=NULL,pairedGroups_matrix=NULL,groupComposite=NULL,clusters=NULL,path,rtrng=100,ppm=50,rtindent=NULL,output=F) #groupComposite=Null
{

  ##generate indicies for groups with suspect peaks, then get EIC's for all of the peaks in each group
  clusters=sort(clusters)
  groups = pairedGroups_matrix[which(pairedGroups_matrix[ ,"clusterNumber"]%in%clusters),"groupname"] #names of all groups in suspect clusters
  EIC = coleGetEIC(xcmsSet,rtrange=rtrng,groupidx=groups,rt="corrected",mzExpMeth="ppm",ppm=ppm)
  if(output)
  {
    cat("\n\nEICs and Groups Output to Global Enviroment\n")
    EIC_out <<- EIC
    group_out <<- groups
  }

  #get names of 12C and 13C samples
  sampleNames = rownames(xcmsSet@phenoData)
  lastCol = ncol(xcmsSet@phenoData)
  C12sampleNames = sampleNames[grepl("12C",xcmsSet@phenoData[ ,lastCol])] #may need to change  "12C" to user defined identifier
  C13sampleNames = sampleNames[grepl("13C",xcmsSet@phenoData[ ,lastCol])] #may need to change  "13C" to user defined identifier
  #If unique folders were not used for the 12C and 13C samples, use the file name
  if(length(C12sampleNames)==0)
  {
    C12sampleNames = sampleNames[grepl("12C",sampnames(xcmsSet))] #may need to change  "12C" to user defined identifier
    C13sampleNames = sampleNames[grepl("13C",sampnames(xcmsSet))] #may need to change  "13C" to user defined identifier
  }

  peaks = peaks(xcmsSet)
  pkindx = groupval(xcmsSet)


  ##write EIC png's for each suspect cluster
  x=1
  for(j in clusters)
  {
    theseGroups = which(pairedGroups_matrix[ ,"clusterNumber"]==j)
    numGroups = length(theseGroups)
    names = pairedGroups_matrix[theseGroups,"groupname"]
#     groupPeaks_idx = which(groupComposite[ ,"group"]%in%names)
    group_idx = which(groupnames(xcmsSet)%in%names)
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


    legend("topright",c("C12","C13","C12_+6","C13_+6","C12_+More","C13_+More"),col = c("black","green","red","blue","purple","brown"),lty=1,ncol=3,xpd=T,seg.len=.5,cex=.75,x.intersp=.45,y.intersp=.75,inset=c(-0.03,-.06)) #lty=1,ncol=3,xpd=T,seg.len=.5,cex=.75,xjust=0.5,text.width=1.6,x.intersp=.75,y.intersp=.75,inset=c(0,-.02)
    dev.off()
    x=x+numGroups
	#cat("Cluster:",j,"\n")
    #cat("New x: ",x,"\n")
  }
  graphics.off()#just in case
}
