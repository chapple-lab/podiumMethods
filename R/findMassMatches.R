
#Note, conversions and massList must be data frames
findMassMatches = function(rawList,refList,ppm=20,rtWind=0,rawmz="[M-H]",rawrt="Rt",refmz="[M-H]",refrt="Rt",short=T,shortList=c(rawmz,rawrt,"default"),shortList2=shortList)
{

  results=NULL
  pb = txtProgressBar(min=0,max=nrow(rawList),style=3,title="Progress")

  if("default"%in%shortList&rtWind==0)
  {
    shortList=shortList[-which(shortList%in%c("default",rawrt))]
  }
  else if("default"%in%shortList)
    shortList=shortList[!(shortList=="default")]


  #test shortList
  if(short&length(shortList)!=length(which(shortList%in%colnames(rawList))))
  {
    stop("Error: Requested return columns not found in input data.")
  }

  for(x in 1:nrow(rawList))
  {
    metab = rawList[x,]
    if(rtWind!=0) # if window is zero, disregard RT matching
    {
      indx = which((refList[,refrt] <= metab[,rawrt]+rtWind)&(refList[,refrt] >= metab[,rawrt]-rtWind))
    }
    else
    {
      indx = 1:nrow(refList)
    }
    massUpper = metab[,rawmz]+metab[,rawmz]*ppm/10^6
    massLower = metab[,rawmz]-metab[,rawmz]*ppm/10^6
    idx = which((refList[indx,refmz] <=massUpper) & (refList[indx,refmz] >= massLower)) #check for Mz match
    if(length(idx)!=0)
    {
      ppm2 = round((metab[,rawmz]-refList[indx[idx],refmz])/refList[indx[idx],refmz]*10^6, digits=2)
      if(short)
        if(length(shortList)==1)
          result=cbind(Peak_Mz=metab[,shortList],match="=>",ppm=ppm2,refList[indx[idx],],row.names=NULL)
          else
            result=cbind(metab[,shortList],match="=>",ppm=ppm2,refList[indx[idx],shortList2],row.names=NULL)
      else
        result=cbind(metab,match="=>",refList[indx[idx],],row.names=NULL)

      results = rbind(results,result)

    }
    setTxtProgressBar(pb,x)
  }

  invisible(results)
}