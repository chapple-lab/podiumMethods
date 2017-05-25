coleGetSpec <-
function(object, sampleidx=sampnames(object), rt=c("raw","corrected"), type=c("multiple","single","average"), ...) {
	## FIXME: unnecessary dependency on profile matrix?
	##generate spectra for each xcmsRaw object specified by sampleidx
	#options(warn=-1)

	files <- filepaths(object)
    #grp <- groups(object)
    #samp <- sampnames(object)
    prof <- profinfo(object)

    rt <- match.arg(rt)

	if (is.numeric(sampleidx))
        sampleidx <- sampnames(object)[sampleidx]
    numsampleidx <- match(sampleidx, sampnames(object))

	specHolder = vector("list",length=length(sampleidx))
	count = 1
	for(samp in seq(along = numsampleidx))
	{
		#cat(sampleidx[samp], "\n")
		lcraw = xcmsRaw(files[numsampleidx[samp]], profmethod = prof$method, profstep = 0)
		if (rt == "corrected")
            lcraw@scantime <- object@rt$corrected[[numsampleidx[samp]]]

		if(length(object@dataCorrection) > 1)
		{
            if(object@dataCorrection[numsampidx[samp]] == 1)
                lcraw<-stitch(lcraw, AutoLockMass(lcraw))
		}

		sel <- profRange(lcraw, ...)

		scans <- list(length(sel$scanidx))
		uniquemz <- numeric()
		for (i in seq(along = sel$scanidx)) {
			scans[[i]] <- getScan(lcraw, sel$scanidx[i], sel$mzrange)
			uniquemz <- unique(c(uniquemz, scans[[i]][,"mz"]))#!!! round(scans[[i]][,"mz"],digits=2)
		}
		uniquemz <- sort(uniquemz)

		intmat <- matrix(nrow = length(uniquemz), ncol = length(sel$scanidx))
		for (i in seq(along = sel$scanidx)) {
			scan <- getScan(lcraw, sel$scanidx[i], sel$mzrange)
			intmat[,i] <- approx(scan, xout = uniquemz)$y
		}

		pts <- cbind(mz = uniquemz, intensity = rowMeans(intmat))
		#make sure all spectra are the same length
		specHolder[[count]] = pts
		count = count+1
		rm(lcraw)
		gc()
	}

	if(type=="multiple")
	{
		return(specHolder)
	}
	else if (type=="single")
	{
		return(specHolder[[1]])
	}
	else if(type=="average")
	{
		##average spectra to create composite spectra
		mzs = mat.or.vec(0,1)
		intens = mat.or.vec(0,1)
		#calculate values for padding specs
		initPts = sapply(specHolder,function(x) x[1 ,"mz"])
		rng = range(initPts)
		#cat("\n\nInit Pts: ",initPts,"\tRange: ",rng,"\n")
		minidx = which(initPts == rng[1])
		#minidx = minidx[1] #!!!!! for rounding
		padNum = findRange(specHolder[[minidx]][ ,"mz"],c(rng[2],rng[2]),NAOK=T) #determine maximum amount of padding needed (first num returned)
		#specHolder_out <<-specHolder

		padHolder=list()
		#pad front of specs (to help  mz vals line up as best as possible for subsequent averaging)
		for(x in 1:length(specHolder))
		{
			temp1 = specHolder[[x]][,"mz"]
			temp2 = specHolder[[x]][,"intensity"]
			pad = findRange(specHolder[[x]][ ,"mz"],c(rng[2],rng[2]),NAOK=T)
			times = abs(padNum[1] - pad[1]) ###WARNING CLUGE!!!!!!!
			#cat("\n\nX: ",x,"\ttimes: ",times,"\n")
			temp1 = c(rep(NA,times), temp1)
			temp2 = c(rep(NA,times), temp2)
			padHolder[[x]] = cbind(mz=temp1,intensity=temp2)
		}
		rm(specHolder)

		expand = max(sapply(padHolder,function(x) length(x)))/2 #determine max length, divide by 2 to get max col length

		#pad end of specs (to prevent recycling)
		for(x in 1:length(padHolder))
		{
			if(length(padHolder[[x]]) < (expand*2))
			{
				temp1 = padHolder[[x]][,"mz"]
				temp2 = padHolder[[x]][,"intensity"]
				length(temp1) = expand
				mzs = cbind(mzs,temp1)
				length(temp2) = expand
				intens = cbind(intens,temp2)
			}
			else
			{
				mzs = cbind(mzs, padHolder[[x]][ ,"mz"])
				intens = cbind(intens, padHolder[[x]][ ,"intensity"])
			}

		}
		mzs = rowMeans(mzs)#,na.rm=T)
		intens = rowMeans(intens)#,na.rm=T)
		matx = cbind(mz=mzs,intensity=intens)
		return(matx)
	}
	else
		stop("Invalid MS type requested.  Valid types are: 'multiple', 'average', or 'single'.")

}
