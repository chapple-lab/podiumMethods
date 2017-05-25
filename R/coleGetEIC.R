coleGetEIC <-
function(object, mzrange, rtrange = 200,
                                        groupidx, sampleidx = sampnames(object),
                                        rt = c("corrected", "raw"), mzExpMeth=c("minMax","ppm"), ppm=50 )
  {

    files <- filepaths(object)
    grp <- groups(object)
    samp <- sampnames(object)
    prof <- profinfo(object)

    rt <- match.arg(rt)

    if (is.numeric(sampleidx))
        sampleidx <- sampnames(object)[sampleidx]
    numsampidx <- match(sampleidx, sampnames(object))


    if (!missing(groupidx)) {
        if (is.numeric(groupidx))
            groupidx <- groupnames(object)[unique(as.integer(groupidx))]
        grpidx <- match(groupidx, groupnames(object, template = groupidx))
    }

    if (missing(mzrange)) {
        if (missing(groupidx))
            stop("No m/z range or groups specified")
        if (any(is.na(groupval(object, value = "mz"))))
            stop('Please use fillPeaks() to fill up NA values !')
        if(mzExpMeth == "minMax") #picks mz limits based on smallest mzmin and largest mzmax across all samples
		{
			mzmin <- -rowMax(-groupval(object, value = "mzmin"))
			mzmax <- rowMax(groupval(object, value = "mzmax"))
			mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
		}
		else if(mzExpMeth == "ppm") #takes average of the mz values from each sample in the group and then does a symmetric ppm expansion around that value
		{
			cat("EIC_ppm= ",ppm,"\n")
			mz = rowMeans(groupval(object, value="mz"))
			mzmin = -(mz*ppm*1e-6)+mz
			mzmax = mz*ppm*1e-6+mz
			mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
		}

    } else if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]
    else if (is.null(dim(mzrange)))
        stop("mzrange must be a matrix")
    colnames(mzrange) <- c("mzmin", "mzmax")

    if (length(rtrange) == 1) {
        if (missing(groupidx))
            rtrange <- matrix(rep(range(object@rt[[rt]][numsampidx]), nrow(mzrange)),
                              ncol = 2, byrow = TRUE)
        else {
            rtrange <- retexp2(grp[grpidx,c("rtmin","rtmax"),drop=FALSE], rtrange)
        }
    } else if (is.null(dim(rtrange)))
        stop("rtrange must be a matrix or single number")
    colnames(rtrange) <- c("rtmin", "rtmax")

	###prevent Rt range from exceeding sample Rt bounds
	##Calculate bounds
	minRt = max(sapply(object@rt[[rt]],function(x) x[1])) #find largest minimum Rt across all samples
	maxRt = min(sapply(object@rt[[rt]],function(x) max(x))) #find smallest maximum Rt across all samples
	##Adjust rtRange matrix
	rtrange[ rtrange[ ,"rtmin"]<minRt, "rtmin"] = minRt
	rtrange[ rtrange[ ,"rtmax"]>maxRt, "rtmax"] = maxRt

    if (missing(groupidx))
        gnames <- character(0)
    else
        gnames <- groupidx

    eic <- vector("list", length(sampleidx))
    names(eic) <- sampleidx

    for (i in seq(along = numsampidx)) { ##changed from along=sampleidx

        cat(sampleidx[i], "\n")
        flush.console()
        lcraw <- xcmsRaw( files[numsampidx[i]], profmethod = prof$method, profstep = 0) ##Changed from files[sampidx[i]]
		#cat("\nProcessing Sample:",lcraw@filepath[1],"\n")
        if(length(object@dataCorrection) > 1){
            if(object@dataCorrection[numsampidx[i]] == 1)
                lcraw<-stitch(lcraw, AutoLockMass(lcraw))
        }
        if (rt == "corrected")
            lcraw@scantime <- object@rt$corrected[[numsampidx[i]]] #changed from object@rt$corrected[[sampidx[i]]]
        if (length(prof) > 2)
            lcraw@profparam <- prof[seq(3, length(prof))]
        currenteic <- getEIC(lcraw, mzrange, rtrange, step = prof$step)
        eic[[i]] <- currenteic@eic[[1]]
        rm(lcraw)
        gc()
    }
    cat("\n")

    invisible(new("xcmsEIC", eic = eic, mzrange = mzrange, rtrange = rtrange,
                  rt = rt, groupnames = gnames))
}
