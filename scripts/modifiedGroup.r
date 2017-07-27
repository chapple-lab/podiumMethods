#setMethod("group.density", "xcmsSet", function(object, bw = 30, minfrac = 0.5, minsamp = 1,
#                                              mzwid = 0.25, max = 50, sleep = 0) {

findEqualGreaterM = coleXcms:::findEqualGreaterM
descendMin = coleXcms:::descendMin
group.density = function(object, bw = 30, minfrac = 0.5, minsamp = 1,mzwid = 0.25, max = 50, sleep = 0) {
    samples <- sampnames(object)
    classlabel <- sampclass(object)
    classnames <- as.character(unique(sampclass(object)))
    classlabel <- as.vector(unclass(classlabel))
    classnum <- table(classlabel)

    peakmat <- peaks(object)
    porder <- order(peakmat[,"mz"])
    peakmat <- peakmat[porder,, drop=FALSE]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])

    mass <- seq(peakmat[1,"mz"], peakmat[nrow(peakmat),"mz"] + mzwid, by = mzwid/2)
    masspos <- findEqualGreaterM(peakmat[,"mz"], mass)

    groupmat <- matrix(nrow = 512, ncol = 7 + length(classnum))
    groupindex <- vector("list", 512)

    endidx <- 0
    num <- 0
    gcount <- integer(length(classnum))
    for (i in seq(length = length(mass)-2)) {
        if (i %% 500 == 0) {
            cat(round(mass[i]), "")
            flush.console()
        }
        startidx <- masspos[i]
        endidx <- masspos[i+2]-1
        if (endidx - startidx < 0)
            next
        speakmat <- peakmat[startidx:endidx,,drop=FALSE]
        den <- density(speakmat[,"rt"], bw, from = retrange[1]-3*bw, to = retrange[2]+3*bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2+length(classnum))
        snum <- 0
        while (deny[maxy <- which.max(deny)] > maxden/20 && snum < max) {
            grange <- descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(speakmat[,"rt"] >= den$x[grange[1]] & speakmat[,"rt"] <= den$x[grange[2]])
            gnum <- classlabel[unique(speakmat[gidx,"sample"])]
            for (j in seq(along = gcount))
                gcount[j] <- sum(gnum == j)
            
			mzmed = median(speakmat[gidx, "mz"])
			if(mzmed>=345.0889&mzmed<=345.0958 )
			{
				cat("\n\nGcount:",gcount,"\nSamps:",as.character(gnum))
			}
			if (! any(gcount >= classnum*minfrac & gcount >= minsamp))
                next
            snum <- snum + 1
            num <- num + 1
### Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            groupmat[num, 1] <- median(speakmat[gidx, "mz"])
            groupmat[num, 2:3] <- range(speakmat[gidx, "mz"])
            groupmat[num, 4] <- median(speakmat[gidx, "rt"])
            groupmat[num, 5:6] <- range(speakmat[gidx, "rt"]) #note: this gives rtmin and rtmax in terms of the central rt for each peak, NOT in terms of the rtmin and rtmax for each peak
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7+seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(porder[(startidx:endidx)[gidx]])
        }
        if (sleep > 0) {
            plot(den, main = paste(round(min(speakmat[,"mz"]), 2), "-", round(max(speakmat[,"mz"]), 2)))
            for (i in seq(along = classnum)) {
                idx <- classlabel[speakmat[,"sample"]] == i
                points(speakmat[idx,"rt"], speakmat[idx,"into"]/max(speakmat[,"into"])*maxden, col = i, pch=20)
            }
            for (i in seq(length = snum))
                abline(v = groupmat[num-snum+i, 5:6], lty = "dashed", col = i)
            Sys.sleep(sleep)
        }
    }
    cat("\n")

    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", classnames)

    groupmat <- groupmat[seq(length = num),,drop=FALSE]
    groupindex <- groupindex[seq(length = num)]

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])

    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)

    groups(object) <- groupmat[uindex,,drop=FALSE]
    groupidx(object) <- groupindex[uindex]

    object
}
#)