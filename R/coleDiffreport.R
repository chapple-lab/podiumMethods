coleDiffreport <-
function(object, groupnames=NULL, class1 = levels(sampclass(object))[1],
                                            class2 = levels(sampclass(object))[2],
                                            filebase = character(), eicmax = 0, eicwidth = 200,
                                            sortpval = TRUE, classeic = c(class1,class2),
                                            value = c("into","maxo","intb"), metlin = FALSE,
                                            h = 480, w = 640, mzdec=2, ppm=50, ...) {

    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

    require(multtest) || stop("Couldn't load multtest")

    value <- match.arg(value)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    samples <- sampnames(object)
    n <- length(samples)
    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]

    values <- groupval(object, "medret", value=value)
    indecies <- groupval(object, "medret", value = "index")

    if (!all(c(class1,class2) %in% classlabel))
        stop("Incorrect Class Labels")

    ## c1 and c2 are column indices of class1 and class2 resp.
    c1 <- which(classlabel %in% class1)
    c2 <- which(classlabel %in% class2)
    ceic <- which(classlabel %in% classeic)
    if (length(intersect(c1, c2)) > 0)
        stop("Intersecting Classes")

    ## Check against missing Values
    if (any(is.na(values[,c(c1,c2)]))) {
        stop("NA values in xcmsSet. Use fillPeaks()")
    }

    mean1 <- rowMeans(values[,c1,drop=FALSE], na.rm = TRUE)
    mean2 <- rowMeans(values[,c2,drop=FALSE], na.rm = TRUE)

    ## Calculate fold change.
    ## For foldchange <1 set fold to 1/fold
    ## See tstat to check which was higher
    fold <- mean2 / mean1
    fold[!is.na(fold) & fold < 1] <- 1/fold[!is.na(fold) & fold < 1]

    testval <- values[,c(c1,c2)]
    testclab <- c(rep(0,length(c1)),rep(1,length(c2)))

    if (min(length(c1), length(c2)) >= 2) {
        tstat <- mt.teststat(testval, testclab, ...)
        pvalue <- pval(testval, testclab, tstat)
    } else {
        message("Too few samples per class, skipping t-test.")
        tstat <- pvalue <- rep(NA,nrow(testval))
    }
    stat <- data.frame(fold = fold, tstat = tstat, pvalue = pvalue)
    if (length(levels(sampclass(object))) >2) {
        pvalAnova<-c()
        for(i in 1:nrow(values)){
            var<-as.numeric(values[i,])
            ano<-summary(aov(var ~ sampclass(object)) )
            pvalAnova<-append(pvalAnova, unlist(ano)["Pr(>F)1"])
        }
        stat<-cbind(stat, anova= pvalAnova)
    }
    if (metlin) {
        neutralmass <- groupmat[,"mzmed"] + ifelse(metlin < 0, 1, -1)
        metlin <- abs(metlin)
        digits <- ceiling(-log10(metlin))+1
        metlinurl <- paste("http://metlin.scripps.edu/metabo_list.php?mass_min=",
                           round(neutralmass - metlin, digits), "&mass_max=",
                           round(neutralmass + metlin, digits), sep="")
        values <- cbind(metlin = metlinurl, values)
    }
    if(length(groupnames)>0)
		twosamp <- cbind(name = groupnames, stat, groupmat, values)
	else
		twosamp <- cbind(name = groupnames(object), stat, groupmat, values)

    if (sortpval) {
        tsidx <- order(twosamp[,"pvalue"])
        twosamp <- twosamp[tsidx,]
        rownames(twosamp) <- 1:nrow(twosamp)
        values<-values[tsidx,]
    } else
        tsidx <- 1:nrow(values)

    #if (length(filebase))
        #write.table(twosamp, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)

    if (eicmax > 0) {
        if (length(unique(peaks(object)[,"rt"])) > 1) {
            ## This looks like "normal" LC data

            eicmax <- min(eicmax, length(tsidx))
            eics <- getEIC(object, rtrange = eicwidth*1.1, sampleidx = ceic,
                           groupidx = tsidx[seq(length = eicmax)],mzExpMeth="ppm",ppm=ppm)

            if (length(filebase)) {
                eicdir <- paste(filebase, "_eic", sep="")
                boxdir <- paste(filebase, "_box", sep="")
                dir.create(eicdir)
                dir.create(boxdir)
                if (capabilities("png")){
                    xcmsBoxPlot(values[seq(length = eicmax),],
                                sampclass(object), dirpath=boxdir, pic="png",  width=w, height=h)
                    png(file.path(eicdir, "%003d.png"), width = w, height = h)
                } else {
                    xcmsBoxPlot(values[seq(length = eicmax),],
                                sampclass(object), dirpath=boxdir, pic="pdf", width=w, height=h)
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                        height = h/72, onefile = FALSE)
                }
            }
            plot(eics, object, rtrange = eicwidth, mzdec=mzdec)

            if (length(filebase))
                dev.off()
        } else {
            ## This looks like a direct-infusion single spectrum
            if (length(filebase)) {
                specdir <- paste(filebase, "_spec", sep="")
                dir.create(specdir)
                if (capabilities("png")){
                    png(file.path(specdir, "%003d.png"), width = w, height = h)
                }else{
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                        height = h/72, onefile = FALSE)
                }
            }

            plotSpecWindow(object, gidxs = tsidx[seq(length = eicmax)], borderwidth=1)

            if (length(filebase))
                dev.off()
        }
    }

    invisible(twosamp)
}
