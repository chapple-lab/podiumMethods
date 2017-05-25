colePrintTIC <-
function(xcmsSet=NULL, pdfname="TICs.pdf",rt=c("raw","corrected"))
{
  files <- filepaths(xcmsSet)
  N <- length(files)
  TIC <- vector("list",N)

  for (i in 1:N)
  {
		cat(files[i],"\n")
		if (!is.null(xcmsSet) && rt == "corrected")
		rtcor <- xcmsSet@rt$corrected[[i]] else
		rtcor <- NULL
		TIC[[i]] <- coleGetTIC(files[i],rtcor=rtcor)
  }

  pdf(pdfname,w=16,h=10)
  cols <- rainbow(N)
  lty = 1:N
  pch = 1:N
  xlim = range(sapply(TIC, function(x) range(x[,1])))
  ylim = range(sapply(TIC, function(x) range(x[,2])))
  plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time", ylab = "TIC")
  for (i in 1:N)
  {
    tic <- TIC[[i]]
    points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
  }
  legend("topright",paste(rownames(xcmsSet@phenoData)), col = cols, lty = lty, pch = pch)
  dev.off()

  invisible(TIC)
}
