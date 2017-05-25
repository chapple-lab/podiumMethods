coleGetTIC <-
function(file,rtcor=NULL)
{
    object = xcmsRaw(file)
	cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity)
}
