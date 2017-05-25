getTIC <-
function(xcmsRaw,rtcor=NULL)
{
  cbind(if (is.null(rtcor)) xcmsRaw@scantime else rtcor, rawEIC(xcmsRaw,mzrange=range(xcmsRaw@env$mz))$intensity)
}
