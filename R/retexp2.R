retexp2 <-
function(peakrange, symExp = 20)
{
    peakrange[,"rtmin"] <- peakrange[,"rtmin"]-symExp/2
    peakrange[,"rtmax"] <- peakrange[,"rtmax"]+symExp/2

    peakrange
}
