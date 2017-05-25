removeGroups <-
function(xcmsSet2,rtrng=NULL)
{
	#Remove all groups and related information from xcmsSet2 that have rtmed within the rt range
	idx = which((xcmsSet2@groups[ ,"rtmed"]<=rtrng[2])&(xcmsSet2@groups[ ,"rtmed"]>=rtrng[1]))
	xcmsSet2@groups = xcmsSet2@groups[-idx, ]
	xcmsSet2@groupidx = xcmsSet2@groupidx[-idx]
	
	return (xcmsSet2)
}
