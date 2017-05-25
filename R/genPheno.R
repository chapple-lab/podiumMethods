genPheno <-
function(xcmsSet)
{
  filePattern = "\\.[Cc][Dd][Ff]$|\\.[Nn][Cc]$|\\.([Mm][Zz])?[Xx][Mm][Ll]$|\\.[Mm][Zz][Dd][Aa][Tt][Aa]$|\\.[Mm][Zz][Mm][Ll]$"
  datePattern = "(19|20)\\d\\d([- /._])(0[1-9]|1[012])\\2(0[1-9]|[12][0-9]|3[01])"
  nameStrings = sampnames(xcmsSet)
  numSamps = length(nameStrings)
  nameStrings = gsub("(12C)|(13C)","",nameStrings,perl=T) #remove 12C and 13C
  nameStrings = gsub(filePattern,"",nameStrings,perl=T) #remove file extensions
  nameStrings = gsub(datePattern,"",nameStrings,perl=T) #remove dates
  # nameStrings = gsub("^|$"," ",nameStrings,perl=T) #prep for next step
  #Old:  delims = strsplit(nameStrings,"\\.+|(?<=[A-EI-QS-Za-ei-qs-z])[0-9]|[0-9](?=[A-EI-QS-Za-ei-qs-z])",perl=T) #remove any numbers that are part of words. protects against destruction of ref, fah, f5h
  delims = strsplit(nameStrings," ",perl=T)
  delims = unlist(delims)
  delims = delims[!(delims%in%"")] #remove "" delimiters
  udelims = unique(delims)
  idx = sapply(udelims,function(x,numSamps) length(delims[delims%in%x])==numSamps,numSamps=numSamps) #idx for removing all ""
  phenotype = paste(udelims[idx],collapse="_")
  return(phenotype)
}
