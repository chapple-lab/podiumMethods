scaleData = function(x,byCol=T,type=c("average","range"))
{
  if(byCol)
  {
    means = colMeans(x)
    for(y in 1:ncol(x))
    {
      if(type=="average")
      {
        x[,y] = (x[,y]-means[y])/means[y] #scale each column by means and center around 0
      }
      else if(type=="range")
      {
        x[,y] = (x[,y]-means[y])/(max(x[,y])-min(x[,y])) #scale each column by the magnitude of the range and center around 0
      }
      else
      {
        stop("Incorrect standardization type")
      }
    }

  }
  else
  {
    means = rowMeans(x)
    for(y in 1:nrow(x))
    {
      if(type=="average")
      {
        x[y, ] =  (x[y, ]-means[y])/means[y] #scale each row by means and center around 0
      }else if(type=="range")
      {
        x[y, ] = (x[y, ]-means[y])/(max(x[y, ])-min(x[y, ])) #scale each row by the magnitude of the range and center around 0
      }else
      {
        stop("Incorrect standardization type")
      }

    }
  }
  # do.call(rbind,x)
  invisible(x)
}