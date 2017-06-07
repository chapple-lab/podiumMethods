getPeakAnnotations = function(object, intval="into") {

  if (!sum(intval == c("into","intb","maxo"))){
    stop("unknown intensity value!")
  }

  #allocate variables for CAMERA output
  adduct   <- vector("character", nrow(object@groupInfo));
  calcMass = vector("double", nrow(object@groupInfo))
  isotopes <- vector("character", nrow(object@groupInfo));
  pcgroup  <- vector("numeric", nrow(object@groupInfo));

  #default polarity set to positive
  polarity <- "+";

  if(length(object@polarity) > 0){
    if(object@polarity == "negative"){
      polarity <- "-";
    }
  }

  #First isotope informationen and adduct informationen
  for(i in seq(along = isotopes)){
    #check if adduct annotation is present for peak i
    if(length(object@derivativeIons) > 0 && !(is.null(object@derivativeIons[[i]]))) {
      #Check if we have more than one annotation for peak i
      if(length(object@derivativeIons[[i]]) > 1) {
        #combine ion species name and rounded mass hypophysis
        names <- paste(object@derivativeIons[[i]][[1]]$name);
        masses = signif(object@derivativeIons[[i]][[1]]$mass,6)
        for(ii in 2:length(object@derivativeIons[[i]])) {
          names <- paste(names, object@derivativeIons[[i]][[ii]]$name);
          masses = paste(masses,signif(object@derivativeIons[[i]][[ii]]$mass, 6), sep=",")
        }
        #save name in vector adduct
        adduct[i] <- names;
        calcMass[i] = masses
      } else {
        #Only one annotation
        adduct[i] <- paste(object@derivativeIons[[i]][[1]]$name);
        calcMass[i] = signif(object@derivativeIons[[i]][[1]]$mass, 6)
      }
    } else {
      #no annotation empty name
      adduct[i] <- "";
      calcMass[i] = NA
    }

    #Check if we have isotope informationen about peak i
    if(length(object@isotopes) > 0&& !is.null(object@isotopes[[i]])) {
      num.iso <- object@isotopes[[i]]$iso;
      #Which isotope peak is peak i?
      if(num.iso == 0){
        str.iso <- "[M]";
      } else {
        str.iso <- paste("[M+", num.iso, "]", sep="")
      }
      #Multiple charged?
      if(object@isotopes[[i]]$charge > 1){
        isotopes[i] <- paste("[", object@isotopes[[i]]$y, "]", str.iso, object@isotopes[[i]]$charge, polarity, sep="");
      }else{
        isotopes[i] <- paste("[", object@isotopes[[i]]$y, "]", str.iso, polarity, sep="");
      }
    } else {
      #No isotope informationen available
      isotopes[i] <- "";
    }
  }

  #Have we more than one pseudospectrum?
  if(length(object@pspectra) < 1){
    pcgroup <- 0;
  } else {
    for(i in seq(along = object@pspectra)){
      index <- object@pspectra[[i]];
      pcgroup[index] <- i;
    }
  }

  return(invisible(data.frame(isotopes, adduct, calcMass, pcgroup, stringsAsFactors=FALSE, row.names=NULL)));
}