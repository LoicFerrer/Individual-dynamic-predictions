require("survival")

NAfix <- function(x, subst=-Inf) {
### Written by Christian Hoffmann; propagate last known non-NA value
### Input:
###     x: numeric vector
###     subst: scalar inidicating which value should replace NA
###         if x starts with a series of NA's
### Output:
###     (numeric) vector, with NA's replaced by last known non-NA value,
###         or 'subst'
    spec <- max(x[!is.na(x)])+1
    x <- c(spec,x)
    while (any(is.na(x))) x[is.na(x)] <- x[(1:length(x))[is.na(x)]-1]
    x[x==spec] <- subst
    x <- x[-1]
    x
}

CumInc <- function(Haz)
{
  ### Calculate cumulative incidence function
  ### Input:
  ###     Haz: dataframe with (at least) columns
  ###         time
  ###         Haz, cumulative Hazard
  ###         cause, the failure cause to which these cumulative hazards belong
  ### Output:
  ###     A dataframe with columns time, CI1 ... CIK 
  ###     containing the estimated cumulative incidences of causes 1 ... K
  ###     and S0 containing the estimated failure-free probability at these
  ###     time-points
  ### Details:
  ###     Haz is assumed to contain all cumulative hazards stacked
  ###     on top of each other. The column cause is used to distinguish
  ###     between the hazards of the different transitions or failure causes
  causes <- unique(Haz$cause)
  K <- length(causes)
  tt <- sort(unique(Haz$time))
  n <- length(tt)
  dfrm <- matrix(NA, n, 3*K+4)
  dfr <- as.data.frame(dfrm)
  names(dfr)[1] <- "time"
  names(dfr)[2:(K+1)] <- paste("Haz",as.character(1:K), sep="")
  names(dfr)[(K+2):(2*K+1)] <- paste("haz",as.character(1:K), sep="")
  names(dfr)[(2*K+2):(3*K+1)] <- paste("CI",as.character(1:K), sep="")
  names(dfr)[3*K+2] <- "hazsum"
  names(dfr)[3*K+3] <- "Hazsum"
  names(dfr)[3*K+4] <- "S0"
  dfr$time <- tt
  for (k in 1:K) { # select the elements for each cause
    wh <- which(Haz$cause==causes[k])
    idx <- match(Haz$time[wh],tt)
    dfr[,k+1][idx] <- Haz$Haz[wh]
    dfr[,k+1] <- NAfix(dfr[,k+1],subst=0)
    dfr[,K+1+k] <- diff(c(0,dfr[,k+1])) # the hazard for cause k
  }
  n <- nrow(dfr)
  for (k in 1:K) # each cumulative hazard is adjusted (cumulative sum of hazard)
    dfr[,k+1] <- cumsum(dfr[,K+1+k])
  ### compute baseline survival S0, need to sum all columns Hazards
  if (K==1) dfr$hazsum <- dfr[,3]
  else dfr$hazsum <- apply(dfr[,((K+2):(2*K+1))],1,sum)
  dfr$S0 <- {
    dfr$hazsum[dfr$hazsum > 1] <- 1 # added
    cumprod(1-dfr$hazsum)
  }
  dfr$Hazsum <- -log(dfr$S0)
  ### compute cumulative incidence function
  for (k in 1:K) {
    dfr[,2*K+k+1] <- cumsum(c(1,dfr$S0[-n])*dfr[,K+1+k])
    dfr[,2*K+k+1][dfr[,2*K+k+1] > 1] <- 1
  }
  #    return(dfr) # this gives complete
  return(dfr[,c(1,(2*K+2):(3*K+1),3*K+4)]) # this gives minimum (CI's, S0)
}



survfitLMCR <- function (...) 
{
  UseMethod("survfitLMCR")
}

print.survfitLMCR <- function (x, ...) 
{
  if (x$simulate) {
    cat("\nPredicted individual cumulative incidences of events\n\n")
    print(x$res)
    cat("\nPredicted individual cumulative incidences of events\n\tbased on", 
        x$M , "Monte Carlo samples\n\n")
    print(x$res.MC)
  }
  else {
    cat("\nPredicted individual cumulative incidences of events\n")
    print(x$res)
  }
  invisible(x)
}

