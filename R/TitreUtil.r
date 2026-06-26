
# https://www.researchgate.net/post/How-I-campute-a-distance-between-a-cumulative-distribution-and-an-empirical-distribution-in-R
# set.seed(112358); x <- rnorm(5); ks.test(x=x, y=pnorm)$statistic

ProbDisagreeGaussian<-function(sigma, dmax=8, verbose=F){
   r <-seq((-1*dmax)-0.5, dmax+0.5, 1)
   d = diff(pnorm(mean=0, sd=sigma, q=r))
   if(verbose){
      print("d=")
      print(d)
   }
   ProbMatrix<-outer(d, d, FUN="*")
   # sum(ProbMatrix)
   diagIndicator <- abs(row(ProbMatrix) - col(ProbMatrix))
   DiagonalSum<-lapply(split(ProbMatrix, diagIndicator), sum)
   return(unlist(DiagonalSum))
}

require(ExtDist)
   
ProbDisagreeLaplace<-function(b, dmax=8, verbose=F){
   r <-seq((-1*dmax)-0.5, dmax+0.5, 1)
   d = diff(pLaplace(mean=0, b=b, q=r))
   if(verbose){
      print("d=")
      print(d)
   }
   ProbMatrix<-outer(d, d, FUN="*")
   # sum(ProbMatrix)
   diagIndicator <- abs(row(ProbMatrix) - col(ProbMatrix))
   DiagonalSum<-lapply(split(ProbMatrix, diagIndicator), sum)
   return(unlist(DiagonalSum))
}

ProbDisagreeLaplace(1)


ProbDisagreeGaussianRepro<-function(sigmaD, meanD, dmax=10, verbose=F){
   # in this case the bivariate distribution is not centred on zero
   r <-seq((-1*dmax)-0.5, dmax+0.5, 1)
   ProbMatrix<-matrix(0, ncol=length(r), nrow=I(length(r)-1))
   for(i in 1:I(length(r)-1)){
      for(j in 1:I(length(r)-1)){
         LowerVector<-c(r[i]  , r[j]  )
         UpperVector<-c(r[i+1], r[j+1])
         MeanVector <-c(meanD, 0)
         sigmaMatrix<-matrix(c(sigmaD, 0, 0, sigmaD), ncol=2)
         ProbMatrix[i,j]<-pmvnorm(
            lower=LowerVector, 
            upper=UpperVector, 
            mean =MeanVector,
            sigma=sigmaMatrix)
      }
   }
   # if(verbose){
   #    print("d=")
   #    print(d)
   # }
   if(verbose){
      print(unlist(list("Sum of probability matrix"=sum(ProbMatrix))))
   }
   # diagIndicator <- abs(row(ProbMatrix) - col(ProbMatrix))
   diagIndicator <- row(ProbMatrix) - col(ProbMatrix)
   DiagonalSum<-lapply(split(ProbMatrix, diagIndicator), sum)
   return(unlist(DiagonalSum))
}


require(mvtnorm)
    round(ProbDisagreeGaussian(sigma=1, dmax=5), 3)
sum(ProbDisagreeGaussian(sigma=1, dmax=5))
    # ProbDisagreeGaussian(sigma=1, dmax=3, verbose=T)

        round(ProbDisagreeGaussianRepro(sigmaD=1, meanD=0, dmax=5), 3)
        round(ProbDisagreeGaussianRepro(sigmaD=1, meanD=3, dmax=5), 3)

    
EFn<-function(location, order){
   M<-matrix(0, nrow=order, ncol=order)
   M[location, location]<-1
   return(M)
}

# EFn(location=1, order=4)
# 
# D<-2
# OrderScalar<-1+(2*D)
# 
# TestMatrix<-outer(1:OrderScalar, 1:OrderScalar, FUN="*")
# TestMatrix
# 
# EFn(1, OrderScalar) %*% TestMatrix %*% EFn(OrderScalar+1-1, OrderScalar)
# 
# # ProbVector<-rep(NA, D+1)
# ProbVector<-rep(NA, OrderScalar)
# 
# 
# SumMatrix<-matrix(0, ncol=OrderScalar, nrow=OrderScalar)
# for(j in 1:OrderScalar){
#    SumMatrix<-SumMatrix + EFn(j, OrderScalar) %*% TestMatrix %*% EFn(OrderScalar+1-j, OrderScalar)
# }
# SumMatrix
# ProbVector[1]<-matrix(1, ncol=OrderScalar, nrow=1) %*% SumMatrix %*% matrix(1, ncol=1, nrow=OrderScalar)
# ProbVector[1]
# 
# 
# for(i in 1:length(ProbVector)){
#    SumMatrix<-matrix(0, ncol=OrderScalar, nrow=OrderScalar)
#    for(j in 1:I(OrderScalar-(i-1))){
#       SumMatrix<-SumMatrix + EFn(j, OrderScalar) %*% TestMatrix %*% EFn((OrderScalar-(i-1))+1-j, OrderScalar)
#    }
#    SumMatrix
#    ProbVector[i]<-matrix(1, ncol=OrderScalar, nrow=1) %*% SumMatrix %*% matrix(1, ncol=1, nrow=OrderScalar)
#    if(i>1){
#       ProbVector[i]<-2*ProbVector[i]
#    }
# }
# ProbVector

ProbDisagreeGaussianMatrix<-function(sigma, dmax=8, verbose=F){
   r <-seq((-1*dmax)-0.5, dmax+0.5, 1)
   d = diff(pnorm(mean=0, sd=sigma, q=r))
   if(verbose){
      print("d=")
      print(d)
   }
   OrderScalar<-length(d)
   ProbMatrix<-outer(d, d, FUN="*")
   if(verbose){
      print(round(ProbMatrix, 2))
      print(unlist(list("sum of probability matrix"=sum(ProbMatrix))))
   }

   ProbVector<-rep(NA, dmax)
   
   for(i in 1:length(ProbVector)){
      SumMatrix<-matrix(0, ncol=OrderScalar, nrow=OrderScalar)
      for(j in 1:I(OrderScalar-(i-1))){
         SumMatrix<-SumMatrix + EFn(j, OrderScalar) %*% ProbMatrix %*% EFn((OrderScalar-(i-1))+1-j, OrderScalar)
      }
      if(verbose){print(round(SumMatrix, 2))}
      ProbVector[i]<-matrix(1, ncol=OrderScalar, nrow=1) %*% SumMatrix %*% matrix(1, ncol=1, nrow=OrderScalar)
      if(i>1){
         ProbVector[i]<-2*ProbVector[i]
      }
   }

   return(ProbVector)
}
# ProbDisagreeGaussianMatrix(sigma=1, dmax=3, verbose=T)


ProbDisagree<-function(p, d){
  ((1-sign(d))*(            ((2*p/((2-p)^2))) *(       (p/2) + (((1-p)^2))/(2-p)) )) +
    ((  sign(d))*(((2*p*((1-p)^d))/((2-p)^2)))  *( ((d-1)*p)   +          (2/(2-p))  ))
}

inverse.logit<-function(x){1/(1+exp(-x))}
logit        <-function(p){log(p/(1-p))}

deviance<-function(logitp, x, verbose=F){
  p<-1/(1+exp(-logitp))
  logL<-log(ProbDisagree(p, x))
  if(verbose){
    print(table(exp(logL), x))
  }
  return(-2*sum(logL))
}

FitGeometric<-function(x){
  temp.nlminb<-nlminb(start=c(0), objective=deviance, x=x)
  
  pFitScalar<-inverse.logit(temp.nlminb$par)
  # print(pFitScalar)
  # print(unlist(list("fitted p"=pFitScalar, "deviance"=temp.nlminb$objective)))
  
  ExpectedxVector<-0:max(x)
  # print(ProbDisagree(pFitScalar, ExpectedxVector))
  # print(sum(ProbDisagree(pFitScalar, ExpectedxVector)))
  
  # print(table(x)/length(x))
  return(list(
    p=pFitScalar,
    FittedProbs=ProbDisagree(pFitScalar, ExpectedxVector),
    deviance=temp.nlminb$objective
  ))
}

devianceGaussianCensored<-function(logsigma, x, censored=rep(0, length(x)), verbose=F){
  sigma<-exp(logsigma)
  
  # assume max(x) is no more than 3, i.e. max number of dilutions of disagreement is no more than 3
  if(max(x)>3){stop("Error: the function is not defined for disagreements of more than 3 dilutions.")}
  
  # allow censoring at one dilution only
  if(any(as.logical(censored))){
    MaxCensoredDifference<-max(x[as.logical(censored)])
    if(MaxCensoredDifference==0){
      stop("Error: censored differences of 0 are not informative and should not be included.")
    }
    if(MaxCensoredDifference> 1){
      stop("Error: censored differences of more than 1 are not covered.")
    }
  }
  if(max(x)>3){stop("Error: the function is not defined for disagreements of more than 3 dilutions.")}
  
  ProbVectorShort<-ProbDisagreeGaussian(sigma=sigma)
  ProbVector<-ifelse(x==0 & !as.logical(censored),   ProbVectorShort["0"], 
                     ifelse(x==1 & !as.logical(censored),   ProbVectorShort["1"],
                            ifelse(x==2 & !as.logical(censored),   ProbVectorShort["2"], 
                                   ifelse(x==3 & !as.logical(censored),   ProbVectorShort["3"], 
                                          ifelse(x==1 &  as.logical(censored), 1-ProbVectorShort["0"], NA)))))
  logL<-log(ProbVector)
  if(verbose){
    print(table(exp(logL), x, censored, exclude=NULL))
  }
  return(-2*sum(logL))
}

devianceLaplaceCensored<-function(logb, x, censored=rep(0, length(x)), verbose=F){
  b<-exp(logb)
  
  # assume max(x) is no more than 3, i.e. max number of dilutions of disagreement is no more than 3
  if(max(x)>3){stop("Error: the function is not defined for disagreements of more than 3 dilutions.")}
  
  # allow censoring at one dilution only
  if(any(as.logical(censored))){
    MaxCensoredDifference<-max(x[as.logical(censored)])
    if(MaxCensoredDifference==0){
      stop("Error: censored differences of 0 are not informative and should not be included.")
    }
    if(MaxCensoredDifference> 1){
      stop("Error: censored differences of more than 1 are not covered.")
    }
  }
  if(max(x)>3){stop("Error: the function is not defined for disagreements of more than 3 dilutions.")}
  
  ProbVectorShort<-ProbDisagreeLaplace(b=b)
  ProbVector<-ifelse(x==0 & !as.logical(censored),   ProbVectorShort["0"], 
                     ifelse(x==1 & !as.logical(censored),   ProbVectorShort["1"],
                            ifelse(x==2 & !as.logical(censored),   ProbVectorShort["2"], 
                                   ifelse(x==3 & !as.logical(censored),   ProbVectorShort["3"], 
                                          ifelse(x==1 &  as.logical(censored), 1-ProbVectorShort["0"], NA)))))
  logL<-log(ProbVector)
  if(verbose){
    print(table(exp(logL), x, censored, exclude=NULL))
  }
  return(-2*sum(logL))
}

devianceSingleError<-function(logitq, x, verbose=F){
  q<-1/(1+exp(-logitq))
  if(max(x)>2){stop("Error: the function is not defined for disagreements of more than 2 dilutions.")}
  ProbVector<-ifelse(x==0, (q^2) + (((1-q)^2)/2),
                     ifelse(x==1, 2*q*(1-q), ((1-q)^2)/2))
  logL<-log(ProbVector)
  if(verbose){
    print(table(exp(logL), x))
  }
  return(-2*sum(logL))
}

devianceSingleErrorCensored<-function(logitq, x, censored=rep(0, length(x)), verbose=F){
  q<-1/(1+exp(-logitq))
  if(max(x)>2){stop("Error: the function is not defined for disagreements of more than 2 dilutions.")}
  
  # allow censoring at one dilution only
  if(any(as.logical(censored))){
    MaxCensoredDifference<-max(x[as.logical(censored)])
    if(MaxCensoredDifference==0){
      stop("Error: censored differences of 0 are not informative and should not be included.")
    }
    if(MaxCensoredDifference> 1){
      stop("Error: censored differences of more than 1 are not covered.")
    }
  }
  
  ProbVectorShort<-unlist(list(
    "0"=(q^2) + (((1-q)^2)/2),
    "1"=2*q*(1-q),
    "2"=((1-q)^2)/2
  ))
  # ProbVector<-ifelse(x==0, (q^2) + (((1-q)^2)/2), 
  #             ifelse(x==1, 2*q*(1-q), ((1-q)^2)/2))
  ProbVector<-ifelse(x==0 & !as.logical(censored),   ProbVectorShort["0"], 
                     ifelse(x==1 & !as.logical(censored),   ProbVectorShort["1"],
                            ifelse(x==2 & !as.logical(censored),   ProbVectorShort["2"], 
                                   ifelse(x==1 &  as.logical(censored), 1-ProbVectorShort["0"], NA))))
  logL<-log(ProbVector)
  if(verbose){
    print(table(exp(logL), x))
  }
  return(-2*sum(logL))
}

