
# for local development only
if(Sys.info()[['user']]=="eidenale"){
  setwd("C:/Users/eidenale/work/project/other/scrub typhus/analysis/TitrationAgreement/examples")
}

# source the workhorse functions
source("../R/TitreUtil.r")

# data reconstructed from Figure 1 in the following:
# Cauchemez S, Horby P, Fox A, Mai le Q, Thanh le T, Thai PQ, et al. Influenza infection rates, measurement errors and the interpretation of paired serology. PLoS Pathog. 2012;8(12):e1003061.
#
cauchemezShortDF<-read.csv(file="../data/cauchemezShort.csv",
                           stringsAsFactors = FALSE)
cauchemezShortDF<-cauchemezShortDF[cauchemezShortDF$count>0,]
cauchemezShortDF$mean<-(as.numeric(substr(cauchemezShortDF$Value, 1, 1))+as.numeric(substr(cauchemezShortDF$Value, 3, 3)))/2
cauchemezShortDF$mean<-ifelse(as.numeric(cauchemezShortDF$censored), NA, cauchemezShortDF$mean)
cauchemezShortDF

CauchemezSummary<-as.data.frame(tapply(cauchemezShortDF$count, 
                                       list(cauchemezShortDF$difference, cauchemezShortDF$censored), 
                                       sum))
# class(CauchemezSummary)
CauchemezSummary[,"1"]<-ifelse(is.na(CauchemezSummary[,"1"]), 0, CauchemezSummary[,"1"])
CauchemezSummary

AbsxCauchemezDF<-data.frame(
  difference=c(
    rep(as.numeric(row.names(CauchemezSummary)), CauchemezSummary[,"0"]), 
    rep(as.numeric(row.names(CauchemezSummary)), CauchemezSummary[,"1"])),
  censored  =c(
    rep( rep(0, nrow(CauchemezSummary)) , CauchemezSummary[,"0"]), 
    rep( rep(1, nrow(CauchemezSummary)) , CauchemezSummary[,"1"]))
)
head(AbsxCauchemezDF)
table(AbsxCauchemezDF$difference, AbsxCauchemezDF$censored)

# data from:
#    White, C., 1973. Serological Epidemiology. Academic Press, New York. Book section. Statistical Methods in Serum Surveys.
# which were re-used in:
#    Thrusfield M. Veterinary Epidemiology. New York: John Wiley & Sons; 2018.

ThrusfieldDF<-read.table("../data/thrusfield.txt",
   stringsAsFactors=FALSE, header=T, sep="\t")

head(ThrusfieldDF)

table(ThrusfieldDF$censored1, ThrusfieldDF$censored2, exclude=NULL)

# tabulate the non-censored ones and use them for heteroscedasticity analysis
ThrusSubset<-!as.logical(ThrusfieldDF$censored1) & !as.logical(ThrusfieldDF$censored2)
table(ThrusSubset)
ThrusTable<-tapply(ThrusfieldDF$count[ThrusSubset], ThrusfieldDF$d[ThrusSubset], sum)
ThrusTable

AbsxThrusfield<-rep(as.numeric(names(ThrusTable)), ThrusTable)
table(AbsxThrusfield)

quantile(rep(ThrusfieldDF$meanDilution, ThrusfieldDF$count), probs=0.5, na.rm = TRUE)
table(rep(ThrusfieldDF$meanDilution, ThrusfieldDF$count)<=3)
table(rep(ThrusfieldDF$meanDilution, ThrusfieldDF$count)<=2.5)

# View(ThrusfieldDF)

ThrusfieldLoMeanSubset<-ThrusfieldDF$meanDilution<=2.5 | is.na(ThrusfieldDF$meanDilution)
table(ThrusfieldLoMeanSubset)

ThrusfieldLoMeanDF<-ThrusfieldDF[ ThrusfieldLoMeanSubset,]
sum(ThrusfieldLoMeanDF$count)
ThrusfieldHiMeanDF<-ThrusfieldDF[!ThrusfieldLoMeanSubset,]
sum(ThrusfieldHiMeanDF$count)

# heteroscedasticity
ThrusLongWithoutCensoredDF<-data.frame(
         mean=rep(ThrusfieldDF$meanDilution[ThrusSubset], ThrusfieldDF$count[ThrusSubset]),
   difference=rep(ThrusfieldDF$d[ThrusSubset]           , ThrusfieldDF$count[ThrusSubset]))
HeteroscadisticityTable<-table(ThrusLongWithoutCensoredDF$mean<=3)
HeteroscadisticityTable

ThrusLongDF<-data.frame(
   difference=rep(ThrusfieldDF$d                                                   , ThrusfieldDF$count),
   censored  =rep(as.numeric(ThrusfieldDF$censored1==1 | ThrusfieldDF$censored2==1), ThrusfieldDF$count))
table(ThrusLongDF$censored)
# length(ThrusLongDF[ThrusLongDF$censored!=1,"difference"])
table(ThrusLongDF[ThrusLongDF$censored!=1,"difference"], exclude=NULL)
table(ThrusLongDF$difference, exclude=NULL)

ThrusLoMeanLongDF<-data.frame(
   difference=rep(ThrusfieldLoMeanDF$d,
                  ThrusfieldLoMeanDF$count),
   censored  =rep(as.numeric(ThrusfieldLoMeanDF$censored1==1 | ThrusfieldLoMeanDF$censored2==1), 
                  ThrusfieldLoMeanDF$count))
dim(ThrusLoMeanLongDF)

ThrusHiMeanLongDF<-data.frame(
   difference=rep(ThrusfieldHiMeanDF$d,
                  ThrusfieldHiMeanDF$count),
   censored  =rep(as.numeric(ThrusfieldHiMeanDF$censored1==1 | ThrusfieldHiMeanDF$censored2==1), 
                  ThrusfieldHiMeanDF$count))
dim(ThrusHiMeanLongDF)

rm(AbsxThrusfield)


# FitGeometric(AbsxThrusfield)
FitGeometric(ThrusLongDF[ThrusLongDF$censored!=1,"difference"])

Hierholzer1<-rep(0:3, c(68, 83, 25, 3))
table(Hierholzer1)
FitGeometric(Hierholzer1)

Hierholzer2<-rep(0:2, c(48, 60, 10))
table(Hierholzer2)
FitGeometric(Hierholzer2)

Hierholzer3<-rep(0:2, c(62, 32, 2))
table(Hierholzer3)
FitGeometric(Hierholzer3)

Hierholzer4<-rep(0:2, c(21, 41, 18))
table(Hierholzer4)
FitGeometric(Hierholzer4)

Hierholzer5<-rep(0:2, c(22, 32, 2))
table(Hierholzer5)
FitGeometric(Hierholzer5)

# Table 6, within day reovirus HA (hemagglutination)
HierholzerReoHAVector<-rep(0:3, c(26, 23, 7, 3))
table(HierholzerReoHAVector)
FitHierholzerReoHAScalar<-FitGeometric(HierholzerReoHAVector)
FitHierholzerReoHAScalar


# https://stackoverflow.com/questions/27935555/get-all-diagonal-vectors-from-matrix

# https://stackoverflow.com/questions/24317929/discretizing-a-continuous-probability-distribution

r <-seq(-2.5, 2.5, 1)
r
d = diff(pnorm(mean=0, sd=.44, q=r))
d
ProbMatrix<-outer(d, d, FUN="*")
sum(ProbMatrix)

diagIndicator <- abs(row(ProbMatrix) - col(ProbMatrix))
diagIndicator
DiagonalSum<-lapply(split(ProbMatrix, diagIndicator), sum)
DiagonalSum
class(DiagonalSum)
sum(unlist(DiagonalSum))

# make this a function of sigma and get ML estimate of it (sigma)

# d is absolute difference

ProbDisagreeGaussian(sigma=1)[1:3]
ProbDisagreeGaussianMatrix(sigma=1)[1:3]

ProbDisagreeGaussian(sigma=1)["0"]

# AbsxCauchemez<-c(rep(0, 74), rep(1, 48), rep(2, 5))
# length(AbsxCauchemez)
# table(AbsxCauchemez)/length(AbsxCauchemez)

 table(AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])/
length(AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])

ProbDisagreeGaussian(sigma=0.5)[1:3]
ProbDisagreeGaussian(sigma=0.45)[1:3]
ProbDisagreeGaussian(sigma=0.44)[1:3]

# estimate sigma by ML
devianceGaussian<-function(logsigma, x, verbose=F){
   sigma<-exp(logsigma)
   # assume max(x) is no more than 3, i.e. max number of dilutions of disagreement is no more than 3
   if(max(x)>3){stop("Error: the function is not defined for disagreements of more than 3 dilutions.")}
   ProbVectorShort<-ProbDisagreeGaussian(sigma=sigma)
   ProbVector<-ifelse(x==0, ProbVectorShort["0"], 
               ifelse(x==1, ProbVectorShort["1"],
               ifelse(x==2, ProbVectorShort["2"], ProbVectorShort["3"])))
   logL<-log(ProbVector)
   if(verbose){
      print(table(exp(logL), x))
   }
   return(-2*sum(logL))
}

devianceLaplace<-function(logb, x, verbose=F){
   b<-exp(logb)
   if(max(x)>3){stop("Error: the function is not defined for disagreements of more than 3 dilutions.")}
   ProbVectorShort<-ProbDisagreeLaplace(b=b)
   ProbVector<-ifelse(x==0, ProbVectorShort["0"], 
               ifelse(x==1, ProbVectorShort["1"],
               ifelse(x==2, ProbVectorShort["2"], ProbVectorShort["3"])))
   logL<-log(ProbVector)
   if(verbose){
      print(table(exp(logL), x))
   }
   return(-2*sum(logL))
}

head(AbsxCauchemezDF)


devianceSingleError(logitq=0, 
   x=AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])

devianceSingleErrorCensored(logitq=0, 
   x       =AbsxCauchemezDF$difference,
   censored=AbsxCauchemezDF$censored)


temp.nlminb<-nlminb(start=c(0), objective=devianceSingleError, 
   x=AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])
inverse.logit(temp.nlminb$par)
# 1 dilution 1 sided error:
(1-inverse.logit(temp.nlminb$par))/2

devianceGaussian(logsigma=log(0.5) , x=AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])
devianceGaussian(logsigma=log(0.45), x=AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])
devianceGaussian(logsigma=log(0.44), x=AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])

devianceLaplace(logb=log(1), x=AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])

devianceGaussianCensored(logsigma=log(0.5) , x=AbsxCauchemezDF$difference, censored=AbsxCauchemezDF$censored, verbose = TRUE)
devianceGaussianCensored(logsigma=log(0.45), x=AbsxCauchemezDF$difference, censored=AbsxCauchemezDF$censored, verbose = TRUE)
devianceGaussianCensored(logsigma=log(0.44), x=AbsxCauchemezDF$difference, censored=AbsxCauchemezDF$censored, verbose = TRUE)

devianceLaplaceCensored(logb=log(1), x=AbsxCauchemezDF$difference, censored=AbsxCauchemezDF$censored, verbose = TRUE)

temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, 
   x=AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])
exp(temp.nlminb$par)

temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, 
   x=AbsxCauchemezDF$difference, censored=AbsxCauchemezDF$censored)
exp(temp.nlminb$par)


 table(AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])/
length(AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])
ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])+1)]
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))

# how likely are 2-dilution errors?
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(1.5,2.5)))
# and 3:
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(2.5,3.5)))

# look at fitted probability of disagreements of up to three dilutions
ProbDisagreeGaussianTwoFold<-ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(AbsxCauchemezDF[!as.logical(AbsxCauchemezDF$censored),"difference"])+2)]
ProbDisagreeGaussianTwoFold



# Now for 1.5 fold dilutions.  
# On the same scale of dilutions as the two-fold ones, the intervals are of 
# width log_2(1.5), or about 0.585.
# And the first interval starts halfway between 0.585 and zero.  

ThreeHalvesDilutionNumber<-(log(1.5, base=2)/2)+(log(1.5, base=2)*(0:4))
ThreeHalvesDilutionNumber

# on the fold change scale we want to get to at least 8
2**ThreeHalvesDilutionNumber

ProbDisagreeGaussianThreeHalves<-function(sigma, r, verbose=T){
   d = diff(pnorm(mean=0, sd=sigma, q=r))
   ProbMatrix<-outer(d, d, FUN="*")
   if(verbose){print(unlist(list("Check by summing probability matrix"=sum(ProbMatrix))))}
   diagIndicator <- abs(row(ProbMatrix) - col(ProbMatrix))
   DiagonalSum<-lapply(split(ProbMatrix, diagIndicator), sum)
   return(unlist(DiagonalSum))
}
DisagreeThreeHalvesCauchemez<-ProbDisagreeGaussianThreeHalves(
   sigma=exp(temp.nlminb$par), 
   r=c(rev(-1*ThreeHalvesDilutionNumber), ThreeHalvesDilutionNumber), 
   verbose=T)
DisagreeThreeHalvesCauchemez[1:5]

# check they sum to 1
sum(DisagreeThreeHalvesCauchemez)
sum(ProbDisagreeGaussianTwoFold)

# CDF
#
cumsum(DisagreeThreeHalvesCauchemez)
rbind(DisagreeThreeHalvesCauchemez[1:4], ThreeHalvesDilutionNumber[1:4]-(log(1.5, base=2)/2))
temp.y<-cumsum(DisagreeThreeHalvesCauchemez[1:4])[c(1:4, 4)]
temp.x<-ThreeHalvesDilutionNumber[c(1:5)]-(log(1.5, base=2)/2)
2**max(temp.x)
rbind(temp.y, temp.x)

png("../figures/ThreeHalvesCDF.png")
yminScalar<--0.05
plot(
   y=temp.y, 
   x=temp.x, type="n", xaxt="n", bty="n",
   ylim=c(yminScalar,1), ylab="cumulative probability  (total both directions)", 
   xlim=c(0,3), xlab="fold change")
lines(
   y=temp.y, 
   x=temp.x, type="s")

ProbDisagreeGaussianTwoFold
lines(
   y=cumsum(ProbDisagreeGaussianTwoFold[1:3])[c(1:3, 3)], 
   x=0:3, type="s", lty=2)

tckScalar<-0.03
axis(1, at=log(1.5, base=2)*I(0:4), labels=1.5**I(0:4), tck=-1*tckScalar)
text(x=0:3, y=rep(yminScalar+1.5*tckScalar,3), labels=2**(0:3))
axis(1, at=0:3, tck=+tckScalar, labels=rep("", length(0:3)))
# https://stackoverflow.com/questions/33343466/put-tick-labels-of-only-x-axis-inside-plotting-area
# axis(1, at=c(0:71), NA, cex.axis=.7, font=1, tck=.01)
dev.off()


# with the 1.5-fold scheme, need a 2-dilution change to get specificity more than 80%, i.e. 2.25-fold
# with the 2-fold scheme, also need a 2-dilution change but that is a 4-fold

# https://stackoverflow.com/questions/33343466/put-tick-labels-of-only-x-axis-inside-plotting-area

# barplot(DisagreeThreeHalvesCauchemez[1:length(ThreeHalvesDilutionNumber)])

# Prob of two readings agreeing if there is 20.2 percent chance of 1-dilution error each way, 
# and no errors of two or more.
p1<-0.202
p0<-1-(2*p1)
p0
(p0*p0)+(2*p1*p1)
# 44%, but the observed value is 58%

# now re-do with the inferred value from the Gaussian method
p1<-diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))
p0<-1-(2*p1)
p0
(p0*p0)+(2*p1*p1)

# length(AbsxThrusfield)
# table(AbsxThrusfield)/length(AbsxThrusfield)

ProbDisagreeGaussian(sigma=0.5)[1:3]
ProbDisagreeGaussian(sigma=0.51)[1:3]
ProbDisagreeGaussian(sigma=0.49)[1:3]

# temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, x=AbsxThrusfield)
# exp(temp.nlminb$par)

temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, 
   x=ThrusLongDF[!as.logical(ThrusLongDF$censored),"difference"])
exp(temp.nlminb$par)


ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(ThrusLongDF[!as.logical(ThrusLongDF$censored),"difference"])+1)]

# table(AbsxThrusfield)/length(AbsxThrusfield)
 table(ThrusLongDF[!as.logical(ThrusLongDF$censored),"difference"])/
length(ThrusLongDF[!as.logical(ThrusLongDF$censored),"difference"])

diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))

# Hierholzer
# should do for the four columns in Table 3 (within day) and the two columns of Table 6

ChisqContingency<-function(O, E){
   return(sum(((O-E)^2)/E))
}

TruncateTableFn<-function(table, threshold=5, force=F, forceLength, verbose=F){
   if(force){
      while(length(table)>forceLength){
         table<-rev(c(sum(rev(table)[1:2]), rev(table)[3:length(table)]))
      }
   }else{
      while(rev(table)[1]<threshold){
         if(verbose){
            print("unforced table")
            print(table)
            print(unlist(list("first element of reversed table"=I(rev(table)[1]))))
            print(unlist(list(threshold=threshold)))
         }
         table<-rev(c(sum(rev(table)[1:2]), rev(table)[3:length(table)]))
      }
   }
   return(table)
}

# Table 3, adenovirus, HA
HierholzerTable3AdenoHA<-rep(0:3, c(53, 62, 5, 1))
table(HierholzerTable3AdenoHA)/length(HierholzerTable3AdenoHA)
temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian        , x=HierholzerTable3AdenoHA)
exp(temp.nlminb$par)
# the following should be the same, because, by default, no censoring is assumed
temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, x=HierholzerTable3AdenoHA)
exp(temp.nlminb$par)
ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(HierholzerTable3AdenoHA)+1)]
FitGeometric(HierholzerTable3AdenoHA)
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))

HierholzerTable3AdenoHI<-rep(0:2, c(65, 28, 5))
table(HierholzerTable3AdenoHI)/length(HierholzerTable3AdenoHI)
# temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, x=HierholzerTable3AdenoHI)
temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, x=HierholzerTable3AdenoHI)
exp(temp.nlminb$par)
ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(HierholzerTable3AdenoHI)+1)]
FitGeometric(HierholzerTable3AdenoHI)
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))

HierholzerTable3MyxHA<-rep(0:2, c(65, 43, 4))
table(HierholzerTable3MyxHA)/length(HierholzerTable3MyxHA)
# temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, x=HierholzerTable3MyxHA)
temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, x=HierholzerTable3MyxHA)
exp(temp.nlminb$par)
ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(HierholzerTable3MyxHA)+1)]
FitGeometric(HierholzerTable3MyxHA)
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))

HierholzerTable3MyxHI<-rep(0:1, c(87, 9))
table(HierholzerTable3MyxHI)/length(HierholzerTable3MyxHI)
# temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, x=HierholzerTable3MyxHI)
temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, x=HierholzerTable3MyxHI)
exp(temp.nlminb$par)
ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(HierholzerTable3MyxHI)+1)]
# table(HierholzerTable3AdenoHA)/length(HierholzerTable3AdenoHA)
FitGeometric(HierholzerTable3MyxHI)
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))

# Table 6, within day reovirus HA (hemagglutination)
HierholzerReoHAVector<-rep(0:3, c(26, 23, 7, 3))
 table(HierholzerReoHAVector)/length(HierholzerReoHAVector)
length(HierholzerReoHAVector)
# temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, x=HierholzerReoHAVector)
temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, x=HierholzerReoHAVector)
exp(temp.nlminb$par)
temp.nlminb$objective
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(1.5,2.5)))

ProbFitGauss<-ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(HierholzerReoHAVector)+1)]
ProbFitGauss
ExpectedNGauss<-ProbFitGauss*length(HierholzerReoHAVector)
ExpectedNGauss
# table(HierholzerReoHAVector)/length(HierholzerReoHAVector)
ProbFitGeometric<-FitGeometric(HierholzerReoHAVector)
ProbFitGeometric
ExpectedNGeometric<-ProbFitGeometric$FittedProbs*length(HierholzerReoHAVector)
ExpectedNGeometric
# the Gaussian fit for this one is not as good as from the geometric
ChisqContingency(O=table(HierholzerReoHAVector), E=ExpectedNGauss)
ChisqContingency(O=table(HierholzerReoHAVector), E=ExpectedNGeometric)

ExpectedNGauss
TruncateTableFn(ExpectedNGauss)
ExpectedNGeometric
TruncateTableFn(ExpectedNGeometric)

TableHierholzerReoHATruncated<-TruncateTableFn(
   table(HierholzerReoHAVector), force=T, forceLength=length(TruncateTableFn(ExpectedNGeometric)))
table(HierholzerReoHAVector)
TableHierholzerReoHATruncated
   
   
HierholzerReoHIVector<-rep(0:2, c(37, 18, 1))
# temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, x=HierholzerReoHIVector)
temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, x=HierholzerReoHIVector)
exp(temp.nlminb$par)
ProbDisagreeGaussian(sigma=exp(temp.nlminb$par))[1:I(max(HierholzerReoHIVector)+1)]
# table(HierholzerReoHAVector)/length(HierholzerReoHAVector)
FitGeometric(HierholzerReoHIVector)
diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))

HierholzerTable3AdenoHADF<-data.frame(difference=HierholzerTable3AdenoHA, censored=rep(0, length(HierholzerTable3AdenoHA)))
HierholzerTable3AdenoHIDF<-data.frame(difference=HierholzerTable3AdenoHI, censored=rep(0, length(HierholzerTable3AdenoHI)))
HierholzerTable3MyxHADF  <-data.frame(difference=HierholzerTable3MyxHA  , censored=rep(0, length(HierholzerTable3MyxHA)))
HierholzerTable3MyxHIDF  <-data.frame(difference=HierholzerTable3MyxHI  , censored=rep(0, length(HierholzerTable3MyxHI)))
HierholzerReoHADF        <-data.frame(difference=HierholzerReoHAVector  , censored=rep(0, length(HierholzerReoHAVector)))
HierholzerReoHIDF        <-data.frame(difference=HierholzerReoHIVector  , censored=rep(0, length(HierholzerReoHIVector)))

# ResultsListOriginal<-list(
#    AbsxCauchemezDF          =AbsxCauchemezDF,
#    ThrusLongDF              =ThrusLongDF,
#    HierholzerTable3AdenoHADF=HierholzerTable3AdenoHADF,
#    HierholzerTable3AdenoHIDF=HierholzerTable3AdenoHIDF,
#    HierholzerTable3MyxHADF  =HierholzerTable3MyxHADF,
#    HierholzerTable3MyxHIDF  =HierholzerTable3MyxHIDF,
#    HierholzerReoHADF        =HierholzerReoHADF,
#    HierholzerReoHIDF        =HierholzerReoHIDF,
#    AbsxCauchemezLoDF        =AbsxCauchemezLoDF,
#    AbsxCauchemezHiDF        =AbsxCauchemezHiDF,
#    ThrusLoMeanLongDF        =ThrusLoMeanLongDF,
#    ThrusHiMeanLongDF        =ThrusHiMeanLongDF
#    )

ResultsList<-list(
  AbsxCauchemezDF          =AbsxCauchemezDF,
  ThrusLongDF              =ThrusLongDF,
  HierholzerTable3AdenoHADF=HierholzerTable3AdenoHADF,
  HierholzerTable3AdenoHIDF=HierholzerTable3AdenoHIDF,
  HierholzerTable3MyxHADF  =HierholzerTable3MyxHADF,
  HierholzerTable3MyxHIDF  =HierholzerTable3MyxHIDF,
  HierholzerReoHADF        =HierholzerReoHADF,
  HierholzerReoHIDF        =HierholzerReoHIDF
)

length(ResultsList)
 names(ResultsList)

# ResultsSigmaVector<-rep(NA, length(ResultsList))
ResultsDF<-data.frame(
   dataset                    =names(ResultsList),
   meanUncensored             =rep(NA, length(ResultsList)),
   varianceUncensored         =rep(NA, length(ResultsList)),
   meanAll                    =rep(NA, length(ResultsList)),
   varianceAll                =rep(NA, length(ResultsList)),
   MMq                        =rep(NA, length(ResultsList)),
   MM1side1dilutionError      =rep(NA, length(ResultsList)),
   MMqchisq                   =rep(NA, length(ResultsList)),
   MMqchisqdf                 =rep(NA, length(ResultsList)),
   MMqchisqp                  =rep(NA, length(ResultsList)),
   MM1side1dilutionDeltaCIlo  =rep(NA, length(ResultsList)),
   MM1side1dilutionDeltaCIhi  =rep(NA, length(ResultsList)),
   MLq                        =rep(NA, length(ResultsList)),
   ML1side1dilutionError      =rep(NA, length(ResultsList)),
   MLqchisq                   =rep(NA, length(ResultsList)),
   MLqchisqdf                 =rep(NA, length(ResultsList)),
   MLqchisqp                  =rep(NA, length(ResultsList)),
   ML1side1dilutionBootCIlo   =rep(NA, length(ResultsList)),
   ML1side1dilutionBootCIhi   =rep(NA, length(ResultsList)),
   MLNdiscreteBootValues      =rep(NA, length(ResultsList)),
   sigma                      =rep(NA, length(ResultsList)),
   b                          =rep(NA, length(ResultsList)),
   Gausschisq                 =rep(NA, length(ResultsList)),
   Gausschisqdf               =rep(NA, length(ResultsList)),
   Gausschisqp                =rep(NA, length(ResultsList)),
   Gauss1side1dilutionError   =rep(NA, length(ResultsList)),
   GaussNdiscreteBootValues   =rep(NA, length(ResultsList)),
   Gauss1side1dilutionBootCIlo=rep(NA, length(ResultsList)),
   Gauss1side1dilutionBootCIhi=rep(NA, length(ResultsList)),
   KSstat                     =rep(NA, length(ResultsList)),
   KSstatLaplace              =rep(NA, length(ResultsList))
   )

par(mfrow=c(2,3))

# ResultsList[[1]]$difference

# for(i in 1:1){
for(i in 1:length(ResultsList)){

   print("##")
   print("##")
   print(names(ResultsList)[i])
   print("##")
   print("##")
   tableResults<-table(ResultsList[[i]]$difference, ResultsList[[i]]$censored)
   print(unlist(list(N=sum(tableResults))))
   print(tableResults)
   print(tableResults/sum(tableResults))
   ResultsDF[i, "meanUncensored"]    <-mean(ResultsList[[i]][!as.logical(ResultsList[[i]]$censored), "difference"])
   ResultsDF[i, "varianceUncensored"]<- var(ResultsList[[i]][!as.logical(ResultsList[[i]]$censored), "difference"])
   print(unlist(list(meanUncensored=ResultsDF[i, "meanUncensored"])))
   ResultsDF[i, "meanAll"]           <-mean(ResultsList[[i]][, "difference"])
   ResultsDF[i, "varianceAll"]       <- var(ResultsList[[i]][, "difference"])
   print("mean with censored values treated as if they were uncensored:")
   print(unlist(list(meanAll=ResultsDF[i, "meanAll"])))
   # print(mean(ResultsList[[i]][, "difference"]))

   ##
   ## model with only 1-dilution errors
   ##
   if(max(ResultsList[[i]][, "difference"])<=2){
      print("Single error model")

      # ResultsDF[i, "MMq"]<-sqrt(1-mean(ResultsList[[i]][!as.logical(ResultsList[[i]]$censored), "difference"]))
      ResultsDF[i, "MMq"]<-sqrt(1-mean(ResultsList[[i]][, "difference"]))
      print("Method of moments value of 1-sided single error (single-error model):")
      ResultsDF[i, "MM1side1dilutionError"]<-(1-ResultsDF[i, "MMq"])/2
      print(ResultsDF[i, "MM1side1dilutionError"])

      # StandardErrorqScalar<-sqrt((ResultsDF[i, "varianceUncensored"])/
      #    (length(ResultsList[[i]][!as.logical(ResultsList[[i]]$censored), "difference"])*4*(ResultsDF[i, "MMq"]**2)))
      StandardErrorqScalar<-sqrt((ResultsDF[i, "varianceAll"])/
         (length(ResultsList[[i]][, "difference"])*4*(ResultsDF[i, "MMq"]**2)))
      
      qHiScalar<-ResultsDF[i, "MMq"]+(qnorm(0.975)*StandardErrorqScalar)
      qLoScalar<-ResultsDF[i, "MMq"]-(qnorm(0.975)*StandardErrorqScalar)
      ResultsDF[i, "MM1side1dilutionDeltaCIhi"]<-(1-qLoScalar)/2
      ResultsDF[i, "MM1side1dilutionDeltaCIlo"]<-(1-qHiScalar)/2
      print("95% CI for 1-sided 1-dilution error, based on delta method SE:")
      print(ResultsDF[i, c("MM1side1dilutionDeltaCIlo", "MM1side1dilutionDeltaCIhi")])

      temp.nlminb<-nlminb(
         start=c(0), 
         objective=devianceSingleErrorCensored,
         x=ResultsList[[i]]$difference, censored=ResultsList[[i]]$censored)
      # temp.nlminb<-nlminb(start=c(0), objective=devianceSingleError, 
      #    x=ResultsList[[i]][, "difference"])
      ResultsDF[i, "MLq"]<-inverse.logit(temp.nlminb$par)
      print("ML value of 1-sided single error (single-error model with censoring):")
      ResultsDF[i, "ML1side1dilutionError"]<-(1-ResultsDF[i, "MLq"])/2
      print((1-ResultsDF[i, "ML1side1dilutionError"])/2)

      # qMM<-ResultsDF[i, "MMq"]
      # ProbFitSingleMM<-c((qMM^2) + (((1-qMM)^2)/2), 2*qMM*(1-qMM), ((1-qMM)^2)/2)
      qML<-ResultsDF[i, "MLq"]
      ProbFitSingleML<-c((qML^2) + (((1-qML)^2)/2), 2*qML*(1-qML), ((1-qML)^2)/2)
      # ExpectedNSingle<-ProbFitSingleMM*length(ResultsList[[i]][, "difference"])
      ExpectedNSingle<-ProbFitSingleML*length(ResultsList[[i]][, "difference"])
      rm(qML)

      tableResults<-table(ResultsList[[i]]$difference)
      forceScalar<-any(as.logical(ResultsList[[i]]$censored))
      ExpectedNSingleTruncated<-TruncateTableFn(
         ExpectedNSingle, force=T, 
         forceLength=ifelse(forceScalar, 2, length(ExpectedNSingleTruncated)))
      TableResultsTruncated   <-TruncateTableFn(
         as.vector(tableResults), force=T, 
         forceLength=ifelse(forceScalar, 2, length(ExpectedNSingleTruncated)))
      print("TableResultsTruncated:")
      print(TableResultsTruncated)
      print("ExpectedNSingleTruncated:")
      print(ExpectedNSingleTruncated)
      print(rbind(TableResultsTruncated, ExpectedNSingleTruncated))
      # ResultsDF[i, "MMqchisq"]  <-ChisqContingency(O=TableResultsTruncated, E=ExpectedNSingleTruncated)
      # ResultsDF[i, "MMqchisqdf"]<-length(TableResultsTruncated)-1
      # ResultsDF[i, "MMqchisqp"] <-1-pchisq(ResultsDF[i, "MMqchisq"], ResultsDF[i, "MMqchisqdf"])
      # print(unlist(list(chisq=ResultsDF[i, "MMqchisq"], df=ResultsDF[i, "MMqchisqdf"], p=ResultsDF[i, "MMqchisqp"])))
      ResultsDF[i, "MLqchisq"]  <-ChisqContingency(O=TableResultsTruncated, E=ExpectedNSingleTruncated)
      ResultsDF[i, "MLqchisqdf"]<-length(TableResultsTruncated)-1
      ResultsDF[i, "MLqchisqp"] <-1-pchisq(ResultsDF[i, "MLqchisq"], ResultsDF[i, "MLqchisqdf"])
      print(unlist(list(chisq=ResultsDF[i, "MLqchisq"], df=ResultsDF[i, "MLqchisqdf"], p=ResultsDF[i, "MLqchisqp"])))

      nBootScalar<-10000
      # nBootScalar<-1000
      if(nBootScalar>1){
         pBootVector<-rep(NA, nBootScalar)

         # DataVector<-ResultsList[[i]][!as.logical(ResultsList[[i]]$censored), "difference"]
         DataVector<-ResultsList[[i]][, "difference"]

         for(j in 1:nBootScalar){
            AbsxBoot<-DataVector[sample(1:length(DataVector), replace=T)]
            
            # MM
            pBootVector[j]<-sqrt(1-mean(AbsxBoot))

         }   
         ResultsDF[i, c("ML1side1dilutionBootCIlo", "ML1side1dilutionBootCIhi")]<-quantile(
            (1-pBootVector)/2, probs=c(0.025, 0.975))
         Prob1ErrorQuantile<-ResultsDF[i, c("ML1side1dilutionBootCIlo", "ML1side1dilutionBootCIhi")]
         hist((1-pBootVector)/2, 
              main=paste0(
                 names(ResultsList)[i],"\nquantiles=",
                 paste(round(Prob1ErrorQuantile, 3), collapse=",")))
         print("Bootstrap distribution of 1-sided 1-dilution (single-error model):")
         print(Prob1ErrorQuantile)
            rm(Prob1ErrorQuantile)
         ResultsDF[i, "MLNdiscreteBootValues"]<-length(table(pBootVector))
         print(unlist(list("number of distinct values"=ResultsDF[i, "MLNdiscreteBootValues"])))
      }

   }

   # Gaussian and Laplace models
   print("Gaussian model")
   # temp.nlminb<-nlminb(start=c(0), objective=devianceGaussian, 
   #    x=ResultsList[[i]][!as.logical(ResultsList[[i]]$censored), "difference"])
   temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored,
      x=ResultsList[[i]]$difference, censored=ResultsList[[i]]$censored)
   print(unlist(list(sigma=exp(temp.nlminb$par))))
   ResultsDF[i, "sigma"]<-exp(temp.nlminb$par)

   temp.nlminbLaplace<-nlminb(start=c(0), objective=devianceLaplaceCensored,
      x=ResultsList[[i]]$difference, censored=ResultsList[[i]]$censored)
   # print(unlist(list(sigma=exp(temp.nlminb$par))))
   ResultsDF[i, "b"]<-exp(temp.nlminbLaplace$par)

   # temp.nlminb$objective
   ResultsDF[i, "Gauss1side1dilutionError"]<-diff(pnorm(mean=0, sd=exp(temp.nlminb$par), q=c(0.5,1.5)))
   print(unlist(list(Prob1Dilution1SidedError=ResultsDF[i, "Gauss1side1dilutionError"])))

   if(nBootScalar>1){
      sigmaBootVector<-rep(NA, nBootScalar)

      # DataVector<-ResultsList[[i]]

      for(j in 1:nBootScalar){
         # AbsxBoot<-DataVector[sample(1:length(DataVector), replace=T)]
         AbsxBootDF<-ResultsList[[i]][sample(1:nrow(ResultsList[[i]]), replace=T),]
         temp.nlminb<-nlminb(start=c(0), objective=devianceGaussianCensored, 
            x=AbsxBootDF$difference, censored=AbsxBootDF$censored)
         sigmaBootVector[j]<-exp(temp.nlminb$par)
      }   
      print("Bootstrap distribution of sigma parameter of Gaussian model:")
      probsVector<-c(0.025, 0.975)
      sigmaQuantile<-quantile(sigmaBootVector, probs=probsVector)
      print(sigmaQuantile)
      ResultsDF[i, "Gauss1side1dilutionBootCIlo"]<-diff(pnorm(mean=0, sd=sigmaQuantile[1], q=c(0.5,1.5)))
      ResultsDF[i, "Gauss1side1dilutionBootCIhi"]<-diff(pnorm(mean=0, sd=sigmaQuantile[2], q=c(0.5,1.5)))
      print("Bootstrap distribution of 1-sided single error (Gaussian model):")
      print(ResultsDF[i, c("Gauss1side1dilutionBootCIlo", "Gauss1side1dilutionBootCIhi")])
      # print(quantile((1-pBootVector)/2, probs=c(0.025, 0.5, 0.975)))
      ResultsDF[i, "GaussNdiscreteBootValues"]<-length(table(pBootVector))
      print(unlist(list("number of distinct values"=ResultsDF[i, "GaussNdiscreteBootValues"])))
   }

   ProbFitGauss<-ProbDisagreeGaussian(sigma=ResultsDF[i, "sigma"])[1:I(max(ResultsList[[i]]["difference"])+1)]
   # print(ProbFitGauss)
   # print(unlist(list(length=length(ResultsList[[i]]))))
   # ExpectedNGauss<-ProbFitGauss*length(ResultsList[[i]])
   ExpectedNGauss<-ProbFitGauss*nrow(ResultsList[[i]])
   print("ExpectedNGauss:")
   print( ExpectedNGauss)

   ProbFitLaplace  <-ProbDisagreeLaplace(b=ResultsDF[i, "b"])[1:I(max(ResultsList[[i]]["difference"])+1)]
   ExpectedNLaplace<-ProbFitLaplace*nrow(ResultsList[[i]])
   print("ExpectedNLaplace:")
   print( ExpectedNLaplace)

      
   # ExpectedNGaussTruncated<-TruncateTableFn(ExpectedNGauss)
   # forceScalar<-any(as.logical(ResultsList[[i]]$censored))
   # print(unlist(list(forceScalar=forceScalar)))
   ExpectedNGaussTruncated<-TruncateTableFn(ExpectedNGauss, 
      force=forceScalar, 
      forceLength=2, verbose=F)

   TableResultsTruncated<-TruncateTableFn(
      tableResults, force=T, forceLength=length(ExpectedNGaussTruncated))
   print("TableResultsTruncated:")
   print( TableResultsTruncated)
   print("ExpectedNGaussTruncated:")
   print( ExpectedNGaussTruncated)
   print(rbind(as.vector(TableResultsTruncated), ExpectedNGaussTruncated))

   # ChisqGauss<-
   ResultsDF[i, "Gausschisq"]  <-ChisqContingency(O=TableResultsTruncated, E=ExpectedNGaussTruncated)
   # dfGauss<-length(TableResultsTruncated)-1
   ResultsDF[i, "Gausschisqdf"]<-length(TableResultsTruncated)-1
   ResultsDF[i, "Gausschisqp"] <-1-pchisq(ResultsDF[i, "Gausschisq"], ResultsDF[i, "Gausschisqdf"])
   # print(unlist(list(chisq=ChisqGauss, df=dfGauss, p=1-pchisq(ChisqGauss, dfGauss))))
   print(unlist(list(chisq=ResultsDF[i, "Gausschisq"], df=ResultsDF[i, "Gausschisqdf"], p=ResultsDF[i, "Gausschisqp"])))

}

format(ResultsDF[, c("MLq", "ML1side1dilutionError", "MLqchisq", "MLqchisqdf", "MLqchisqp")], scientific=F)

summary(ResultsDF$b)

par(mfrow=c(1,1))

write.csv(ResultsDF, file="../figures/Results.csv", row.names = FALSE)

# omit heteroscedasticity results
dim(ResultsDF)
ResultsDF<-ResultsDF[!ResultsDF$dataset %in% c("AbsxCauchemezLoDF", "AbsxCauchemezHiDF","ThrusLoMeanLongDF","ThrusHiMeanLongDF"),]
dim(ResultsDF)


#  View(ResultsDF)
# names(ResultsDF)
# min(ResultsDF$MMqchisqp, na.rm=T)
min(ResultsDF$MLqchisqp, na.rm=T)
min(ResultsDF$Gausschisqp, na.rm=T)
max(abs(ResultsDF$MMq-ResultsDF$MLq), na.rm=T)
rbind( ResultsDF$MM1side1dilutionError,ResultsDF$ML1side1dilutionError,
   abs(ResultsDF$MM1side1dilutionError-ResultsDF$ML1side1dilutionError))
SubsetVector<-as.logical(c(0, 0, 0, 1, 1, 1, 0, 1))
max(abs(ResultsDF$MM1side1dilutionError[SubsetVector]-ResultsDF$ML1side1dilutionError[SubsetVector]), na.rm=T)
# plot(x=ResultsDF$variance, y=(1-ResultsDF$MMq)**2, type="n")
FittedVariance<-(1-ResultsDF$MMq)*(2-((1+ResultsDF$MMq)*(1-(ResultsDF$MMq**2))))

# the "variance" variable doesn't exist (but there is "VarianceAll")
# plot(x=ResultsDF$variance[!is.na(ResultsDF$MMq)], 
#      y=FittedVariance[!is.na(ResultsDF$MMq)], 
#      type="n")
# text(x=ResultsDF$variance, y=FittedVariance, labels=letters[1:8])
# abline(a=0, b=1, lty=2)

#  make the plots for Gaussian as done previously for geometric in "error.png"
SigmaLength<-400
sigmaVector<-seq(0.001, 2, length=SigmaLength)
    bVector<-seq(0.001, 1, length=SigmaLength)
# ProbDisagreeGaussian(sigmaVector)

ErrorDF<-AgreeDF<-data.frame(
   sigma            =rep(NA, SigmaLength),
   p0               =rep(NA, SigmaLength),
   p1               =rep(NA, SigmaLength),
   p2               =rep(NA, SigmaLength),
   p3plus           =rep(NA, SigmaLength),
   p0ThreeHalves    =rep(NA, SigmaLength),
   p1ThreeHalves    =rep(NA, SigmaLength),
   p2ThreeHalves    =rep(NA, SigmaLength),
   p3ThreeHalves    =rep(NA, SigmaLength),
   p4plusThreeHalves=rep(NA, SigmaLength),
   b                =rep(NA, SigmaLength),
   p0Laplace        =rep(NA, SigmaLength),
   p1Laplace        =rep(NA, SigmaLength),
   p2Laplace        =rep(NA, SigmaLength),
   p3plusLaplace    =rep(NA, SigmaLength)
)

for(i in 1:length(sigmaVector)){
   AgreeDF[i,c("sigma", "p0", "p1", "p2")]<-c(sigmaVector[i], ProbDisagreeGaussian(sigmaVector[i])[c("0", "1", "2")])
   ErrorDF[i,c("sigma", "p0", "p1", "p2")]<-c(sigmaVector[i], 
        diff(pnorm(mean=0, sd=sigmaVector[i], q=c(-0.5,0.5))),
      2*diff(pnorm(mean=0, sd=sigmaVector[i], q=c( 0.5,1.5))),
      2*diff(pnorm(mean=0, sd=sigmaVector[i], q=c( 1.5,2.5))) )
   ErrorDF[i,c("p0ThreeHalves", "p1ThreeHalves", "p2ThreeHalves", "p3ThreeHalves")]<-c(
        diff(pnorm(mean=0, sd=sigmaVector[i], q=c( -log(1.5, base=2)/2,   log(1.5, base=2)/2))),
      2*diff(pnorm(mean=0, sd=sigmaVector[i], q=c(  log(1.5, base=2)/2, 3*log(1.5, base=2)/2))),
      2*diff(pnorm(mean=0, sd=sigmaVector[i], q=c(3*log(1.5, base=2)/2, 5*log(1.5, base=2)/2))),
      2*diff(pnorm(mean=0, sd=sigmaVector[i], q=c(5*log(1.5, base=2)/2, 7*log(1.5, base=2)/2))) )

      AgreeDF[i,c("b", "p0Laplace", "p1Laplace", "p2Laplace")]<-
         c(bVector[i], ProbDisagreeLaplace(bVector[i])[c("0", "1", "2")])

}


AgreeDF$p3plus           <-1-(AgreeDF$p0            + AgreeDF$p1            + AgreeDF$p2)
AgreeDF$p4plusThreeHalves<-1-(AgreeDF$p0ThreeHalves + AgreeDF$p1ThreeHalves + AgreeDF$p2ThreeHalves + AgreeDF$p3ThreeHalves)
AgreeDF$p3plusLaplace    <-1-(AgreeDF$p0Laplace     + AgreeDF$p1Laplace     + AgreeDF$p2Laplace)
# View(AgreeDF)
# View(ErrorDF)

# where does the one dilution error reach its maximum
plot(x=ErrorDF$sigma, y=(ErrorDF$p1)/2)
ErrorDF[which.max(ErrorDF$p1),]

png("../figures/errorGaussian.png", width = 480*16/9, height = 480)
   par(mfrow=c(1,2))
   # par(mgp=c(2.2, 1, 0))
   par(mgp=c(2, 0.75, 0))
   # par(oma=c(0, 1, 0, 1))

   LegendHeadDF<-data.frame(
      yAgree=c(0.10),
      yError=c(0.07),
      textAgree=c("# of 2-fold\ndilutions\nby which\ndisagree"),
      textError=c("error in # of\n2-fold dilutions\n(solid lines) or\n1.5-fold dilutions\n(dashed lines)")
      )
   # LegendHeadDF

   LegendDF<-data.frame(
      yAgree   =c(0.92,       0.70, 0.25),
      yError   =c(0.98,       0.86, 0.36),
      textAgree=c("\u22642", "\u22641", "0"),
      textError=c("\u22642", "\u22641", "0"))
      # textAgree=c("disagree\nby", "\u22642", "\u22641", "0\n(agree)"),
      # textError=c("error of"    , "\u22642", "\u22641", "0\n(correct)"))
   # LegendXScalar<-0.89

   LegendThreeHalvesDF<-data.frame(
      yError   =c(0.98,       0.86, 0.36, 0.20),
      textError=c("\u22643", "\u22642", "\u22641", "0"))

   LegendXAgreeScalar<-LegendXErrorScalar<-1.09
   # LegendXErrorScalar           <-1.06 
   # LegendXErrorThreeHalvesScalar<-1.08
   LegendXErrorJitterScalar<-0.03

   # ymax<-1.03
   # xmax<-1.11
   # ymax<-1.20
   ymax<- 1.03
   ymin<--0.04
   xmax<-1.20

   sigmaMin<-0.2
   sigmaMax<-1
   sigmaSubsetVector<-AgreeDF$sigma>sigmaMin & AgreeDF$sigma<=sigmaMax

   plot(x=AgreeDF$sigma[sigmaSubsetVector], y=AgreeDF$p0[sigmaSubsetVector], type="l", 
        ylim=c(ymin,ymax), xlim=c(sigmaMin, xmax),
        xlab=expression(paste("standard deviation (", sigma, ") of discretized Gaussian model")),
          ylab="cumulative probability of observed difference in between-replicate dilutions")
   lines(x=AgreeDF$sigma[sigmaSubsetVector], y=AgreeDF$p0[sigmaSubsetVector]+AgreeDF$p1[sigmaSubsetVector])
   lines(x=AgreeDF$sigma[sigmaSubsetVector], y=AgreeDF$p0[sigmaSubsetVector]+AgreeDF$p1[sigmaSubsetVector]+AgreeDF$p2[sigmaSubsetVector])
   text(x=LegendXAgreeScalar, y=LegendHeadDF$yAgree, labels=LegendHeadDF$textAgree, adj=0.5)

   y0scalar<-min(AgreeDF$p0[sigmaSubsetVector])
   y1scalar<-min(AgreeDF$p0[sigmaSubsetVector]+AgreeDF$p1[sigmaSubsetVector])
   y2scalar<-min(AgreeDF$p0[sigmaSubsetVector]+AgreeDF$p1[sigmaSubsetVector]+AgreeDF$p2[sigmaSubsetVector])

   # text(x=LegendXAgreeScalar, y=    LegendDF$yAgree, labels=    LegendDF$textAgree, adj=0.5)
   text(x=LegendXAgreeScalar, y=c(y2scalar, y1scalar, y0scalar), labels=    LegendDF$textAgree, adj=0.5)

   print("plotting data values for each dataset")
   for(i in 1:length(ResultsList)){
      print("")
      print(names(ResultsList)[i])
      letterScalar<-letters[i]
      print(letterScalar)
   
      xScalar<-ResultsDF[i, "sigma"]
      bScalar<-ResultsDF[i, "b"]
      
      # pick out closest value of sigma 
      print(head(sigmaSubsetVector))
      sigmaClosest<-AgreeDF[which.min(abs(AgreeDF$sigma-xScalar)),"sigma"]
          bClosest<-AgreeDF[which.min(abs(AgreeDF$b    -bScalar)),"b"]
      
      print(unlist(list(sigma=xScalar, sigmaClosest=sigmaClosest, bClosest=bClosest)))
      
      if(any(as.logical(ResultsList[[i]]$censored))){
         # there are censored values
         # plot only the proportion of agreement
         CumulativeProportion<-mean(as.numeric(ResultsList[[i]]$difference==0))
         
      if(i<=nrow(ResultsDF)){
         ResultsDF[i,"KSstat"]<-abs(
            AgreeDF[AgreeDF$sigma==sigmaClosest, "p0"]
            -CumulativeProportion)
         ResultsDF[i,"KSstatLaplace"]<-abs(
            AgreeDF[AgreeDF$b    ==bClosest    , "p0"]
            -CumulativeProportion)
      }
      print(CumulativeProportion)

      }else{
         # no censored values
         CumulativeProportion <-(cumsum(table(ResultsList[[i]]$difference))/length(ResultsList[[i]]$difference))
         # drop last element
         CumulativeProportion<-rev(rev(CumulativeProportion)[2:length(CumulativeProportion)])
         
         Proportion0Fitted<-AgreeDF[AgreeDF$sigma==sigmaClosest,"p0"]
         Proportion1Fitted<-AgreeDF[AgreeDF$sigma==sigmaClosest,"p1"]
         Proportion2Fitted<-AgreeDF[AgreeDF$sigma==sigmaClosest,"p2"]
         Proportion3Fitted<-AgreeDF[AgreeDF$sigma==sigmaClosest,"p3plus"]

         Proportion0FittedLaplace<-AgreeDF[AgreeDF$b==bClosest,"p0Laplace"]
         Proportion1FittedLaplace<-AgreeDF[AgreeDF$b==bClosest,"p1Laplace"]
         Proportion2FittedLaplace<-AgreeDF[AgreeDF$b==bClosest,"p2Laplace"]
         Proportion3FittedLaplace<-AgreeDF[AgreeDF$b==bClosest,"p3plusLaplace"]
         print(unlist(list(
            Proportion0FittedLaplace=Proportion0FittedLaplace,
            Proportion1FittedLaplace=Proportion1FittedLaplace,
            Proportion2FittedLaplace=Proportion2FittedLaplace,
            Proportion3FittedLaplace=Proportion3FittedLaplace
            )))
         if(length(CumulativeProportion)==1){
            CumulativeProportionFitted       <-as.vector(Proportion0Fitted)
            CumulativeProportionFittedLaplace<-as.vector(Proportion0FittedLaplace)
         }else{
         if(length(CumulativeProportion)==2){
            CumulativeProportionFitted<-c(
               Proportion0Fitted, 
               Proportion0Fitted+Proportion1Fitted)
            CumulativeProportionFittedLaplace<-c(
               Proportion0FittedLaplace, 
               Proportion0FittedLaplace+Proportion1FittedLaplace)
         }else{
         if(length(CumulativeProportion)==3){
            CumulativeProportionFitted<-c(
               Proportion0Fitted, 
               Proportion0Fitted+Proportion1Fitted,
               Proportion0Fitted+Proportion1Fitted+Proportion2Fitted)
            CumulativeProportionFittedLaplace<-c(
               Proportion0FittedLaplace, 
               Proportion0FittedLaplace+Proportion1FittedLaplace,
               Proportion0FittedLaplace+Proportion1FittedLaplace+Proportion2FittedLaplace)
         }else{
            CumulativeProportionFitted<-c(
               Proportion0Fitted, 
               Proportion0Fitted+Proportion1Fitted,
               Proportion0Fitted+Proportion1Fitted+Proportion2Fitted,
               Proportion0Fitted+Proportion1Fitted+Proportion2Fitted+Proportion3Fitted)
            CumulativeProportionFittedLaplace<-c(
               Proportion0FittedLaplace, 
               Proportion0FittedLaplace+Proportion1FittedLaplace,
               Proportion0FittedLaplace+Proportion1FittedLaplace+Proportion2FittedLaplace,
               Proportion0FittedLaplace+Proportion1FittedLaplace+Proportion2FittedLaplace+Proportion3FittedLaplace)
         }
         }}
         # print(round(CumulativeProportion, 5))
         # print(round(CumulativeProportionFitted, 5))
         print(round(rbind(
            CumulativeProportion, 
            CumulativeProportionFitted,
            CumulativeProportionFittedLaplace
            ), 5))
         if(i<=nrow(ResultsDF)){
            ResultsDF[i,"KSstat"       ]<-max(abs(CumulativeProportion-CumulativeProportionFitted))
            ResultsDF[i,"KSstatLaplace"]<-max(abs(CumulativeProportion-CumulativeProportionFittedLaplace))
         }
      }
   
      # Jitterx      <-0
      JitterVectorx<-JitterVectory<-rep(0, length(CumulativeProportion))
      # if(letterScalar=="a"){JitterVectory[1:2]<-rep( 0.015, 2)}
      # if(letterScalar=="e"){JitterVectory[1:2]<-rep(-0.015, 2)}
      print(JitterVectorx)
      xPlot<-rep(xScalar, length(CumulativeProportion)) + JitterVectorx
      print(xPlot)
      text(x=xPlot, 
           y=CumulativeProportion + JitterVectory,
           labels=rep(letterScalar, length(CumulativeProportion)))
   }

   #
   # second plot
   #
   sigmaSubsetVector<-AgreeDF$sigma>sigmaMin & AgreeDF$sigma<=sigmaMax
    plot(x=ErrorDF$sigma[sigmaSubsetVector], y=ErrorDF$p0[sigmaSubsetVector], type="l", 
         ylim=c(ymin,ymax), xlim=c(sigmaMin, xmax),
         xlab=expression(paste("standard deviation (", sigma, ") of discretized Gaussian model")),
         ylab="cumulative probability of dilutions different from unobserved target value")
   lines(x=ErrorDF$sigma[sigmaSubsetVector], y=ErrorDF$p0[sigmaSubsetVector]+ErrorDF$p1[sigmaSubsetVector])
   lines(x=ErrorDF$sigma[sigmaSubsetVector], y=ErrorDF$p0[sigmaSubsetVector]+ErrorDF$p1[sigmaSubsetVector]+ErrorDF$p2[sigmaSubsetVector])
   # text(x=LegendXScalar, y=LegendDF$yError, labels=LegendDF$textError, adj=0.5)
   text(x=LegendXErrorScalar, 
        y=LegendHeadDF$yError, labels=LegendHeadDF$textError, adj=0.5)

   # for 2-fold dilutions
   y0scalar<-min(ErrorDF$p0[sigmaSubsetVector])
   y1scalar<-min(ErrorDF$p0[sigmaSubsetVector]+ErrorDF$p1[sigmaSubsetVector])
   y2scalar<-min(ErrorDF$p0[sigmaSubsetVector]+ErrorDF$p1[sigmaSubsetVector]+ErrorDF$p2[sigmaSubsetVector])
   text(x=LegendXErrorScalar-LegendXErrorJitterScalar, 
        y=c(y2scalar, y1scalar, y0scalar), 
        labels=    LegendDF$textError, adj=0.5)

   # for 1.5 fold dilutions
   y0scalar<-min(ErrorDF$p0ThreeHalves[sigmaSubsetVector])
   y1scalar<-min(ErrorDF$p0ThreeHalves[sigmaSubsetVector]+ErrorDF$p1ThreeHalves[sigmaSubsetVector])
   y2scalar<-min(ErrorDF$p0ThreeHalves[sigmaSubsetVector]+ErrorDF$p1ThreeHalves[sigmaSubsetVector]+ErrorDF$p2ThreeHalves[sigmaSubsetVector])
   y3scalar<-min(ErrorDF$p0ThreeHalves[sigmaSubsetVector]+ErrorDF$p1ThreeHalves[sigmaSubsetVector]+ErrorDF$p2ThreeHalves[sigmaSubsetVector]+ErrorDF$p3ThreeHalves[sigmaSubsetVector])
   text(x=LegendXErrorScalar+LegendXErrorJitterScalar, 
        y=c(y3scalar, y2scalar, y1scalar, y0scalar), 
        labels=    LegendThreeHalvesDF$textError, adj=0.5)

   # lines for 1.5 fold dilutions
   lines(
      x=ErrorDF$sigma[sigmaSubsetVector], 
      y=ErrorDF$p0ThreeHalves[sigmaSubsetVector], lty=4)
   lines(
      x=ErrorDF$sigma[sigmaSubsetVector], 
      y=ErrorDF$p0ThreeHalves[sigmaSubsetVector]+ErrorDF$p1ThreeHalves[sigmaSubsetVector], lty=4)
   lines(
      x=ErrorDF$sigma[sigmaSubsetVector], 
      y=ErrorDF$p0ThreeHalves[sigmaSubsetVector]+ErrorDF$p1ThreeHalves[sigmaSubsetVector]+ErrorDF$p2ThreeHalves[sigmaSubsetVector], lty=4)
   lines(
      x=ErrorDF$sigma[sigmaSubsetVector], 
      y=ErrorDF$p0ThreeHalves[sigmaSubsetVector]+ErrorDF$p1ThreeHalves[sigmaSubsetVector]+ErrorDF$p2ThreeHalves[sigmaSubsetVector]+ErrorDF$p3ThreeHalves[sigmaSubsetVector], lty=4)

   # # at.y<-seq(0, 1, 0.2)
   # # axis(4, at=at.y,  labels=as.character(1-at.y))
   # # mtext("specificity of using one more than this number as a threshold", side = 4, line=2)
   # 

   par(mfrow=c(1,1))
   par(mgp=c(3, 1, 0))
   # par(oma=c(0, 0, 0, 0))
   
dev.off()

cbind(
   ResultsDF[,c("dataset", "KSstat", "KSstatLaplace")], 
   ResultsDF$KSstat - ResultsDF$KSstatLaplace)
# View(ResultsDF)
# View(ErrorDF)

png("../figures/errorInterval.png", width = 480*16/9, height = 480)

xJitter<-0.25

plot(
   x=c(0.5, 10.5), 
   y=c(0, max(
      c(ResultsDF$MM1side1dilutionDeltaCIhi),
      c(ResultsDF$ML1side1dilutionBootCIhi),
      c(ResultsDF$Gauss1side1dilutionBootCIhi), na.rm=T
      )),
   type="n",
   xlab="", ylab="probability of a 1-sided 1-dilution error",
   xaxt="n"
)
mtext("dataset", side=1, at=mean(1:8)-0.025, line=2.2)

axis(1, at=1:8, labels=letters[1:8])

pchVector<-c(1, 2, 4)

UncensoredSubsetVector<-as.logical(c(0, 0, 1, 1, 1, 1, 1, 1))
points(x=(1:8)[UncensoredSubsetVector]-xJitter, 
       y=ResultsDF$MM1side1dilutionError[UncensoredSubsetVector], pch=pchVector[1])
segments(x0=(1:8)[UncensoredSubsetVector]-xJitter, 
         y0=ResultsDF$MM1side1dilutionDeltaCIlo[UncensoredSubsetVector], 
         y1=ResultsDF$MM1side1dilutionDeltaCIhi[UncensoredSubsetVector])

points(x=(1:8)-0, y=ResultsDF$ML1side1dilutionError, pch=pchVector[2])
segments(x0=(1:8)-0, 
         y0=ResultsDF$ML1side1dilutionBootCIlo, 
         y1=ResultsDF$ML1side1dilutionBootCIhi)

points(x=(1:8)+xJitter, y=ResultsDF$Gauss1side1dilutionError, pch=pchVector[3])
segments(x0=(1:8)+xJitter, 
         y0=ResultsDF$Gauss1side1dilutionBootCIlo, 
         y1=ResultsDF$Gauss1side1dilutionBootCIhi)

legend("topright", pch=pchVector, 
       legend=c(
          "Single-error model,\nmethod of moments,\ndelta method.\n", 
          "Single-error model,\nmaximum likelihood,\nbootstrap.\n",
          "Gaussian model,\nmaximum likelihood,\nbootstrap.\n"
          ), bty="n", cex=0.9)

dev.off()
