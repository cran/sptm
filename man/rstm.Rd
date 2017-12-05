\name{rstm}
\alias{rstm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate failure time from a semiparametric transformation model
}
\description{
Simulate failure time from a semiparametric transformation model
}
\usage{
rstm(n, family = c("PH", "PO", "P2"), linear.predictors, baseline.hazard = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
integer. Sample size
}
  \item{family}{
string. 
}
  \item{linear.predictors}{
vector. It can also be a matrix of 1 column, the dimension will be dropped
}
  \item{baseline.hazard}{
numeric.
}
}
\details{
Called by sim.fong
}
\examples{

n=100
beta= c(log(.5), log(.7), log(1.2)) 
t0=2.9999
init = c(log(0.0373*t0),beta)        
ft=rstm (n, family="PH", runif(n,1,2), baseline.hazard=0.032) 


}
