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
