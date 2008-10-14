\name{stm}
\alias{stm}
\alias{getFixedEf.stm}
\title{Fit a semiparametric transformation model}
\description{Fit a semiparametric transformation model}

\usage{
stm (formula, dat, strata.formula, phase2.ind=NULL, imputation.formula=NULL, 
    family=c("PH","PO","P2"), ee=c("fine2","fine1","kong"), var.est.type=c("1","2"),
    t0, init=NULL, maxit=1000,
    intermediate=FALSE, verbose=FALSE, show.time.elapsed=TRUE)


\method{getFixedEf}{stm} (object, ...)


}

\arguments{
  \item{formula}{formula. Regression model of interest}
  \item{dat}{data frame.}
  \item{strata.formula}{formula.}
  \item{phase2.ind}{Boolean vector. If TRUE, phase II samples; if FALSE, phase I samples. If NULL, will try to infer from which subjects have phase II variables. Should not be 0/1}
  \item{imputation.formula}{formula. If not NULL, calibration weighting is done}

  \item{family}{string.}
  \item{ee}{string. Type of design matrix used in estimating equation}
  \item{var.est.type}{string. 1: one-stage estimator, 2: two-stage estimator}

  \item{t0}{numeric. Should be close to the end of study time}
  \item{init}{numerical vector.}
  \item{maxit}{integer. Maximum number of iterations in the optimization process}

  \item{intermediate}{Boolean.}
  \item{verbose}{Boolean.}
  \item{show.time.elapsed}{Boolean.}

  \item{object}{an object of type stm}
  \item{...}{additional arguments}

}
\details{
Fit stm both with and without calibration. Calls stm.internal.
}
\value{
An object of type stm
}
\examples{

n=100
beta= c(log(.5), log(.7), log(1.2)) 
t0=2.9999
init = c(log(0.0373*t0),beta)        
dat = sim.fong(n, family="PH", beta, random.censoring="0", design="CC", auxiliary="weak", seed=1)

est = stm(formula=Surv(X,d) ~ z + s + z:s, dat, strata.formula=~d, family="PH", t0=t0, init=init,
    var.est.type="1", verbose=3)

}
