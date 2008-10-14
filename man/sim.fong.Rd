\name{sim.fong}
\alias{sim.fong}
\title{
Data Simulation as in Fong and Gilbert (2014)
}
\description{
Simulate data as in Fong and Gilbert (2014).
}
\usage{
sim.fong (n, family=c("PH","PO","P2"), beta, 
    random.censoring=c("0\%","20\%","60\%"), prevalence=0.1, non.adherence.ratio=0,
    design=c("FULL","CC"), auxiliary=c("weak","good","excellent","none"), 
    seed=NULL, var.S=1, var.W=1) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
integer. Sample size
}
  \item{family}{
string. Link functions in the semiparametric transformation model
}
  \item{beta}{
numerical vector. Coefficients of the linear model
}
  \item{random.censoring}{
string. Random censoring in addition to administrative censoring
}
  \item{prevalence}{
numerical. Proportion of cases among z==0 when there is no random censoring and non-adherence ratio is 0
}
  \item{design}{
string. Full cohort or case-cohort (finite population sampling)
}
  \item{auxiliary}{
string. }
  \item{seed}{
integer. Random generator seed}
  \item{var.S}{
numeric. Variance of the phase II covariate s}
  \item{var.W}{
numeric. Variance of the baseline covariate w
}
  \item{non.adherence.ratio}{ratio of non-adherent}
}
\details{
The number of rows is the size of the full cohort.
Adherence ratio works as a Bernoulli variable.
Prevalence is used to compute baseline hazard function based on some empirical evidence.
}
\value{
If design is FULL, returns a data frame of:
    \item{ft}{failure time}
    \item{C}{censoring time}
    \item{X}{smaller of the ft and C}
    \item{d}{event indicator}
    \item{z}{baseline covariate z}
    \item{s}{phase II covariate s}

If design is CC, returns a data frame of:
    \item{ft}{failure time}
    \item{C}{censoring time}
    \item{X}{smaller of the ft and C}
    \item{d}{event indicator}
    \item{z}{baseline covariate z}
    \item{s}{phase II covariate s}
    \item{w}{baseline auxiliary covariate w}
}
\examples{

dat = sim.fong(n=10000, family="PH", beta=c(log(.5), log(.7), log(1.2)), design="CC", 
    auxiliary="weak", seed=1, prevalence=0.1, non.adherence.ratio=0, random.censoring="0")
mean(dat$d[dat$z==0])

dat = sim.fong(n=10000, family="PH", beta=c(log(.5), log(.7), log(1.2)), design="CC", 
    auxiliary="weak", seed=1, prevalence=0.1, non.adherence.ratio=0.15, random.censoring="0")
sum(dat$d & !is.na(dat$s))
sum(!dat$d & !is.na(dat$s)) / sum(dat$d & !is.na(dat$s))

dat = sim.fong(n=10000, family="PH", beta=c(log(.5), log(.7), log(1.2)), design="CC", 
    auxiliary="weak", seed=1, prevalence=0.1, non.adherence.ratio=0.15, random.censoring="20")
sum(dat$d & !is.na(dat$s))
sum(!dat$d & !is.na(dat$s)) / sum(dat$d & !is.na(dat$s))

}
