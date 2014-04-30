## ======================================================================
## Copyright 2009--2014, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================
## F_distribution.R
## R functions involving the noncentral F distribution
##
## Supported by National Science Foundation grants DMS-0604963 and
## DMS-0906864.
##
## Reference:
## L. Wei, P. F. Craigmile, and W. M. King (2012). Spectral-based
## noncentral F mixed effect models, with application to otoacoustic
## emissions. Journal of Time Series Analysis, 33, 850-862.
## ======================================================================

F.mean <- function (df1, df2, ncp) {
  ## Calculate the mean of the F distribution with degrees of
  ## freedom, df1 and df2, and noncentrality parameter 'ncp'.

  ##  b <- df2 / (df2 - 2)
  ##  a <- b / df1
  ##  a * ncp + b

  df2 * (ncp + df1) / (df1 * (df2 - 2))
}


F.var.from.mean <- function (df1, df2, mean) {
  ## Calculate the variance of the F distribution with degrees of
  ## freedom, df1 and df2, at a given 'mean' value.

  rs <- F.roots(df1, df2)
  2 * (mean - rs[1]) * (mean - rs[2]) / (df2 - 4)
}


F.ncp.from.mean <- function (df1, df2, mean) {
  ## Calculate the noncentrality parameter of the F distribution
  ## with degrees of freedom, df1 and df2, at a given 'mean' value.
  
  ##  b <- df2 / (df2 - 2)
  ##  a <- b / df1
  ##  (mean - b) / a

   df1 * (df2-2) * mean / df2 - df1
}



F.c.values <- function (df1, df2) {
  ## Calculate the coefficients of the polynomial characterizing the variance
  ## of a noncentral F distribution.

  df2.minus.4 <- df2 - 4

  c0 <- -2 * df2^2 / (df1 * (df2-2) * df2.minus.4)
  c1 <- 4 * df2 / (df1 * df2.minus.4)
  c2 <- 2 / df2.minus.4

  c(c0, c1, c2)
}


F.roots <- function (df1, df2, use.polyroot=F) {
  ## Calculate the roots of the polynomial characterizing the variance
  ## of a noncentral F distribution.

  if (use.polyroot) {
    sort(Re(polyroot(F.c.values(df1, df2))))
  }
  else {
    pm <- sqrt( (df1+df2-2) / (df2-2) )
    (df2/df1) * c(-1 - pm, -1 + pm)
  }
}



F.sim <- function (X, beta, df1, df2) {
  ## simulate from an non-central F distribution on df1 and df2
  ## degrees of freedom with noncentrality parameter exp(X %*% beta).
   
  ncp <- drop(exp(X %*% beta))

  rf(nrow(X), df1, df2, ncp)
}


F.sim.mean <- function (X, beta, df1, df2) {
  ## simulate from an non-central F distribution on df1 and df2
  ## degrees of freedom with mean exp(X %*% beta).
   
  mu  <- drop(exp(X %*% beta))
  ncp <- F.ncp.from.mean(df1, df2, mu)

  rf(nrow(X), df1, df2, ncp)
}



F.ml.beta <- function (y, X, df1, df2, beta0=rep(0, ncol(X)), extras=FALSE) {
  ## Calculate the maximum likelihood estimation of beta for y
  ## assuming the data follows a noncentral F distribution on df1 and
  ## df2 degrees of freedom with noncentrality parameter X %*% beta.
  ## Start the optimizer at 'beta0'.
  ## Include extra quantities in the output if 'extras=TRUE'.

  F.minus.ll.beta <- function (beta, df1, df2, X=X, y=y) {

    -sum(df(y, df1, df2, exp(drop(X %*% drop(beta))), log=TRUE))
  }

  opt <- optim(beta0, F.minus.ll.beta, df1=df1, df2=df2, X=X, y=y)
  beta.hat <- opt$par

  if (!extras) {
    
    list(beta.hat = beta.hat)
  } else {
    
    eta.hat <- drop(X %*% beta.hat)
    mu.hat  <- F.mean(df1, df2, exp(eta.hat))
    var.hat <- F.var.from.mean(df1, df2, mu.hat)
    
    list(beta.hat = beta.hat,
         eta.hat  = eta.hat,
         mu.hat   = mu.hat,
         var.hat  = var.hat)
  }
}






## Next two functions were used for quasi-likelihood estimation

F.link.constraint <- function (df2) {

  log(df2/(df2-2))
}


F.check.Xbeta.constraint <- function (X, beta, df2) {
  ## check whether the constraint holds for beta (TRUE for yes!)

  all(X %*% beta >= F.link.constraint(df2))
}


