source('hmm.R')
library(mvtnorm)

MNormHmm <- setRefClass("MNormHmm",
  contains = "Hmm",
  fields = list(
  ),
  methods = list(
    initialize = function(m, X, A, priors, pdf.params, do.plot=FALSE,
                          nstart=10, iter.max=20, ...) {
      if (missing(A) && missing(priors) && missing(pdf.params)) {
        callSuper(m, X, ...)
        mean.init <- list()
        sigma.init <- list()
        km <- kmeans(X, m, iter.max=iter.max, nstart=nstart)
        if (ncol(X) > 2) {
          do.plot <- FALSE
        }
        if (do.plot) {
          xlim <- c(min(X[,1]), max(X[,1]))
          ylim <- c(min(X[,2]), max(X[,2]))
          axes <- TRUE
        }
        for (i in 1:m) {
          X.i <- X[km$cluster==i,]
          mean.init[[i]] <- colMeans(X.i)
          sigma.init[[i]] <- cov(X.i)
          if (do.plot) {
            plot(X.i, xlim=xlim, ylim=ylim, col=i+1, axes=axes)
            points(km$centers[i,1], km$centers[i,2], col=i+1, pch=9, cex=5)
            par(new=TRUE)
            axes <- FALSE
          }
        }
        par(new=FALSE)
        pdf.params <<- list(
          mean = mean.init,
          sigma = sigma.init
        )
      } else if (missing(m) && missing(X)) {
        callSuper(A=A, priors=priors)
        pdf.params <<- pdf.params
      } else {
        stop("Either m, X or A, priors and pdf.params should be specified as arguments.")
      }
      pdf <<- dmvnorm
      rng <<- rmvnorm
    },
    getD = function() {
      length(pdf.params$mean[[1]])
    },
    getDf = function() {
      m <- getM()
      d <- getD()
      # row-stochastic A + priors + vector of means + symmetric Sigma
      m*(m - 1) + (m - 1) + (d + d*(d + 1)/2)*m
    },
    toUpperTri = function(v) {
      k <- length(v)
      d <- (-1 + sqrt(1 + 8*k)) / 2
      M <- matrix(0, d, d)
      for (i in 0:(d-1)) {
        s <- i*(i+1)/2
        M[1:(i+1),(i+1)] <- v[(s+1):(s+i+1)]
      }
      M
    },
    getWrkParams = function() {
      d <- getD()
      m <- getM()
      # sigma.wrk <- numeric(m*d*(d + 1)/2)
      sigma.wrk <- numeric()
      for (i in 1:m) {
        L <- chol(pdf.params$sigma[[i]])
        diag.idcs <- as.logical(diag(d))
        L[diag.idcs] <- log(L[diag.idcs])
        sigma.wrk <- c(sigma.wrk, L[L != 0])
      }
      as.numeric(
        c(unlist(pdf.params$mean),
          sigma.wrk,
          callSuper())
      )
    },
    setWrkParams = function(params) {
      # mean, sigma (log-Cholesky transformed), A (by columns excluding diag), priors:
      #               m*d + m*d*(d + 1)/2 + m*(m - 1) + (m - 1)
      m <- getM()
      d <- getD()
      callSuper(params, m*d + m*d*(d + 1)/2)
      mean.nat <- list()
      sigma.nat <- list()
      diag.idcs <- as.logical(diag(d))
      for (i in 1:m) {
        mean.nat[[i]] <- params[((i - 1)*d + 1):(i*d)]
        L <- params[(m*d + (i-1)*d*(d + 1)/2 + 1):(m*d + i*d*(d + 1)/2)]
        L <- toUpperTri(L)
        L[diag.idcs] <- exp(L[diag.idcs])
        sigma.nat[[i]] <- t(L) %*% L
      }
      pdf.params <<- list(mean = mean.nat,
                          sigma = sigma.nat)
    }
  )
)
