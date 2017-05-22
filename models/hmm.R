Hmm <- setRefClass("Hmm",
  fields = list(
    A = "matrix",
    priors = "numeric",
    pdf = "function",
    pdf.params = "list",
    rng = "function",
    n.iter = "numeric"
  ),
  methods = list(
    init.from.data = function(m, r = 10, A.init = 'main.diag') {
      x <- r / (r - 1 + m)
      if (A.init == 'main.diag') {
        A <<- diag(x, m, m)
      } else if (A.init == 'sec.diag') {
        if (m == 1) {
          # Otherwise `apply` simplifies to vector
          A <<- matrix(1)
        } else {
          A <<- apply(diag(x, m, m), 2, rev)
        }
      } else {
        stop('Unknown initialization method.')
      }
      A[!A] <<- x / r
      priors <<- rep(1, m) / m + runif(m, -1/r*1/m, 1/r*1/m)
      priors <<- priors / sum(priors)
    },
    init.from.params = function(A, priors) {
      A <<- A
      priors <<- priors
    },
    initialize = function(m, X, A, priors, ...) {
      if (missing(A) && missing(priors)) {
        init.from.data(m, ...)
      } else if (missing(m) && missing(X)) {
        init.from.params(A, priors)
      } else {
        stop("Either m, X or A, priors should be given as arguments.")
      }
    },
    show = function() {
      cat("A:\n")
      methods::show(round(A, 2))
      cat("\npriors:\n")
      methods::show(round(priors, 2))
      cat("\npdf.params:\n")
      methods::show(rapply(pdf.params, function(x) round(x, 2), how='list'))
    },
    getM = function() {
      nrow(A)
    },
    getDf = function() {
      m <- getM()
      # row-stochastic A + priors + pdf.params
      m*(m-1) + (m-1) + length(pdf.params) * m
    },
    getWrkParams = function() {
      # TODO: stationary case
      log(c((A / diag(A))[!diag(nrow(A))],
            priors[-1] / priors[1]))

    },
    setWrkParams = function(params, A.offset) {
      # TODO: stationary case
      m <- getM()
      A <<- diag(m)
      A[!A] <<-  exp(params[(A.offset + 1):((A.offset + 1) + m * (m - 1) - 1)])
      A <<- A / rowSums(A)
      priors <<- 1
      if (m > 1) {
        priors <<- c(priors, exp(params[(A.offset + m * (m - 1) + 1):length(params)]))
      }
      priors <<- priors / sum(priors)
    },
    evalPdfScalar = function(x) {
      do.call(pdf, c(x, pdf.params))
    },
    evalPdfVect = function(xx) {
      par.zip <- c(do.call(rbind, pdf.params))
      m <- getM()
      out <- numeric(m)
      n.par <- length(pdf.params)
      for (j in 1:m) {
        par.cur <- par.zip[((j-1)*n.par+1):(j*n.par)]
        out[j] <- do.call(pdf, c(list(x=xx), par.cur))
      }
      out
    },
    alphaNorm = function(X, compute.L = FALSE) {
      if (is.vector(X)) {
        T <- length(X)
        alpha <- evalPdfScalar(X[1]) * priors
      } else { # matrix
        T <- nrow(X)
        alpha <- evalPdfVect(X[1,]) * priors
      }
      if (compute.L) {
        L.log <- log(sum(alpha))
      }
      alpha.norm <- alpha / sum(alpha)
      if (T > 1) {
        for (t in 2:T) {
          if (is.vector(X)) {
            bb.t <- evalPdfScalar(X[t])
          } else {
            bb.t <- evalPdfVect(X[t,])
          }
          z <- alpha.norm %*% A * bb.t
          z.sum <- sum(z)
          if (compute.L) {
            L.log <- L.log + log(z.sum)
          }
          alpha.norm <- z / z.sum
        }
      }
      if (compute.L) {
        return(L.log)
      } else {
        return(alpha.norm)
      }
    },
    logL = function(X) {
      alphaNorm(X, compute.L = TRUE)
    },
    mLogL = function(par.wrk, X) {
      setWrkParams(par.wrk)
      mllk <- -logL(X)
      if (!is.finite(mllk)) {
        cat('Oops! mllk is not finite.\n')
        # FIXME: NaNs produced
        # cat('Inf / NA / NaN: ', mllk, '\n')
        # methods::show(.self)
        # mllk <- .Machine$integer.max * sign(mllk)
        mllk <- .Machine$integer.max
      }
      return(mllk)
    },
    fit = function(X, iterlim=1000, method='nlm', ...) {
      if (method == 'nlm') {
        nlm.out <- nlm(.self$mLogL, getWrkParams(), X, iterlim=iterlim, ...)
        setWrkParams(nlm.out$estimate)
        n.iter <<- nlm.out$iterations
        return(nlm.out)
      } else {
        # browser()
        # TODO: upgrade to 3.2.3 to fix the error
        # TODO: try different methods and optimizers
        # TODO: control
        # TODO: proper parnames (for confints)
        # TODO: remove n.iter?
        # TODO: process error codes
        ### HACK
        par.init <- getWrkParams()
        n.par <- length(par.init)
        parnames(mLogL) <- as.character(1:n.par)
        names(par.init) <- as.character(1:n.par)
        ### /HACK
        nlm.out <- mle2(mLogL,
                        start=par.init,
                        data=list(X=X),
                        control=list(maxit=iterlim)) # iterlim!
        # nlm.out <- mle2(.self$mLogL,
        #                 vecpar=TRUE,
        #                 start=list(par.wrk=getWrkParams()),
        #                 data=list(xx=xx),
        #                 optimizer='nlm',
        #                 iterlim=iterlim)
        browse()
        # r <- profile(nlm.out) # can take up to 30 minutes
        # confint(r)
        #     2.5 %     97.5 %
        # 1   2.471940  2.6729458
        # 2   2.700080  3.3890329
        # 3   2.767046  3.7076717
        # 4   3.280958  3.5433704
        # 5  -5.778032  6.4607278
        # 6         NA         NA
        # 7         NA         NA
        # 8         NA         NA
        # 9  -6.803150         NA
        # 10        NA         NA
        # 11        NA -1.6854847
        # 12        NA         NA
        # 13 -6.676508 -0.2812028
        # 14        NA -1.7843541
        # 15 -5.255174 14.4263307
        # 16        NA         NA
        # 17        NA         NA
        # 18        NA         NA
        # 19        NA         NA
        # TODO: what do NA's mean?
        # ci <- confint(r)
        # plot(r, which=rownames(na.exclude(ci)))
        # setWrkParams(nlm.out$estimate)
        # n.iter <<- nlm.out$iterations
        # return(nlm.out)
      }
    },
    aic = function(X) {
      -2 * logL(X) + 2*getDf()
    },
    bic = function(X) {
      -2 * logL(X) + getDf() * log(if (is.vector(X)) length(X) else nrow(X))
    },
    walkChain = function(T) {
      components <- 1:getM()
      state <- numeric(T)
      state[1] <- sample(components, 1, prob=priors)
      for (t in 2:T) {
        state[t] <- sample(components, 1, prob=A[state[t-1],])
      }
      return(state)
    },
    genSample = function(T) {
      state <- walkChain(T)
      do.call(rng, c(T, rapply(pdf.params,
                               function(params) params[state],
                               how='list')))
    },
    forecastDist = function(X, x.probes, h.max) {
      if (is.vector(x.probes)) {
        n <- length(x.probes)
        # B.probe <- t(outer(x.probes, pdf.params$lambda, pdf))  # 4 times faster
        B.probe <- sapply(x.probes, evalPdfScalar)
      } else {
        n <- nrow(x.probes)
        B.probe <- apply(x.probes, 1, evalPdfVect)
      }
      F <- matrix(0, h.max, n)
      cc <- alphaNorm(X, compute.L = FALSE)
      for (h in 1:h.max) {
        cc <- cc %*% A
        F[h,] <- cc %*% B.probe
      }
      # return(F)
      # FIXME: is normalization necessary? Works fine w/o it in the case of PoisHmm
      return(F / rowSums(F))
    },
    forecast = function(X, x.probes, h.max, method='mode') {
      T <- if (is.vector(X)) length(X) else nrow(X)
      F <- forecastDist(X, x.probes, h.max)
      fc.idcs <- numeric(h.max)
      if (method == 'mode') {
        fc.idcs <- apply(F, 1, which.max)
      } else if (method == 'mean') {
        if (is.vector(x.probes)) {
          fc.idcs <- as.numeric(F %*% x.probes)
          for (h in 1:length(fc.idcs)) {
            fc.idcs[h] <- which.min(abs(x.probes - fc.idcs[h]))
          }
        } else {
          # FIXME: implement
          stop("Not Implemented")
        }
      } else if (method == 'median') {
        for (h in 1:length(fc.idcs)) {
          fc.idcs[h] <- which.min(abs(cumsum(F[h,]) - 0.5))
        }
      } else {
        stop(paste('Unknown method', method))
      }
      fc <- list()
      fc$'forecast' <- if (is.vector(x.probes)) x.probes[fc.idcs] else
                                                fc.idcs
      conf.levels <- c(0.8, 0.95)
      fc$'level' <- conf.levels * 100
      for (b in c('lower', 'upper')) {
        fc[[b]] <- matrix(0, h.max, length(conf.levels))
        colnames(fc[[b]]) <- fc$level
      }
      for (lvl.idx in 1:length(conf.levels)) {
        conf.level <- conf.levels[lvl.idx]
        ui <- numeric(h.max)
        li <- numeric(h.max)
        for (h in 1:h.max) {
          # print(h)
          pivot <- fc.idcs[h]
          accum <- F[h, pivot]
          ui[h] <- pivot
          li[h] <- pivot
          while ((li[h] > 1) || (ui[h] < length(F[h,]))) {
            if (li[h] > 1)  {
              li[h] <- li[h] - 1
              if ((F[h, li[h]] + accum) > conf.level) {
                break
              } else {
                accum <- accum + F[h, li[h]]
              }
            }
            if (ui[h] < length(F[h,])) {
              ui[h] <- ui[h] + 1
              if ((F[h, ui[h]] + accum) > conf.level) {
                break
              }
              accum <- accum + F[h, ui[h]]
            }
            # print(c(li[h], ui[h], accum))
          }
          # print(paste('lower', h, conf.level, '<-', x.probes[li[h]]))
          # print(paste('upper', h, conf.level, '<-', x.probes[ui[h]]))
          fc$lower[h, lvl.idx] <- if (is.vector(x.probes)) x.probes[li[h]] else
                                                           li[h]
          fc$upper[h, lvl.idx] <- if (is.vector(x.probes)) x.probes[ui[h]] else
                                                           ui[h]
        }
      }
      return(fc)
    },
    betaNorm = function(xx) {
      if (!is.vector(xx)) {
        # TODO: implement
        stop("Not Implemented")
      }
      # dim(t(t(v))) = 3 1, so R works with column-vector by default
      m <- getM()
      beta.T <- rep(1/m, m)
      if (length(xx) > 1) {
        for (t in (length(xx)-1):1) {
          bb.next <- do.call(pdf, c(xx[t+1], pdf.params))
          zz <- (A * bb.next) %*% beta.T
          beta.T <- zz / sum(zz)
        }
      }
      return(t(beta.T))
    },
    interpolateDist = function(xx, x.probes, t.probes) {
      if (!is.vector(xx)) {
        # TODO: implement
        stop("Not Implemented")
      }
      m <- getM()
      T <- length(xx)
      # browser()
      D <- matrix(0, length(t.probes), length(x.probes))
      for (j in 1:length(t.probes)) {
        t.0 <- t.probes[j]
        if ((t.0 < 1) || (t.0 > T)) {
          stop(paste('Invalid time', t.0))
        }
        beta.norm <- betaNorm(xx[t.0:T])
        # B.xx is an l x m matrix, l = |x.probes|
        B.xx <- t(sapply(x.probes, function(x) do.call(pdf, c(x, pdf.params))))
        if (t.0 == 1) {
          zz <- t(priors)
        } else {
          alpha.norm <- alphaNorm(xx[1:(t.0-1)])
          zz <- alpha.norm %*% A
        }
        # Matrix by vector scalar multiplication works by columns, hence t(B.xx)
        num <- t(as.vector(zz) * t(B.xx)) %*% t(beta.norm)
        denom <- as.vector(zz %*% t(beta.norm))
        D[j,] <- t(num / denom)
      }
      return(D)
    },
    interpolate = function(xx, x.probes, t.probes, method='mode') {
      if (!is.vector(xx)) {
        # TODO: implement
        stop("Not Implemented")
      }
      D <- interpolateDist(xx, x.probes, t.probes)
      if (method == 'mode') {
        ii <- x.probes[apply(D, 1, which.max)]
      } else if (method == 'mean') {
        M <- matrix(rep(D %*% x.probes, length(x.probes)),
                    length(x.probes), length(t.probes), byrow=TRUE)
        ii <- apply(abs(x.probes - M), 2, which.min)
      } else if (method == 'median') {
        ii <- apply(abs(t(apply(D, 1, cumsum)) - 0.5), 1, which.min)
      } else {
        stop(paste('Unknown method', method))
      }
      return(ii)
    }
  )
)