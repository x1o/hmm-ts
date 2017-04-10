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
      priors <<- rep(1, m) / m
      priors <<- priors + runif(m, -1/r*1/m, 1/r*1/m)
    },
    init.from.params = function(A, priors) {
      A <<- A
      priors <<- priors
    },
    initialize = function(m, xx, A, priors, ...) {
      if (missing(A) && missing(priors)) {
        init.from.data(m, ...)
      } else if (missing(m) && missing(xx)) {
        init.from.params(A, priors)
      } else {
        stop("Either m, xx or A, priors should be given as arguments.")
      }
    },
    show = function() {
      cat("A:\n")
      methods::show(round(A, 2))
      cat("\npriors:\n")
      methods::show(round(priors, 2))
      cat("\npdf.params:\n")
      methods::show(lapply(pdf.params, function(x) round(x, 2)))
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
    alphaNorm = function(xx, compute.L = FALSE) {
      T <- length(xx)
      alpha <- do.call(pdf, c(xx[1], pdf.params)) * priors
      if (compute.L) {
        L.log <- log(sum(alpha))
      }
      alpha.norm <- alpha / sum(alpha)
      if (T > 1) {
        for (t in 2:T) {
          bb.t <- do.call(pdf, c(xx[t], pdf.params))
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
    logL = function(xx) {
      alphaNorm(xx, compute.L = TRUE)
    },
    mLogL = function(par.wrk, xx) {
      setWrkParams(par.wrk)
      mllk <- -logL(xx)
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
    fit = function(xx, iterlim=1000, method='nlm', ...) {
      if (method == 'nlm') {
        nlm.out <- nlm(.self$mLogL, getWrkParams(), xx, iterlim=iterlim, ...)
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
                        data=list(xx=xx),
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
    forecastDist = function(xx, x.probes, h.max) {
      F <- matrix(0, h.max, length(x.probes))
      cc <- alphaNorm(xx, compute.L = FALSE)
      # B.probe <- t(outer(x.probes, pdf.params$lambda, pdf))  # 4 times faster
      B.probe <- sapply(x.probes, function(x) do.call(pdf, c(x, pdf.params)))
      for (h in 1:h.max) {
        cc <- cc %*% A
        F[h,] <- cc %*% B.probe
      }
      # return(F)
      # FIXME: is normalization necessary? Works fine w/o it in the case of PoisHmm
      return(F / rowSums(F))
    },
    forecast = function(xx, x.probes, h.max, method='mode') {
      T <- length(xx)
      F <- forecastDist(xx, x.probes, h.max)
      if (method == 'mode') {
        fc.idcs <- apply(F, 1, which.max)
      } else if (method == 'mean') {
        fc.idcs <- as.numeric(F %*% x.probes)
        # FIXME: rewrite loop
        for (h in 1:length(fc.idcs)) {
          fc.idcs[h] <- which.min(abs(x.probes - fc.idcs[h]))
        }
      } else if (method == 'median') {
        # FIXME: rewrite loop
        fc.idcs <- numeric(h.max)
        for (h in 1:length(fc.idcs)) {
          fc.idcs[h] <- which.min(abs(cumsum(F[h,]) - 0.5))
        }
      } else {
        stop(paste('Unknown method', method))
      }
      fc <- list()
      fc$'forecast' <- x.probes[fc.idcs]
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
        # FIXME: rewrite loop
        for (h in 1:h.max) {
          # print(h)
          pivot <- fc.idcs[h]
          accum <- F[h, pivot]
          ui[h] <- pivot
          li[h] <- pivot
          #FIXME: broken logic: check before adding
          while ((li[h] > 1) || (ui[h] < length(F[h,]))) {
            if (li[h] > 2)  {
              if ((F[h, li[h] - 1] + accum) > conf.level) {
                break
              }
              li[h] <- li[h] - 1
              accum <- accum + F[h, li[h]]
            }
            if (ui[h] < length(F[h,])-1) {
              if ((F[1, ui[h] + 1] + accum) > conf.level) {
                break
              }
              ui[h] <- ui[h] + 1
              accum <- accum + F[h, ui[h]]
            }
            # print(c(li[h], ui[h], accum))
          }
          # print(paste('lower', h, conf.level, '<-', x.probes[li[h]]))
          # print(paste('upper', h, conf.level, '<-', x.probes[ui[h]]))
          fc$lower[h, lvl.idx] <- x.probes[li[h]]
          fc$upper[h, lvl.idx] <- x.probes[ui[h]]
        }
      }
      return(fc)
    },
    betaNorm = function(xx) {
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
    interpolateDist = function(xx, x.probes, t.0) {
      m <- getM()
      T <- length(xx)
      if ((t.0 < 1) || (t.0 > T)) {
        stop(paste('Invalid time', t.0))
      }
      beta.norm <- betaNorm(xx[t.0:T])
      # B.xx is a s x m matrix, s = |x.probes|
      B.xx <- t(sapply(x.probes, function(x) do.call(pdf, c(x, pdf.params))))
      if (t.0 == 1) {
        zz <- t(priors)
      } else {
        alpha.norm <- alphaNorm(xx[1:(t.0-1)])
        zz <- alpha.norm %*% A
      }
      denom <- as.vector(zz %*% t(beta.norm))
      # Matrix by vector scalar multiplication works by columns, hence t(B.xx)
      dd <- (t(as.vector(zz) * t(B.xx)) %*% t(beta.norm)) / denom
      return(t(dd))
    },
    aic = function(xx) {
      -2 * logL(xx) + 2*getDf()
    },
    bic = function(xx) {
      -2 * logL(xx) + getDf() * log(length(xx))
    }
  )
)