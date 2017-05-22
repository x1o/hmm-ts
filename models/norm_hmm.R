source('hmm.R')

NormHmm <- setRefClass("NormHmm",
  contains = "Hmm",
  fields = list(
  ),
  methods = list(
    initialize = function(m, xx, A, priors, pdf.params, ...) {
      if (missing(A) && missing(priors) && missing(pdf.params)) {
        callSuper(m, xx, ...)
        pdf.params <<- list(
          mean = seq(min(xx), max(xx), length.out = m + 2)[2:(m + 1)],
          sd = seq(0, sd(xx), length.out = m + 2)[2:(m + 1)]
        )
      } else if (missing(m) && missing(xx)) {
        callSuper(A=A, priors=priors)
        pdf.params <<- pdf.params
      } else {
        stop("Either m, xx or A, priors and pdf.params should be specified as arguments.")
      }
      pdf <<- dnorm
      rng <<- rnorm
    },
    getWrkParams = function() {
      as.numeric(
        c(pdf.params$mean,
          log(pdf.params$sd),
          callSuper())
      )
    },
    setWrkParams = function(params) {
      # mean, sd, A (by columns excluding diag), priors : m + m + m(m-1) + (m-1)
      # m <- -1 + sqrt(2 + length(params))
      m <- getM()
      callSuper(params, 2*m)
      pdf.params <<- list(mean = params[1:m],
                          sd = exp(params[(m + 1):(2 * m)]))
    }
  )
)