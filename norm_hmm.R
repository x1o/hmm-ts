source('hmm.R')

NormHmm <- setRefClass("NormHmm",
  contains = "Hmm",
  fields = list(
  ),
  methods = list(
    initialize = function(m, xx) {
      callSuper(m, xx)
      pdf <<- dnorm
      pdf.params <<- list(
        mean = seq(min(xx), max(xx), length.out = m + 2)[2:(m + 1)],
        sd = seq(0, sd(xx), length.out = m + 2)[2:(m + 1)]
      )
      rng <<- rnorm
    },
    getWrkParams = function() {
      c(pdf.params$mean,
        log(pdf.params$sd),
        callSuper())
    },
    setWrkParams = function(params) {
      # mean, sd, A (by columns excluding diag), priors : m + m + m(m-1) + (m-1)
      # m <- -1 + sqrt(2 + length(params))
      m <- getM()
      callSuper(params, 2)
      pdf.params <<- list(mean = params[1:m],
                          sd = exp(params[(m + 1):(2 * m)]))
    }
  )
)