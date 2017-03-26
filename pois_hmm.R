source('hmm.R')

PoisHmm <- setRefClass("PoisHmm",
  contains = "Hmm",
  fields = list(
  ),
  methods = list(
    initialize = function(m, xx) {
      callSuper(m, xx)
      pdf <<- dpois
      pdf.params <<- list(
        lambda = quantile(xx, seq(0, 1, 1/(m + 1)))[2:(m + 1)]
      )
      rng <<- rpois
    },
    getWrkParams = function() {
      as.numeric(
        c(log(pdf.params$lambda),
          callSuper())
      )
    },
    setWrkParams = function(params) {
      # lambda, A (by columns excluding diag), priors : m + m(m-1) + (m-1)
      # m <- (-1 + sqrt(5 + 4 * length(params))) / 2
      m <- getM()
      callSuper(params, 1)
      pdf.params <<- list(lambda = exp(params[1:m]))
    }
  )
)