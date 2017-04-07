source('hmm.R')

PoisHmm <- setRefClass("PoisHmm",
  contains = "Hmm",
  fields = list(
  ),
  methods = list(
    initialize = function(m, xx, A, priors, pdf.params) {
      if (missing(A) && missing(priors) && missing(pdf.params)) {
        callSuper(m, xx)
        pdf.params <<- list(
          lambda = quantile(xx, seq(0, 1, 1/(m + 1)))[2:(m + 1)]
        )
      } else if (missing(m) && missing(xx)) {
        callSuper(A=A, priors=priors)
        pdf.params <<- pdf.params
      } else {
        stop("Either m, xx or A, priors and pdf.params should be specified as arguments.")
      }
      pdf <<- dpois
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
      callSuper(params, m)
      pdf.params <<- list(lambda = exp(params[1:m]))
    }
  )
)