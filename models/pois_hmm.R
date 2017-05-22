source('hmm.R')

PoisHmm <- setRefClass("PoisHmm",
  contains = "Hmm",
  fields = list(
  ),
  methods = list(
    initialize = function(m, xx, A, priors, pdf.params, ...) {
      if (missing(A) && missing(priors) && missing(pdf.params)) {
        callSuper(m, xx, ...)
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
      pp <- as.numeric(
        c(log(pdf.params$lambda),
          callSuper())
      )
      non.f.idcs <- !is.finite(pp)
      if (sum(non.f.idcs) > 0) {
        cat('Oops! Some of the parameters aren\'t finite.\n')
        pp[non.f.idcs] <- .Machine$integer.max * sign(pp[non.f.idcs])
      }
      return(pp)
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