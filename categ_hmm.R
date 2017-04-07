source('hmm.R')

ddiscr <- function(x, prob) {
  prob[,x]
}

rdiscr <- function(n, prob) {
  out <- numeric(n)
  if (!is.list(prob)) {
    out <- sample(1:length(prob), n, replace=TRUE, prob=prob)
  } else {
    for (i in 1:n) {
      j <- (i - 1) %% length(prob) + 1
      out[i] <- sample(1:length(prob[[j]]), 1, replace=TRUE, prob=prob[[j]])
    }
  }
  return(out)
}

str.to.idcs <- function(str, alphabet) {
  match(unlist(strsplit(str, '')), alphabet)
}

idcs.to.str <- function(idcs, alphabet) {
  paste(alphabet[idcs], collapse='')
}

CategHmm <- setRefClass("CategHmm",
  contains = "Hmm",
  fields = list(
  ),
  methods = list(
    initialize = function(m, xx, A, priors, pdf.params) {
      if (missing(A) && missing(priors) && missing(pdf.params)) {
        callSuper(m, xx)
        # FIXME: better initialization
        k <- length(unique(xx))
        m <- getM()
        pdf.params <<- list(prob=matrix(1/k, m, k, byrow=TRUE))
      } else if (missing(m) && missing(xx)) {
        callSuper(A=A, priors=priors)
        pdf.params <<- pdf.params
      } else {
        stop("Either m, xx or A, priors and pdf.params should be specified as arguments.")
      }
      pdf <<- ddiscr
      rng <<- rdiscr
    },
    getWrkParams = function() {
      as.numeric(
        c(log(as.numeric(t(pdf.params$prob[,-1] / pdf.params$prob[,1]))),
          callSuper())
      )
    },
    setWrkParams = function(params) {
      # prob, A (by columns excluding diag), priors : m*(k-1) + m(m-1) + (m-1)
      # m <- (-1 + sqrt(5 + 4 * length(params))) / 2
      m <- getM()
      # k <- (length(params) - m*(m-1) - (m-1))/m + 1
      k <- getK()
      callSuper(params, m*(k-1))
      prob <- cbind(
        matrix(1, m, 1),
        matrix(exp(params[1:(m*(k-1))]), m, k-1, byrow=TRUE)
      )
      pdf.params <<- list(prob=prob/rowSums(prob))
    },
    genSample = function(T) {
      state <- walkChain(T)
      probs <- split(pdf.params$prob, row(pdf.params$prob))
      rng(T, probs[state])
    },
    getK = function() {
      ncol(pdf.params$prob)
    },
    getDf = function() {
      m <- getM()
      # row-stochastic A + priors + pdf.params
      m*(m-1) + (m-1) + (getK() - 1) * m
    }
  )
)