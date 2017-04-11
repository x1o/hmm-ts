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

str.to.idcs <- function(s, alphabet) {
  match(unlist(strsplit(s, '')), alphabet)
}

idcs.to.str <- function(idcs, alphabet) {
  paste(alphabet[idcs], collapse='')
}

is.string <- function(s) {
  is.character(s) && (length(s) == 1)
}

str.to.chars <- function(s) {
  unlist(strsplit(s, ''))
}

# TODO: this should probably be split into CategHmm and StringHmm

CategHmm <- setRefClass("CategHmm",
  contains = "Hmm",
  fields = list(
    qq = "character"
  ),
  methods = list(
    initialize = function(m, xx, A, priors, pdf.params, qq=character(0), ...) {
      if (is.string(qq)) {
        qq <- str.to.chars(qq)
      }
      if (missing(A) && missing(priors) && missing(pdf.params)) {
        # xx is either a "string", or a character or numeric vector
        if (is.character(xx)) {
          if (length(nchar(xx)) == 1) {
            # "string"
            xx <- str.to.chars(xx)
          }
          # char vector
          qq <<- sort(unique(xx))
          xx <- match(xx, .self$qq)
          k <- length(.self$qq)
        } else {
          k <- length(unique(xx))
          if (length(qq) > 0) {
            if (k != length(qq)) {
              stop("Observation alphabet length is not equal to the length of the sample alphabet.")
            }
            qq <<- qq
          }
        }
        # now xx is numeric
        callSuper(m, xx, ...)
        # FIXME: fair initialization?
        m <- getM()
        # h <- hist(xx, breaks=0:k, plot=FALSE)
        # pdf.params <<- list(prob=matrix(rep(h$density, m), m, k, byrow=TRUE))
        pdf.params <<- list(prob=matrix(1/k, m, k, byrow=TRUE))
      } else if (missing(m) && missing(xx)) {
        callSuper(A=A, priors=priors)
        pdf.params <<- pdf.params
        if (length(qq) > 0) {
          qq <<- qq
        }
      } else {
        stop("Either m, xx or A, priors and pdf.params should be specified as arguments.")
      }
      pdf <<- ddiscr
      rng <<- rdiscr
      if (length(.self$qq) > 0) {
        colnames(pdf.params$prob) <<- .self$qq
      }
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
    xx.to.numeric = function(xx) {
      if (is.character(xx)) {
        if (length(qq) == 0) {
          stop('No observation alphabet specified.')
        }
        if (length(nchar(xx)) == 1) {
          xx <- str.to.chars(xx)
        }
        xx <- match(xx, .self$qq)
      }
      return(xx)
    },
    fit = function(xx, ...) {
      callSuper(xx.to.numeric(xx), ...)
    },
    logL = function(xx, ...) {
      callSuper(xx.to.numeric(xx), ...)
    },
    genSample = function(T) {
      state <- walkChain(T)
      probs <- split(pdf.params$prob, row(pdf.params$prob))
      xx <- rng(T, probs[state])
      if (length(qq) > 0) {
        xx <- idcs.to.str(xx, qq)
      }
      return(xx)
    },
    interpolateDist = function(xx, x.probes, t.probes) {
      callSuper(xx.to.numeric(xx), xx.to.numeric(x.probes), xx.to.numeric(t.probes))
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