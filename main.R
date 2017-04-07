source('pois_hmm.R')
source('import_earthquakes.R')

# xx <- earthquakes.usgs.annual$count
#
# hmm <- PoisHmm(m=2, xx=xx)
# hmm$fit(xx)
# print(hmm)
#
# A <- matrix(c(0.6, 0.4,
#               0.3, 0.7), 2, 2, byrow=TRUE)
# priors <- c(0.1, 0.9)
# pdf.params <- list(lambda = c(1, 2))
# hmm <- PoisHmm(A=A, priors=priors, pdf.params=pdf.params)
# print(hmm)
# hmm$fit(xx)

TestClass <- setRefClass("TestClass",
  fields = list(
    a = "numeric",
    b = "numeric"
  ),
  methods = list(
    # initialize = function(a=NULL, b=NULL) {
    #   if (is.null(a)) {
    #     a <<- 3
    #   } else {
    #     a <<- a
    #   }
    #   if (is.null(b)) {
    #     b <<- 4
    #   } else {
    #     b <<- b
    #   }
    # }
    initialize = function(a, b) {
      if (missing(a)) {
        a <<- 3
      } else {
        a <<- a
      }
      if (missing(b)) {
        b <<- 4
      } else {
        b <<- b
      }
    },
    show = function() {
      print(' [*] a =')
      cat(' [v] b =', a, '\n')
    }
  )
)
#
t <- TestClass(1)
print(t)
# t <- TestClass(a=1)
# print(t)
# t <- TestClass(b=2)
# print(t)
#
# t <- TestClass$new(a = 8, b = 9)
# print(t)