nat2wrk.norm.mix <- function(par.nat) {
  # 3m "natural" to 3m - 1 "working" paramteters
  (c(par.nat$mu,
     log(par.nat$sigma),
     log(par.nat$delta[-1] / par.nat$delta[1])))
}

wrk2nat.norm.mix <- function(par.wrk) {
  # 3m - 1 "working" to 3m "natural" parameters
  m <- (length(par.wrk) + 1) / 3
  delta <- 1
  if (m > 1) {
    delta <- c(delta, exp(par.wrk[(2*m+1):(3*m-1)]))
  }
  list(mu=par.wrk[1:m],
       sigma=exp(par.wrk[(m+1):(2*m)]),
       delta=delta / sum(delta))
}

mllk.norm.mix <- function(par.wrk, x) {
  par.nat <- wrk2nat(par.wrk)
  pdf <- function(param) {
    dnorm(x, mean=param[1], sd=param[2])
  }
  -sum(log(apply(rbind(par.nat$mu, par.nat$sigma), 2, pdf) %*% par.nat$delta))
}

dnorm.mix <- function(params, weights) {
  function(x) {
    as.numeric(apply(rbind(params$mu, params$sigma),
                     2,
                     function(p) dnorm(x, p[1], p[2])) %*% weights)
  }
}

# p.23
# norm.mix.mean <- function(pdf.params, weights) {
#   as.numeric(pdf.params %*% weights)
# }

# norm.mix.var <- function(pdf.params, weights) {
#   as.numeric(weights %*% (pdf.params + pdf.params^2) - mixture.mean(pdf.params, weights)^2)
# }