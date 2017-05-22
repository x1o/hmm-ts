nat2wrk.pois.mix <- function(par.nat) {
  # 2m "natural" to 2m - 1 "working" paramteters
  log(c(par.nat$lambda,
        par.nat$delta[-1] / par.nat$delta[1]))
}

wrk2nat.pois.mix <- function(par.wrk) {
  # 2m - 1 "working" to 2m "natural" parameters
  m <- (length(par.wrk) + 1) / 2
  delta <- 1
  if (m > 1) {
    delta <- c(delta, exp(par.wrk[(m+1):(2*m-1)]))
  }
  list(lambda=exp(par.wrk[1:m]),
       delta=delta / sum(delta))
}

mllk.pois.mix <- function(par.wrk, x) {
  par.nat <- wrk2nat(par.wrk)
  # f <- function(lambda) dpois(x, lambda)
  # -sum(log(sapply(par.nat$lambda, f) %*% par.nat$delta))
  -sum(log(outer(x, par.nat$lambda, dpois) %*% par.nat$delta))  # 20% faster
}

pois.mix <- function(pdf, pdf.params, weights) {
  function(x) {
    as.numeric(sapply(pdf.params, function(p) pdf(x, p)) %*% weights)
  }
}

pois.mix.mean <- function(pdf.params, weights) {
  as.numeric(pdf.params %*% weights)
}

pois.mix.var <- function(pdf.params, weights) {
  as.numeric(weights %*% (pdf.params + pdf.params^2) - pois.mix.mean(pdf.params, weights)^2)
}