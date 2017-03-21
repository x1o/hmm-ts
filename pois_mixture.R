source('import_earthquakes.R')

# xx <- earthquakes.zuc.annual$count
xx <- earthquakes.usgs.annual$count

# plot(earthquakes.usgs.annual$year, xs, type="o", pch=20)
# # h <- hist(xs, breaks=20, plot=FALSE)
# # h$density <- h$density / length(h$breaks)
# # plot(h, freq=FALSE, col="grey")
# lambda=mean(xs)
# hist(xs, breaks=30, freq=FALSE, ylab='', xlab='', col='grey',
#      xlim=c(0, max(xs)), ylim=c(0, dpois(floor(lambda), lambda)))
# par(new=TRUE)
# plot(dpois(0:max(xs), lambda), pch=16, axes=FALSE,
#      xlim=c(0, max(xs)), ylim=c(0, dpois(floor(lambda), lambda)))

nat2wrk <- function(par.nat) {
  # 2m "natural" to 2m - 1 "working" paramteters
  log(c(par.nat$lambda,
        par.nat$delta[-1] / par.nat$delta[1]))
}

wrk2nat <- function(par.wrk) {
  # 2m - 1 "working" to 2m "natural" parameters
  m <- (length(par.wrk) + 1) / 2
  delta <- 1
  if (m > 1) {
    delta <- c(delta, exp(par.wrk[(m+1):(2*m-1)]))
  }
  list(lambda=exp(par.wrk[1:m]),
       delta=delta / sum(delta))
}

mllk <- function(par.wrk, x) {
  par.nat <- wrk2nat(par.wrk)
  # f <- function(lambda) dpois(x, lambda)
  # -sum(log(sapply(par.nat$lambda, f) %*% par.nat$delta))
  -sum(log(outer(x, par.nat$lambda, dpois) %*% par.nat$delta))  # 20% faster
}

mixture <- function(pdf, pdf.params, weights) {
  function(x) {
    as.numeric(sapply(pdf.params, function(p) pdf(x, p)) %*% weights)
  }
}

mixture.mean <- function(pdf.params, weights) {
  as.numeric(pdf.params %*% weights)
}

mixture.var <- function(pdf.params, weights) {
  as.numeric(weights %*% (pdf.params + pdf.params^2) - mixture.mean(pdf.params, weights)^2)
}

test.pois.mixture <- function() {
  mfrow.old <- par()$mfrow
  par(mfrow=c(2, 2))
  # list(lambda=c(5, 10, 15, 25, 30), delta=c(1, 1, 1, 1, 1)/5),
  for (par.init in list(
    list(lambda=c(10), delta=c(1)),
    list(lambda=c(6, 18), delta=c(1, 1)/3),
    list(lambda=c(6, 12, 18), delta=c(1, 1, 1)/3),
    list(lambda=c(10, 20, 25, 30), delta=c(1, 1, 1, 1)/4)
  )) {
    m <- length(par.init$lambda)
    cat("m =", m, "\n")
    par.opt.wrk <- nlm(mllk, nat2wrk(par.init), xs)
    par.opt.nat <- wrk2nat(par.opt.wrk$estimate)
    cat("mllk =", par.opt.wrk$minimum, "\n")
    par.opt.nat
    eta <- mixture(dpois, par.opt.nat$lambda, par.opt.nat$delta)
    
    h <- hist(xs, breaks=20, plot=FALSE)
    ymax <- max(dpois(floor(mean(xs)), mean(xs)), max(h$density))
    # h$density <- h$density / length(h$breaks)
    plot(h, freq=FALSE, col="grey", main=paste("m = ", m), ylab='', xlab='', axes=FALSE,
         xlim=c(0, max(xs)), ylim=c(0, ymax))
    par(new=TRUE)
    plot(eta(0:max(xs)), type="l", 
         xlim=c(0, max(xs)), ylim=c(0, ymax))
    
    cat("mean =", mixture.mean(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    cat("variance =", mixture.var(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    cat("\n")
  }
  
  cat("Data mean: ", mean(xs), "\n")
  cat("Data variance: ", var(xs), "\n")
  par(mfrow=mfrow.old)
}

# test.pois.mixture()