test.pois.mixt <- function() {
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
    par.opt.wrk <- nlm(mllk.pois.mix, nat2wrk.pois.mix(par.init), xx)
    par.opt.nat <- wrk2nat.pois.mix(par.opt.wrk$estimate)
    cat("mllk =", par.opt.wrk$minimum, "\n")
    par.opt.nat
    eta <- pois.mix(dpois, par.opt.nat$lambda, par.opt.nat$delta)

    h <- hist(xx, breaks=20, plot=FALSE)
    ymax <- max(dpois(floor(mean(xx)), mean(xx)), max(h$density))
    # h$density <- h$density / length(h$breaks)
    plot(h, freq=FALSE, col="grey", main=paste("m = ", m), ylab='', xlab='', axes=FALSE,
         xlim=c(0, max(xx)), ylim=c(0, ymax))
    par(new=TRUE)
    plot(eta(0:max(xx)), type="l",
         xlim=c(0, max(xx)), ylim=c(0, ymax))

    cat("mean =", pois.mix.mean(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    cat("variance =", pois.mix.var(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    cat("\n")
  }

  cat("Data mean: ", mean(xx), "\n")
  cat("Data variance: ", var(xx), "\n")
  par(mfrow=mfrow.old)
}

test.norm.mix <- function(xx, m) {
  mfrow.old <- par()$mfrow
  par(mfrow=c(floor(sqrt(m)), ceiling(m/floor(sqrt(m)))))
  par.init.list <- list()
  sdmax <- 1.5*sd(xx)
  for (k in seq(1:m)) {
    par.init.list[[k]] <- list(mu=seq(min(xx), max(xx), length.out=k+2)[2:(k+1)],
                               sigma=seq(0, sdmax, length.out=k+2)[2:(k+1)],
                               delta=rep(1, k)/k)
  }
  # print(par.init.list)
  for (par.init in par.init.list) {
    m <- length(par.init$mu)
    cat("m =", m, "\n")
    par.opt.wrk <- nlm(mllk.norm.mix, nat2wrk.norm.mix(par.init), xx)
    par.opt.nat <- wrk2nat.norm.mix(par.opt.wrk$estimate)
    # par.opt.wrk
    cat("mllk =", par.opt.wrk$minimum, "\n")
    # print(par.opt.nat)
    p.eta <- dnorm.mix(par.opt.nat[1:2], par.opt.nat$delta)
    hist(xx, breaks=30, freq=FALSE, col="grey",
         main=paste("m = ", m, "; -log(L) = ", round(par.opt.wrk$minimum, 4)),
         xlab="Mean magnitude")

    # h <- hist(eq.annual$mag.mean, breaks=30, plot=FALSE)
    # h$density <- h$density / length(h$breaks)
    # plot(h, freq=FALSE, col="grey", main=paste("m = ", m))

    curve(p.eta, from=min(xx), to=max(xx),col="black", add=TRUE)
    # cat("mean =", mixture.mean(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    # cat("variance =", mixture.var(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    # cat("\n")
  }

  cat("Data mean: ", mean(xx), "\n")
  cat("Data variance: ", var(xx), "\n")

  par(mfrow=mfrow.old)
}