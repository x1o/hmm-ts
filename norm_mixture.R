
nat2wrk <- function(par.nat) {
  # 3m "natural" to 3m - 1 "working" paramteters
  (c(par.nat$mu,
     log(par.nat$sigma),
     log(par.nat$delta[-1] / par.nat$delta[1])))
}

wrk2nat <- function(par.wrk) {
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

mllk.norm.mixture <- function(par.wrk, x) {
  par.nat <- wrk2nat(par.wrk)
  pdf <- function(param) {
    dnorm(x, mean=param[1], sd=param[2])
  }
  -sum(log(apply(rbind(par.nat$mu, par.nat$sigma), 2, pdf) %*% par.nat$delta))
}

dnorm.mixture <- function(params, weights) {
  function(x) {
    as.numeric(apply(rbind(params$mu, params$sigma),
                     2,
                     function(p) dnorm(x, p[1], p[2])) %*% weights)
  }
}

# p.23
# mixture.mean <- function(pdf.params, weights) {
#   as.numeric(pdf.params %*% weights)
# }
# 
# mixture.var <- function(pdf.params, weights) {
#   as.numeric(weights %*% (pdf.params + pdf.params^2) - mixture.mean(pdf.params, weights)^2)
# }


test.mixtures <- function(xs, m) {
  mfrow.old <- par()$mfrow
  par(mfrow=c(floor(sqrt(m)), ceiling(m/floor(sqrt(m)))))
  par.init.list <- list()
  sdmax <- 1.5*sd(xs)
  for (k in seq(1:m)) {
    par.init.list[[k]] <- list(mu=seq(min(xs), max(xs), length.out=k+2)[2:(k+1)],
                            sigma=seq(0, sdmax, length.out=k+2)[2:(k+1)],
                            delta=rep(1, k)/k)
  }
  # print(par.init.list)
  for (par.init in par.init.list) {
    m <- length(par.init$mu)
    cat("m =", m, "\n")
    par.opt.wrk <- nlm(mllk.norm.mixture, nat2wrk(par.init), xs)
    par.opt.nat <- wrk2nat(par.opt.wrk$estimate)
    # par.opt.wrk
    cat("mllk =", par.opt.wrk$minimum, "\n")
    # print(par.opt.nat)
    p.eta <- dnorm.mixture(par.opt.nat[1:2], par.opt.nat$delta)
    hist(xs, breaks=30, freq=FALSE, col="grey",
         main=paste("m = ", m, "; -log(L) = ", round(par.opt.wrk$minimum, 4)),
         xlab="Mean magnitude")
    
    # h <- hist(eq.annual$mag.mean, breaks=30, plot=FALSE)
    # h$density <- h$density / length(h$breaks)
    # plot(h, freq=FALSE, col="grey", main=paste("m = ", m))
    
    curve(p.eta, from=min(xs), to=max(xs),col="black", add=TRUE)
    # cat("mean =", mixture.mean(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    # cat("variance =", mixture.var(par.opt.nat$lambda, par.opt.nat$delta), "\n")
    # cat("\n")
  }
  
  cat("Data mean: ", mean(xs), "\n")
  cat("Data variance: ", var(xs), "\n")
  
  par(mfrow=mfrow.old)
}

test.mixtures(xs, 9)
