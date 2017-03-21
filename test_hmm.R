testHmmFit <- function(HmmClass, xx, m.max, is.discrete.hmm) {
  mfrow.old <- par()$mfrow
  par(mfrow=c(floor(sqrt(m.max)), ceiling(m.max / floor(sqrt(m.max)))))
  hmm.test.list <- list()
  for (m in seq(1:m.max)) {
    hmm <- HmmClass(m, xx)
    hmm.test.list[[m]] <- hmm
  }
  for (hmm in hmm.test.list) {
    print(hmm)
    m <- hmm$getM()
    cat("--> m =", m, "<--\n")
    nlm.out <- hmm$fit(xx)
    # print(nlm.out)
    cat("mllk =", nlm.out$minimum,
        "\nn.iter =", nlm.out$iterations, "\n")
    h <- hist(xx, breaks=nclass.FD, plot=FALSE)
    # h$density <- h$density / length(h$breaks)
    xmax <- max(xx)
    xmin <- min(xx)
    # FIXME: not only for Pois
    if (is.discrete.hmm) {
      ymax <- ymax <- max(h$density)
      for (j in 1:m) {
        ymax <- max(ymax, hmm$pdf(xmin:xmax, hmm$pdf.params$lambda[j]))
      }
    } else {
      ymax <- max(h$density)
    }
    ymax <- ymax * 1.2
    plot(h, col='grey', freq=FALSE,
         main=paste("m = ", m, "; -log(L) = ", round(nlm.out$minimum, 4)),
         ylim=c(0, ymax), xlim=c(xmin, xmax))
    for (j in 1:m) {
      if (is.discrete.hmm) {
        par(new=TRUE)
        # FIXME: kludge please fix me
        plot(hmm$pdf(xmin:xmax, hmm$pdf.params$lambda[j]),
             type='l', axes=FALSE,
             ylab='', xlab='', col='red',
             ylim=c(0, ymax), xlim=c(0, xmax))
      } else {
        # FIXME: kludge please fix me
        # curve(do.call(pdf, c(x=x, ps)), from=xmin, to=xmax) +
        # + do.call(rbind, pdf.params) doesn't seem to work
        curve(hmm$pdf(x, hmm$pdf.params$mean[j], hmm$pdf.params$sd[j]),
              add=TRUE,
              from=xmin, to=xmax)
      }
    }
    # FIXME: plot mixture density
    # par(new=TRUE)
    # plot(dpois.hmm(min(xs):max(xs), par.opt), type='l', col='blue', axes=FALSE,
    #      ylab='', xlab='',
    #      ylim=c(0, 0.10), xlim=c(0, max(xs)))
  }
  cat("Data mean: ", mean(xx), "\n")
  cat("Data variance: ", var(xx), "\n")
  par(mfrow=mfrow.old)
}

testHmmForecast <- function(HmmClass, xx, h.max, is.discrete.hmm, discr.rk=50) {
  mfrow.old <- par()$mfrow
  par(mfrow=c(floor(sqrt(h.max)), ceiling(h.max/floor(sqrt(h.max)))))
  if (is.discrete.hmm) {
    x.probes <- min(xx):max(xx)
  } else {
    x.probes <- seq(min(xx), max(xx), length.out = discr.rk)
  }
  hmm <- HmmClass(3, xx)
  hmm$fit(xx)
  F <- hmm$forecastDist(xx, x.probes, h = h.max)
  for (h in 1:h.max) {
    plot(x.probes, F[h,], type="h", main=paste("h =", h),
         xlim=range(xx), ylim=c(0, max(F)))
    par(new=TRUE)
    argmax = which.max(F[h,])
    plot(x.probes[argmax], F[h, argmax], type="h", lwd=3, axes=FALSE,
         xlim=range(xx), ylim=c(0, max(F)))
    axis(3, at=round(x.probes[argmax], 2))
    cat('h =', h, ':', c(argmax, F[h, argmax]), '\n')
  }
  par(mfrow=mfrow.old)
  plot(xx, type='o', pch=20,
       xlim=c(0, length(xx) + h.max), ylim=range(xx))
  par(new=TRUE)
  plot((length(xx)+1):(length(xx)+h.max), x.probes[apply(F, 1, which.max)],
       type='o', pch=20, col='red', axes=FALSE,
       xlim=c(0, length(xx) + h.max), ylim=range(xx))
  abline(v=length(xx), lty=2)
}
