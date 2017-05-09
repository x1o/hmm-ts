plotDataSummary <- function(xx) {
  par(mfrow=c(2,2))
  plot(xx, type='o', pch=20, xlab='')
  hist(xx, main='', xlab='')
  acf(xx, length(xx), main='', xlab='')
  par(mfrow=c(1,1))
}

plotForecast <- function(fc, xx, h.max, xx.true=numeric(0), ylim, ...) {
  T <- length(xx)
  if (missing(ylim)) {
    ylim <- range(c(xx, xx.true))
  }
  plot(xx, type='o', pch=20, xlab='', ylab='', bty='n',
       xlim=c(1, T+h.max), ylim=ylim, ...)
  tt.fc <- (T+1):(T+h.max)
  band.colors <- c('#B1B5CE', '#DBDBDF')
  for (j in rev(1:ncol(fc$upper))) {
    polygon(c(tt.fc, rev(tt.fc)),
            c(fc$lower[,j], rev(fc$upper[,j])),
            col=band.colors[j], border=FALSE)
  }
  par(new=TRUE)
  plot(tt.fc, fc$forecast, col='blue',
       type='o', pch=20, axes=FALSE, xlab='', ylab='',
       xlim=c(1, T+h.max), ylim=ylim, ...)
  abline(v=T, lty=2)
  if (length(xx.true) != 0) {
    par(new=TRUE)
    plot(tt.fc, xx.true,
         xlim=c(1, T+h.max), ylim=ylim,
         xlab='', ylab='',
         type='o', pch=21, lty='dashed', axes=FALSE, ...)
  }
  box()
  # close-up
  # plot((T.train+1):T, xx[(T.train+1):T], col='green', type='o', pch=20, ylim=ylim)
  # par(new=TRUE)
  # plot((T.train+1):T, x.probes[apply(F, 1, which.max)],
  #      type='o', pch=1, axes=FALSE, ylim=ylim)
  # plotCI((T.train+1):T, x.probes[apply(F, 1, which.max)],
  #        ui=x.probes[ui], li=x.probes[li], ylim=ylim,
  #        add=TRUE, pch=NA)
}

plotForecastDist <- function(D, xx, x.probes, h.plot=c(), xx.true=c()) {
  plot.idcs <- 1:nrow(D)
  if (length(h.plot) != 0) {
    plot.idcs <- h.plot
  }
  n.plots <- length(plot.idcs)
  par(mfrow=c(floor(sqrt(n.plots)), ceiling(n.plots/floor(sqrt(n.plots)))))
  for (h in plot.idcs) {
    xlim <- c(min(x.probes), max(x.probes))
    ylim <- c(0, max(D[h,]))
    # plot(x.probes, D[h,], type="h", main='',
    #      xlab='Observations', ylab='Probabilities',
    #      xlim=range(xx), ylim=c(0, max(D)))
    plot(x.probes, D[h,], type="h", main='',
         xlab='Observations', ylab='Probabilities',
         xlim=xlim, ylim=ylim)
    title(paste("h =", h), line = 2)
    par(new=TRUE)
    if (length(xx.true) > 0) {
      # if (is.character(xx) && (length(xx) == 1)) {
      #   # a string
      #   x.true <- substr(xx.true, h, h)
      # } else {
      #   x.true <- xx.true[h]
      # }
      x.true <- xx.true[h]
      x.idx.true <- which(x.probes == x.true)
      plot(x.true, D[h, x.idx.true],
           type="h", lwd=3, axes=FALSE,
           xlab='', ylab='',
           xlim=xlim, ylim=ylim)
      # axis(3, at=round(x.probes[argmax], 2))
      # axis(1, at=1:33, labels=hmm$qq)
      # axis(2, at=seq(0, 0.14, by=0.02))
      # box()
    } else {
      argmax = which.max(D[h,])
      # plot(x.probes[argmax], D[h, argmax], type="h", lwd=3, axes=FALSE,
      #      xlab='', ylab='',
      #      xlim=range(xx), ylim=c(0, max(D)))
      plot(x.probes[argmax], D[h, argmax], type="h", lwd=3, axes=FALSE,
           xlab='', ylab='',
           xlim=xlim, ylim=ylim)
      axis(3, at=round(x.probes[argmax], 2))
      # cat('h = ', h, ': F[', argmax, '] = ', D[h, argmax], '\n', sep='')
    }
  }
  par(mfrow=c(1,1))
}

plotParams <- function(pp) {
  library('corrplot')
  if (!is.matrix(pp)) {
    pp <- matrix(pp,  1, length(pp))
  }
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                             "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                             "#4393C3", "#2166AC", "#053061"))
  corrplot(pp,
           method='circle', is.corr=FALSE, addCoef.col = 'black',
           col=col2(100)[30:100])
}

plotModel <- function(hmm) {
  # FIXME: relative width
  # def.par <- par(no.readonly = TRUE)
  # layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  # plotParams(hmm$A)
  # plotParams(hmm$priors)
  # plotParams(hmm$pdf.params[[1]])
  # for (k in 1:length(hmm$pdf.params)) {
  #
  # }
  # par(def.par)
}