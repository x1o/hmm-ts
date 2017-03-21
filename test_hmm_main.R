source('import_earthquakes.R')
source('pois_hmm.R')
source('norm_hmm.R')
source('test_hmm.R')

tested.hmms <- list(
  list(PoisHmm, earthquakes.zuc.annual$count, TRUE)
  # list(NormHmm, earthquakes.usgs.annual$mag.mean, FALSE)
)

for (test.case in tested.hmms) {
  HmmClass <- test.case[[1]]
  xx <- test.case[[2]]
  is.discrete.hmm <- test.case[[3]]
  
  ##### Data 'analysis'
  plot(xx, type="o", pch=20)
  plot(earthquakes.usgs.annual$year, earthquakes.usgs.annual$mag.mean,
       type="o", pch=20, ylab="Magnitude", xlab="Year")
  hist(xx, breaks=20, col='grey')
  curve(dnorm(x, mean(xx), sd(xx)),
        from=min(xx), to=max(xx), add=TRUE)
  acf(xx)
  
  ##### Tests
  ### Init
  cat('---------- Init ----------\n')
  hmm <- HmmClass(3, xx)
  print(hmm)

  ### Parameter conversions
  cat('---------- Parameter conversions ----------\n')
  par.wrk <- hmm$getWrkParams()
  hmm$setWrkParams(par.wrk)
  print(hmm)
  # FIXME: hmm$setWrkParams(hmm$getWrkParams()) not working (promise...)

  ### alphaNorm
  cat('---------- alphaNorm ----------\n')
  print(hmm$alphaNorm(xx))
  print(hmm$alphaNorm(xx, compute.L = TRUE))

  ### fit
  cat('---------- fit ----------\n')
  hmm$fit(xx)
  print(hmm)
  testHmmFit(HmmClass, xx, 4, is.discrete.hmm)

  ### Forecasting
  cat('---------- Forecasting ----------\n')
  F <- hmm$forecastDist(xx, 1:50, 4)
  plot(F[4,], type="h", xlim=c(0,50), ylim=c(0, 0.1), main=paste("h =", 4))
  testHmmForecast(HmmClass, xx, 16, is.discrete.hmm)

  ### Sampling
  cat('---------- Sampling ----------\n')
  hist(hmm$genSample(100000), breaks=30, col='grey', freq=FALSE)
  #
  # FIXME
  # xx <- earthquakes.usgs.annual$mag.mean
  # h.max <- 32
  # T <- length(xx)
  # T.train <- T-h.max
  # xx.train <- xx[1:(T.train)]
  # hmm <- NormHmm(3, xx.train)
  # hmm$fit(xx.train)
  # # F <- hmm$forecastDist(xx, x.probes, 16)
  # # hmm$plotForecastDist(F, xx, x.probes)
  # plot(xx.train, type='o', pch=20,
  #      xlim=c(1, T), ylim=range(xx))
  # par(new=TRUE)
  # plot((T.train+1):T, hmm$genSample(h.max),
  #      type='o', pch=1, axes=FALSE,
  #      xlim=c(1, T), ylim=range(xx))
  # abline(v=T.train, lty=2)
  # par(new=TRUE)
  # plot((T.train+1):T, xx[(T.train+1):T], col='green',
  #      type='o', pch=20,
  #      xlim=c(1, T), ylim=range(xx))
  
  ### AIC & BIC
  cat('---------- AIC & BIC ----------\n')
  selectModelIc(HmmClass, xx, 10, do.plot = TRUE)
}
