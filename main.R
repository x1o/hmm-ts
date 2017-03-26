# TODO:
# - E, D of norm
# - MC ACF
# - Fit MC
# - Find MC stationary distribution
# - stationary mc distribution -- when important? one extra multiplication is
#   not a big deal...
# - NULL => stationary
# - hmm pdf in the case of a stationary MC
# - HMM for m = 1 (make sure it's working)
# - sd or sd^2?
# - hist of gen.sample -- unimodal??
# - test on model solution
# - confidence intervals
# - compare pois hmm component densities with zuc
# - proper initialize, missing(), ...
# - optim, constrOptim
# - forecast for continuous
# - while fitting, check the code; if not 1, act appropriately

# - В качестве прогноза попробовать среднее (может будет прыгать)
# - Почитать Hasan и другие статьи по hmm+forecast
# - Понять, почему важно иметь предсказание распределения (для геологов) и т.д.
# - Задачи, где СММ выбирают модель
# - Максимум hist и curve / pdf в графиках

# fix plotForecast --- xlim via x.probes
# plotForecast -> plotForecastDist
# plotForecast --- mode, mean, median

# - prediction on smoothed ts?
# - is hmm necessary? if acf != 0, just work with differences

source('import_earthquakes.R')
source('pois_hmm.R')
source('norm_hmm.R')

xx <- earthquakes.zuc.annual$count
T <- 100
tt <- seq(0, 4*pi, length.out=T)
xx <- 1 + 2*tt
xx <- xx + sin(tt)
xx <- runif(100, -1, 1)
plot(xx)
hmm <- PoisHmm(4, xx)
hmm$fit(xx)

# x.probes <- 0:50
# # hmm$plotForecastDist(hmm$forecastDist(xx, x.probes, 16), xx, x.probes)
# fc <- hmm$forecast(xx, x.probes, 16, 'median')

# x <- 1:5
# set.seed(1001)
# y <- rbinom(5,prob=x/(1+x),size=10)
# mfun <- function(p, lol) {
#   # a <- p[1]
#   # b <- p[2]
#   # -sum(dbinom(y,prob=a*x/(b+x),size=10,log=TRUE))
#   if (lol != 'hi') return(0)
#   -sum(dbinom(y,prob=p[1]*x/(p[2]+x),size=10,log=TRUE))
# }
# # optim(fn=mfun,par=c(1,1))
# # parnames(mfun) <- c("a","b")
# parnames(mfun) <- as.character(1:2)
# v <- c(1,1)
# names(v) <- as.character(1:2)
# mle2(minuslogl=mfun,start=v,data=list(lol='hi'),method="Nelder-Mead")
