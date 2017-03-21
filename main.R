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
x.probes <- 0:50
hmm <- PoisHmm(4, xx)
hmm$fit(xx)
# hmm$plotForecastDist(hmm$forecastDist(xx, x.probes, 16), xx, x.probes)
fc <- hmm$forecast(xx, x.probes, 16, 'median')
