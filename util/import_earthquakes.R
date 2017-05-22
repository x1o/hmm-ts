require(dplyr)

earthquakes.usgs <- read.csv('../../data/earthquakes.csv')
earthquakes.usgs.annual <- earthquakes.usgs %>%
  select(year = time, mag) %>%
  transform(year = as.numeric(format(as.Date(year), "%Y"))) %>%
  group_by(year) %>%
  summarise(mag.mean = mean(mag), count = n())

earthquakes.zuc.annual <- read.delim('../../data/earthquakes-zucchini.txt', FALSE)
names(earthquakes.zuc.annual)[1] <- 'year'
names(earthquakes.zuc.annual)[2] <- 'count'
