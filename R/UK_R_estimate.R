library(tidyverse)
library(EpiEstim)
library(here) 

setwd(here()) # should set working directory to project folder; if not, specify manually

# over-write the EpiEstim package functions with these, which are 
# simply edited so that they also output the 10th and 90th quantiles
source("R/estimate_r_edited_for_quantile_output.R")

# also need to re-source "utilities.R" from EpiEstim:
source("R/utilities.R")

# load population data for LTLAs of UK (from ONS):
population.by.ltla <- read_csv("data/population.by.ltla.csv")

ltlas.alphabetical <- sort(population.by.ltla$Area)
ltla.nums <- 1:length(ltlas.alphabetical)

##### prepare backcalculation functions

# zeta is the distribution for the delay from infection to case confirmation
zeta.max <- 35 # hard-code that no delays are ever longer than that
delays.considered <- seq(0, zeta.max)

# incubation period
s.meanlog <- 1.58
s.sdlog   <- 0.47
s <- function(tau) {dlnorm(tau, sdlog = s.sdlog, meanlog = s.meanlog)}
s.discrete <- vapply(delays.considered, s, numeric(1))
s.discrete <- s.discrete / sum(s.discrete)

symptom.to.swab.mean <- 1.2
symptom.to.swab.sd   <- 1
symptom.to.swab.shape <- symptom.to.swab.mean^2 / symptom.to.swab.sd^2
symptom.to.swab.rate <- symptom.to.swab.shape / symptom.to.swab.mean
symptom.to.swab.pdf <- function(t) {
  dgamma(t, shape = symptom.to.swab.shape, rate = symptom.to.swab.rate)
}
symptom.to.swab.discrete <- vapply(delays.considered, symptom.to.swab.pdf, numeric(1))
symptom.to.swab.discrete <- symptom.to.swab.discrete / sum(symptom.to.swab.discrete)



# Zeta is the convolution of s and symptom to swab time.
# Careful: we include day zero but vector indexing is 1-based
zeta <- rep(NA, zeta.max + 1)
for (t in delays.considered) {
  convolution.sum.range <- seq(0, t)
  zeta[[t + 1]] <- sum(s.discrete[convolution.sum.range + 1] *
                         symptom.to.swab.discrete[t - convolution.sum.range + 1])
}
zeta <- zeta / sum(zeta) # unit normalise



### download cases by specimen date for all available UK LTLAs:

ltla.spec.date <- read_csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumCasesBySpecimenDate&format=csv")

ltla.spec.date <- ltla.spec.date %>% arrange(date)    # sort into ascending date order                          

# filter by date: as we only want a point estimate at the most recent date, to save time we don't need to go all the way back through the epidemic
start.date <- as.Date("2020-12-01")
last.date <- as.Date(max(ltla.spec.date$date))
all.dates <- seq(start.date, last.date, by="day")
dat.UK.ltla <- ltla.spec.date %>% filter(date %in% all.dates)


### compute
ltlas.incidence <- lapply(ltlas.alphabetical, function(area) {
  print(area)
  dat.area <- dat.UK.ltla %>% filter(areaName == area)
  
  # # the most recent date is sometimes missing, e.g. when only one of England and Wales has updated.
  # # when this is the case, replicate the last entry for those missing one
  # if (!(last.date %in% dat.area$date)) {
  #   dat.area[nrow(dat.area)+1,] <- dat.area[nrow(dat.area),]
  #   dat.area[nrow(dat.area),"date"] <- last.date
  # }
  # stopifnot(!any(setdiff(all.dates,dat.area$date) > min(dat.area$date)) ) # check that missing dates are all at the start, so we are safe to fill in cumulative cases as zeroes

  
  # the most recent date or two are sometimes missing, e.g. when only one of England and Wales has updated. 
  # NB this is particularly true as of 17th April, when Wales moved to six day reporting.
  # When this is the case, replicate the last entry or two for those missing them
  dat.area <- dat.area %>%
    complete(date = as.Date(union(dat.area$date, c(last.date - 1, last.date)), origin="1970-01-01"), fill=list(areaName = area, areaType= unique(dat.area$areaType), areaCode = unique(dat.area$areaCode), cumCasesBySpecimenDate=max(dat.area$cumCasesBySpecimenDate)))
  
  # then fill in any necessary zeroes at the start
  dat.area <- dat.area %>%
    complete(date = all.dates, fill=list(areaName = area, areaType= unique(dat.area$areaType), areaCode = unique(dat.area$areaCode), cumCasesBySpecimenDate=0))
  
  dat.area$Incidence <- rep(0,nrow(dat.area))
  dat.area$Incidence[1] <- dat.area$cumCasesBySpecimenDate[1]
  for(row in 2:nrow(dat.area)) dat.area$Incidence[row] <- dat.area$cumCasesBySpecimenDate[row] - dat.area$cumCasesBySpecimenDate[row-1]
  
  dat.area
})

# backcalculate the daily new infections
ltlas.incidence.backcalculation <- lapply(ltla.nums, function(x) {
  # Get the case counts time series, fill in any missing dates, and add in dates
  # beforehand reaching back to the maximum possible delay (which induces NAs,
  # all of which we set to zero).
  df <- cbind.data.frame(
    "dates" = ltlas.incidence[[x]]$date,
    "counts" = ltlas.incidence[[x]]$Incidence
  )
  df <- df %>%
    complete(dates = seq.Date(min(dates) - zeta.max, max(dates), by="day")) %>%
    replace_na(list(counts = 0))
  
  df$infections <- 0
  for (days.plus.1.since.start in seq(1, nrow(df))) {
    integration.range <- seq(days.plus.1.since.start,
                             min(nrow(df), days.plus.1.since.start + zeta.max))
    df$infections[[days.plus.1.since.start]] <-
      sum(df$counts[integration.range] *
            zeta[integration.range - days.plus.1.since.start + 1])
  }
  df
})

# setup start and end times for the windows over which to estimate R. 
# We get EpiEstim warnings about "estimating R too early in the epidemic"
# but this is ok as we are not interested in the accuracy of our early estimates.
t_start <- seq(2,as.numeric(last.date - start.date) + zeta.max - 5) 
t_end <- t_start + 7 - 1

ltlas.R.backcalculated <- lapply(ltla.nums, function(i) {
  print(ltlas.alphabetical[[i]])
  df <- ltlas.incidence.backcalculation[[i]]
  
  # df contains data for the ltla: date, swab count, inferred new infections
  # use inferred new infections to estimate R
  
  ltla.incidence <- cbind.data.frame(
    "dates" = df$dates,
    "I" = df$infections
  )
  
  ltla.R.backcalculated <- estimate_R(ltla.incidence, 
                                      method="parametric_si",
                                      config = make_config(list(
                                        t_start = t_start,
                                        t_end = t_end,
                                        mean_si = 5.5,  # NB now using the generation time distribution because we're using inferred times of infection, not cases
                                        std_si = 2.14,
                                        mean_prior = 1,
                                        std_prior = 1))
  )
  
  # change to mode
  ltla.R.backcalculated$R$`Mean(R)` <- ltla.R.backcalculated$R$`Mean(R)` - (ltla.R.backcalculated$R$`Std(R)`)^2 / ltla.R.backcalculated$R$`Mean(R)`
  
  ltla.R.backcalculated
})

# collect results of R estimates by area and date
df.R.ltlas <- cbind.data.frame(
  "Dates" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$dates[t_end] - 4)), # R estimates are labelled by the end of the week over which they were calculated; shift it to the middle
  "R" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$R$`Mean(R)`)),
  "Area" = unlist(lapply(ltla.nums, function(area) rep(ltlas.alphabetical[[area]], nrow(ltlas.R.backcalculated[[area]]$R)))),
  #"AreaCode" = unlist(lapply(1:nrow(ltla.codes), function(area) rep(ltla.codes$Code[[area]], nrow(ltlas.R.backcalculated[[area]]$R)))),
  "lower" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$R$`Quantile.0.025(R)`)),
  "tenth" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$R$`Quantile.0.10(R)`)),
  "twentyfifth" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$R$`Quantile.0.25(R)`)),
  "seventyfifth" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$R$`Quantile.0.75(R)`)),
  "ninetieth" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$R$`Quantile.0.90(R)`)),
  "upper" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$R$`Quantile.0.975(R)`))
)

# recover date formatting
df.R.ltlas$Dates <- as.Date(df.R.ltlas$Dates, origin="1970-01-01") 


df.incidence.ltlas <- cbind.data.frame(
  "Dates" = unlist(lapply(ltla.nums, function(area) ltlas.incidence.backcalculation[[area]]$dates)),
  "Incidence" = unlist(lapply(ltla.nums, function(area) ltlas.incidence.backcalculation[[area]]$infections)),
  "Area"= unlist(lapply(ltla.nums, function(area) rep(ltlas.alphabetical[[area]], length(ltlas.incidence.backcalculation[[area]]$dates))))#,
  #"AreaCode" = unlist(lapply(1:nrow(ltla.codes), function(area) rep(ltla.codes$Code[[area]], length(ltlas.incidence.backcalculation[[area]]$dates))))
)

# combine with population data
df.incidence.ltlas <- left_join(df.incidence.ltlas, population.by.ltla)

# scale incidence per 100,000 population
df.incidence.ltlas <- df.incidence.ltlas %>%
  mutate("scaled_per_capita" = Incidence / population * 100000)

# recover date formatting
df.incidence.ltlas$Dates <- as.Date(df.incidence.ltlas$Dates,  origin = "1970-01-01")


# calculate nowcast projections
projected.cases.ltlas <- cbind.data.frame(
  "Dates" = unlist(lapply(ltla.nums, function(area) ltlas.R.backcalculated[[area]]$dates[t_end] - 4)),
  "Projection" = unlist(lapply(ltla.nums, function(area) {
    sapply(ltlas.R.backcalculated[[area]]$R$t_end, function(i) { # the first t_end is 8
      mean(ltlas.R.backcalculated[[area]]$I[(i-6):i]) * # average incidence over that week
        ltlas.R.backcalculated[[area]]$R$`Mean(R)`[[which(ltlas.R.backcalculated[[area]]$R$t_end == i)]] # last R value
    })
  }
  )),
  "Area"= unlist(lapply(ltla.nums, function(area) rep(ltlas.alphabetical[[area]], nrow(ltlas.R.backcalculated[[area]]$R))))#,
  #"AreaCode" = unlist(lapply(1:nrow(ltla.codes), function(area) rep(ltla.codes$Code[[area]], nrow(ltlas.R.backcalculated[[area]]$R))))
)

# recover date formatting
projected.cases.ltlas$Dates <- as.Date(projected.cases.ltlas$Dates,  origin = as.Date("1970-01-01"))

# combine with population data
projected.cases.ltlas <- left_join(projected.cases.ltlas, population.by.ltla)

# scale nowcast per 100,000 population
projected.cases.ltlas <- projected.cases.ltlas %>%
  mutate("scaled_per_capita" = Projection / population * 100000)

R.trim <- 9 # we censor by 9 days to account for incomplete data and the need to "see ahead" how infections increased or decreased beyond the date in question

# get the estimate of R and nowcast value for the last date before we censor
last.date.before.censored <- max(df.R.ltlas$Dates) - R.trim
ltlas.latest.R <- df.R.ltlas %>% 
  filter(Dates == last.date.before.censored)
ltlas.latest.nowcast <- projected.cases.ltlas %>% 
  filter(Dates == last.date.before.censored)  

# combine R and nowcast to get a weighted estimate per ltla
combined.R.and.nowcast <- left_join(ltlas.latest.R, ltlas.latest.nowcast)
combined.R.and.nowcast <- combined.R.and.nowcast %>%
  mutate("R times N" = R * scaled_per_capita) %>%
  mutate("Rlower times N" = lower * scaled_per_capita) %>%
  mutate("Rupper times N" = upper * scaled_per_capita) %>%
  mutate("R25th times N" = twentyfifth * scaled_per_capita) %>%
  mutate("R75th times N" = seventyfifth * scaled_per_capita) %>%
  mutate("R10th times N" = tenth * scaled_per_capita) %>%
  mutate("R90th times N" = ninetieth * scaled_per_capita) 

uplift <- 0.06 # recent R estimates have been consistently underestimated by 0.06 during 2021, in comparison to the final estimate obtained with a week or more's further data. Keeping this under review.
# NB at the time of writing, setting R.trim=12 requires an uplift of 0.02; setting R.trim=15 requires no uplift

combined.estimate <- sum(combined.R.and.nowcast$`R times N`) / sum(combined.R.and.nowcast$scaled_per_capita) + uplift
combined.estimate.lower <- sum(combined.R.and.nowcast$`Rlower times N`) / sum(combined.R.and.nowcast$scaled_per_capita) + uplift
combined.estimate.upper <- sum(combined.R.and.nowcast$`Rupper times N`) / sum(combined.R.and.nowcast$scaled_per_capita) + uplift
combined.estimate.25th <- sum(combined.R.and.nowcast$`R25th times N`) / sum(combined.R.and.nowcast$scaled_per_capita) + uplift
combined.estimate.75th <- sum(combined.R.and.nowcast$`R75th times N`) / sum(combined.R.and.nowcast$scaled_per_capita) + uplift
combined.estimate.10th <- sum(combined.R.and.nowcast$`R10th times N`) / sum(combined.R.and.nowcast$scaled_per_capita) + uplift
combined.estimate.90th <- sum(combined.R.and.nowcast$`R90th times N`) / sum(combined.R.and.nowcast$scaled_per_capita) + uplift

results <- cbind.data.frame("Central estimate" = combined.estimate,
                         "2.5th quantile" = combined.estimate.lower,
                         "10th quantile" = combined.estimate.10th,
                         "25th quantile" = combined.estimate.25th,
                         "75th quantile" = combined.estimate.75th,
                         "90th quantile" = combined.estimate.90th,
                         "97.5th quantile" = combined.estimate.upper)

### Convert our estimates of R to a growth rate r
# Assume generation time is gamma distributed with
alpha <- 6.6
beta <- 1.2
# then using Wallinga and Lipsitch: https://royalsocietypublishing.org/doi/10.1098/rspb.2006.3754 we have
# R <- function(r) ((beta + r)/beta)^alpha   
# which gives 

r <- function(R) R^(1/alpha) * beta - beta
full.results <- rbind.data.frame(
  "R" = results,
  "r" = r(results)
)

paste0("Estimates for ",format(last.date.before.censored, "%d %B %Y"),":")
full.results




