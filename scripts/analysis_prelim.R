# analyzing bee and flower data for phenological skew

library(sn)

fit.sn <- function(abs, times, plot=FALSE, ...){
  abs <- abs*10
  if(sum(abs) == 0 | length(abs) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999))
  #if(length(unique(abs)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999))
  if(length(unique(times)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999))
  if(sum(abs) < 10) return(list(mean=-8888, sd=-8888, skew=-8888, skew_se=-8888))
  
  abs_expand <- rep(times, times=round(abs))
  #model <- selm(abs_expand ~ 1)
  model <- tryCatch(selm(abs_expand ~ 1), error=function(e) "error") # catching errors so it doesn't abort whole analysis
  if(class(model)[1] != "selm") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777))
  summ <- summary(model)
  pars <- summ@param.table
  mean <- pars[1,1]
  sd <- pars[2,1]
  skew <- pars[3,1]
  skew_se <- pars[3,2]
  if(plot) plot(model, ...)
  # I probably want to report total pop size, unique time points, and some goodness of fit?
  return(list(mean=mean, sd=sd, skew=skew, skew_se=skew_se))
}

bee_skews[mean == -7777]
bee_skews <- bees[, fit.sn(bowl_adj, DOY), by=.(species, site, year)]
bee_skews <- bee_skews[mean %!in% c(-7777, -8888, -9999),]
bee_skews <- bee_skews[!is.na(skew_se),]
bee_skews[skew > 0.9]
bee_skews[skew < 0.5 & skew > -0.5]
bees[site == "Beaver" & species == "Lasioglossum ruidosense", fit.sn(bowl_adj, DOY, plot=T)]
bees[site == "Gothic" & species == "Panurginus cressoniellus", fit.sn(bowl_adj, DOY, plot=T)]
bees[site == "Seans" & species == "Panurginus cressoniellus" & year == 2016, fit.sn(bowl_adj, DOY, plot=T)]


hist(bee_skews$mean, 20)
hist(bee_skews$sd, 20)
hist(bee_skews$skew, 20)
hist(bee_skews$skew_se, 20)

bees[1:6, fit.sn(bowl_adj, DOY)]
bees[site == "Willey" & species == "Andrena algida", fit.sn(bowl_adj, DOY)]
bees[site == "Willey" & species == "Andrena algida",]


flower_skews <- flowers[, fit.sn(floralcount, doy), by=.(species, plot, year)] # waiting on this to run
flower_skews <- flower_skews[mean %!in% c(-7777, -8888, -9999),]
flower_skews <- flower_skews[!is.na(skew_se),]
flower_skews[skew > 0.9]
flowers[plot == "RM5" & species == "Erigeron flagellaris", fit.sn(floralcount, doy, plot=T, which=c(2))]

hist(flower_skews$mean, 20)
hist(flower_skews$sd, 20)
hist(flower_skews$skew, 20)
hist(flower_skews$skew_se, 20)


# prelim figs for Becky
par(mfrow=c(2,1))
bees[site == "Gothic" & species == "Panurginus cressoniellus", fit.sn(bowl_adj, DOY, plot=T, which=c(2), main="Panurginus cressoniellus at Gothic")]
flowers[plot == "RM5" & species == "Erigeron flagellaris" & year == 2015, fit.sn(floralcount, doy, plot=T, which=c(2), main="Erigeron falgellaris at RM5")]

hist(bee_skews$skew, 20, main="Bee skews")
hist(flower_skews$skew, 20, main="Flower skews")
spar(mfrow=c(2,1))

# working example

ab <- c(0,0,10,20,35,50,20,5,0)
ab <- c(0,0,10,120,65,50,20,5,1)
time <- seq(110,190, by=10)
ab_expand <- rep(time, times=ab)
plot(ab ~ time)
hist(ab_expand)

model <- selm(ab_expand ~ 1)
summary(model)
plot(model)
summ <- summary(model)
summ@param.table

predict_ab <- predict(model)
predict_ab <- predict(model, data.frame(time = 110:190), param.type="CP")

profile(model, param.name="cp")
     


test <- try(sum(c(2,NA,"a")), silent=T)

tryCatch(sum(c(2,4,"a")), error=function(e) "hello")

tryCatch(1, finally = print("Hello"))
e <- simpleError("test error")
## Not run: 
stop(e)
tryCatch(stop(e), finally = print("Hello"))
tryCatch(stop("fred"), finally = print("Hello"))

## End(Not run)
tryCatch(stop(e), error = function(e) e, finally = print("Hello"))

### does population size bias skew estimates?
library(sn)
big_pop <- rsn(1000, alpha=5)
little_pop <- rsn(100, alpha=5)
hist(big_pop, 100, col=rgb(1,0,0,0.5))
hist(little_pop, 100, col=rgb(0,1,0,0.5), add=T)
model_big <- selm(big_pop ~ 1)
model_little <- selm(little_pop ~ 1)
summary(model_big)
summary(model_little)



# looks like no, but it smaller samples decrease accuracy of estimates, not surprisingly
