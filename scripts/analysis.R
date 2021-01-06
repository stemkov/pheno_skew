# analyzing bee and flower data for phenological skew

library(sn)
library(data.table)
library(tools)
library(mgsub)
library(lubridate)
library(visreg)

`%!in%` = Negate(`%in%`)

setwd("/home/michael/Documents/Grad School/Research Projects/pheno_skew")

if(!(exists("bees") | exists("flowers"))){
  bees <- fread("clean_data/bees.csv", stringsAsFactors = FALSE)
  flowers <- fread("clean_data/flowers.csv", stringsAsFactors = FALSE)
}

fit.sn <- function(abs, times, plot=FALSE, ...){
  # abs <- c(0,0,1,10,4,3,2,0,1)
  # times <- seq(110,190, by=10)
  #abs <- abs*10
  if(sum(abs) == 0 | length(abs) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999, skew_p=-9999))
  #if(length(unique(abs)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999))
  if(length(unique(times)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999, skew_p=-9999))
  if(sum(abs) < 10) return(list(mean=-8888, sd=-8888, skew=-8888, skew_se=-8888, skew_p=-8888))
  
  abs_expand <- rep(times, times=round(abs))
  #model <- selm(abs_expand ~ 1)
  model <- tryCatch(selm(abs_expand ~ 1), error=function(e) "error") # catching errors so it doesn't abort whole analysis
  if(class(model)[1] != "selm") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777, skew_p=-7777))
  summ <- summary(model)
  pars <- summ@param.table
  mean <- pars[1,1]
  sd <- pars[2,1]
  skew <- pars[3,1]
  skew_se <- pars[3,2]
  skew_p <- pars[3,4]
  if(plot) plot(model, ...)
  # I probably want to report total pop size, unique time points, and some goodness of fit?
  return(list(mean=mean, sd=sd, skew=skew, skew_se=skew_se, skew_p=skew_p))
}

### I want to do the following:
# - skews for individual time-series. for each year/species/site timeseries
# - aggredate across genera, whole community
# - aggredate across years
# - try a full model with multiple predictors - this seems to have problems running with more than just year and species, and only estimates one gamma parameter
# - does skew vary between years? ... maybe predicted by drought index?
# - are skews predicted by strength of interaction?
# - show that population size doesn't matter *

### setting up color palettes
bee_pal <- c(point = "#f09c00", line = "#c48000")
flr_pal <- c(point = "#a655cf", line = "#7b04b8")

### full model with multiple predictors - this takes really long time to run

# bees_expand <- bees[rep(1:.N,ab)][,index:=1:.N,by=doy]
# bees_expand$year <- as.factor(bees_expand$year)
# test <- selm(doy ~ year + species, data=bees_expand)
# summary(test)
# plot(test)
# 
# hist(bees_expand$doy)

### accessory skew distribution figure function
skew.fig <- function(skew = 0, percent = NA, x1=0, x2=0.5, y1=0, y2=0.5, line_col="red", fill_col = "pink", text_cex=2, shade_percent=0){
  # line_col <- "red"
  # fill_col <- "pink"
  # skew <- 0.8
  # percent <- 75
  # text_cex <- 2
  # shade_percent <- 50
  
  #par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
  par(fig = c(x1,x2,y1,y2), new = T, mgp=c(3,0.7,0)) 
  
  #  xi=0, omega=1, alpha=0, tau=0,
  pars <- cp2dp(c(0,0,1,skew), "SN")
  vals <- seq(-3, 3, by=0.1)
  pdf <- dsn(vals, dp=pars[c(1,3,4)])
  plot(pdf ~ vals, type="l", lwd=2, col=line_col, axes=F, xlab="", ylab="")
  poly_x <- c(vals, max(vals), min(vals))
  poly_y <- c(pdf, 0, 0)
  polygon(poly_x, poly_y,
          col=adjustcolor(fill_col,0.5), border=NA)
  if(shade_percent > 0){
    y_thresh <- max(pdf)*(shade_percent/100)
    poly_y[which(poly_y > y_thresh)] <- y_thresh
    polygon(poly_x, poly_y,
            col="white", border=NA,
            density=25, angle=45)
  }
  lines(pdf ~ vals, lwd=2, col=line_col)
  if(!is.na(percent)) text((vals[which(pdf == max(pdf))]+median(vals))/2, max(pdf)/2.5, paste0(percent,"%"), cex=text_cex, col=adjustcolor(line_col, 1,.5,.5,.5))
  
  par(mar=c(5,4,4,2))
}

### BEES each species/site/year time-series individually - meta analysis method
bees_female <- bees[sex == "F",] # getting rid of male timeseries

n_combs <- uniqueN(bees_female[,c("site", "species", "year", "dataset")])
bee_skews <- bees_female[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy)}, by=.(species, site, year, dataset)]
print(paste0(round(nrow(bee_skews[mean == -9999,])/nrow(bee_skews)*100, 1), "% of time-series too short"))
print(paste0(round(nrow(bee_skews[mean == -8888,])/nrow(bee_skews)*100, 1), "% too little catch"))
print(paste0(round(nrow(bee_skews[mean == -7777,])/nrow(bee_skews)*100, 1), "% error in fitting model"))
print(paste0(round(nrow(bee_skews[is.na(skew_se),])/nrow(bee_skews)*100, 1), "% likely stopped sampling too early/late"))

bee_skews <- bee_skews[mean %!in% c(-7777, -8888, -9999),]
bee_skews <- bee_skews[!is.na(skew_se),]
#merge(bee_skews, bees[,.()], by)

par(mfrow=c(2,1), mar=c(4.1,4,2,2))
hist(bee_skews$mean, seq(min(bee_skews$mean), max(bee_skews$mean), length.out = 20),
     main="Mean phenology", xlab="Day of year")
hist(bee_skews$sd, seq(min(bee_skews$sd), max(bee_skews$sd), length.out = 20),
     main="Phenology standard deviation", xlab="Days")

png("figures/bee_skews_lowest.png", width=700, height=500)
par(mfrow=c(1,1), mar=c(5,4,4,2))
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
breaks <- seq(min(bee_skews$skew), max(bee_skews$skew), length.out = 20)
perc_breaks <- c(nrow(bee_skews[skew < breaks[length(breaks)/2]]),
                 nrow(bee_skews[skew > breaks[length(breaks)/2] & skew < breaks[length(breaks)/2+1]]),
                 nrow(bee_skews[skew > breaks[length(breaks)/2+1]]))
perc_breaks <- round((perc_breaks/nrow(bee_skews))*100,1)
hist(bee_skews$skew, breaks,
     main="Bee Skews", xlab="Marginal skewness", border=FALSE,
     col=c(rep(bee_pal["point"],9), bee_pal["line"], rep(bee_pal["point"],9)))
skew.fig(skew=-0.8, percent=perc_breaks[1], line_col = bee_pal["point"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.05,x2=0.35,y1=0.35,y2=0.7)
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0.8, percent=perc_breaks[3], line_col = bee_pal["point"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.65,x2=0.95,y1=0.4,y2=0.75)
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0, percent=perc_breaks[2], line_col = bee_pal["line"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.25,x2=0.54,y1=0.6,y2=0.95)
dev.off()

png("figures/bee_mean_sd_lowest.png", width=900, height=450)
par(mfrow=c(1,2), mar=c(4.1,4,2,1))
bee_skews_just <- bee_skews[abs(skew) > 0.05,] # just skews that were more than 5% different from 0
model <- lm(skew ~ mean, data = bee_skews_just)
#plot(skew ~ mean, data = bee_skews_just)
summary(model)
visreg(model, xlab="Mean phenology", ylab="Marginal skewness", ylim=c(-1,1),
       points = list(col=bee_pal["point"]), line.par=list(col=bee_pal["line"]))

model <- lm(abs(skew) ~ log(sd), data = bee_skews_just)
#plot(abs(skew) ~ sd, data = bee_skews_just)
summary(model)
visreg(model, xlab="Phenology standard deviation", ylab="Absolute marginal skewness",
       points = list(col=bee_pal["point"]), line.par=list(col=bee_pal["line"]))
par(mfrow=c(1,1), mar=c(5,4,4,2))
dev.off()


### FLOWERS each species/site/year time-series individually
n_combs <- uniqueN(flowers[,c("site", "species", "year")])
flr_skews <- flowers[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy)}, by=.(species, site, year)]
print(paste0(round(nrow(flr_skews[mean == -9999,])/nrow(flr_skews)*100, 1), "% of time-series too short"))
print(paste0(round(nrow(flr_skews[mean == -8888,])/nrow(flr_skews)*100, 1), "% too little catch"))
print(paste0(round(nrow(flr_skews[mean == -7777,])/nrow(flr_skews)*100, 1), "% error in fitting model"))
print(paste0(round(nrow(flr_skews[is.na(skew_se),])/nrow(flr_skews)*100, 1), "% likely stopped sampling too early/late"))

flr_skews <- flr_skews[mean %!in% c(-7777, -8888, -9999),]
flr_skews <- flr_skews[!is.na(skew_se),]
#merge(flr_skews, bees[,.()], by)

par(mfrow=c(2,1), mar=c(4.1,4,2,2))
hist(flr_skews$mean, seq(min(flr_skews$mean), max(flr_skews$mean), length.out = 20),
     main="Mean phenology", xlab="Day of year")
hist(flr_skews$sd, seq(min(flr_skews$sd), max(flr_skews$sd), length.out = 20),
     main="Phenology standard deviation", xlab="Days")

png("figures/flower_skews_lowest.png", width=700, height=500)
par(mfrow=c(1,1), mar=c(5,4,4,2))
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
breaks <- seq(min(flr_skews$skew), max(flr_skews$skew), length.out = 20)
perc_breaks <- c(nrow(flr_skews[skew < breaks[length(breaks)/2]]),
                 nrow(flr_skews[skew > breaks[length(breaks)/2] & skew < breaks[length(breaks)/2+1]]),
                 nrow(flr_skews[skew > breaks[length(breaks)/2+1]]))
perc_breaks <- round((perc_breaks/nrow(flr_skews))*100,1)
hist(flr_skews$skew, breaks,
     main="Flower Skews", xlab="Marginal skewness", border=FALSE,
     col=c(rep(flr_pal["point"],9), flr_pal["line"], rep(flr_pal["point"],9)))
skew.fig(skew=-0.85, percent=perc_breaks[1], line_col = flr_pal["point"], fill_col = flr_pal["point"], text_cex = 1,
         x1=0.05,x2=0.35,y1=0.17,y2=0.52)
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0.85, percent=perc_breaks[3], line_col = flr_pal["point"], fill_col = flr_pal["point"], text_cex = 1,
         x1=0.65,x2=0.95,y1=0.4,y2=0.75)
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0, percent=perc_breaks[2], line_col = flr_pal["line"], fill_col = flr_pal["point"], text_cex = 1,
         x1=0.25,x2=0.54,y1=0.6,y2=0.95)
dev.off()

png("figures/flower_mean_sd_lowest.png", width=900, height=450)
par(mfrow=c(1,2), mar=c(4.1,4,2,1))
flr_skews_just <- flr_skews[abs(skew) > 0.05,] # just skews that were more than 5% different from 0
#flr_skews_just <- flr_skews
model <- lm(skew ~ mean, data = flr_skews_just)
#plot(skew ~ mean, data = flr_skews_just)
summary(model)
visreg(model, xlab="Mean phenology", ylab="Marginal skewness", ylim=c(-1,1),
       points = list(col=flr_pal["point"]), line.par=list(col=flr_pal["line"]))

model <- lm(abs(skew) ~ log(sd), data = flr_skews_just)
#plot(abs(skew) ~ sd, data = flr_skews_just)
summary(model)
visreg(model, xlab="Phenology standard deviation", ylab="Absolute marginal skewness",
       points = list(col=flr_pal["point"]), line.par=list(col=flr_pal["line"]))
par(mfrow=c(1,1), mar=c(5,4,4,2))
dev.off()


# denoting Bombus and others - I expect skews to be different
bee_skews[, group := ifelse(grepl("Bombus", species), "Bombus", "other"), by=species]
# placeholder in flower_skews just so I can rbind later
flr_skews$group <- "NA"
flr_skews$dataset <- "rmbl"

### Bombus vs others
bombus_model <- lm(skew ~ group, data=bee_skews)
summary(bombus_model) # There is a pretty major difference in skew between bombus and others - as expected
par(mfrow=c(1,1), mar=c(5,4,4,2))
visreg(bombus_model)

png("figures/bee_skews_bombus.png", width=700, height=500)
par(mfrow=c(1,1), mar=c(5,4,4,2))
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
breaks <- seq(min(bee_skews$skew), max(bee_skews$skew), length.out = 20)
perc_breaks <- c(nrow(bee_skews[skew < breaks[length(breaks)/2]]),
                 nrow(bee_skews[skew > breaks[length(breaks)/2] & skew < breaks[length(breaks)/2+1]]),
                 nrow(bee_skews[skew > breaks[length(breaks)/2+1]]))
perc_bombus <- c(nrow(bee_skews[skew < breaks[length(breaks)/2] & group == "Bombus"]),
                 nrow(bee_skews[skew > breaks[length(breaks)/2] & skew < breaks[length(breaks)/2+1] & group == "Bombus"]),
                 nrow(bee_skews[skew > breaks[length(breaks)/2+1] & group == "Bombus"]))
perc_bombus <- round((perc_bombus/perc_breaks)*100,1)
perc_breaks <- round((perc_breaks/nrow(bee_skews))*100,1)
hist(bee_skews$skew, breaks,
     main="Bee Skews", xlab="Marginal skewness", border=FALSE,
     col=c(rep(bee_pal["point"],9), bee_pal["line"], rep(bee_pal["point"],9)))
hist(bee_skews[group=="Bombus", skew], breaks,
     main="Bee Skews", xlab="Marginal skewness", border=FALSE,
     col="white",
     density=25, angle=46, add=T)
legend("topleft", c("Bombus", "Other bees", "Flowers"), fill=c(bee_pal["point"], bee_pal["point"], flr_pal["line"]),
       density=c(NA,NA,NA), border=c("white","white","white"), cex=1, inset=0.05)
legend("topleft", "", fill="white", bty="n",
       density=25, border="white", cex=1, inset=0.05)
skew.fig(skew=-0.85, percent=perc_breaks[1], line_col = bee_pal["point"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.05,x2=0.35,y1=0.35,y2=0.7,
         shade_percent = perc_bombus[1])
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0.85, percent=perc_breaks[3], line_col = bee_pal["point"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.65,x2=0.95,y1=0.4,y2=0.75,
         shade_percent = perc_bombus[3])
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0, percent=perc_breaks[2], line_col = bee_pal["line"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.25,x2=0.54,y1=0.6,y2=0.95,
         shade_percent = perc_bombus[2])
dev.off()




### joint figure
png("figures/mean_sd_lowest.png", width=1000, height=550)

bee_skews$guild <- "bee"
flr_skews$guild <- "flower"
skews <- rbind(flr_skews, bee_skews)
skews$guild <- relevel(as.factor(skews$guild), "flower") # to put bee points on top

par(mfrow=c(1,2), mar=c(4.1,4,2,1))
skews_just <- skews[abs(skew) > 0.05,] # just skews that were more than 5% different from 0
#flr_skews_just <- flr_skews
mean_model <- lm(skew ~ mean*guild, data = skews_just)
#plot(skew ~ mean, data = flr_skews_just)
summary(mean_model)
visreg(mean_model, "mean", by="guild", xlab="Mean phenology (day of year)", ylab="Marginal skewness", ylim=c(-1,1), overlay=T,
       points = list(col=c(flr_pal["point"], bee_pal["line"]),
                     cex=c(0.2,0.7),
                     pch=c(19,21),
                     bg=c(1,bee_pal["point"])),
       line.par = list(col=c(flr_pal["line"], bee_pal["line"])),
       fill.par = list(col=adjustcolor(c(flr_pal["point"], bee_pal["point"]), 0.5)),
       legend=F)

sd_model <- lm(abs(skew) ~ log(sd)*guild, data = skews_just)
#plot(abs(skew) ~ sd, data = bee_skews_just)
summary(sd_model)
visreg(sd_model, "sd", by="guild", xlab="Phenological breadth (days)", ylab="Absolute marginal skewness", overlay=T, ylim=c(0,1),
       points = list(col=c(flr_pal["point"], bee_pal["line"]),
                     cex=c(0.2,0.7),
                     pch=c(19,21),
                     bg=c(1,bee_pal["point"])),
       line.par = list(col=c(flr_pal["line"], bee_pal["line"])),
       fill.par = list(col=adjustcolor(c(flr_pal["point"], bee_pal["point"]), 0.5)),
       legend=F)
legend("right", c("Flowers", "Bees"), col=c(flr_pal["line"], bee_pal["line"]),
       lwd=2, cex=1, inset=0.05)
par(mfrow=c(1,1), mar=c(5,4,4,2))
dev.off()



### aggregating by all species/sites across years
# it shows the typical distributions for whole communities at sites
# it DOES NOT show the total community distributions - those would have more complex shapes that sn isn't appropriate for
# dashed lines are those years in which we only have bombus, so they look really left skewed
# 2012 is a drought year, 2019 is a relatively wet year (following a drought year of 2018)
bee_years <- bees_female[, .(ab = sum(ab),
                             date = unique(date),
                             n_species = length(unique(species))),
                          by = .(year, dataset, doy, site)]
flr_years <- flowers[, .(ab = sum(ab),
                             date = unique(date),
                             n_species = length(unique(species))),
                         by = .(year, doy, site)]

n_combs <- uniqueN(bee_years[,c("year", "dataset", "site")])
bee_yr_skews <- bee_years[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy, plot=F, which=2, main=paste(year, dataset, site))}, by=.(year, dataset, site)]
bee_yr_skews <- bee_yr_skews[mean %!in% c(-7777, -8888, -9999),]
bee_yr_skews <- bee_yr_skews[!is.na(skew_se),]

n_combs <- uniqueN(flr_years[,c("year", "site")])
flr_yr_skews <- flr_years[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy, plot=F, which=2, main=paste(year, site))}, by=.(year, site)]
flr_yr_skews <- flr_yr_skews[mean %!in% c(-7777, -8888, -9999),]
flr_yr_skews <- flr_yr_skews[!is.na(skew_se),]

dist.plot <- function(mean=0, sd=1, skew=0, year=NA, year_sub=c(2012,2019), xlim=c(100,250), col="red", ...){
  # col <- "red"
  # mean <- 100
  # sd <- 20
  # skew <- 0.8
  # xlim <- c(0,200)
  print(paste(year, "skew:", skew))
  pars <- cp2dp(c(mean,0,sd,skew), "SN")
  vals <- seq(xlim[1], xlim[2], by=1)
  pdf <- dsn(vals, dp=pars[c(1,3,4)])
  lines(pdf ~ vals, col=col, ...)
  if(!is.na(year) & year %in% year_sub) text(vals[which(pdf == max(pdf))], max(pdf)*1.05, year, col=adjustcolor(col, 2))
  #plot(pdf ~ vals, type="l", lwd=2, col=line_col, axes=F, xlab="", ylab="")
}

png("figures/year_dists.png", width=800, height=600)
par(mfrow=c(1,1), mar=c(5,4,4,2), mgp=c(3, 0.7, 0))
xlim <- c(100,275)
plot(NA, xlim=xlim, ylim=c(0,.034),
     yaxt="n", ylab="",
     xlab="Day of year")
title(ylab="Abundance probability", line=1, cex.lab=1, family="Calibri Light")
flr_yr_skews[year %in% c(2009,2010,2011,2013,2014,2015,2016,2017,2018), dist.plot(mean(mean), mean(sd), mean(skew),
                                              year = unique(year),  xlim=xlim,
                                              col = adjustcolor(flr_pal["line"], 0.5),
                                              lwd=2),
             by=year]


bee_yr_skews[dataset=="rmbl" & year %in% c(2009,2010,2011,2013,2014,2015,2016,2017,2018) , dist.plot(mean(mean), mean(sd), mean(skew),
                                         year = unique(year),  xlim=xlim,
                                         col = adjustcolor(bee_pal["line"], 0.6),
                                         lwd=2),
             by=year]

bee_yr_skews[dataset=="rmbl" & year == 2020 , dist.plot(mean(mean), mean(sd), mean(skew),
                                                                  year = unique(year),  xlim=xlim,
                                                                  col = adjustcolor(bee_pal["line"], 0.6),
                                                                  lwd=2, lty=2),
             by=year] # these years only have bombus

# years to highlight

flr_yr_skews[year %in% c(2012, 2019), dist.plot(mean(mean), mean(sd), mean(skew),
                                                year = unique(year),  xlim=xlim,
                                                col = "white",
                                                lwd=5),
             by=year]
flr_yr_skews[year %in% c(2012, 2019), dist.plot(mean(mean), mean(sd), mean(skew),
                                                                                  year = unique(year),  xlim=xlim,
                                                                                  col = adjustcolor(flr_pal["line"], 0.8),
                                                                                  lwd = 3),
             by=year]

bee_yr_skews[dataset=="rmbl" & year %in% c(2012) , dist.plot(mean(mean), mean(sd), mean(skew),
                                                                  year = unique(year),  xlim=xlim,
                                                                  col = "white",
                                                                  lwd=5),
             by=year]
bee_yr_skews[dataset=="rmbl" & year %in% c(2012) , dist.plot(mean(mean), mean(sd), mean(skew),
                                                             year = unique(year),  xlim=xlim,
                                                             col = adjustcolor(bee_pal["line"], 0.8),
                                                             lwd=3),
             by=year]

bee_yr_skews[dataset=="rmbl" & year == 2019 , dist.plot(mean(mean), mean(sd), mean(skew),
                                                       year = unique(year),  xlim=xlim,
                                                       col = "white",
                                                       lwd=5, lty=2),
             by=year]
bee_yr_skews[dataset=="rmbl" & year == 2019 , dist.plot(mean(mean), mean(sd), mean(skew),
                                                       year = unique(year),  xlim=xlim,
                                                       col = adjustcolor(bee_pal["line"], 0.8),
                                                       lwd=3, lty=2),
             by=year]
legend("topleft", c("Flowers", "All bees", "Bombus only"), col=c(flr_pal["line"], bee_pal["line"], bee_pal["line"]),
       lwd=2, cex=1, inset=0.05, lty=c(1,1,2))
par(mfrow=c(1,1), mar=c(5,4,4,2), mgp=c(3,0.7,0))
dev.off()


# I still want to check whether interacting species have similar skews:
# for each bee species (at RMBL) calculate average skewness of associated flower species, weighted by relative frequency of visits
# and run bee_skew ~ associated_plant_skew ... that's it



### do interacting species have more similar skews?

bee_assoc <- bee_skews[dataset == "rmbl", .(n = .N,
                                            n_sites = length(unique(site)),
                                            n_years = length(unique(year)),
                                            mean = mean(mean),
                                            sd = mean(sd),
                                            skew = mean(skew),
                                            skew_se = mean(skew_se)),
                       by=species] # there are pretty few (28) bee species at rmbl for which I was able to get skew estimates
flr_assoc <- flr_skews[, .(n = .N,
                           n_sites = length(unique(site)),
                           n_years = length(unique(year)),
                           mean = mean(mean),
                           sd = mean(sd),
                           skew = mean(skew),
                           skew_se = mean(skew_se)),
                       by=species] # 142 flower species species

# reloading these because in data cleaning I removed netting data from rmbl_bees
rmbl_bees <- fread("raw_data/rmbl_bees/bees_2020-12-16.csv", stringsAsFactors = FALSE)
rmbl_bombus <- fread("raw_data/rmbl_bees/bombus_2020-12-16.csv", stringsAsFactors = FALSE)
rmbl_bombus$genus_species <- paste(rmbl_bombus$genus, rmbl_bombus$genus_species)

get.assoc <- function(bee_species, parameter="skew"){
  #bee_species <- "Bombus flavifrons"
  # bee_species <- "Lasioglossum sedi"
  # bee_species <- "Lasioglossum prasinogaster"
  # bee_species <- "Calliopsis teucrii"
  if(grepl("Bombus", bee_species)){
    flrs <- rmbl_bombus[genus_species == bee_species, plant_species]
  } else{
    flrs <- rmbl_bees[genus_species == bee_species & method %in% c("Net", "Net, PM", "Net, AM", "Net, Am"), details]
  }
  
  flr_table <- table(flrs)
  #flr_data <- flr_assoc[species %in% names(flr_table), ]
  flr_data <- flr_assoc[species %in% names(flr_table) | gsub(" .*", "", species) %in% gsub(" sp\\.", "", names(flr_table)), ] # second logic term mainly to add in "Crisium sp."
  if(nrow(flr_data) == 0) return(as.numeric(NA))
  flr_data[, weight := length(which(flrs == species | gsub(" sp.", "", flrs) == gsub(" .*", "", species))), by=species]
  flr_param <- weighted.mean(flr_data$skew, flr_data$weight, na.rm=T)
  return(flr_param)
}

bee_assoc[, flr_skew := get.assoc(species, "skew"), by=species]

png("figures/assoc_skew.png", width=700, height=600)
assoc_model <- lm(skew ~ flr_skew, data=bee_assoc)
summary(assoc_model)
plot(skew ~ flr_skew, data=bee_assoc, xlim=c(-0.1,0.47), ylim=c(-0.5,0.75),
     col=bee_pal["line"], pch=20, cex=1.2,
     xlab="Associated floral skew", ylab="Bee skew")
abline(h=0, v=0, lty=2)
abline(assoc_model, lwd=2, col="red")
label_y <- bee_assoc$skew; label_y[21] <- label_y[21]-0.02; label_y[4] <- label_y[4]+0.02; label_y[22] <- label_y[22]+0.02; label_y[27] <- label_y[27]+0.02
text(bee_assoc$flr_skew-0.005, label_y,
     labels = bee_assoc$species, cex=0.7,
     adj=1)
dev.off()


### Bombus and flower example distributions
bee_skews[group == "Bombus" & skew < - 0.7 & skew_se < .15, ]

# bees_female[species == "Bombus bifarius" & year == 2017 & site == "Lypps", ]
# 
# bee_ex <- bees_female[species == "Bombus bifarius" & year == 2017 & site == "Lypps", ]
# bee_ex <- bees_female[species == "Bombus bifarius" & year == 2018 & site == "Seans", ]
# bee_ex <- bees_female[species == "Bombus appositus" & year == 2015 & site == "Little", ]
bee_ex <- bees_female[species == "Bombus bifarius" & year == 2018 & site == "Seans", ]
fit.sn(bee_ex$ab, bee_ex$doy, plot=T, which=2)

plot(ab ~ doy, data=bee_ex)

bee_skews

### summary stats
table(bee_skews[dataset == "rmbl", site])
table(bee_skews[dataset == "mary", site])
table(flr_skews$site)

length(table(bee_skews[dataset == "rmbl", site]))
length(table(bee_skews[dataset == "mary", site]))
length(table(flr_skews$site))

# test <- bees[dataset == "mary", .(ab, doy, year), by = .(site, date, species)]
# table(test[site == "MD877a7fb7" & year == 2009, date])

print(paste(nrow(bee_skews), "bee time-series.",
            nrow(bee_skews[dataset == "rmbl",]), "from RMBL.",
            nrow(bee_skews[dataset == "mary",]), "from Mid-Atlantic"))

print(paste(nrow(flr_skews), "flower time-series.",
            length(unique(flr_skews$species)), "flower species"))

print(paste(length(unique(bee_skews$species)), "species.",
            length(unique(bee_skews[dataset == "rmbl", species])), "from RMBL.",
            length(unique(bee_skews[dataset == "mary", species])), "from Mid-Atlantic"))

print(paste(round((nrow(bee_skews[group == "Bombus",])/nrow(bee_skews))*100, 1), "% bombus total. ",
            perc_bombus[1], "% bombus on left. "))














### preliminary analysis

# bee_skews[mean == -7777]
# bee_skews <- bees[, fit.sn(bowl_adj, DOY), by=.(species, site, year)]
# bee_skews <- bee_skews[mean %!in% c(-7777, -8888, -9999),]
# bee_skews <- bee_skews[!is.na(skew_se),]
# bee_skews[skew > 0.9]
# bee_skews[skew < 0.5 & skew > -0.5]
# bees[site == "Beaver" & species == "Lasioglossum ruidosense", fit.sn(bowl_adj, DOY, plot=T)]
# bees[site == "Gothic" & species == "Panurginus cressoniellus", fit.sn(bowl_adj, DOY, plot=T)]
# bees[site == "Seans" & species == "Panurginus cressoniellus" & year == 2016, fit.sn(bowl_adj, DOY, plot=T)]
# 
# 
# hist(bee_skews$mean, 20)
# hist(bee_skews$sd, 20)
# hist(bee_skews$skew, 20)
# hist(bee_skews$skew_se, 20)
# 
# bees[1:6, fit.sn(bowl_adj, DOY)]
# bees[site == "Willey" & species == "Andrena algida", fit.sn(bowl_adj, DOY)]
# bees[site == "Willey" & species == "Andrena algida",]
# 
# 
# flower_skews <- flowers[, fit.sn(floralcount, doy), by=.(species, plot, year)] # waiting on this to run
# flower_skews <- flower_skews[mean %!in% c(-7777, -8888, -9999),]
# flower_skews <- flower_skews[!is.na(skew_se),]
# flower_skews[skew > 0.9]
# flowers[plot == "RM5" & species == "Erigeron flagellaris", fit.sn(floralcount, doy, plot=T, which=c(2))]
# 
# hist(flower_skews$mean, 20)
# hist(flower_skews$sd, 20)
# hist(flower_skews$skew, 20)
# hist(flower_skews$skew_se, 20)
# 
# 
# # prelim figs for Becky
# par(mfrow=c(2,1))
# bees[site == "Gothic" & species == "Panurginus cressoniellus", fit.sn(bowl_adj, DOY, plot=T, which=c(2), main="Panurginus cressoniellus at Gothic")]
# flowers[plot == "RM5" & species == "Erigeron flagellaris" & year == 2015, fit.sn(floralcount, doy, plot=T, which=c(2), main="Erigeron falgellaris at RM5")]
# 
# hist(bee_skews$skew, 20, main="Bee skews")
# hist(flower_skews$skew, 20, main="Flower skews")
# spar(mfrow=c(2,1))
# 
# # working example
# 
# ab <- c(0,0,10,20,35,50,20,5,0)
# ab <- c(0,0,10,120,65,50,20,5,1)
# time <- seq(110,190, by=10)
# ab_expand <- rep(time, times=ab)
# plot(ab ~ time)
# hist(ab_expand)
# 
# model <- selm(ab_expand ~ 1)
# summary(model)
# plot(model)
# summ <- summary(model)
# summ@param.table
# 
# predict_ab <- predict(model)
# predict_ab <- predict(model, data.frame(time = 110:190), param.type="CP", param.name="gamma1")
# 
# profile(model, param.name="cp")
# profile(model)     
# 
# 
# test <- try(sum(c(2,NA,"a")), silent=T)
# 
# tryCatch(sum(c(2,4,"a")), error=function(e) "hello")
# 
# tryCatch(1, finally = print("Hello"))
# e <- simpleError("test error")
# ## Not run: 
# stop(e)
# tryCatch(stop(e), finally = print("Hello"))
# tryCatch(stop("fred"), finally = print("Hello"))
# 
# ## End(Not run)
# tryCatch(stop(e), error = function(e) e, finally = print("Hello"))
# 
# ### does population size bias skew estimates?
# library(sn)
# big_pop <- rsn(1000, alpha=5)
# little_pop <- rsn(100, alpha=5)
# hist(big_pop, 100, col=rgb(1,0,0,0.5))
# hist(little_pop, 100, col=rgb(0,1,0,0.5), add=T)
# model_big <- selm(big_pop ~ 1)
# model_little <- selm(little_pop ~ 1)
# summary(model_big)
# summary(model_little)
# 
# 
# 
# # looks like no, but it smaller samples decrease accuracy of estimates, not surprisingly
# 
# 
# 
# 
# 
# ### trying to understand omega, alpha, and gamma
# 
# data <- rsn(1000,100,10,10)
# hist(data, 100)
# 
# # omega is just SD
# # alpha is skew
# # gamma is also skew...?
# 
# # I think alpha and gamma have to do with whether CP or DP is used
# # yep, you can look at both the CP and DP estimates with"
# # summary(model, "dp")
# # summary(model, "cp")
# # I don't think it matters which I use for this project
# 
