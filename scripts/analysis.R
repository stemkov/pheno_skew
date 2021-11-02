# analyzing bee and flower data for phenological skew

library(sn)
library(data.table)
library(tools)
library(mgsub)
library(lubridate)
library(visreg)
library(plotly)
library(pbapply)
library(moments)
library(sfsmisc)
library(fields)
library(climateR)
library(ncdf4)

`%!in%` = Negate(`%in%`)

setwd("/home/michael/Documents/Grad School/Research Projects/pheno_skew")

if(!(exists("bees") | exists("flowers"))){
  bees <- fread("clean_data/bees.csv", stringsAsFactors = FALSE)
  flowers <- fread("clean_data/flowers.csv", stringsAsFactors = FALSE)
}

fit.sn <- function(abs, times, method="sn", family="SN", plot=FALSE, ...){
  # abs <- c(0,0,1,10,4,3,2,0,1)
  # times <- seq(110,190, by=10)
  #abs <- abs*10
  if(method %!in% c("sn", "fGarch", "moments")) stop("unknown method")
  if(sum(abs) == 0 | length(abs) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999, skew_t=-7777, skew_p=-9999))
  #if(length(unique(abs)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999))
  if(length(unique(times)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999, skew_t=-7777, skew_p=-9999))
  if(sum(abs) < 10) return(list(mean=-8888, sd=-8888, skew=-8888, skew_se=-8888, skew_t=-7777, skew_p=-8888))
  
  abs_expand <- rep(times, times=round(abs))
  if(method == "sn"){
    #model <- selm(abs_expand ~ 1)
    model <- tryCatch(selm(abs_expand ~ 1, family=family), error=function(e) "error") # catching errors so it doesn't abort whole analysis
    if(class(model)[1] != "selm") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777, skew_t=-7777, skew_p=-7777))
    summ <- summary(model)
    pars <- summ@param.table
    mean <- pars[1,1]
    sd <- pars[2,1]
    skew <- pars[3,1]
    skew_se <- pars[3,2]
    skew_t <- pars[3,3]
    skew_p <- pars[3,4]
    if(plot) plot(model, ...)
  }
  if(method == "fGarch"){
    #model <- snormFit(abs_expand)
    model <- tryCatch(snormFit(abs_expand), error=function(e) "error") # catching errors so it doesn't abort whole analysis
    if(class(model)[1] != "list") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777, skew_t=-7777, skew_p=-7777))
    pars <- model$par
    mean <- pars["mean"]
    sd <- pars["sd"]
    skew <- pars["xi"]
    skew_se <- -6666
    skew_t <- -6666
    skew_p <- -6666
    if(plot){
      doy_seq <- seq(range(times)[1], range(times)[2], length.out=100)
      d <- dsnorm(doy_seq, mean, sd, skew)
      hist(abs_expand, xlim=range(doy_seq))
      par(new=T)
      plot(d ~ doy_seq, type="l", axes=F, bty="n", xlab="", ylab="")
    }
  }
  if(method == "moments"){
    #skew <- moments::skewness(abs_expand)
    #skew_test <- agostino.test(abs_expand)
    skew_test <- tryCatch(agostino.test(abs_expand), error=function(e) "error") # catching errors so it doesn't abort whole analysis
    if(class(skew_test) != "htest") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777, skew_t=-7777, skew_p=-7777))
    skew <- skew_test$statistic[1]
    mean <- mean(abs_expand)
    sd <- sd(abs_expand)
    skew_se <- -6666
    skew_t <- -6666
    skew_p <- skew_test$p.value
    if(skew > 4 | skew < -4) return(list(mean=-5555, sd=-5555, skew=-5555, skew_se=-5555, skew_t=-5555, skew_p=-5555))
    if(plot) hist(abs_expand)
  }
  # I probably want to report total pop size, unique time points, and some goodness of fit?
  return(list(mean=mean, sd=sd, skew=skew, skew_se=skew_se, skew_t=skew_t, skew_p=skew_p))
}

get.overlap <- function(dist1, dist2, x){
  #sum(pmin(dist1, dist2))/ sum(dist1) # only works when I don't scale
  dist1_area <- integrate.xy(x, dist1)
  dist2_area <- integrate.xy(x, dist2)
  total <- dist1_area + dist2_area
  intersection <- integrate.xy(x, pmin(dist1, dist2))
  overlap <- (2*intersection) / total
  return(overlap)
} 

maximize.overlap <- function(skew1, skew2=0, x = seq(60,140,by=0.2), mean = 100, sd = 10, plot=F){
  if(FALSE){
    x <- seq(60,140,by=0.2)
    mean <- 100
    sd <- 10
    skew1 <- 10
    skew2 <- -0.2
    skew1 <- -4
    skew2 <- 4
    plot=T
  }
  sn1 <- dsn(x, mean, sd, skew1, 0)
  sn2 <- dsn(x, mean, sd, skew2, 0)
  if(plot){plot(sn1 ~ x, type="l", lwd=3); lines(sn2 ~ x, lwd=3, col="blue")} 
  
  middle_x <- seq(round(length(x)/2/2), round(length(x)/2*1.5))
  shifts <- x[middle_x] # only need this range for skew -5 to 5
  scales <- seq(0.4,1,length.out=10) # only need this range for skew -5 to 5
  transforms <- expand.grid(shifts, scales)
  transforms <- data.table(shift = transforms$Var1,
                           scale = transforms$Var2)
  
  if(abs(skew1) > abs(skew2)){
    transforms[, overlap := get.overlap(sn1, dsn(x,shift,sd*scale,skew2), x), by=.(shift,scale)]
    which_maximized <- transforms[which.max(overlap),]
    if(plot) lines(dsn(x,which_maximized$shift, sd*which_maximized$scale, skew2) ~ x, lty=2, col="blue", lwd=2)
  } else if (abs(skew2) > abs(skew1)){
    transforms[, overlap := get.overlap(dsn(x,shift,sd*scale,skew1), sn2, x), by=.(shift,scale)]
    which_maximized <- transforms[which.max(overlap),]
    if(plot) lines(dsn(x,which_maximized$shift, sd*which_maximized$scale, skew1) ~ x, lty=2, col="black", lwd=2)
  } else if(skew1 == skew2){
    return(1)
  } else{
    transforms[, overlap := get.overlap(dsn(x,shift,sd,skew1), sn2, x), by=.(shift)]
    which_maximized <- transforms[which.max(overlap),]
    if(plot) lines(dsn(x,which_maximized$shift, sd, skew1) ~ x, lty=2, col="black", lwd=2)
  }
  
  return(max(transforms$overlap))
}

# to plot distributions and their overlap
plot.overlap <- function(x, mean1=100, mean2=100, sd1=10, sd2=10, skew1=0, skew2=0,
                         col1 = flr_pal["line"], col2 = bee_pal["line"], col1_shade = flr_pal["point"], col2_shade = bee_pal["point"]){
  if(FALSE){
    x <- c(70:130)
    mean1 <- 100; mean2 <- 110
    sd1 <- 10; sd2 <- 10
    skew1 <- 0; skew2 <- 0.5
  }
  dist1 <- dsn(x, mean1, sd1, skew1, 0)
  dist2 <- dsn(x, mean2, sd2, skew2, 0)
  intersection <- pmin(dist1, dist2)
  
  plot(dist1 ~ x, type="l", lwd=3, col=col1, xaxt="n", yaxt="n", ylab="", xlab="")
  lines(dist2 ~ x, lwd=3, col=col2)
  #lines(intersection ~ x, lwd=4, col=col_inter_line)
  polygon(x = c(x, min(x)),
          y = c(intersection, intersection[1]),
          col=col1_shade, border=NA, density=10, angle=45)
  polygon(x = c(x, min(x)),
          y = c(intersection, intersection[1]),
          col=col2_shade, border=NA, density=10, angle=-45)
}

### Structure for the analysis:
# - 1. show that skewness is important for phenological overlap
# - 2. skews for individual time-series. for each year/species/site timeseries
# - -  Demonstrate average skewness and how it varies with mean and SD
# - 3. aggregate across years
# - -  does skew vary between years? ... maybe predicted by drought index?
# - 4. how much loss of overlap does skewness account for?
# - -  compared to if all distributions were assumed to be normal
# - 5ish. are skews predicted by strength of interaction? (might not do this)


### setting up color palettes
bee_pal <- c(point = "#f09c00", line = "#c48000")
flr_pal <- c(point = "#a655cf", line = "#7b04b8")

### full model with multiple predictors - this takes really long time to run

# bees[, species_site_year := paste(species,site,year, sep="_")] # making species/site/year "dummy" variable for full model
# 
# bees_summary <- bees[, .(total_ab = sum(ab),
#                          n_days = length(ab)), by=.(species_site_year)]
# hist(bees_summary$total_ab, 100)
# hist(bees_summary$n_days)
# bees_to_cut <- bees_summary[total_ab < 100 | n_days < 3, species_site_year]
# bees_sub <- bees[sex == "F" & species_site_year %!in% bees_to_cut,]
# bees_expand <- bees_sub[rep(1:.N,ab)][,index:=1:.N,by=doy]
# bees_expand$year <- as.factor(bees_expand$year)
# bees_expand$species <- as.factor(bees_expand$species)
# bees_expand$species_site_year <- as.factor(bees_expand$species_site_year)
# 
# full_lm <- lm(doy ~ species_site_year, data = bees_expand)
# summary(full_lm)
# 
# full_sn <- selm(doy ~ species_site_year, data=bees_expand)
# summary(full_sn)

# this isn't very satisfying - it just provides one skewness estimate

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
bees_sub <- bees[sex == "F",] # getting rid of male and ambiguous timeseries

n_combs <- uniqueN(bees_sub[,c("site", "species", "year", "dataset")])
bee_skews <- bees_sub[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy, method="moments")}, by=.(species, site, year, dataset)]
print(paste0(round(nrow(bee_skews[mean == -9999,])/nrow(bee_skews)*100, 1), "% of time-series too short"))
print(paste0(round(nrow(bee_skews[mean == -8888,])/nrow(bee_skews)*100, 1), "% too little catch"))
print(paste0(round(nrow(bee_skews[mean == -7777,])/nrow(bee_skews)*100, 1), "% error in fitting model"))
print(paste0(round(nrow(bee_skews[is.na(skew_se),])/nrow(bee_skews)*100, 1), "% likely stopped sampling too early/late"))
print(paste0(round(nrow(bee_skews[mean == -5555,])/nrow(bee_skews)*100, 1), "% extreme skew value"))

bee_skews <- bee_skews[mean %!in% c(-5555, -7777, -8888, -9999),]
bee_skews <- bee_skews[!is.na(skew_se),]


### FLOWERS each species/site/year time-series individually
n_combs <- uniqueN(flowers[,c("site", "species", "year")])
flr_skews <- flowers[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy, method="moments")}, by=.(species, site, year)]
print(paste0(round(nrow(flr_skews[mean == -9999,])/nrow(flr_skews)*100, 1), "% of time-series too short"))
print(paste0(round(nrow(flr_skews[mean == -8888,])/nrow(flr_skews)*100, 1), "% too little catch"))
print(paste0(round(nrow(flr_skews[mean == -7777,])/nrow(flr_skews)*100, 1), "% error in fitting model"))
print(paste0(round(nrow(flr_skews[is.na(skew_se),])/nrow(flr_skews)*100, 1), "% likely stopped sampling too early/late"))
print(paste0(round(nrow(flr_skews[mean == -5555,])/nrow(flr_skews)*100, 1), "% extreme skew value"))

flr_skews <- flr_skews[mean %!in% c(-5555, -7777, -8888, -9999),]
flr_skews <- flr_skews[!is.na(skew_se),]


### conceptual demonstration of mean, SD, and skew effects on overlap


png("figures/conceptual.png", width=900, height=300)

layout(matrix(c(1,2,3), nrow=1))
par(mar=c(4,3.5,2,1))

x_vals <- c(60:140)
plot.overlap(x_vals, mean1 = 100, mean2 = 110)
mtext("Abundance", 2, 2)
mtext("  Mean change",3, -2, adj=0)

par(mar=c(4,2,2,1))
plot.overlap(x_vals, sd1 = 10, sd2 = 15)
mtext("Time", 1, 2)
mtext("  Breadth change",3, -2, adj=0)

par(mar=c(4,2,2,2))
plot.overlap(x_vals, skew1=0, skew2=-5, mean2 = 106, sd2 = 18)
mtext("  Skewness \n  change",3, -3.5, adj=0)

legend("topright", c("Flower", "Bee"), col=c(flr_pal["line"], bee_pal["line"]),
       lwd=3, cex=1.5, inset=0.05)

dev.off()


# combining bee and flower skews
flr_skews$dataset <- "rmbl"
bee_skews$guild <- "bee"
flr_skews$guild <- "flower"
skews <- rbind(flr_skews, bee_skews)
skews$guild <- relevel(as.factor(skews$guild), "flower") # to put bee points on top


### Skewness histograms
png("figures/skewness.png", width=600, height=800)
par(mfrow=c(1,1), mar=c(5,4,1,3.5))
par(fig = c(0,1,0,1))

hist(skews[guild=="flower", skew], 100,
     main="", xlab="Skewness", ylab="Frequency",
     border=F, col=flr_pal["line"],
     xlim=c(-4,4), ylim=c(-500,700),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,700, by=100), labels=seq(0,700, by=100), tick=T, las=2)
hist(skews[guild=="flower" & skew_p > 0.05, skew], 30,
     main="", xlab="Skewness", ylab="Frequency",
     border=F, col="white",
     yaxt="n", xaxt="n", add=T,
     density=15, angle=45)

par(new=T)
hist(skews[guild=="bee", skew], 100,
     main="", xlab="", border=F,
     col=bee_pal["point"],
     ylim=c(0,100), xlim=c(-4,4),
     yaxt="n", ylab=NA, xaxt="n")
axis(2, at=seq(0,30, by=10), labels=seq(0,30, by=10), tick=T, las=2)
hist(skews[guild=="bee" & skew_p > 0.05, skew], 30,
     main="", xlab="Skewness", ylab="Frequency",
     border=F, col="white",
     yaxt="n", xaxt="n", add=T,
     density=15, angle=45)

axis(1, at=seq(-4,4,by=1), labels=seq(-4,4,by=1), tick=T)
abline(v=0)

legend("topright", c("Flowers", "Bees"), fill=c(flr_pal["line"], bee_pal["point"]),
       density=c(NA,NA), border=c("white","white","white"), cex=1, inset=0.05)

flr_prec_left <- (skews[guild=="flower" & skew_p < 0.05 & skew < 0, .N] / skews[guild=="flower", .N]) * 100
flr_prec_right <- (skews[guild=="flower" & skew_p < 0.05 & skew > 0, .N] / skews[guild=="flower", .N]) * 100
bee_prec_left <- (skews[guild=="bee" & skew_p < 0.05 & skew < 0, .N] / skews[guild=="bee", .N]) * 100
bee_prec_right <- (skews[guild=="bee" & skew_p < 0.05 & skew > 0, .N] / skews[guild=="bee", .N]) * 100

par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=-0.85, percent=round(flr_prec_left, 1), line_col = flr_pal["line"], fill_col = flr_pal["point"], text_cex = 1,
         x1=0.05,x2=0.40,y1=0.5,y2=0.7)
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0.85, percent=round(flr_prec_right, 1), line_col = flr_pal["line"], fill_col = flr_pal["point"], text_cex = 1,
         x1=0.6,x2=0.95,y1=0.5,y2=0.7)
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=-0.85, percent=round(bee_prec_left, 1), line_col = bee_pal["line"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.05,x2=0.40,y1=0.20,y2=0.4)
par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
skew.fig(skew=0.85, percent=round(bee_prec_right, 1), line_col = bee_pal["line"], fill_col = bee_pal["point"], text_cex = 1,
         x1=0.6,x2=0.95,y1=0.20,y2=0.4)

dev.off()

### Joint mean and SD figure
png("figures/mean_sd.png", width=1000, height=550)

par(mfrow=c(1,2), mar=c(4.1,4,2,1))
mean_model <- lm(skew ~ mean*guild, data = skews)
#plot(skew ~ mean, data = flr_skews_just)
summary(mean_model)
visreg(mean_model, "mean", by="guild", xlab="Mean phenology (day of year)", ylab="Skewness", overlay=T,
       points = list(col=c(flr_pal["point"], bee_pal["line"]),
                     cex=c(0.2,0.7),
                     pch=c(19,21),
                     bg=c(1,bee_pal["point"])),
       line.par = list(col=c(flr_pal["line"], bee_pal["line"])),
       fill.par = list(col=adjustcolor(c(flr_pal["point"], bee_pal["point"]), 0.5)),
       legend=F)
abline(h=0)
box()

sd_model <- lm(abs(skew) ~ sd*guild, data = skews)
summary(sd_model)
visreg(sd_model, "sd", by="guild", xlab="Phenological breadth (days)", ylab="Absolute skewness", overlay=T,
       points = list(col=c(flr_pal["point"], bee_pal["line"]),
                     cex=c(0.2,0.7),
                     pch=c(19,21),
                     bg=c(1,bee_pal["point"])),
       line.par = list(col=c(NA, bee_pal["line"])),
       fill.par = list(col=adjustcolor(c(NA, bee_pal["point"]), 0.5)),
       legend=F)
flr_line <- as.data.table(visreg(sd_model, "sd", by="guild", plot=F)$fit)
flr_line <- flr_line[guild == "flower" & sd < max(skews[guild == "flower", sd]),] # cutting out flower line beyond data range 
lines(flr_line$visregFit ~ flr_line$sd, lwd=3, col=flr_pal["line"])
polygon(x = c(flr_line$sd, rev(flr_line$sd)),
        y = c(flr_line$visregLwr, rev(flr_line$visregUpr)),
        col=adjustcolor(flr_pal["point"], 0.5), border = NA)
legend("topright", c("Flowers", "Bees"), col=c(flr_pal["line"], bee_pal["line"]),
       lwd=3, cex=1, inset=0.05)
box()
par(mfrow=c(1,1), mar=c(5,4,4,2))
dev.off()


### demonstration of skewness effect on overlap
n <- 50
skew_vals <- seq(-5,5, length.out=n)
skew_scan <- as.data.table(expand.grid(skew_vals, skew_vals))
colnames(skew_scan) <- c("skew1", "skew2")

if(!file.exists("clean_data/skew_scan.RDS")){
  overlaps_vals <- pbmapply(maximize.overlap, skew1 = skew_scan$skew1, skew2 = skew_scan$skew2)
  skew_scan$overlap <- overlaps_vals
  saveRDS(skew_scan,"clean_data/skew_scan.RDS")
} else{
  skew_scan <- readRDS("clean_data/skew_scan.RDS")
}

# getting flr and bee quantiles for plotting
flr_low <- skews[guild == "flower", quantile(skew, 0.025)]
flr_high <- skews[guild == "flower", quantile(skew, 0.975)]
bee_low <- skews[guild == "bee", quantile(skew, 0.025)]
bee_high <- skews[guild == "bee", quantile(skew, 0.975)]

# heatmap
png("figures/heatmap.png", width=700, height=600)

skew_scan_mat <- matrix(skew_scan$overlap * 100, nrow = n, ncol = n)

par(mar=c(4.5,4.2,2,3))
image.plot(x = skew_vals, y = skew_vals, z = skew_scan_mat,
      col=colorRampPalette(c("red","orange","white"))(100),
      xlab = "Flower skewness", ylab="Bee skewness",
      legend.shrink = 0.75)

abline(v=flr_low, col=flr_pal["line"], lwd=2, lty=3)
abline(v=flr_high, col=flr_pal["line"], lwd=2, lty=3)
abline(h=bee_low, col=bee_pal["line"], lwd=2, lty=3)
abline(h=bee_high, col=bee_pal["line"], lwd=2, lty=3)
lines(x = c(flr_low, flr_low),
      y = c(bee_low, bee_high),
      col = flr_pal["line"], lwd=3)
lines(x = c(flr_high, flr_high),
      y = c(bee_low, bee_high),
      col = flr_pal["line"], lwd=3)
lines(x = c(flr_low, flr_high),
      y = c(bee_low, bee_low),
      col = bee_pal["line"], lwd=3)
lines(x = c(flr_low, flr_high),
      y = c(bee_high, bee_high),
      col = bee_pal["line"], lwd=3)

mtext("Maximum possible overlap (%)",
      side = 4, line = 0.75)

dev.off()


### Does drought predict skewness?

# getting drought index data
nc <- nc_open("raw_data/terraclim.nc")
pdsi_raw <- ncvar_get(nc, varid="PDSI")
pdsi_time <- ncvar_get(nc, varid="time")
nc_close(nc)

pdsi <- data.table(drought = pdsi_raw, days_since_1900 = pdsi_time)

pdsi[, date := as.Date(days_since_1900, origin = "1900-01-01")]
pdsi[, year := year(date)]
pdsi[, month := month(date)]

drought_data <- pdsi[, .(yearly_drought = mean(drought),
                         summer_drought = mean(.SD[month < 9 & month > 4,drought])),
                     by=.(year)]
plot(summer_drought ~ yearly_drought, data=drought_data)

test <- merge(skews, drought_data, by="year")
test <- test[dataset != "mary",]

test_lm <- lm(skew ~ summer_drought*guild + site, data = test)
summary(test_lm)
visreg(test_lm, "summer_drought", by="guild")

