### Analysis of skewness in bee and flower phenological distributions
# Michael Stemkovski
# m.stemkovski@gmail.com
#
# run cleaning.R before this script

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
library(flextable)
library(taxize)

`%!in%` = Negate(`%in%`)

setwd("/home/michael/Documents/Grad School/Research Projects/pheno_skew")

if(!(exists("bees") | exists("flowers"))){
  bees <- fread("clean_data/bees.csv", stringsAsFactors = FALSE)
  flowers <- fread("clean_data/flowers.csv", stringsAsFactors = FALSE)
}

# core function for calculating skewness
# multiple methods are supported, but I just used moments in the final analysis
# the name fit.sn is outdated because I started off fitting skew-normal distributions
fit.sn <- function(abs, times, method="sn", family="SN", plot=FALSE, min_ab=10, ...){
  # abs <- c(0,0,1,10,4,3,2,0,1)
  # times <- seq(110,190, by=10)
  #abs <- abs*10
  if(method %!in% c("sn", "fGarch", "moments")) stop("unknown method")
  if(sum(abs) == 0 | length(abs) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999, skew_t=-7777, skew_p=-9999, n=-9999))
  #if(length(unique(abs)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999))
  if(length(unique(times)) < 3) return(list(mean=-9999, sd=-9999, skew=-9999, skew_se=-9999, skew_t=-7777, skew_p=-9999, n=-9999))
  if(sum(abs) < min_ab) return(list(mean=-8888, sd=-8888, skew=-8888, skew_se=-8888, skew_t=-7777, skew_p=-8888, n=as.double(sum(abs))))
  
  abs_expand <- rep(times, times=round(abs))
  if(method == "sn"){
    #model <- selm(abs_expand ~ 1)
    model <- tryCatch(selm(abs_expand ~ 1, family=family), error=function(e) "error") # catching errors so it doesn't abort whole analysis
    if(class(model)[1] != "selm") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777, skew_t=-7777, skew_p=-7777, n=-7777))
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
    if(class(model)[1] != "list") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777, skew_t=-7777, skew_p=-7777, n=-7777))
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
    if(class(skew_test) != "htest") return(list(mean=-7777, sd=-7777, skew=-7777, skew_se=-7777, skew_t=-7777, skew_p=-7777, n=-7777))
    skew <- skew_test$statistic[1]
    mean <- mean(abs_expand)
    sd <- sd(abs_expand)
    skew_se <- -6666
    skew_t <- -6666
    skew_p <- skew_test$p.value
    if(skew > 5 | skew < -5) return(list(mean=-5555, sd=-5555, skew=-5555, skew_se=-5555, skew_t=-5555, skew_p=-5555, n=-5555))
    if(plot) hist(abs_expand)
  }
  
  n <- as.double(sum(abs))
  
  return(list(mean=mean, sd=sd, skew=skew, skew_se=skew_se, skew_t=skew_t, skew_p=skew_p, n=n))
}

if(!exists("rmbl_bees_w_zeros")){rmbl_bees_w_zeros <- fread("clean_data/rmbl_bees_w_zeros.csv", stringsAsFactors = FALSE)}
if(!exists("flowers_w_zeros")){flowers_w_zeros <- fread("clean_data/flowers_w_zeros.csv", stringsAsFactors = FALSE)}
# check for whether there is censoring (lack of zeros before or after non-zeros)
# returns either "before", "after", "both", or "none" for where the censoring is occurring
is.ab.censored <- function(abs){
  zero_positions <- which(abs == 0)
  #non_zero_positions <- which(abs != 0)
  
  if(length(zero_positions) == 0) return("both")
  if(1 %!in% zero_positions & length(abs) %!in% zero_positions) return("both")
  if(1 %!in% zero_positions) return("before")
  if(length(abs) %!in% zero_positions) return("after")
  
  return("none")
}

# checks if TS for site/year/species is censored - uses is.ab.censored
check.if.censored <- function(site_v, year_v, species_v, sex_v="F", dataset_v="bees"){
  if(dataset_v=="flowers") rel_abs <- flowers_w_zeros[site == site_v & year == year_v & species == species_v, ab]
  if(dataset_v=="bees") rel_abs <- rmbl_bees_w_zeros[site == site_v & year == year_v & species == species_v & sex == sex_v, ab]
  if(length(rel_abs) == 0) stop("no records found - wrong inputs")
  return(is.ab.censored(rel_abs))
}
# check.if.censored("Beaver", 2015, "Dufourea harveyi")
# check.if.censored("Beaver", 2015, "Bombus bifarius Q")
# check.if.censored("Beaver", 2015, "Lasioglossum sedi")
# plot(flowers_w_zeros[species == "Claytonia lanceolata" & year == 2015 & site == "RM_a", ab])
# check.if.censored("RM_a", 2015, "Claytonia lanceolata", dataset_v = "flowers")
# plot(flowers_w_zeros[species == "Taraxacum officinale" & year == 2015 & site == "RM_a", ab])
# check.if.censored("RM_a", 2015, "Taraxacum officinale", dataset_v = "flowers")

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

maximize.overlap(-5,5, plot=T)

# to plot distributions and their overlap
plot.overlap <- function(x, mean1=100, mean2=100, sd1=10, sd2=10, skew1=0, skew2=0,
                         col1 = flr_pal["line"], col2 = bee_pal["line"], col1_shade = flr_pal["point"], col2_shade = bee_pal["point"]){
  if(FALSE){
    x <- c(70:130)
    mean1 <- 100; mean2 <- 100
    sd1 <- 10; sd2 <- 10
    skew1 <- 0; skew2 <- -5
    col1 = flr_pal["line"]; col2 = bee_pal["line"]; col1_shade = flr_pal["point"]; col2_shade = bee_pal["point"]
  }
  dist1 <- dsn(x, mean1, sd1, skew1, 0)
  dist2 <- dsn(x, mean2, sd2, skew2, 0)
  intersection <- pmin(dist1, dist2)
  
  plot(dist1 ~ x, type="l", lwd=3, col=col1, xaxt="n", yaxt="n", ylab="", xlab="")
  lines(dist2 ~ x, lwd=3, col=col2)
  #lines(intersection ~ x, lwd=4, col=col_inter_line)
  polygon(x = c(x, min(x)),
          y = c(intersection, intersection[1]),
          col=col1_shade, border=NA, density=12, angle=45)
  polygon(x = c(x, min(x)),
          y = c(intersection, intersection[1]),
          col=col2_shade, border=NA, density=12, angle=-45)
  
  print(paste(mean(rep(x,round(dist1*10000))),",", mean(rep(x,round(dist2*10000)))))
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
skew.fig <- function(skew = 0, percent = NA, x1=0, x2=0.5, y1=0, y2=0.5, line_col="red", fill_col = "pink", text_cex=2, shade_percent=0, reset_window=FALSE){
  # line_col <- "red"
  # fill_col <- "pink"
  # skew <- 0.8
  # percent <- 75
  # text_cex <- 2
  # shade_percent <- 50
  
  #par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
  if(reset_window) par(fig = c(0,1,0,1), mar=c(5,4,4,2)); print("nice")
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
bees_sub <- bees_sub[dataset == "rmbl",] # cutting out Mid-Atlantic bees for sake of comparability and simplicity

n_combs <- uniqueN(bees_sub[,c("site", "species", "year", "dataset")])
bee_skews <- bees_sub[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy, method="moments")}, by=.(species, site, year, dataset)]
bees_total_n <- bee_skews[mean >= 0, sum(n)]
bees_below_minimum_n <- bee_skews[mean == -8888, sum(n)]
print(paste0(round(nrow(bee_skews[mean == -9999,])/nrow(bee_skews)*100, 1), "% of time-series too short ", bee_skews[mean == -9999,.N]))
print(paste0(round(nrow(bee_skews[mean == -8888,])/nrow(bee_skews)*100, 1), "% too little catch ", bee_skews[mean == -8888,.N]))
print(paste0(round(nrow(bee_skews[mean == -7777,])/nrow(bee_skews)*100, 1), "% error in fitting model ", bee_skews[mean == -7777,.N]))
print(paste0(round(nrow(bee_skews[is.na(skew_se),])/nrow(bee_skews)*100, 1), "% likely stopped sampling too early/late ", bee_skews[is.na(skew_se),.N]))
print(paste0(round(nrow(bee_skews[mean == -5555,])/nrow(bee_skews)*100, 1), "% extreme skew value ", bee_skews[mean == -5555,.N]))

bee_skews <- bee_skews[mean %!in% c(-5555, -7777, -8888, -9999),]
bee_skews <- bee_skews[!is.na(skew_se),]

#bee_skews_rmbl_mary <- bee_skews # saving mid-atlantic bees for supplement figure
#bee_skews <- bee_skews[dataset == "rmbl",] # cutting out mid-atlantic bees for main analysis


### FLOWERS each species/site/year time-series individually
n_combs <- uniqueN(flowers[,c("site", "species", "year")])
#flr_skews <- flowers[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy, method="moments")}, by=.(species, site, year)]
flr_skews <- flowers[, {cat(.GRP/n_combs*100,"%\n"); fit.sn(ab, doy, method="moments", min_ab = 100)}, by=.(species, site, year)]

flrs_total_n <- flr_skews[mean >= 0, sum(n)]
flrs_below_minimum_n <- flr_skews[mean == -8888, sum(n)]
print(paste0(round(nrow(flr_skews[mean == -9999,])/nrow(flr_skews)*100, 1), "% of time-series too short ", flr_skews[mean == -9999,.N]))
print(paste0(round(nrow(flr_skews[mean == -8888,])/nrow(flr_skews)*100, 1), "% too little catch ", flr_skews[mean == -8888,.N]))
print(paste0(round(nrow(flr_skews[mean == -7777,])/nrow(flr_skews)*100, 1), "% error in fitting model ", flr_skews[mean == -7777,.N]))
print(paste0(round(nrow(flr_skews[is.na(skew_se),])/nrow(flr_skews)*100, 1), "% likely stopped sampling too early/late ", bee_skews[is.na(skew_se),.N]))
print(paste0(round(nrow(flr_skews[mean == -5555,])/nrow(flr_skews)*100, 1), "% extreme skew value ", flr_skews[mean == -5555,.N]))

flr_skews <- flr_skews[mean %!in% c(-5555, -7777, -8888, -9999),]
flr_skews <- flr_skews[!is.na(skew_se),]



### checking for truncation
bee_skews[, censored := check.if.censored(site, year, species, dataset_v="bees"), by=.(site, year, species)]
table(bee_skews$censored)
sum(table(bee_skews$censored))

n_combs <- uniqueN(flr_skews[,c("site", "species", "year")])
flr_skews[, censored := {cat(.GRP/n_combs*100,"%\n"); check.if.censored(site, year, species, dataset_v="flowers")}, by=.(site, year, species)]
table(flr_skews$censored)
sum(table(flr_skews$censored))


### combining bee and flower skews
flr_skews$dataset <- "rmbl"
bee_skews$guild <- "bee"
flr_skews$guild <- "flower"
skews <- rbind(flr_skews, bee_skews)
skews$guild <- relevel(as.factor(skews$guild), "flower") # to put bee points on top
skews[, genus := gsub(" .*", "", species)]

### Summary statistics
table(skews$guild)
table(skews[,c("guild", "dataset")])
paste("# flower species:",length(unique(skews[guild == "flower",species])), "from", length(unique(skews[guild == "flower",genus])), "genera") 
paste("# bee species:",length(unique(skews[guild == "bee",species])), "from", length(unique(skews[guild == "bee",genus])), "genera")

overall_model <- lm(skew ~ guild, data=skews)
summary(overall_model)

paste("mean date flowers:", gsub("2000-","",as.Date(round(mean(skews[guild == "flower", mean])), origin="2000-01-01")))
paste("mean date bees:", gsub("2000-","",as.Date(round(mean(skews[guild == "bee", mean])), origin="2000-01-01")))

paste("mean breadth within 1SD flowers:", round(mean(skews[guild == "flower", sd])*2))
paste("mean breadth within 1SD bees:", round(mean(skews[guild == "bee", sd])*2))

paste("mean skew flowers:",round(mean(skews[guild == "flower", skew]),3))
paste("mean skew bees:",round(mean(skews[guild == "bee", skew]),3))

paste(round((skews[guild=="flower" & skew_p < 0.05 & skew < 0, .N] / skews[guild=="flower", .N]) * 100,2), "% flowers left-skewed")
paste(round((skews[guild=="flower" & skew_p < 0.05 & skew > 0, .N] / skews[guild=="flower", .N]) * 100,2), "% flowers right-skewed")
paste(round((skews[guild=="flower" & skew_p >= 0.05, .N] / skews[guild=="flower", .N]) * 100,2), "% flowers symmetrical")
paste(round((skews[guild=="bee" & skew_p < 0.05 & skew < 0, .N] / skews[guild=="bee", .N]) * 100,2), "% bees left-skewed")
paste(round((skews[guild=="bee" & skew_p < 0.05 & skew > 0, .N] / skews[guild=="bee", .N]) * 100,2), "% bees right-skewed")
paste(round((skews[guild=="bee" & skew_p >= 0.05, .N] / skews[guild=="bee", .N]) * 100,2), "% bees symmetrical")

paste(bees_below_minimum_n, "of", bees_total_n, "bee records removed because total N below 10.", round(bees_below_minimum_n/bees_total_n,4)*100, "%")
paste(flrs_below_minimum_n, "of", flrs_total_n, "flowers removed because total N below 100.", round(flrs_below_minimum_n/flrs_total_n,4)*100, "%")


### species-specific examples

bee_skews[, .(mean =  gsub("2000-","",as.Date(round(mean(mean)), origin="2000-01-01")),
              skew = mean(skew),
              n = .N,
              p = median(skew_p) < 0.05), by=.(species)][order(mean)]


print(flr_skews[, .(mean = gsub("2000-","",as.Date(round(mean(mean)), origin="2000-01-01")),
                    skew = mean(skew),
                    n = .N,
                    p = median(skew_p) < 0.05), by=.(species)][order(mean)], nrows=135)

### comparison with Thomson 1980 and Rabinowitz 1981
thomson <- fread("raw_data/thomson_skews.csv")
rabinowitz <- fread("raw_data/rabinowitz_skews.csv")
mean(thomson$skewness)
mean(rabinowitz$skewness)
mean(flr_skews$skew)


### conceptual demonstration of mean, SD, and skew effects on overlap
#png("figures/conceptual.png", width=3, height=7, units = "in", res=300, pointsize = 9)
svg("figures/conceptual.svg", width=3, height=7)

  layout(matrix(c(1,2,3), nrow=3),
         heights=c(0.9, 0.95, 1.1))
  par(mar=c(0,3,1,1))
  
  x_vals <- c(60:140)
  plot.overlap(x_vals, mean1 = 100, mean2 = 110)
  
  mtext("  Mean \n  difference",3, -3.5, adj=0)
  
  par(mar=c(0,3,1,1))
  plot.overlap(x_vals, sd1 = 10, sd2 = 15)
  mtext("  Breadth \n  difference",3, -3.5, adj=0)
  mtext("Abundance", 2, 1.5)
  
  par(mar=c(3,3,1,1))
  # peaks match up, but means are different
  # plot.overlap(x_vals, skew1=0, skew2=-5, mean2 = 106, sd2 = 18)
  # mtext("  Skewness \n  difference",3, -3.5, adj=0)
  # mtext("Time", 1, 1.5)
  
  plot.overlap(x_vals, skew1=0, skew2=-10, mean1=100, mean2 = 114.93, sd2 = 19.2) # the means args here are different because they are dsn means, not actual distribution means. the print results show the distribution means, and I manually lined them up to 100
  mtext("  Skewness \n  difference",3, -3.5, adj=0)
  mtext("Time", 1, 1.5)
  
dev.off()


### Skewness histograms
#png("figures/skewness.png", width=600, height=800, res=80)
svg("figures/skewness.svg", width=7, height=10)
  par(mfrow=c(1,1), mar=c(5,4.2,1,3.5))
  par(fig = c(0,1,0,1))
  
  breaks <- seq(-5, 5, by=0.1)
  
  hist(skews[guild=="flower", skew], breaks=breaks,
       main="", xlab="Skewness", ylab="Frequency",
       border=F, col=flr_pal["line"],
       xlim=c(-4,5), ylim=c(-200,300),
       yaxt="n", xaxt="n",
       cex.lab = 1.4, cex.axis=1.4)
  axis(2, at=seq(0,300, by=50), labels=seq(0,300, by=50), tick=T, las=2, cex.axis=1.4)
  hist(skews[guild=="flower" & skew_p > 0.05, skew], breaks = breaks,
       main="", xlab="Skewness", ylab="Frequency",
       border=F, col="white",
       yaxt="n", xaxt="n", add=T,
       density=15, angle=45)
  
  par(new=T)
  hist(skews[guild=="bee", skew], breaks=breaks,
       main="", xlab="", border=F,
       col=bee_pal["point"],
       ylim=c(0,60), xlim=c(-4,5),
       yaxt="n", ylab=NA, xaxt="n")
  axis(2, at=seq(0,20, by=5), labels=seq(0,20, by=5), tick=T, las=2, cex.axis=1.4)
  hist(skews[guild=="bee" & skew_p > 0.05, skew], breaks=breaks,
       main="", xlab="Skewness", ylab="Frequency",
       border=F, col="white",
       yaxt="n", xaxt="n", add=T,
       density=15, angle=45)
  
  axis(1, at=seq(-4,5,by=1), labels=seq(-4,5,by=1), tick=T, cex.axis=1.4)
  abline(v=0 ,lwd=5, col="white"); abline(v=0)
  
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
           x1=0.55,x2=0.90,y1=0.5,y2=0.7)
  par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
  skew.fig(skew=-0.85, percent=round(bee_prec_left, 1), line_col = bee_pal["line"], fill_col = bee_pal["point"], text_cex = 1,
           x1=0.05,x2=0.40,y1=0.20,y2=0.4)
  par(fig = c(0,1,0,1), mar=c(4.1,3.9,2,1))
  skew.fig(skew=0.85, percent=round(bee_prec_right, 1), line_col = bee_pal["line"], fill_col = bee_pal["point"], text_cex = 1,
           x1=0.55,x2=0.90,y1=0.20,y2=0.4)

dev.off()

### Joint mean and SD figure
#png("figures/mean_sd.png", width=450, height=1000)
#png("figures/mean_sd.png", width=3, height=6, units = "in", pointsize = 7, res=200)
svg("figures/mean_sd.svg", width=4, height=7.5)
  par(mfrow=c(2,1), mar=c(4.1,4,0.5,0.5))
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
  legend("topright", c("Flowers", "Bees"), col=c(flr_pal["line"], bee_pal["line"]),
         lwd=3, cex=0.8, inset=0.05)
  box()
  
  sd_model <- lm(abs(skew) ~ sd*guild, data = skews)
  summary(sd_model)
  visreg(sd_model, "sd", by="guild", xlab="Phenological breadth (days)", ylab="Absolute skewness", overlay=T,
         points = list(col=c(flr_pal["point"], bee_pal["line"]),
                       cex=c(0.2,0.7),
                       pch=c(19,21),
                       bg=c(1,bee_pal["point"])),
         line.par = list(col=c(NA, NA)),
         fill.par = list(col=adjustcolor(c(NA, NA), 0.5)),
         legend=F)
  # flower line
  flr_line <- as.data.table(visreg(sd_model, "sd", by="guild", plot=F)$fit)
  flr_line <- flr_line[guild == "flower" & sd < max(skews[guild == "flower", sd]),] # cutting out flower line beyond data range 
  lines(flr_line$visregFit ~ flr_line$sd, lwd=3, col=flr_pal["line"])
  polygon(x = c(flr_line$sd, rev(flr_line$sd)),
          y = c(flr_line$visregLwr, rev(flr_line$visregUpr)),
          col=adjustcolor(flr_pal["point"], 0.5), border = NA)
  # bee line
  bee_line <- as.data.table(visreg(sd_model, "sd", by="guild", plot=F)$fit)
  bee_line <- bee_line[guild == "bee" & sd > min(skews[guild == "bee", sd]),]
  lines(bee_line$visregFit ~ bee_line$sd, lwd=3, col=bee_pal["line"])
  polygon(x = c(bee_line$sd, rev(bee_line$sd)),
          y = c(bee_line$visregLwr, rev(bee_line$visregUpr)),
          col=adjustcolor(bee_pal["point"], 0.5), border = NA)
  
  box()
  
  par(mfrow=c(1,1), mar=c(5,4,4,2))
dev.off()


### demonstration of skewness effect on overlap
n <- 50
skew_vals <- seq(-5,5, length.out=n)
skew_scan <- as.data.table(expand.grid(skew_vals, skew_vals))
colnames(skew_scan) <- c("skew1", "skew2")

if(!file.exists("clean_data/skew_scan5.RDS")){
  overlaps_vals <- pbmapply(maximize.overlap, skew1 = skew_scan$skew1, skew2 = skew_scan$skew2)
  skew_scan$overlap <- overlaps_vals
  saveRDS(skew_scan,"clean_data/skew_scan5.RDS")
} else{
  skew_scan <- readRDS("clean_data/skew_scan5.RDS")
}

# getting flr and bee quantiles for plotting
flr_low <- skews[guild == "flower", quantile(skew, 0.025)]
flr_high <- skews[guild == "flower", quantile(skew, 0.975)]
bee_low <- skews[guild == "bee" & dataset == "rmbl", quantile(skew, 0.025)]
bee_high <- skews[guild == "bee" & dataset == "rmbl", quantile(skew, 0.975)]

# heatmap
#png("figures/heatmap.png", width=700, height=600)
svg("figures/heatmap.svg", width=7.5, height=6)

  skew_scan_mat <- matrix(skew_scan$overlap * 100, nrow = n, ncol = n)
  
  par(mar=c(4.5,4.2,2,3))
  image.plot(x = skew_vals, y = skew_vals, z = skew_scan_mat,
        col=colorRampPalette(c("red","orange","white"))(100),
        xlab = "Flower skewness", ylab="Bee skewness",
        legend.shrink = 0.75)
  
  #white
  abline(v=flr_low, col="white", lwd=4, lty=1)
  abline(v=flr_high, col="white", lwd=4, lty=1)
  abline(h=bee_low, col="white", lwd=4, lty=1)
  abline(h=bee_high, col="white", lwd=4, lty=1)
  #color
  abline(v=flr_low, col=flr_pal["line"], lwd=2, lty=3)
  abline(v=flr_high, col=flr_pal["line"], lwd=2, lty=3)
  abline(h=bee_low, col=bee_pal["line"], lwd=2, lty=3)
  abline(h=bee_high, col=bee_pal["line"], lwd=2, lty=3)
  #white
  lines(x = c(flr_low, flr_low), y = c(bee_low, bee_high), col = "white", lwd=9)
  lines(x = c(flr_high, flr_high), y = c(bee_low, bee_high), col = "white", lwd=9)
  lines(x = c(flr_low, flr_high), y = c(bee_low, bee_low), col = "white", lwd=9)
  lines(x = c(flr_low, flr_high), y = c(bee_high, bee_high), col = "white", lwd=9)
  #color
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
  
  # adding a sequqnce of skew figs along y and x axes
  
  #skew.fig(-0.9, x1=0, x2=0.5, y1=0, y2=0.5, line_col="red", fill_col = "pink", reset_window=T)

dev.off()

# What's the distribution of overlap within the bee/flower box?
relevant_overlap <- skew_scan[skew1 >= flr_low & skew1 <= flr_high &
                                skew2 >= bee_low & skew2 <= skew2, overlap]
hist(relevant_overlap,40)
paste0("Skewness could account for up to ", round(1 - min(relevant_overlap),2)*100, "% of overlap losses based on the data") 

# what's the distribution of overlap in the whole simulation?
hist(skew_scan$overlap,40)
paste0("The maximum loss in overlap due to skewness from the simulation is ", round(1 - min(skew_scan$overlap),2)*100, "%") 


### Summary tables
#readRDS("clean_data/bees_summary.RDS")
bees_summary <- bee_skews[, .(mean_date =  gsub("2000-","",as.Date(round(mean(mean)), origin="2000-01-01")),
                              mean_skew = round(mean(skew),2),
                              mean_sd = round(mean(sd),2),
                              mean_abundance = round(mean(n)),
                              n_timeseries = .N,
                              n_sites = length(unique(site)),
                              n_years = length(unique(year))), by=.(species)]

bees_summary[, genus := gsub(" .*","",species)]
bees_summary[, family := tax_name(genus, get="family", db="ncbi")$family, by=.(genus)]
bees_summary <- bees_summary[order(family, species)]
setcolorder(bees_summary, c("family", "species", "mean_date", "mean_skew", "mean_sd", "mean_abundance", "n_timeseries", "n_sites", "n_years"))
bees_summary[, genus := NULL]
saveRDS(bees_summary, "clean_data/bees_summary.RDS")

flrs_summary <- flr_skews[, .(mean_date =  gsub("2000-","",as.Date(round(mean(mean)), origin="2000-01-01")),
                              mean_skew = round(mean(skew),2),
                              mean_sd = round(mean(sd),2),
                              mean_abundance = round(mean(n)),
                              n_timeseries = .N,
                              n_sites = length(unique(site)),
                              n_years = length(unique(year))), by=.(species)]

flrs_summary[, genus := gsub(" .*","",species)]
flrs_summary[, family := tax_name(genus, get="family", db="ncbi")$family, by=.(genus)]
flrs_summary <- flrs_summary[order(family, species)]
setcolorder(flrs_summary, c("family", "species", "mean_date", "mean_skew", "mean_sd", "mean_abundance", "n_timeseries", "n_sites", "n_years"))
flrs_summary[, genus := NULL]
flrs_summary[is.na(family), species] # these weren't found automatically, so I'm filling by hand
flrs_summary[is.na(family), "family"] <- c("Crassulaceae", "Caprifoliaceae", "Orchidaceae",
                                         "Asteraceae", "Rosaceae", "Gentianaceae", "Crassulaceae") # this is rickety, sorry
flrs_summary <- flrs_summary[order(family, species)]
saveRDS(flrs_summary, "clean_data/flrs_summary.RDS")

bees_table <- flextable(bees_summary)
bees_table <- fontsize(bees_table, size=9)
bees_table <- fontsize(bees_table, size=9, part="header")
bees_table <- theme_vanilla(bees_table)
bees_table <- width(bees_table, j=2, 0.3)
save_as_html(bees_table, path="manuscript/bees_table.html")
bees_table

flrs_table <- flextable(flrs_summary)
flrs_table <- fontsize(flrs_table, size=9)
flrs_table <- fontsize(flrs_table, size=9, part="header")
flrs_table <- theme_vanilla(flrs_table)
flrs_table <- width(flrs_table, j=2, 0.3)
save_as_html(flrs_table, path="manuscript/flrs_table.html")
flrs_table


### Supplementary figures

# skewness by dataset
png("figures/skewness_mid_atlantic.png", width=800, height=500)
  par(mfrow=c(1,2))
  hist(bee_skews_rmbl_mary[dataset == "rmbl", skew], 100, xlim=c(-4,5), main="Rocky Mountain", xlab="Skewness")
  hist(bee_skews_rmbl_mary[dataset == "mary", skew], 100, xlim=c(-4,5), main="Mid-Atlantic", xlab="Skewness", ylab="")
  dataset_model <- lm(skew ~ dataset, data = bee_skews_rmbl_mary)
  summary(dataset_model)
dev.off()


# skew by truncation group
png("figures/truncation.png", width=5, height=8, units="in", res=300, pointsize = 9)
  
  par(mfrow=c(4,2))
  hist(bee_skews[censored == "none",skew], 20, xlab="Skewness", main="not truncated (bees)"); abline(v=mean(bee_skews[censored == "none",skew]), col="red")
  hist(bee_skews[censored == "before",skew], 20, xlab="Skewness", main="left-truncated (bees)"); abline(v=mean(bee_skews[censored == "before",skew]), col="red")
  hist(bee_skews[censored == "after",skew], 20, xlab="Skewness", main="right-truncated (bees)"); abline(v=mean(bee_skews[censored == "after",skew]), col="red")
  hist(bee_skews[censored == "both",skew], 20, xlab="Skewness", main="doubly-truncated (bees)"); abline(v=mean(bee_skews[censored == "both",skew]), col="red")
  
  hist(flr_skews[censored == "none",skew], 40, xlab="Skewness", main="not truncated (flowers)"); abline(v=mean(flr_skews[censored == "none",skew]), col="red")
  hist(flr_skews[censored == "before",skew], 40, xlab="Skewness", main="left-truncated (flowers)"); abline(v=mean(flr_skews[censored == "before",skew]), col="red")
  hist(flr_skews[censored == "after",skew], 40, xlab="Skewness", main="right-truncated (flowers)"); abline(v=mean(flr_skews[censored == "after",skew]), col="red")
  hist(flr_skews[censored == "both",skew], 40, xlab="Skewness", main="doubly-truncated (flowers)"); abline(v=mean(flr_skews[censored == "both",skew]), col="red")
  par(mfrow=c(1,1))

dev.off()

# mean and SD correlations by truncation group
png("figures/truncation_mean_sd.png", width=8, height=9, units="in", res=300, pointsize = 9)
  
  #bee_skews[, censored := mgsub(censored, c("after", "before", "both"), c("right-", "left-", "doubly-")),] # changing censoring names for the legend

  par(mfrow=c(2,2))
  
  mean_model_trunc_bees <- lm(skew ~ mean*censored, data=bee_skews)
  summary(mean_model_trunc_bees)
  line_pal <- brewer.pal(4, "Dark2")
  visreg(mean_model_trunc_bees, "mean", by="censored", overlay=T,
         points.par = list(col=line_pal), line.par=list(col=line_pal), fill.par=list(col=adjustcolor(line_pal,0.2)))
  mtext("Bees    ", side=3, line=-2, adj=1, cex=1.5)
  
  sd_model_trunc_bees <- lm(skew ~ sd*censored, data=bee_skews)
  summary(sd_model_trunc_bees)
  line_pal <- brewer.pal(4, "Dark2")
  visreg(sd_model_trunc_bees, "sd", by="censored", overlay=T,
         points.par = list(col=line_pal), line.par=list(col=line_pal), fill.par=list(col=adjustcolor(line_pal,0.2)))
  mtext("Bees    ", side=3, line=-2, adj=1, cex=1.5)
  
  mean_model_trunc_flrs <- lm(skew ~ mean*censored, data=flr_skews)
  summary(mean_model_trunc_flrs)
  # making the "before" "after" and "both" colors more transparents because they're based on just a few points each in flowers
  line_pal[1] <- adjustcolor(line_pal[1], 0.3)
  line_pal[2] <- adjustcolor(line_pal[2], 0.3)
  line_pal[3] <- adjustcolor(line_pal[3], 0.3)
  visreg(mean_model_trunc_flrs, "mean", by="censored", overlay=T,
         ylim=c(-2,5), points = list(col=line_pal, cex=0.15), line=list(col=line_pal), fill=list(col=adjustcolor(line_pal,0.2)))
  mtext("    Flowers", side=3, line=-2, adj=0, cex=1.5)
  
  sd_model_trunc_flrs <- lm(skew ~ sd*censored, data=flr_skews)
  summary(sd_model_trunc_flrs)
  visreg(sd_model_trunc_flrs, "sd", by="censored", overlay=T,
         ylim=c(-2,5), points.par = list(col=line_pal, cex=0.15), line.par=list(col=line_pal), fill.par=list(col=adjustcolor(line_pal,0.2)))
  mtext("    Flowers", side=3, line=-2, adj=0, cex=1.5)
  
  par(mfrow=c(1,1))

dev.off()

# skewness correlated with population size
png("figures/pop_size.png", width=6, height=3, units="in", res=300, pointsize = 9)
  par(mfrow=c(1,2))
  # log-transforming just to visualize the data better
  hist(log(skews[guild == "bee",n]), 50, main="Bees", xlab="Log total number of bees caught")
  hist(log(skews[guild == "flower",n]), 50, main="Flowers", xlab="Log total number of flowers counted") # note this isn't # of unique flowers because there are recounts
  par(mfrow=c(1,1))
dev.off()

summary(bee_skews$n)
summary(flr_skews$n)

#pop_skew_model <- lm(skew ~ n*guild, data = skews) # not using this full model because "n" is a different measure b/w bees and flowers b/c recounting of flowers
pop_skew_bees <- lm(skew ~ n, data = bee_skews)
summary(pop_skew_bees)
hist(resid(pop_skew_bees), 100)

pop_skew_flowers <- lm(skew ~ n, data = flr_skews)
summary(pop_skew_flowers)
hist(resid(pop_skew_flowers), 100)

png("figures/pop_vs_skew.png", width=6, height=3, units="in", res=300, pointsize = 9)
  par(mfrow=c(1,2))
  visreg(pop_skew_bees, main="Bees", xlab="Total number of bees caught")
  visreg(pop_skew_flowers, main="Flowers", xlab="Total number of flowers counted")
  par(mfrow=c(1,1))
dev.off()

# is variance related to pop size?
bee_test <- data.table(resids = abs(resid(pop_skew_bees)), n = bee_skews$n)
bee_var_mod <- rq(resids ~ n, data=bee_test, tau=0.68)
summary(bee_var_mod)
visreg(bee_var_mod)
bptest(pop_skew_bees)

flr_test <- data.table(resids = abs(resid(pop_skew_flowers)), n = flr_skews$n)
flr_var_mod <- rq(resids ~ n, data=flr_test, tau=0.68)
summary(flr_var_mod)
visreg(flr_var_mod, ylab="Absolute residuals")

bptest(pop_skew_flowers)

var(flr_skews[, skew])
var(flr_skews[n < 100, skew])
var(flr_skews[n < 10000, skew])
var(flr_skews[n < 10000, skew])
var(flr_skews[n > 10000, skew])




### Does drought predict skewness? - probably won't include this in the main text

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

# plopping in snowmelt date
bare_ground <- fread("raw_data/billy_barr_bare_ground.csv")
bare_ground[, date_bare_ground := mdy(date_bare_ground)]
bare_ground[, snowmelt_doy := yday(date_bare_ground)]

skews_weather <- merge(skews, drought_data, by="year")
skews_weather <- merge(skews_weather, bare_ground, by="year")
skews_weather <- skews_weather[dataset != "mary",]

plot(skew ~ snowmelt_doy, data = skews_weather[guild=="flower",])
plot(skew ~ snowmelt_doy, data = skews_weather[guild=="bee",])

weather_model <- lm(skew ~ summer_drought*guild + snowmelt_doy*guild + site, data = skews_weather)
summary(weather_model)
visreg(weather_model, "summer_drought", by="guild")
visreg(weather_model, "snowmelt_doy", by="guild")

species_weather_model_bees <- with(skews_weather[guild=="bee",], lmer(skew ~ snowmelt_doy + (1|species)))
summary(species_weather_model_bees)
visreg(species_weather_model_bees)
visreg(species_weather_model_bees, "snowmelt_doy", by="species")



