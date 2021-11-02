# cleaning bee and flower data for phenological skew analysis

library(data.table)
library(tools)
library(mgsub)
library(lubridate)

`%!in%` = Negate(`%in%`)

setwd("/home/michael/Documents/Grad School/Research Projects/pheno_skew")

# old flower data
# load("raw_data/rmbl_flowers/rmbl.RData")
# flowers <- data.table(all.data)
# rm(all.data)

load("raw_data/rmbl_flowers/rmbl_flowers_2019.RData")
flowers <- data.table(data)
rm(data)

flowers <- flowers[!is.na(doy),]

#bees <- fread("raw_data/sp_time_series_2019_10_13_no_sings.csv")
#bees <- fread("raw_data/data_summary_2019_10_14_4.csv")
rmbl_bees <- fread("raw_data/rmbl_bees/bees_2020-12-16.csv", stringsAsFactors = FALSE)
rmbl_bombus <- fread("raw_data/rmbl_bees/bombus_2020-12-16.csv", stringsAsFactors = FALSE)
rmbl_bombus[, genus_species := paste(genus, genus_species)]
mary_bees <- fread("raw_data/maryland_bees/2OccurrenceLevel_WithTrapInfo.csv", stringsAsFactors = FALSE)

unwanted_taxa <- fread("raw_data/unwanted_taxa.csv", stringsAsFactors = FALSE)
rmbl_sites <- fread("raw_data/rmbl_sites.csv", stringsAsFactors = FALSE)
rmbl_effort <- fread("raw_data/rmbl_effort.csv", stringsAsFactors = FALSE)

# cleaning unwanted taxa
rmbl_bees <- rmbl_bees[genus_species != "" & genus_species %!in% unwanted_taxa$rmbl_bees_species & genus %!in% unwanted_taxa$rmbl_bees_genus & genus != "Bombus",]
rmbl_bombus <- rmbl_bombus[genus_species != "" & genus_species %!in% unwanted_taxa$rmbl_bees_species,]
mary_bees <- mary_bees[grouped_name != "" & grouped_name %!in% unwanted_taxa$mary_bees_species & Genus %!in% unwanted_taxa$mary_bees_genus,]
flowers <- flowers[species != "" & species %!in% unwanted_taxa$rmbl_flrs_species,]

# I want to separate queens and workers in analysis
rmbl_bombus[, genus_species := paste(genus_species, caste)]

# weird little cleaning
mary_bees <- mary_bees[trapdays == 1,] # removing when traps were left out for more than 1 day
rmbl_bees <- rmbl_bees[method %in% c("Bowl","bowl"),] # removing netted bees
rmbl_bees$sex <- toupper(rmbl_bees$sex)
rmbl_bombus$caste <- mgsub(rmbl_bombus$caste, c("Q", "W"), c("F", "F")) # collapsing queens and workers into females - following Maryland dataset
mary_bees$sex <- mgsub(mary_bees$sex, c("female", "male"), c("F", "M"))
mary_bees <- mary_bees[sex %in% c("F", "M"),]

# removing timeseries with uneven numbers of traps over time
bowl_disparities <- mary_bees[, sd(as.numeric(NTraps), na.rm=T), by=.(SiteID_Year)]
bowl_disparities <- bowl_disparities[V1 > 0, ]
mary_bees <- mary_bees[SiteID_Year %!in% bowl_disparities$SiteID_Year, ] 

# adding lat and lon
rmbl_bees[, site := toTitleCase(site), by=site]
rmbl_bombus[, site := toTitleCase(site), by=site]
setkey(rmbl_sites, "site")
setkey(rmbl_bees, "site")
rmbl_bees <- merge(rmbl_bees, rmbl_sites, by="site")
rmbl_bombus <- rmbl_bombus[site %in% rmbl_sites$site,]
setkey(rmbl_bombus, "site")
rmbl_bombus <- merge(rmbl_bombus, rmbl_sites, by="site")

# cutting out low sampling effort days - assuming that the rest accurately captured abundance
mary_bees[, `:=`(time1 = as_datetime(gsub("x", "0", time1)), time2 = as_datetime(gsub("x", "0", time2)))]
#as.numeric(mary_bees[150, time2] - mary_bees[150, time1])
#mary_bees[150, .(time1, time2)]
mary_bees[, trap_time := as.numeric(time2-time1)/3600]
hist(mary_bees$trap_time)
abline(v=3, col="red")
mary_bees <- mary_bees[trap_time >= 3,] # removing days when traps were out for less than an hour

get.effort <- function(site_v, year_v, date_v, method_v){
  record <- rmbl_effort[site == site_v & year == year_v & date_sampled == date_v,]
  if(method_v == "bowl") time <- as.numeric(hm(record$bowl_time))/3600 #time <- as.numeric(gsub("\\:","\\.",record$bowl_time))
  if(method_v == "net") time <- as.numeric(record$total_time)
  return(time)
}
rmbl_bees[, trap_time := get.effort(site, year, date_sampled, "bowl"), by=.(site, year, date_sampled)]
hist(rmbl_bees$trap_time,100, xlab="Bowl hours", main="")
abline(v=3, col="red")
rmbl_bees <- rmbl_bees[trap_time >= 3,]

rmbl_bombus[, trap_time := get.effort(site, year, date_sampled, "net"), by=.(site, year, date_sampled)]
hist(rmbl_bombus$trap_time,100, xlab="Bowl hours", main="")
abline(v=0.95, col="red")
rmbl_bombus <- rmbl_bombus[trap_time >= 1,]

# combining data
bee_data <- data.table(year = c(rmbl_bees$year, rmbl_bombus$year, mary_bees$year),
                       date = c(dmy(paste(rmbl_bees$date_sampled, rmbl_bees$year, sep="-")), dmy(paste(rmbl_bombus$date_sampled, rmbl_bombus$year, sep="-")), ymd(mary_bees$startdate)),
                       site = c(rmbl_bees$site, rmbl_bombus$site, mary_bees$SiteID),
                       dataset = c(rep("rmbl", nrow(rmbl_bees)+nrow(rmbl_bombus)), rep("mary", nrow(mary_bees))),
                       lat = c(rmbl_bees$lat, rmbl_bombus$lat, mary_bees$latitude),
                       lon = c(rmbl_bees$lon, rmbl_bombus$lon, mary_bees$longitude),
                       species = c(rmbl_bees$genus_species, rmbl_bombus$genus_species, mary_bees$grouped_name),
                       sex = c(rmbl_bees$sex, rmbl_bombus$caste, mary_bees$sex),
                       effort = c(rmbl_bees$trap_time, rmbl_bombus$trap_time, mary_bees$trap_time))
bee_data
bee_data[dataset == "rmbl" & grepl("Bombus", species),]

# cleaning "sp." or uncertain groups
# I'm getting rid of species containing sp. spp. ssp. aff. cf. af. w- nr
# I'm removing bombus with "unknown" caste
# keeping species with "/" because they seem not to duplicate
ambiguities <- paste(c("sp\\.", "spp\\.", "ssp\\.", "aff\\.", "cf\\.", "af\\.", "w\\-", " nr ", "unknown"),collapse="|")
bee_data <- bee_data[!grepl(ambiguities, species), ]
flowers <- flowers[!grepl(ambiguities, species),]

# totaling catch by species/sex/site/date, making it analogous to flower data
rm(bees)
bees <- bee_data[, .(year = unique(year),
                     dataset = unique(dataset),
                     lat = unique(lat),
                     lon = unique(lon),
                     effort = mean(effort),
                     ab = .N),
                 by = .(species, sex, site, date)]
bees$doy <- yday(bees$date)

# aggregating flower counts by site
flowers[, site := gsub("[[:digit:]]", "", plot), by=plot]
flowers <- flowers[, .(year = unique(year),
                       ab = sum(floralcount),
                       doy = unique(doy)),
                    by = .(date, species, site)]

#colnames(flowers) <- c("year", "site", "species", "date", "doy", "ab")

write.csv(bees, "clean_data/bees.csv", row.names = FALSE)
write.csv(flowers, "clean_data/flowers.csv", row.names = FALSE)

