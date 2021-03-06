# ---------------------------------------------------------
# Script Information 
# ---------------------------------------------------------
# Calculating PDSI on the Regional MIP output.
# This will take the monthly met output from the models
# and combine them with the environmental drivers we provided
# models
#
# ------------------
# Workflow
# ------------------
# 0. Define file paths, etc.
# 1. Extract Data: Monthly T & P, soil texture & depth
# 2. Reformat Data, calculate AWC
# 3. Calculate PDSI
# 4. Format Data, send it out
# ------------------
# ---------------------------------------------------------


# ---------------------------------------------------------
# 0. Define file paths, etc.
# ---------------------------------------------------------

# Path to where the MIP output is
# - this will give you monthly met vars
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where the model environmental drivers are; these can be downloaded
# from Cyverse
path.soil <- "~/Dropbox/PalEON_CR/env_regional/phase2_env_drivers_v2/soil/"

# The Model Met outputs aren't matching, so lets pull the raw data
path.met <- "~/Desktop/Research/PalEON_MIP_Region/phase2_met_regional_v2_monthly/"

# This is the path to my pdsi code
# This is currently part of my met ensemble code and code 
# and can be found here: https://github.com/PalEON-Project/modeling_met_ensemble.git
path.pdsi <- "~/Desktop/Research/PalEON_CR/met_ensemble/scripts/"

# Path where you want to save the PDSI time series  
path.out <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/"
# ---------------------------------------------------------

# ---------------------------------------------------------
# 1. Extract Data: Monthly T & P, soil texture & depth
# ---------------------------------------------------------
# load in the paleon key
load(file.path(path.data, "PalEON_siteInfo_all.RData"))
paleon <- paleon[,2:4]
paleon$latlon <- as.factor(paleon$latlon)
summary(paleon)

# Load in the monthly temp & precip time series for each site
tair1    <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.tair.rds"   ))
precip1  <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.precipf.rds"))
tair2   <-  readRDS(file.path(path.data, "ED2/ED2.tair.rds"   ))
precip2 <-  readRDS(file.path(path.data, "ED2/ED2.precipf.rds"   ))
sw1 <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.swdown.rds"   ))
sw2 <- readRDS(file.path(path.data, "ED2/ED2.swdown.rds"   ))

# plot(tair[(nrow(tair1)-100*12):nrow(tair1),1] ~ tair2[(nrow(tair1)-100*12):nrow(tair1),1]); abline(a=0, b=1, col="red")
# plot(precip[(nrow(tair1)-10*12):nrow(tair1),2] ~ precip2[(nrow(tair1)-10*12):nrow(tair1),2]); abline(a=0, b=1, col="red")

# plot(sw1[,1] ~ sw2[,1]); abline(a=0, b=1, col="red")

# Extracting soil data for the appropriate sites 
# Note: We need this for the "upper" and "lower" soil layers
library(ncdf4)
tair.nc <- nc_open(file.path(path.met, "tair.nc"))
precipf.nc <- nc_open(file.path(path.met, "precipf.nc"))
sand.t <- nc_open(file.path(path.soil, "paleon_soil_t_sand.nc"))
sand.s <- nc_open(file.path(path.soil, "paleon_soil_s_sand.nc"))
clay.t <- nc_open(file.path(path.soil, "paleon_soil_t_clay.nc"))
clay.s <- nc_open(file.path(path.soil, "paleon_soil_s_clay.nc"))
depth  <- nc_open(file.path(path.soil, "paleon_soil_soil_depth.nc"))

# Getting the lat/lon index for each point to extract soil
lon <- ncvar_get(sand.t, "longitude")
lat <- ncvar_get(sand.t, "latitude")

lon2 <- ncvar_get(tair.nc, "lon")
lat2 <- ncvar_get(tair.nc, "lat")

paleon$x.ind <- apply(data.frame(paleon$lon), 1, FUN=function(x){which(lon == x)})
paleon$y.ind <- apply(data.frame(paleon$lat), 1, FUN=function(x){which(lat == x)})

tair.raw <- matrix(NA, nrow=nrow(tair1), ncol=ncol(tair1)) # A place holder matrix
precipf.raw <- matrix(NA, nrow=nrow(tair1), ncol=ncol(tair1)) # A place holder matrix
for(i in 1:nrow(paleon)){
  x.ind  <- which(lon == paleon[i,"lon"])
  y.ind  <- which(lat == paleon[i,"lat"])
  x.ind2 <- which(lon2 == paleon[i,"lon"])
  y.ind2 <- which(lat2 == paleon[i,"lat"])
  
  tair.raw[,i]       <- ncvar_get(tair.nc, "tair", c(x.ind2, y.ind2, 1), c(1,1,13932))
  precipf.raw[,i]    <- ncvar_get(precipf.nc, "precipf", c(x.ind2, y.ind2, 1), c(1,1,13932))
  paleon[i,"sand.t"] <- ncvar_get(sand.t, "t_sand", c(x.ind, y.ind), c(1,1))
  paleon[i,"sand.s"] <- ncvar_get(sand.s, "s_sand", c(x.ind, y.ind), c(1,1))
  paleon[i,"clay.t"] <- ncvar_get(clay.t, "t_clay", c(x.ind, y.ind), c(1,1))
  paleon[i,"clay.s"] <- ncvar_get(clay.s, "s_clay", c(x.ind, y.ind), c(1,1))
  paleon[i,"depth"]  <- ncvar_get(depth , "soil_depth", c(x.ind, y.ind), c(1,1)) # cm
}



# Comparing my monthly averages with the models
plot(tair1[,1] ~ tair.raw[,1], main="LPJ-GUESS vs. Actual"); abline(a=0, b=1, col="red")
plot(tair2[,1] ~ tair.raw[,1], main="ED2 vs. Actual"); abline(a=0, b=1, col="red")
plot(tair2[,1] ~ tair1[,1], main="LPJ-GUESS vs. ED2"); abline(a=0, b=1, col="red")

plot(precip1[,1] ~ precipf.raw[,1], main="LPJ-GUESS vs. Actual"); abline(a=0, b=1, col="red")
plot(precip2[,1] ~ precipf.raw[,1], main="ED2 vs. Actual"); abline(a=0, b=1, col="red")
plot(precip2[,1] ~ precip1[,1], main="LPJ-GUESS vs. ED2"); abline(a=0, b=1, col="red")

# ---------------------------------------------------------


# ---------------------------------------------------------
# 2. Reformat Data, calculate AWC, PDSI
# ---------------------------------------------------------
# Calculate AWC for each site
# Assuming 30 cm for top soil; this is safe for all of our sites
source(file.path(path.pdsi, "calc.awc.R"))
paleon$awc.t <- calc.awc(paleon$sand.t, paleon$clay.t)
paleon$awc.s <- calc.awc(paleon$sand.s, paleon$clay.s)

# Calculating water capacity **in inches**
paleon$wcap.t <- paleon$awc.t * 30 * 1/2.54 # 30 cm depth * 1 in / 2.54 cm
paleon$wcap.s <- paleon$awc.s * (paleon$depth-30) * 1/2.54 # 30 cm depth * 1 in / 2.54 cm
summary(paleon)

library(R.matlab)
dayz <- readMat(file.path(path.pdsi, "PDSI_fromBenCook/PDSICODE/daylennh.mat"))$dayz

# Formatting the climate data and calculating PDSI
source(file.path(path.pdsi, "pdsi1.R"))
pdsi.final <- matrix(NA, nrow=nrow(tair.raw), ncol=ncol(tair.raw)) # A place holder matrix
pdsi.all <- array(dim=c(1161, 12, ncol(tair.raw))) # Save it all so we can make some easier graphs
pdsiM.all <- array(dim=c(1161, 12, ncol(tair.raw))) # Save it all so we can make some easier graphs
pdsiH.all <- array(dim=c(1161, 12, ncol(tair.raw))) # Save it all so we can make some easier graphs

pb <- txtProgressBar(min = 0, max = ncol(tair.raw), style = 3)
# problem site: #196

sw.corner <- which(paleon$lat==min(paleon$lat) & paleon$lon==min(paleon$lon))
problem <- which(paleon$lat==36.75 & paleon$lon==min(paleon$lon))
problem2 <- 127
site.min2 <- 120
paleon[problem,]

problem3 <- 249
# yrs.prob3 <- 587 # Starts in year 286 though
for(i in 1:ncol(tair.raw)){
  if(is.na(mean(tair.raw[,i]))) next
  # Format the climate data into monthly columns
  TEMP1   <- matrix(tair.raw[,i], nrow=nrow(tair.raw)/12, ncol=12, byrow=T)
  PRECIP1 <- matrix(precipf.raw[,i], nrow=nrow(precipf.raw)/12, ncol=12, byrow=T)
  row.names(TEMP1) <- 850:2010
  colnames (TEMP1) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  row.names(PRECIP1) <- 850:2010
  colnames (PRECIP1) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  # Convert Temp
  # TEMP native units: K
  C2F <- function(x){x*9/5 + 32}
  TEMP1 <- C2F(TEMP1-273.15)
  
  # Convert Precip: 
  # PRECIP native units: kg/m2/s = mm/s 
  # mm/s -> in /mo:  mm/s*s/day*dpm * in/mm
  library(lubridate)
  sec2day = 1/(60*60*24)
  dpm <- days_in_month(1:12)
  rows.leap <- which(leap_year(850:2010))
  PRECIP1 <- PRECIP1*1/sec2day # mm/s to mm/day
  PRECIP1 <- PRECIP1/25.4 # mm to in
  PRECIP1 <- t(apply(PRECIP1, 1, function(x){x*dpm}))
  PRECIP1[rows.leap,2] <- PRECIP1[rows.leap,2]*29/28
  
  datmet <- list(Temp=TEMP1, Precip=PRECIP1)
  
  # Adding in everything else we need
  datother <- list()
  datother$pdsi.fun <- path.pdsi
  datother$metric <- F
  datother$lat <- paleon$lat[i]
  datother$watcap <- list(awcs=paleon$wcap.t[i], awcu=paleon$wcap.s[i])
  datother$yrs.calib <- c(1931, 1990) # copied from B. Cook
  datother$dayz      <- dayz
  datother$daylength <- NULL
  
  siteID <- paleon$latlon[i]
  
  # Run the actual PDSI calculation
  pdsi.out <- pdsi1(datmet, datother, metric=F, siteID, method.PE="Thornthwaite", snow=NULL, snowopts=NULL, penopts=NULL, datpen=NULL)
  
  setTxtProgressBar(pb, i)
  
  pdsi.all[,,i] <- pdsi.out$X
  pdsiM.all[,,i] <- pdsi.out$XM
  pdsiH.all[,,i] <- pdsi.out$XH
  pdsi.final[,i] <- as.vector(t(pdsi.out$X))
}


saveRDS(pdsi.final, file.path(path.out, "PalEON_Regional_Extract/Met/pdsi.rds"))
saveRDS(tair.raw, file.path(path.out, "PalEON_Regional_Extract/Met/tair.rds"))
saveRDS(precipf.raw, file.path(path.out, "PalEON_Regional_Extract/Met/precipf.rds"))


# yrs=850:2010
# rows.cal <- which(yrs %in% 1931:1990)
# par(mfrow=c(3,1))
# plot(rowMeans(pdsi.out$T), type="l", col="darkgoldenrod2"); abline(h=mean(rowMeans(pdsi.out$T[rows.cal,])), lty="dashed")
# plot(rowMeans(pdsi.out$P), type="l", col="blue"); abline(h=mean(rowMeans(pdsi.out$P[rows.cal,])), lty="dashed")
# plot(rowMeans(pdsi.out$X), type="l", col="red"); abline(h=mean(rowMeans(pdsi.out$X[rows.cal,])), lty="dashed")
# par(mfrow=c(1,1))
# plot(rowMeans(pdsi.out$PE), type="l", col="red"); abline(h=mean(rowMeans(pdsi.out$PE[rows.cal,])), lty="dashed")
# plot(rowMeans(pdsi.out$P), type="l", col="red"); abline(h=mean(rowMeans(pdsi.out$P[rows.cal,])), lty="dashed")
# plot(rowMeans(pdsi.out$W), type="l", col="red"); abline(h=mean(rowMeans(pdsi.out$W[rows.cal,])), lty="dashed")
# plot(rowMeans(pdsi.out$RO), type="l", col="red"); abline(h=mean(rowMeans(pdsi.out$RO[rows.cal,])), lty="dashed")
# plot(rowMeans(pdsi.out$R), type="l", col="red"); abline(h=mean(rowMeans(pdsi.out$R[rows.cal,])), lty="dashed")
# 

pdsi.ann <- apply(pdsi.all, c(1,3), mean, na.rm=T)

# Finding out which sites have unreasonably values
pdsi.min <- apply(pdsi.ann, 2, min)
length(pdsi.min)
sites.low <- which(pdsi.min < -20)
length(sites.low)
summary(paleon[sites.low,])

site.min <- which(pdsi.min == min(pdsi.min, na.rm=T))
paleon[site.min,]

pdsi <- data.frame(year=850:2010, 
                   mean=apply(pdsi.ann, 1, mean, na.rm=T),
                   lwr =apply(pdsi.ann, 1, quantile, 0.025, na.rm=T),
                   upr =apply(pdsi.ann, 1, quantile, 0.975, na.rm=T))


library(ggplot2)
ggplot(data=pdsi) +
  geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), alpha=0.5) +
  geom_line(aes(x=year, y=mean), size=1)

# 
# 

# ---------------------------------------------------------


