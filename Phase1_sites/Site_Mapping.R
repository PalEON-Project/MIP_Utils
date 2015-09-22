# Plotting PalEON MIP Site Locations
library(raster)
library(rgdal)
library(maps)

setwd("../")

spatdat.dir <- "SpatialData"
figures.dir <- "Figures"

sites <- data.frame(Site=c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL", "PMB" ), 
					Lon =c(-72.18, -79.73, -89.53, -94.58, -95.17, -82.83) + 360,
					Lat =c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61)
					)
summary(sites)
coordinates(sites) <- c("Lon", "Lat")
# plot(sites, pch=19)

# paleon.domain <- raster(file.path(spatdat.dir, "paleon_full_ll_v0.1.tif"))
#plot(paleon.domain)

tavg <- raster(file.path(spatdat.dir, "PRISM_Climate", "PRISM_Tavg_Year_PalEON_Mask"))

pdf(file.path(figures.dir, "PalEON_MIP_SiteLocations.pdf"), width=10, height=6)
par(mar=c(5,5,1,0))
plot(tavg, col=bpy.colors(250), xlab="Longitude", ylab="Latitude", font.lab=2, cex.lab=1.5, legend=F, cex.axis=1.25)
plot(tavg, legend.only=T, legend.args=list(text=expression(bold(paste("Temperature ("^"o","C)"))), side=1, line=2.7, font=2, cex=1.5, cex.lab= 1.25), col=bpy.colors(250), smallplot=c(0.45, 0.8, 0.88, 0.91), axis.args=list(cex.axis=1.25), horizontal=T)
map("state", plot=T, add=T, lty="solid", col="gray30", lwd=1.5)
plot(sites, add=T, pch=19, cex=2, col="green3")
dev.off()

setveg.biom <- read.csv(file.path(spatdat.dir, "plss_biomass_v0.9-1.csv"))
setveg.biom$Total <- rowSums(setveg.biom[,5:ncol(setveg.biom)])
summary(setveg.biom)
