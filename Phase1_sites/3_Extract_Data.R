# ----------------------------------------------------
# Graphs & Analysis for Poster at 2015 NACP Meeting
# Christine R. Rollinson
# crollinson@gmail.com 
# 18 Jan 2015
# ----------------------------------------------------

library(ncdf4)
library(ggplot2)
library(grid)
library(nlme)

setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/")
outputs <- "phase1a_output_variables"
years <- 850:2010

large.axes <- theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=20), axis.text.y=element_text(color="black", size=20), axis.title.x=element_text(face="bold", size=24, vjust=-1),  axis.title.y=element_text(face="bold", size=24, vjust=2.5), plot.margin=unit(c(2,2,2,2), "lines"))

# -------------------------------------------------------------------
# Reading in .nc files of variables created using script in MIP_formatting_ModelLoop_Yr.R
# -------------------------------------------------------------------
models.updated <- c("lpj.guess", "lpj.wsl", "ed2", "ed2.lu", "jules.stat", "jules.triffid", "", "sibcasa", "clm.cn", "clm.bgc")

# Open the .nc files I saved from the two "MIP CODE" scripts
gpp <- nc_open(file.path(outputs, "GPP.annual.nc"))
agb <- nc_open(file.path(outputs, "AGB.annual.nc"))
lai <- nc_open(file.path(outputs, "LAI.annual.nc"))
npp <- nc_open(file.path(outputs, "NPP.annual.nc"))
nee <- nc_open(file.path(outputs, "NEE.annual.nc"))
auto.resp   <- nc_open(file.path(outputs, "RespirationAuto.annual.nc"))
hetero.resp <- nc_open(file.path(outputs, "RespirationHetero.annual.nc"))
fire        <- nc_open(file.path(outputs, "Fire.annual.nc"))
soilmoist   <- nc_open(file.path(outputs, "SoilMoist.annual.nc"))
soilcarb    <- nc_open(file.path(outputs, "TotSoilCarb.annual.nc"))
evap        <- nc_open(file.path(outputs, "Evap.annual.nc"))
# transp <- nc_open(file.path(outputs, "Transp.annual.nc"))
evergreen   <- nc_open(file.path(outputs, "FracEvergreen.annual.nc"))
deciduous   <- nc_open(file.path(outputs, "FracDeciduous.annual.nc"))
grass       <- nc_open(file.path(outputs, "FracGrass.annual.nc"))

tair    <- nc_open(file.path(outputs, "tair.annual.nc"))
precipf <- nc_open(file.path(outputs, "precipf.annual.nc"))
swdown  <- nc_open(file.path(outputs, "swdown.annual.nc"))
lwdown  <- nc_open(file.path(outputs, "lwdown.annual.nc"))
wind    <- nc_open(file.path(outputs, "wind.annual.nc"))
psurf   <- nc_open(file.path(outputs, "psurf.annual.nc"))
qair    <- nc_open(file.path(outputs, "qair.annual.nc"))


model.names <- c(ncvar_get(gpp, "ModelNames"))
years <- ncvar_get(gpp, "Year")-1000

sites <- names(gpp$var)[1:(length(names(gpp$var))-1)]

# extract the data into lists that are easy to use
#sites <- "PHA"
GPP <- AGB <- LAI <- NPP <- NEE <- RA <- RH <- FIRE <- SOILMOIST <- SOILCARB <- EVAP <- TRANSP <- EVERGREEN <- DECIDUOUS <- GRASS <- list()
TAIR <- PRECIPF <- SWDOWN <- LWDOWN <- WIND <- PSURF <- QAIR <- list()
for(i in 1:length(sites)){
	GPP[[i]]       <- data.frame(t(ncvar_get(gpp        , sites[i])))
	AGB[[i]]       <- data.frame(t(ncvar_get(agb        , sites[i])))
	LAI[[i]]       <- data.frame(t(ncvar_get(lai        , sites[i])))
	NPP[[i]]       <- data.frame(t(ncvar_get(npp        , sites[i])))
	NEE[[i]]       <- data.frame(t(ncvar_get(nee        , sites[i])))
	RA[[i]]        <- data.frame(t(ncvar_get(auto.resp  , sites[i])))
	RH[[i]]        <- data.frame(t(ncvar_get(hetero.resp, sites[i])))
	FIRE[[i]]      <- data.frame(t(ncvar_get(fire       , sites[i])))
	SOILCARB[[i]]  <- data.frame(t(ncvar_get(soilcarb, sites[i])))
	SOILMOIST[[i]] <- data.frame(t(ncvar_get(soilmoist, sites[i])))
	# TRANSP[[i]]    <- data.frame(t(ncvar_get(transp, sites[i])))
	EVAP[[i]]      <- data.frame(t(ncvar_get(evap, sites[i])))
	EVERGREEN[[i]] <- data.frame(t(ncvar_get(evergreen, sites[i])))
	DECIDUOUS[[i]] <- data.frame(t(ncvar_get(deciduous, sites[i])))
	GRASS[[i]]     <- data.frame(t(ncvar_get(grass    , sites[i])))

	TAIR[[i]]    <- data.frame(t(ncvar_get(tair   , sites[i])))
	PRECIPF[[i]] <- data.frame(t(ncvar_get(precipf, sites[i])))
	SWDOWN[[i]]  <- data.frame(t(ncvar_get(swdown , sites[i])))
	LWDOWN[[i]]  <- data.frame(t(ncvar_get(lwdown , sites[i])))
	WIND[[i]]    <- data.frame(t(ncvar_get(wind   , sites[i])))
	PSURF[[i]]   <- data.frame(t(ncvar_get(psurf  , sites[i])))
	QAIR[[i]]    <- data.frame(t(ncvar_get(qair   , sites[i])))


	names(GPP[[i]]) <- names(AGB[[i]]) <- names(LAI[[i]]) <- names(NPP[[i]]) <- names(NEE[[i]]) <- names(RA[[i]]) <- names(RH[[i]]) <- names(FIRE[[i]]) <- names(SOILCARB[[i]]) <- names(SOILMOIST[[i]]) <- names(EVAP[[i]]) <- names(TAIR[[i]]) <- names(PRECIPF[[i]]) <- names(SWDOWN[[i]]) <- names(LWDOWN[[i]]) <- names(WIND[[i]]) <- names(PSURF[[i]]) <- names(QAIR[[i]]) <- names(EVERGREEN[[i]]) <- names(DECIDUOUS[[i]]) <- names(GRASS[[i]]) <- model.names
	
	row.names(GPP[[i]]) <- row.names(AGB[[i]]) <- row.names(LAI[[i]]) <- row.names(NPP[[i]]) <- row.names(NEE[[i]]) <- row.names(RA[[i]]) <- row.names(RH[[i]]) <- row.names(SOILCARB[[i]]) <- row.names(SOILMOIST[[i]]) <- row.names(EVAP[[i]]) <- row.names(TAIR[[i]]) <- row.names(PRECIPF[[i]]) <- row.names(SWDOWN[[i]]) <- row.names(LWDOWN[[i]]) <- row.names(WIND[[i]]) <- row.names(PSURF[[i]]) <- row.names(QAIR[[i]]) <- row.names(EVERGREEN[[i]]) <- row.names(DECIDUOUS[[i]]) <- row.names(GRASS[[i]]) <- years

}
names(GPP) <- names(AGB) <- names(LAI) <- names(NPP) <- names(NEE) <- names(RA) <- names(RH) <- names(FIRE) <- names(SOILCARB) <- names(SOILMOIST) <- names(EVAP) <- names(TAIR) <- names(PRECIPF) <- names(SWDOWN) <- names(LWDOWN) <- names(WIND) <- names(PSURF) <- names(QAIR) <- names(EVERGREEN) <- names(DECIDUOUS) <- names(GRASS) <- sites

nc_close(gpp); nc_close(agb); nc_close(lai); nc_close(npp); nc_close(nee);nc_close(auto.resp); nc_close(hetero.resp); nc_close(fire); nc_close(soilcarb); nc_close(soilmoist); nc_close(evap); nc_close(tair);  nc_close(precipf); nc_close(swdown); nc_close(lwdown); nc_close(wind); nc_close(psurf); nc_close(qair); nc_close(evergreen); nc_close(deciduous); nc_close(grass)

# Bind everything together into a single data frame to write as a csv
for(i in 1:length(sites)){
	gpp.df <- stack(GPP[[i]])
	names(gpp.df) <- c("GPP", "Model")
	gpp.df$Site <- as.factor(names(GPP)[i])
	gpp.df$Year <- years
	gpp.df <- gpp.df[,c("Model", "Site", "Year", "GPP")]

	agb.df <- stack(AGB[[i]])
	names(agb.df) <- c("AGB", "Model")
	agb.df$Site <- as.factor(names(AGB)[i])

	lai.df <- stack(LAI[[i]])
	names(lai.df) <- c("LAI", "Model")
	lai.df$Site <- as.factor(names(LAI)[i])

	npp.df <- stack(NPP[[i]])
	names(npp.df) <- c("NPP", "Model")
	npp.df$Site <- as.factor(names(NPP)[i])

	nee.df <- stack(NEE[[i]])
	names(nee.df) <- c("NEE", "Model")
	nee.df$Site <- as.factor(names(NEE)[i])

	ra.df <- stack(RA[[i]])
	names(ra.df) <- c("AutoResp", "Model")
	ra.df$Site <- as.factor(names(RA)[i])

	rh.df <- stack(RH[[i]])
	names(rh.df) <- c("HeteroResp", "Model")
	rh.df$Site <- as.factor(names(RH)[i])

	fire.df <- stack(FIRE[[i]])
	names(fire.df) <- c("Fire", "Model")
	fire.df$Site <- as.factor(names(FIRE)[i])

	soilcarb.df <- stack(SOILCARB[[i]])
	names(soilcarb.df) <- c("SoilCarb", "Model")
	soilcarb.df$Site <- as.factor(names(SOILCARB)[i])

	soilmoist.df <- stack(SOILMOIST[[i]])
	names(soilmoist.df) <- c("SoilMoist", "Model")
	soilmoist.df$Site <- as.factor(names(SOILMOIST)[i])

	evap.df <- stack(EVAP[[i]])
	names(evap.df) <- c("Evap", "Model")
	evap.df$Site <- as.factor(names(EVAP)[i])

	# transp.df <- stack(TRANSP[[i]])
	# names(transp.df) <- c("Transp", "Model")
	# transp.df$Site <- as.factor(names(TRANSP)[i])
	evergreen.df <- stack(EVERGREEN[[i]])
	names(evergreen.df) <- c("Evergreen", "Model")
	evergreen.df$Site <- as.factor(names(EVERGREEN)[i])

	deciduous.df <- stack(DECIDUOUS[[i]])
	names(deciduous.df) <- c("Deciduous", "Model")
	deciduous.df$Site <- as.factor(names(DECIDUOUS)[i])

	grass.df <- stack(GRASS[[i]])
	names(grass.df) <- c("Grass", "Model")
	grass.df$Site <- as.factor(names(GRASS)[i])
	
	tair.df <- stack(TAIR[[i]])
	names(tair.df) <- c("tair", "Model")
	tair.df$Site <- as.factor(names(TAIR)[i])
	tair.df$Year <- as.numeric(row.names(TAIR[[i]]))

	precipf.df <- stack(PRECIPF[[i]])
	names(precipf.df) <- c("precipf", "Model")
	precipf.df$Site <- as.factor(names(PRECIPF)[i])
	precipf.df$Year <- as.numeric(row.names(PRECIPF[[i]]))

	swdown.df <- stack(SWDOWN[[i]])
	names(swdown.df) <- c("swdown", "Model")
	swdown.df$Site <- as.factor(names(SWDOWN)[i])
	swdown.df$Year <- as.numeric(row.names(SWDOWN[[i]]))

	lwdown.df <- stack(LWDOWN[[i]])
	names(lwdown.df) <- c("lwdown", "Model")
	lwdown.df$Site <- as.factor(names(LWDOWN)[i])
	lwdown.df$Year <- as.numeric(row.names(LWDOWN[[i]]))

	wind.df <- stack(WIND[[i]])
	names(wind.df) <- c("wind", "Model")
	wind.df$Site <- as.factor(names(WIND)[i])
	wind.df$Year <- as.numeric(row.names(WIND[[i]]))

	psurf.df <- stack(PSURF[[i]])
	names(psurf.df) <- c("psurf", "Model")
	psurf.df$Site <- as.factor(names(PSURF)[i])
	psurf.df$Year <- as.numeric(row.names(PSURF[[i]]))

	qair.df <- stack(QAIR[[i]])
	names(qair.df) <- c("qair", "Model")
	qair.df$Site <- as.factor(names(QAIR)[i])
	qair.df$Year <- as.numeric(row.names(QAIR[[i]]))



	# write in the temp & precip data for models that didn't give it back and just assume it's fine
	for(y in years){
		tair.df[(tair.df$Model=="linkages" | tair.df$Model=="lpj.wsl") & tair.df$Year==y, "tair"] <- tair.df[tair.df$Model=="lpj.guess" & tair.df$Year==y, "tair"]
		precipf.df[(precipf.df$Model=="linkages" | precipf.df$Model=="lpj.wsl") & precipf.df$Year==y, "precipf"] <- precipf.df[precipf.df$Model=="lpj.guess" & precipf.df$Year==y, "precipf"]
		swdown.df[(swdown.df$Model=="lpj.wsl") & swdown.df$Year==y, "swdown"] <- swdown.df[swdown.df$Model=="sibcasa" & swdown.df$Year==y, "swdown"]
		lwdown.df[(lwdown.df$Model=="lpj.wsl") & lwdown.df$Year==y, "lwdown"] <- lwdown.df[lwdown.df$Model=="sibcasa" & lwdown.df$Year==y, "lwdown"]
	}	

	if(i == 1) {
		df1 <- cbind(gpp.df, agb.df$AGB, lai.df$LAI, npp.df$NPP, nee.df$NEE, ra.df$AutoResp, rh.df$HeteroResp, fire.df$Fire, soilcarb.df$SoilCarb, soilmoist.df$SoilMoist, evap.df$Evap, evergreen.df$Evergreen, deciduous.df$Deciduous, grass.df$Grass, tair.df$tair, precipf.df$precipf, swdown.df$swdown, lwdown.df$lwdown, wind.df$wind, psurf.df$psurf, qair.df$qair)
	} else {
		df1 <- rbind(df1, cbind(gpp.df, agb.df$AGB, lai.df$LAI, npp.df$NPP, nee.df$NEE, ra.df$AutoResp, rh.df$HeteroResp, fire.df$Fire, soilcarb.df$SoilCarb, soilmoist.df$SoilMoist, evap.df$Evap, evergreen.df$Evergreen, deciduous.df$Deciduous, grass.df$Grass, tair.df$tair, precipf.df$precipf, swdown.df$swdown, lwdown.df$lwdown, wind.df$wind, psurf.df$psurf, qair.df$qair))
	}
}
names(df1) <- c("Model", "Site", "Year", "GPP", "AGB", "LAI", "NPP", "NEE", "AutoResp", "HeteroResp", "Fire", "SoilCarb", "SoilMoist", "Evap", "Evergreen", "Deciduous", "Grass", "tair", "precipf", "swdown", "lwdown", "wind", "psurf", "qair")
summary(df1)


df1$Updated <- as.factor(ifelse(df1$Model %in% models.updated, "Yes", "No"))
summary(df1)

write.csv(df1, file.path(outputs, "PalEON_MIP_Yearly.csv"), row.names=F)
# -------------------------------------------------------------------


