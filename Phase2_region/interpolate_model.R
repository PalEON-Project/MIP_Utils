# Function to do a very simple interpolation of ED for the QA/QC figures
# Note: This has very narrow functionality at the moment and was designed purely to make
#       the QA/QC code cleaner.  This is a very dumb, linear interpolation, but it works
#       & has very few built-in assumptions

# Arguments:
#  1. dat.all = a dataframe with all models & data added together
#  2. mod.interp = the model that needs to be interpolated
#  3. mod.mask = the model that you want to use as a mask for filling in the interpolated data
#  4. var.interp = the variable we're interpolating
interp.mod <- function(dat.all, mod.missing, mod.mask, var.interp){
  library(akima)
  
  lat <- unique(dat.all$lat)
  lon <- unique(dat.all$lon)
  
  for(x in unique(dat.all[dat.all$Model==mod.missing,"Time"])){
    dat.interp <- SpatialPointsDataFrame(coords=dat.all[!is.na(dat.all[,var.interp]) & dat.all$Model==mod.missing & dat.all$Time==x ,c("lon", "lat")], data=dat.all[!is.na(dat.all[,var.interp]) & dat.all$Model==mod.missing & dat.all$Time==x,])
    
    # Somewheres along the way, this started returning a list rather than a data frame
    interp.list <- interp(x=dat.interp, z=var.interp,
                          xo=lon, nx=length(lon),
                          yo=lat, ny=length(lat),
                          linear=T, extrap=F)
    
    interp.out <- data.frame(lon=interp.list$x, 
                             lat=rep(interp.list$y, each=length(interp.list$x)), 
                             z=stack(data.frame(interp.list$z))[,1])
    # interp.out <- data.frame(interp(x=dat.interp, z=var.interp,
    #                                 xo=lon, nx=length(lon),
    #                                 yo=lat, ny=length(lat),
    #                                 linear=T, extrap=F))

    for(i in unique(dat.all[!is.na(dat.all[,var.interp]) & dat.all$Model==mod.mask,"lat"])){
      for(j in unique(dat.all[dat.all$lat==i & !is.na(dat.all[,var.interp]) & dat.all$Model==mod.mask,"lon"])){
        if(!is.na(dat.all[dat.all$lat==i & dat.all$lon==j & dat.all$Model==mod.missing & dat.all$Time==x,var.interp])) next
        # find closest interpolated grid cell(s)
        use <- which(sqrt((j-interp.out$lon)^2 + (i-interp.out$lat)^2)==min(sqrt((j-interp.out$lon)^2 + (i-interp.out$lat)^2))) # Entering the 
        dat.all[dat.all$lat==i & dat.all$lon==j & dat.all$Model==mod.missing & dat.all$Time==x,"flag"] <- "interpolated"
        dat.all[dat.all$lat==i & dat.all$lon==j & dat.all$Model==mod.missing & dat.all$Time==x,var.interp] <- mean(interp.out[use,"z"]) # Mean helps if there's >1 cell of equal distance
      }
    }
  }
 return(dat.all) 
}