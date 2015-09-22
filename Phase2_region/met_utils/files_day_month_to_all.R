# -------------------------------------------------------------------------
# Aggregate sub-daily (6 hrly) data into daily means for those models that 
# need it (LPJ)
# Original: Christy Rollinson, crollinson@gmail.com, 21 Sept, 2015
#
# Note: This script will stick to writing to monthly files so that we have 
#       a consistent met driver format, but we'll distribute a helper script 
#       that will concatenate these files into a single .nc
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Set up libraries, file paths, useful variables, etc
# -------------------------------------------------------------------------
library(ncdf4)
library(ncdf4)
basedir <- "/projectnb/dietzelab/paleon/met_regional/bias_corr/corr_timestamp_v2/"
outpath <- "/projectnb/dietzelab/paleon/met_regional/bias_corr/corr_timestamp_v2/daily_singlefile"

if(!dir.exists(outpath)) dir.create(outpath)

vars  <- c("lwdown","precipf_corr","psurf","qair","swdown","tair","wind")
dpm   <- c(31,28,31,30,31,30,31,31,30,31,30,31) #days per month
dpm.l <- c(31,29,31,30,31,30,31,31,30,31,30,31) #leap year days per month
mv    <- 1e30    # Missing value
fillv   <- 1e+30
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Loop through each variable to do aggregation
# -------------------------------------------------------------------------
for(v in 1:length(vars)){
  print(paste("--------------------", vars[v], "--------------------", sep=" "))
  # Get a list of files (should be 1 file per month)
  files <- dir(file.path(basedir,vars[v]), ".nc")
    
  for(f in 1:length(files)){
    #format time as days since 850-01-01 midnight
    tmp  <- strsplit(files[f],"_")
    year <- as.numeric(tmp[[1]][2])
    mon  <- as.numeric(substring(tmp[[1]][3],1,2))
    print(year)

    nc.file <- nc_open(paste(basedir,vars[v],"/",files[f],sep=""))
    data.df <- ncvar_get(nc.file, vars[v])
    time.df <- ncvar_get(nc.file, "time" )    
    lat     <- ncvar_get(nc.file, "lat" )
    lon     <- ncvar_get(nc.file, "lon" )
    nc_close(nc.file)
    
    # Either create a new object or append the current data to what we alrady have
    # This should make for a more flexible framework for if our files aren't 1:length(files)

    if (exists("dat.new")) { # if we already have created the new objects to store everything, append the new info
      dat.new <- abind(dat.new, data.df, along=3) # Just adds
      nc.time <- c(nc.time, time.df)
    } else { # if the new objects don't exist, put the current data there
      dat.new <- data.df
      nc.time <- time.df
    }
   } # close file loop
    # convert decimal days to just day ID & setting up dim defs
    nc_time_units <- paste('days since 0850-01-01 00:00:00', sep='')
    dim.time      <- ncdim_def("time",nc_time_units,nc.time,unlim=TRUE)
    
    # Looping through days to find the mean
    # NOTE: na.rm=F so that if there's a problem with missing data, we'll find it
 
    # Print correct units 
    if (vars[v] == 'lwdown') {
      nc_variable_long_name=paste('Incident (downwelling) longwave ',
                                  'radiation averaged over the time step of the forcing data', sep='')
      nc_variable_units='W m-2'
    } else if (vars[v] == 'precipf') {
      nc_variable_long_name=paste('The per unit area and time ',
                                  'precipitation representing the sum of convective rainfall, ',
                                  'stratiform rainfall, and snowfall', sep='')
      nc_variable_units='kg m-2 s-1'
    } else if (vars[v] == 'psurf') {
      nc_variable_long_name='Pressure at the surface'
      nc_variable_units='Pa'
    } else if (vars[v] == 'qair') {
      nc_variable_long_name=
        'Specific humidity measured at the lowest level of the atmosphere'
      nc_variable_units='kg kg-1'
    } else if (vars[v] == 'swdown') {
      nc_variable_long_name=paste('Incident (downwelling) radiation in ',
                                  'the shortwave part of the spectrum averaged over the time ',
                                  'step of the forcing data', sep='')
      nc_variable_units='W m-2'
    } else if (vars[v] == 'tair') {
      nc_variable_long_name='2 meter air temperature'
      nc_variable_units='K'
    } else if (vars[v] == 'wind') {
      nc_variable_long_name='Wind speed'
      nc_variable_units='m s-1'
    }
    
    # Make a few dimensions we can use
    dimY <- ncdim_def( "lat", "longitude: degrees", lat )
    dimX <- ncdim_def( "lon", "latitude: degrees", lon )
     
    nc_var  <- ncvar_def(vars[v],nc_variable_units,
                            list(dimX,dimY,dim.time), fillv, longname=nc_variable_long_name,prec="double")
    
    ofname  <- paste(outpath,vars[v],"/",vars[v],'.nc',sep="")
    newfile <- nc_create( ofname, nc_var ) # Initialize file 
    
    ncatt_put( newfile, nc_var, 'days since 850', nc.time)
    ncatt_put( newfile, 0, 'description',"PalEON formatted Phase 1 met driver")
 
    ncvar_put(newfile, nc_var, dat.new) # Write netCDF file
    
    nc_close(newfile)  

} # Close variable loop
# -------------------------------------------------------------------------


