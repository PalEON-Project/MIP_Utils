# Purpose: Utility script to easily extract variable(s) of interest from the raw output for a given model
# This function takes the raw output and returns it in the form of a list.  
#    Within the list there is an array of lon x lat x Time (x PFT/nsoil if necessary)  for each variable
# NOTE: This will return the raw timestep of the output, which means you may get a mix of monthly and
#       yearly output from different models
# 
# This function takes the following arguments
#    1. model -- the name of the model you're interested in; this should correspond to the PREFIX of the model directory
#    2. model.dir -- the file path to where you extracted the compressed model directory; DO NOT MOVE FILES AROUND INSIDE THAT DIRECTORY!
#    3. vars  -- a vector of the variables you wish to extract for each site
#    4. xmin  -- the minimum (western0most) longitude of interest
#    5. xmax  -- the maximum (eastern-most) longitude of interest
#    6. ymin  -- the minimum (southern-most) latitude of interest
#    7. ymax  -- the maximum (northern-most) latitude of interest
#    8. yrmin -- the first year of interest; defaults to 850 (first year of simulations)
#    9. yrmax -- the last year of interest; defaults to 2010 (last year of simulations)
#
# Note: This will by default return the entire time series at the raw time step provided by each model


extract.paleon.site <- function(model, model.dir, vars, xmin=-100, xmax=-60, ymin=35, ymax=50, yrmin=850, yrmax=2010){
  library(ncdf4); library(abind); library(car)

  
  if(!model %in% c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL")) stop("Model hasn't been checked yet.  Missing PFT info")
    
  # Defining the PFT table internally just to keep everything together
  # This is ugly, but it keeps everybody on the same page
  # NOTE: These came from this google doc: https://docs.google.com/spreadsheets/d/1f0LXVnHLglppztsFbYzDdaC5TMSZYUGg4OM__MnmaDM/edit?usp=sharing
  if(model %in% c("ED2", "ED2-LU")){
    pft.lab=c("grass.c4", "tropic.early", "tropic.mid", "tropic.late", "grass.c3.temp", "pine.north", "pine.south", "conifer.late", "temp.decid.early", "temp.decid.mid", "temp.decid.late","ag1", "ag2", "ag3", "ag4","grass.c3.subtrop","Araucaria")
    soil=c(4.00, 3.00, 2.17, 1.50, 1.10, 0.80, 0.60, 0.45, 0.30, 0.20, 0.12, 0.06)
  }
  if(model == "LPJ-GUESS"){
    pft.lab=c("BNE", "BINE", "BNS", "BIBS", "TeBS", "BeIBS", "TeBE", "TrBE", "TrIBE", "TrBR", "C3G", "C4G", "Total")
    soil=c("A", "B") # Replace with depths when possible
  }
  if(model == "LPJ-WSL"){
    pft.lab=c("TrBE", "TrBR", "TeNE", "TeBE", "TeBS", "BNE", "BBS", "C3G", "C4G")
    soil=c("A", "B") # Replace with depths when possible
  }
  
  
  
  # If we're setting mins &/or maxes outside of our bounding box, reset them and send a message
  if(xmin  < -100){ warning("xmin outside PalEON MIP spatial domain; resetting to -100"); xmin=-100 }
  if(xmax  >  -60){ warning("xmax outside PalEON MIP spatial domain; resetting to -60"); ymax=-60 }
  if(ymin  <   35){ warning("ymin outside PalEON MIP spatial domain; resetting to 35"); ymin=35 }
  if(ymax  >   50){ warning("ymax outside PalEON MIP spatial domain; resetting to 50"); ymax=50 }
  if(yrmin <  850){ warning("yrmin outside PalEON MIP temporal domain; resetting to 850"); yrmin=850 }
  if(yrmax > 2010){ warning("yrmax outside PalEON MIP temporal domain; resetting to 2010"); yrmax=2010 }

  # Translating yrmin and yrmax into the rounded file names we asked for in paleon
  file.min <- ifelse(substr(yrmin,1,1)==8, 850, ifelse(substr(yrmin, 1,1)==9, 900, as.numeric(paste0(substr(yrmin,1,2),"00"))))
  file.max <- ifelse(substr(yrmax,1,1)==8, 850, ifelse(substr(yrmax, 1,1)==9, 900, as.numeric(paste0(substr(yrmax,1,2),"00"))))
  
  # Figure out which files we need for the years we're interested in
  files.use <- vector()
  if(model=="LPJ-GUESS"){ # LPJ-GUESS did not conform to the convention we wanted, so it requires a special case
    dir.month <- dir(model.dir, "month")
    dir.month <- dir.month[!substr(dir.month, nchar(dir.month)-2, nchar(dir.month))==".gz"] # exclude any compressed files
    files.all <- dir(file.path(model.dir, dir.month), ".nc")
    
    # figuring out which files to use
    files.use <- vector()
    for(i in 1:length(files.all)){
      if(as.numeric(strsplit(files.all[i], "_")[[1]][2])>=file.min & as.numeric(strsplit(files.all[i], "_")[[1]][2])<=file.max ){
        files.use <- c(files.use, i) 
      }
    }
  } else { # For everything else, the format shoudl be [MODEL]_[start year].nc
    # get list of available .nc files
    files.all <- dir(model.dir, ".nc")
    
    # Subset the files based on the years we're interested in
    files.use <- which(as.numeric(substr(files.all,nchar(files.all)-6, nchar(files.all)-3))>=file.min & as.numeric(substr(files.all,nchar(files.all)-6, nchar(files.all)-3))<=file.max)
  } # End file listing 
  
  if(!length(files.use)>0){ stop("No files meet extraction criteria") }
  
  # Making a blank list to store things in
  mod.out <- list()

  # Looping through the files!
  # NOTE: This loop works for models with 1 nc file per time step or LPJ-GUESS monthly
  for(i in files.use){
    
    if(model=="LPJ-GUESS"){
      ncT <- nc_open(file.path(model.dir, dir.month, files.all[i]))
    } else {
      ncT <- nc_open(file.path(model.dir, files.all[i]))
    } # End opening .nc file
    
    # finding out the frist year in the file
    if(model=="LPJ-GUESS"){ # LPJ-GUESS isn't following our convention so it requries a special case
      yr.file <- as.numeric(strsplit(files.all[i], "_")[[1]][2])
    } else {
      yr.file <- as.numeric(substr(files.all[i],nchar(files.all[i])-6, nchar(files.all[i])-3))
    }
    
    
    # Figure out our time indices
    lat <- ncvar_get(ncT, "lat")
    lon <- ncvar_get(ncT, "lon")
    ind.lat <- which(lat>=ymin & lat<=ymax)
    ind.lon <- which(lon>=xmin & lon<=xmax)
    
    vars.skip <- vector()
    for(v in vars){
      
      # Skip the variable if we can't find it
      if(!(v %in% names(ncT$var))) {
        if(!v %in% vars.skip){ # If we haven't already noted that this one is missing, do so now
          if(!model=="LPJ-GUESS") warning("Warning: ", v, " found!  skipped")
          vars.skip <- c(vars.skip, v) # stash it for LPJ-GUESS later
        }
        next
      } # end skip variable case
      
      dat.all <- ncvar_get(ncT, v)
      
      # Figuring out some of the dimension & stats of the data
      # NOTE: there should be 30 latitudes & 80 longitudes
      # NOTE: we have to do this by variable because pft location throws things off
      dims.dat <- dim(dat.all)
      if(length(dims.dat)<3 | length(dims.dat)>4){stop("Dimensions of output are weird! Can't continue")}
      pft <- ifelse(length(dims.dat)==4, TRUE, FALSE) # Note: I'm using PFT as a general wrapper here for pft or soil
      dim.lat <- which(dims.dat==30)
      dim.lon <- which(dims.dat==80)
      
      if(length(dim.lat)>1 | length(dim.lon)>1){ stop("2 dimensions meet lat or lon criteria! Need to stop and fix things!")}
      
      # Finding our time dimension
      if(pft == FALSE){
        # If there are no pfts, time is our remaining dimension
        dim.time <- c(1:length(dims.dat))[!c(1:length(dims.dat)) %in% c(dim.lat, dim.lon)]
      } else {
        # The PFT dims should match the length of PFT names I have
        dim.pft <- which(dims.dat==length(pft.lab) | dims.dat==length(soil))  
        if(!length(dim.pft)==1) stop("PFT/Soil dimension doesn't line up with our PFT/soils lists.  Stop and clarify (and update function)")
        
        dim.time <- c(1:length(dims.dat))[!c(1:length(dims.dat)) %in% c(dim.lat, dim.lon, dim.pft)]
      } # End Time dimension define
      
      # Reshape to a common dimensions format
      # Note: I have latitude in rows and longitude in columns so that if you take 1 slice it should look right spatially
      if(pft==FALSE){
        dat.all <- aperm(dat.all, c(dim.lat, dim.lon, dim.time)) 
      } else {
        dat.all <- aperm(dat.all, c(dim.lat, dim.lon, dim.time, dim.pft)) 
      } # end data rearranging
      
      # Figure out which times we want
      # NOTE: This will require figuring out if the time step is annual or monthly
      #       and we're going to do this a clunky way because the time step variables have not
      #       been reliable
      
      # Start index
      if(dim(dat.all)[3]%%12==0){ # if the dimension is evenly divisible by 12, we have monthly data!
        ind.min <- (max(yrmin, yr.file) - yr.file)*12+1 
      } else {
        ind.min <- (max(yrmin, yr.file) - yr.file)+1
      } # end find start index

      # End index
      if(i < files.use[length(files.use)]){ # if we're going to need another file after this, take everything
        ind.max <- dim(dat.all)[3]
      } else {
        if(dim(dat.all)[3]%%12==0){ # if the dimension is evenly divisible by 12, we have monthly data!
          ind.max <- (yrmax - yr.file + 1)*12
        } else {
          ind.max <- (yrmax - yr.file + 1)
        }
      } # end find end index
      
      
      # We now have all of our indices!! Let's extract the data! (Finally!)
      # Note, the syntax here is designed to make sure we don't drop dimensions as we go along
      if(pft==FALSE){
        dat.temp <- array(dat.all[ind.lat,ind.lon,ind.min:ind.max], dim=c(length(ind.lat), length(ind.lon), length(ind.min:ind.max)))
      } else {
        dat.temp <- array(dat.all[ind.lat,ind.lon,ind.min:ind.max,], dim=c(length(ind.lat), length(ind.lon), length(ind.min:ind.max), dim(dat.all)[4]))
      }
      
      # Write the data into our output list
      if(i == files.use[1]){ # if this is our first time through start a new layer in the list
        mod.out[[v]] <- dat.temp
      } else { # If we already have a list going, add the output to what we already have
        mod.out[[v]] <- abind(mod.out[[v]], dat.temp, along=3)
      }
    rm(dat.temp, dat.all)  
    } # end vars loop
    nc_close(ncT)
    
  } # end files loop

  
  # If this is LPJ-GUESS, lets check the annual file for leftovers
  if(model=="LPJ-GUESS"){
    files.ann <- dir(model.dir, c("ann", ".nc"))
    files.ann <- files.ann[!substr(files.ann, nchar(files.ann)-2, nchar(files.ann))==".gz"]
    
    # At the moment, only 1 annual file
    for(i in 1:length(files.ann)){
      ncT <- nc_open(file.path(model.dir, files.ann[i]))
      
      # Just to make life easier
      if(length(files.ann) == 1){ 
        yr.file <- 850 
      } else { 
        yr.file <- as.numeric(strsplit(files.all[i], "_")[[1]][2])
      }
  
      # Figure out our time indices
      lat <- ncvar_get(ncT, "lat")
      lon <- ncvar_get(ncT, "lon")
      ind.lat <- which(lat>=ymin & lat<=ymax)
      ind.lon <- which(lon>=xmin & lon<=xmax)
      
      for(v in vars.skip){
        # Skip the variable if we can't find it
        if(!(v %in% names(ncT$var))) {
          warning("Warning: ", v, " not in additional output!  will not be included in output. ")
          next
        } # end skip variable case
        
        dat.all <- ncvar_get(ncT, v)
        
        # Figuring out some of the dimension & stats of the data
        # NOTE: there should be 30 latitudes & 80 longitudes
        # NOTE: we have to do this by variable because pft location throws things off
        dims.dat <- dim(dat.all)
        if(length(dims.dat)<3 | length(dims.dat)>4){stop("Dimensions of output are weird! Can't continue")}
        pft <- ifelse(length(dims.dat)==4, TRUE, FALSE) # Is the output by pft?
        dim.lat <- which(dims.dat==30)
        dim.lon <- which(dims.dat==80)
        
        if(length(dim.lat)>1 | length(dim.lon)>1){ stop("2 dimensions meet lat or lon criteria! Need to stop and fix things!")}
        
        # Finding our time dimension
        if(pft == FALSE){
          # If there are no pfts, time is our remaining dimension
          dim.time <- c(1:length(dims.dat))[!c(1:length(dims.dat)) %in% c(dim.lat, dim.lon)]
        } else {
          # If we have PFTs, that should be our smallest dimension unless we're working in the 2000-2010 file and 
          # the output is an annual timestep (safe assumption unless we start running with 50 PFTs, which would be insane)
          if(!yr.file==2000){
            dim.pft <- which(dims.dat==min(dims.dat))  
            dim.time <- c(1:length(dims.dat))[!c(1:length(dims.dat)) %in% c(dim.lat, dim.lon, dim.pft)]
          } else {
            # if we are looking at the 2000-2010 file, the dimensions of time should be 11 for annual output, 11*12 for monthly
            # if for some reason the model is missing 2010, it should be 10 or 10*12; if two dimensions == 9, need to stop and create a special case
            dim.time <- which(dims.dat %in% c(11, 11*12, 10, 10*12))
            if(length(dim.time)>1){ stop("2 dimensions meet time criteria! Need to stop and fix things!")}
            dim.pft <- c(1:length(dims.dat))[!c(1:length(dims.dat)) %in% c(dim.lat, dim.lon, dim.time)]
          } # end 2000 file check
        } # End Time dimension define
        
        # Reshape to a common dimensions format
        # Note: I have latitude in rows and longitude in columns so that if you take 1 slice it should look right spatially
        if(pft==FALSE){
          dat.all <- aperm(dat.all, c(dim.lat, dim.lon, dim.time)) 
        } else {
          dat.all <- aperm(dat.all, c(dim.lat, dim.lon, dim.time, dim.pft)) 
        } # end data rearranging
        
        # Figure out which times we want
        # NOTE: This will require figuring out if the time step is annual or monthly
        #       and we're going to do this a clunky way because the time step variables have not
        #       been reliable
        
        # Start index
        if(dim(dat.all)[3]%%12==0){ # if the dimension is evenly divisible by 12, we have monthly data!
          ind.min <- (max(yrmin, yr.file) - yr.file)*12+1 
        } else {
          ind.min <- (max(yrmin, yr.file) - yr.file)+1
        } # end find start index
        
        # End index
        if(i < length(files.ann)){ # if we're going to need another file after this, take everything
          ind.max <- dim(dat.all)[3]
        } else {
          if(dim(dat.all)[3]%%12==0){ # if the dimension is evenly divisible by 12, we have monthly data!
            ind.max <- (yrmax - yr.file + 1)*12
          } else {
            ind.max <- (yrmax - yr.file + 1)
          }
        } # end find end index
        
        
        # We now have all of our indices!! Let's extract the data! (Finally!)
        # Note, the syntax here is designed to make sure we don't drop dimensions as we go along
        if(pft==FALSE){
          dat.temp <- array(dat.all[ind.lat,ind.lon,ind.min:ind.max], dim=c(length(ind.lat), length(ind.lon), length(ind.min:ind.max)))
        } else {
          dat.temp <- array(dat.all[ind.lat,ind.lon,ind.min:ind.max,], dim=c(length(ind.lat), length(ind.lon), length(ind.min:ind.max), dim(dat.all)[4]))
        }
        
        # Write the data into our output list
        # Only 1 annual file so this is commented out for right now
        if(i == files.use[1]){ # if this is our first time through start a new layer in the list
          mod.out[[v]] <- dat.temp
        } else { # If we already have a list going, add the output to what we already have
          mod.out[[v]] <- abind(mod.out[[v]], dat.temp, along=3)
        }
        rm(dat.temp, dat.all)  
      } # end vars loop
      nc_close(ncT)
    } # End annual loop for LPJ-GUESS
  } # End LPJ-GUESS special case
  
  
  # Add some dimension names to make life easier
  for(v in names(mod.out)){
    # Find the time stamps (in decimal years)
    if(dim(mod.out[[v]])[3]%%12==0){ 
      lab.time <- round(yrmin + (1:dim(mod.out[[v]])[3] - 1)/12, 3)
    } else {
      lab.time <- yrmin + (1:dim(mod.out[[v]])[3] - 1)
    }

    # Write in the dimnames
    if(length(dim(mod.out[[v]]))==4){
      if(dim(mod.out[[v]])[4]==length(pft.lab)){
        dimnames(mod.out[[v]]) <- list(lat=lat[ind.lat], lon=lon[ind.lon], time=lab.time, pft=pft.lab)
      } else {
        dimnames(mod.out[[v]]) <- list(lat=lat[ind.lat], lon=lon[ind.lon], time=lab.time, soil.layer=soils)
      }
      
    } else {
      dimnames(mod.out[[v]]) <- list(lat=lat[ind.lat], lon=lon[ind.lon], time=lab.time)
    }
  }
  
  
  return(mod.out)
} # End Function
