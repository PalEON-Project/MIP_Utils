# Purpose: Utility script to easily extract variable(s) of interest from the raw output for a given model
# This function takes the raw output and returns it in the form of a list.  
#    Within the list there is an array of Years (x PFT/nsoil if necessary) x Sites for each variable
# 
# This function takes the following arguments
#    1. model -- the name of the model you're interested in; this should correspond to the PREFIX of the model directory
#    2. model.dir -- the file path to the model you're working with; you must specify this in case there are multiple versions floating around
#    3. sites -- a vector of the sites you wish to extract data for
#    4. vars  -- a vector of the variables you wish to extract for each site
#
# Note: This will by default return the entire time series at the raw time step provided by each model


extract.paleon.site <- function(model, model.dir, sites, vars){
  library(ncdf4); library(abind); library(car)
  model.out <- list()
  
  # recoding some variable names that have slight typos (that weren't worth emailing about)
  if(model %in% c("JULES", "JULES-TRIFFID")){
    vars <- recode(vars, "'TotLivBiom'='TotLivBio'")
  }
  if(model %in% c("SiBCASA")){
    vars <- recode(vars, "'Transp'='Tranp'")
  }
  
  sites.data <- vector()
  for(i in 1:length(sites)){
    out.file <- list() # create a new outfile for each site

    # In case we have the compressed files or soemthign else in the folder, only grab whatever is the directory
    if(model=="LPJ-WSL") {
      dir.site <- ""
    } else {
      dir.site <- dir(model.dir, sites[i])[which(file.info(file.path(model.dir, dir(model.dir, sites[i])))$isdir)]
    }
    
    if(length(dir.site)>1){ print(paste0("Site: ", sites[i], " Multiple directories for this site! check your file structures")); next }

    # Note: LPJ-GUESS gives us separate monthly & annual files; we want to loop through the monthly files first to extract what we can
    if(model=="LPJ-GUESS"){
      files.mo <- dir(file.path(model.dir, dir.site), c("month", ".nc"))
      files.ann <- dir(file.path(model.dir, dir.site), c("annual", ".nc"))
      
      # for LPJ-GUESS, we'll need to track which variables do not have monthly output
      vars.ann <- vector()
      for(f in 1:length(files.mo)){
        ncT <- nc_open(file.path(model.dir, dir.site, files.mo[f]))        
        for(v in vars){
          # If this variable doesn't have monthly output; save it to check annual
          if(!(v %in% names(ncT$var))){ vars.ann <- c(vars.ann, v); next}

          # Extract the data
          var.temp <- ncvar_get(ncT, v)
          # if the years are columns instead of rows, transpose the data frame
          if(length(dim(var.temp))==1) var.temp <- array(var.temp, dim=c(length(var.temp), 1))
          if(dim(var.temp)[2]>dim(var.temp)[1]) var.temp <- t(var.temp)
          
          if(f==1){ # If this is our first time through for this variable, add the layer to the list
            out.file[[v]] <- var.temp
          } else {
            out.file[[v]] <- rbind(out.file[[v]], var.temp)
          }
          
        }
        nc_close(ncT)
      }
      vars.ann <- unique(vars.ann)
      if(length(vars.ann)>0){
        for(f in 1:length(files.ann)){
          ncT <- nc_open(file.path(model.dir, dir.site, files.ann[f]))
          for(v in vars.ann){
            # If this variable doesn't have monthly output; save it to check annual
            if(!(v %in% names(ncT$var))) next # if we don't have this variable, just skip over it
            
            # Extract the data
            var.temp <- ncvar_get(ncT, v)
            # if the years are columns instead of rows, transpose the data frame
            if(length(dim(var.temp))==1) var.temp <- array(var.temp, dim=c(length(var.temp), 1))
            if(dim(var.temp)[2]>dim(var.temp)[1]) var.temp <- t(var.temp)
            
            if(f==1){ # If this is our first time through for this variable, add the layer to the list
              out.file[[v]] <- var.temp
            } else {
              out.file[[v]] <- rbind(out.file[[v]], var.temp)
            }
            
          }
          nc_close(ncT)
        }
      } # End annual loop
    # End LPJ-GUESS Special Case
    } else { 
      if(model=="LPJ-WSL"){
        files.site <- dir(file.path(model.dir, dir.site), c(sites[i], ".nc"))  
      } else {
        files.site <- dir(file.path(model.dir, dir.site), ".nc")
      }
      
      files.site <- files.site[!substr(files.site, nchar(files.site)-2, nchar(files.site))=="var"] # Make sure we don't get the weird linkages files
    
      if(!length(files.site)>0){ print(paste0("Site: ", sites[i], " There are no files for this site!")); next }
    
      for(f in 1:length(files.site)){
         ncT <- nc_open(file.path(model.dir, dir.site, files.site[f]))
         for(v in vars){
           # If this variable doesn't have monthly output; save it to check annual
           if(!(v %in% names(ncT$var))) next # if we don't have this variable, just skip over it
           
           # Extract the data
           var.temp <- ncvar_get(ncT, v)
           # if the years are columns instead of rows, transpose the data frame
           if(length(dim(var.temp))==1) var.temp <- array(var.temp, dim=c(length(var.temp), 1))
           if(dim(var.temp)[2]>dim(var.temp)[1] | (f==length(files.site) & model.name=="LINKAGES" & v %in% c("Fcomp", "AGB.pft"))) var.temp <- t(var.temp)
           
           if(f==1){ # If this is our first time through for this variable, add the layer to the list
             out.file[[v]] <- var.temp
           } else {
             out.file[[v]] <- rbind(out.file[[v]], var.temp)
           }
           
         }
         nc_close(ncT)
      }
    } # End getting the data for 1 site

    # Bind our sites together if we have multiple sites
    if(i == 1){
      for(v in names(out.file)){
        if(ncol(out.file[[v]])>1){
          model.out[[v]] <- array(out.file[[v]], dim=c(dim(out.file[[v]]),1))
        } else {
          model.out[[v]] <- out.file[[v]]
        }
      }
    } else {
      for(v in names(model.out)){
        model.out[[v]] <- abind(model.out[[v]], out.file[[v]], along=length(dim(model.out[[v]])))
      }
    } # End binding sites together
    sites.data <- c(sites.data, sites[i])
  } # End Sites
  
  for(v in names(model.out)){
    names(dimnames(model.out[[v]]))[[1]] <- "Time"
    names(dimnames(model.out[[v]]))[[length(dim(model.out[[v]]))]] <- "Site"
    dimnames(model.out[[v]])[[length(dim(model.out[[v]]))]] <- sites.data
  }
  
  # Fix any wonky names in specific models (undoing lines 19-24)
  names(model.out) <- recode(names(model.out), "'TotLivBio'='TotLivBiom'; 'Tranp'='Transp'")
  return(model.out)
} # End Function
