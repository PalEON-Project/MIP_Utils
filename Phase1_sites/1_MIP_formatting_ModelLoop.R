# Doing some EDA to make sure the ED runs are at least somewhat on par with others
library(ncdf4)
library(car)
setwd("~/Desktop/Research/PalEON_CR/PalEON_MIP_Site/")

# ------------------------------------------------
# Setting up to compare the inital data from the models
# ------------------------------------------------
model.dir <- "phase1a_model_output"
#model.dir <- "phase1a_model_output/"

#~/Desktop/PalEON_CR/PalEON_MIP_Site/phase1a_model_output
# Models for which we have data
model.list <- dir(model.dir)
model.list

# Sites
site.list <- c("PHA", "PHO", "PUN", "PBL", "PDL", "PMB")
#site.list <- c("PHA", "PBL", "PDL", "PMB")
# useful numbers
yr2sec <- 1/(365*24*60*60)
mo2sec <- 1/(12*24*60*60)

# ------------------------------------------------------------------------
# Extracting Variables names to make life easier
# ------------------------------------------------------------------------
# Setting up directories to pull an example file
dir.ed <- file.path(model.dir, "ED2.v7", site.list[1])
files.ed <- dir(dir.ed)

dir.ed.lu <- file.path(model.dir, "ED2-LU.v8", site.list[1])
files.ed.lu <- dir(dir.ed.lu)

dir.clm.bgc <- file.path(model.dir, "CLM45BGC.v5.1", paste0(site.list[1], ".CLM45BGC"))
dir.clm.cn <- file.path(model.dir, "CLM45CN.v3.1", paste0(site.list[1], ".CLM45CN"))
files.clm.bgc <- dir(dir.clm.bgc)
files.clm.cn <- dir(dir.clm.cn)

dir.lpj.g <- file.path(model.dir, "LPJ-GUESS.v6", paste(site.list[1], "LPJ-GUESS", sep="_"))
files.lpj.g <- dir(dir.lpj.g)
index <- gregexpr("month",files.lpj.g[2])[[1]][1] # LPJ-GUESS has separate annual and monthly files & we just want the monthly
files.lpj.g.m <- files.lpj.g[substr(files.lpj.g, index, index+4)=="month"]
files.lpj.g.y <- files.lpj.g[substr(files.lpj.g, index, index+5)=="annual"]

dir.lpj.w <- file.path(model.dir, "LPJ-WSL.v5")
files.lpj.w <- dir(dir.lpj.w)

dir.jules.s <- file.path(model.dir, "JULES.v2", paste(site.list[1], "JULES_STATIC", sep="_"))
files.jules.s <- dir(dir.jules.s)

dir.jules.triff <- file.path(model.dir, "JULES_TRIFFID.v1", paste(site.list[1], "JULES_TRIFFID", sep="_"))
files.jules.triff <- dir(dir.jules.triff)

dir.linkages <- file.path(model.dir, "LINKAGES.v1.3", paste(site.list[1], "LINKAGES", sep="_"))
files.linkages <- dir(dir.linkages, ".nc")

dir.sib <- file.path(model.dir, "SiBCASA.v1", paste(site.list[1], "SiBCASA", sep="_"))
files.sib <- dir(dir.sib, ".nc")

# Opening an example file from each model
ed          <- nc_open(file.path(dir.ed, files.ed[1]))
ed.lu       <- nc_open(file.path(dir.ed.lu, files.ed.lu[1]))
clm.bgc     <- nc_open(file.path(dir.clm.bgc, files.clm.bgc[1]))
clm.cn      <- nc_open(file.path(dir.clm.cn, files.clm.cn[1]))
lpj.g.m     <- nc_open(file.path(dir.lpj.g, files.lpj.g.m[1]))
lpj.g.y     <- nc_open(file.path(dir.lpj.g, files.lpj.g.y[1]))
lpj.w       <- nc_open(file.path(dir.lpj.w, paste(site.list[1], "LPJ-wsl.850.nc", sep=".")))
jules.s     <- nc_open(file.path(dir.jules.s, files.jules.s[1]))
jules.triff <- nc_open(file.path(dir.jules.triff, files.jules.triff[1]))
linkages    <- nc_open(file.path(dir.linkages, files.linkages[1]))
sib         <- nc_open(file.path(dir.sib, files.sib[1]))

# extracting variable names
ed.var <- names(ed$var)
ed.lu.var <- names(ed.lu$var)
clm.bgc.var <- names(clm.bgc$var)
clm.cn.var <- names(clm.bgc$var)
lpj.g.var.m <- names(lpj.g.m$var)
lpj.g.var.y <- names(lpj.g.y$var)[!(names(lpj.g.y$var) %in% lpj.g.var.m)] # only take yearly what we can't get monthly
lpj.g.var <- c(lpj.g.var.m, lpj.g.var.y)
lpj.w.var <- names(lpj.w$var)
jules.s.var <- names(jules.s$var)[4:length(jules.s$var)]
jules.s.var2 <- c("TotLivBiom", jules.s.var[2:length(jules.s.var)])
jules.triff.var <- names(jules.triff$var)[4:length(jules.triff$var)]
jules.triff.var2 <- recode(jules.triff.var, "'TotLivBio'='TotLivBiom'")
linkages.var <- names(linkages$var)
sib.var <- names(sib$var)
sib.var2 <- recode(sib.var, "'Tranp'='Transp'")

# PFT-level variables need to be dealt with slightly differently than single-string variables
var.diversity <- c("BA", "Dens", "Fcomp", "PFT", "fpc", "pft-vegc", "pft-lai", "pft-npp", "pft-diam", "pft-height", "nind", "estrate")

# summary(clm$var)
# summary(ncvar_get(clm, "Fcomp"))
#ncvar_get(clm, "pft")

# -----------------------------------
# Soil variables have different layers and need to be indexed accordingly
# -----------------------------------
#   These indices lump things into the 0.5 and 1.5 m depths of lpj.g
soil.var <- c("SoilDepth", "SoilMoist", "SoilTemp")

soil.ed <- ncvar_get(ed, "SoilDepth")
soil.ed.5 <- which(abs(soil.ed)<=0.5); vol.ed <- vector(length=length(soil.ed))
for(i in 1:(length(soil.ed)-1)){
  vol.ed[length(soil.ed)] <- abs(soil.ed[length(soil.ed)])
  vol.ed[i] <- abs(abs(soil.ed[i]) - abs(soil.ed[i+1]))
  }

soil.ed.lu <- ncvar_get(ed.lu, "SoilDepth")
soil.ed.lu.5 <- which(abs(soil.ed.lu)<=0.5); vol.ed.lu <- vector(length=length(soil.ed.lu))
for(i in 1:(length(soil.ed.lu)-1)){
  vol.ed.lu[length(soil.ed.lu)] <- abs(soil.ed.lu[length(soil.ed.lu)])
  vol.ed.lu[i] <- abs(abs(soil.ed.lu[i]) - abs(soil.ed.lu[i+1]))
  }


soil.clm.bgc <- ncvar_get(clm.bgc, "SoilDepth")
soil.clm.bgc.5 <- which(abs(soil.clm.bgc)<=0.5); vol.clm.bgc <- vector(length=length(soil.clm.bgc))
for(i in 1:(length(soil.clm.bgc)-1)){
  vol.clm.bgc[length(soil.clm.bgc)] <- abs(soil.clm.bgc[length(soil.clm.bgc)])
  vol.clm.bgc[i] <- abs(abs(soil.clm.bgc[i]) - abs(soil.clm.bgc[i+1]))
}

soil.clm.cn <- ncvar_get(clm.cn, "SoilDepth")
soil.clm.cn.5 <- which(abs(soil.clm.cn)<=0.5); vol.clm.cn <- vector(length=length(soil.clm.cn))
for(i in 1:(length(soil.clm.cn)-1)){
  vol.clm.cn[length(soil.clm.cn)] <- abs(soil.clm.cn[length(soil.clm.cn)])
  vol.clm.cn[i] <- abs(abs(soil.clm.cn[i]) - abs(soil.clm.cn[i+1]))
}

soil.lpj.g <- ncvar_get(lpj.g.m, "SoilDepth")
soil.lpj.g.5 <- which(abs(soil.lpj.g)<=0.5); vol.lpj.g <- vector(length=length(soil.lpj.g))
for(i in 1:(length(soil.lpj.g)-1)){
  vol.lpj.g[length(soil.lpj.g)] <- abs(soil.lpj.g[length(soil.lpj.g)])
  vol.lpj.g[i] <- abs(abs(soil.lpj.g[i]) - abs(soil.lpj.g[i+1]))
}


soil.lpj.w <- ncvar_get(lpj.w, "soil.depths")
soil.lpj.w.5 <- which(abs(soil.lpj.w)<=0.5); vol.lpj.w <- vector(length=length(soil.lpj.w))
for(i in 1:(length(soil.lpj.w)-1)){
  vol.lpj.w[length(soil.lpj.w)] <- abs(soil.lpj.w[length(soil.lpj.w)])
  vol.lpj.w[i] <- abs(abs(soil.lpj.w[i]) - abs(soil.lpj.w[i+1]))
}

# soil.jules.s <- ncvar_get(jules.s, "soil.depths")
vol.jules.s <- c(0.1, 0.25, 0.65, 2)
soil.jules.s <- vol.jules.s[1]
for(i in 2:length(vol.jules.s)){
	soil.jules.s[i] <- soil.jules.s[i-1] + vol.jules.s[i]
} 
soil.jules.s.5 <- which(abs(soil.jules.s)<=0.5)


# soil.jules.triff <- ncvar_get(jules.triff, "soil.depths")
vol.jules.triff <- c(0.1, 0.25, 0.65, 2)
soil.jules.triff <- vol.jules.triff[1]
for(i in 2:length(vol.jules.triff)){
	soil.jules.triff[i] <- soil.jules.triff[i-1] + vol.jules.triff[i]
} 
soil.jules.triff.5 <- which(abs(soil.jules.triff)<=0.5)

soil.sib <- ncvar_get(sib, "SoilDepth")
soil.sib.5 <- which(abs(soil.sib)<=0.5); vol.sib <- vector(length=length(soil.sib))
for(i in 1:(length(soil.sib)-1)){
  vol.sib[length(soil.sib)] <- abs(soil.sib[length(soil.sib)])
  vol.sib[i] <- abs(abs(soil.sib[i]) - abs(soil.sib[i+1]))
}


# Closing files
nc_close(ed); 
nc_close(ed.lu); 
nc_close(clm.bgc);
nc_close(clm.cn) 
nc_close(lpj.g.m); 
nc_close(lpj.g.y); 
nc_close(lpj.w); 
nc_close(jules.s)
nc_close(jules.triff)
nc_close(linkages)
nc_close(sib)
# ------------------------------------------------------------------------
# EXTRACTING MODEL OUTPUTS
# ------------------------------------------------------------------------
# -----------------------------------
# ED 2.1
# -----------------------------------
ed <- list()
ed.diversity <- list()
for(s in 1:length(site.list)){
  dir.ed <- file.path(model.dir, "ED2.v7", site.list[s])
  files.ed <- dir(dir.ed)
  
  #  nee.temp <- npp.temp <- rh.temp <- ah.temp <- gpp.temp <- vector()
  ed.var.list <- list()
  div.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.ed)){
    ncMT <- nc_open(file.path(dir.ed, files.ed[i]))
    for(v in 1:length(ed.var)){
      if(i == 1) temp <- vector() else temp <- ed.var.list[[v]]
      if(ed.var[v] %in% var.diversity[1:3]){
      	temp <- c(temp, colSums(ncvar_get(ncMT, ed.var[[v]])))
      } else if(ed.var[v] %in% soil.var[2:3]){
        soil.temp <- t(ncvar_get(ncMT, ed.var[v]))[,soil.ed.5]
        for(q in 1:ncol(soil.temp)){
        	soil.temp[,q] <- soil.temp[,q]* vol.ed[soil.ed.5[q]]/sum(vol.ed[soil.ed.5])
        }
        temp <- c(temp, rowSums(soil.temp))
      } else {      
      temp <- c(temp, ncvar_get(ncMT, ed.var[v])) }
      ed.var.list[[v]] <- temp
    }
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# ----------------------
    if(i == 1){ evg <- decid <- grass <- vector() 
    } else { 
	  evg   <- ed.var.list[["Evergreen"]]
	  decid <- ed.var.list[["Deciduous"]]
	  grass <- ed.var.list[["Grass"]]
    }

    ed.var.list[["Evergreen"]] <- c(evg  , colSums(ncvar_get(ncMT, "Fcomp")[6:8 ,]))
    ed.var.list[["Deciduous"]] <- c(decid, colSums(ncvar_get(ncMT, "Fcomp")[9:11,]))
    ed.var.list[["Grass"    ]] <- c(grass, colSums(ncvar_get(ncMT, "Fcomp")[c(1,5,12:16),]))
    # ----------------------
    nc_close(ncMT)      
  }
  names(ed.var.list) <- c(ed.var, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:(length(ed.var)+3)){
    if(s == 1){
      ed[[v]] <- data.frame(ed.var.list[[v]]) 
    } else {
      ed[[v]][,s] <- ed.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(ed) <- c(ed.var, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(ed)){
  names(ed[[i]]) <- site.list
}
# -----------------------------------


# -----------------------------------
# ed.lu 2.1 - Land Use
# -----------------------------------
ed.lu <- list()
ed.lu.diversity <- list()
for(s in 1:length(site.list)){
  dir.ed.lu <- file.path(model.dir, "ED2-LU.v8", site.list[s])
  files.ed.lu <- dir(dir.ed.lu)
  
  #  nee.temp <- npp.temp <- rh.temp <- ah.temp <- gpp.temp <- vector()
  ed.lu.var.list <- list()
  div.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.ed.lu)){
    ncMT <- nc_open(file.path(dir.ed.lu, files.ed.lu[i]))
    for(v in 1:length(ed.lu.var)){
      if(i == 1) temp <- vector() else temp <- ed.lu.var.list[[v]]
      if(ed.lu.var[v] %in% var.diversity[1:3]){
      	temp <- c(temp, colSums(ncvar_get(ncMT, ed.lu.var[[v]])))
      } else if(ed.lu.var[v] %in% soil.var[2:3]){
        soil.temp <- t(ncvar_get(ncMT, ed.lu.var[v]))[,soil.ed.lu.5]
        for(q in 1:ncol(soil.temp)){
        	soil.temp[,q] <- soil.temp[,q]* vol.ed.lu[soil.ed.lu.5[q]]/sum(vol.ed.lu[soil.ed.lu.5])
        }
        temp <- c(temp, rowSums(soil.temp))
      } else {      
      temp <- c(temp, ncvar_get(ncMT, ed.lu.var[v])) }
      ed.lu.var.list[[v]] <- temp
    }
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# ----------------------
    if(i == 1){ evg <- decid <- grass <- vector() 
    } else { 
	  evg   <- ed.lu.var.list[["Evergreen"]]
	  decid <- ed.lu.var.list[["Deciduous"]]
	  grass <- ed.lu.var.list[["Grass"]]
    }

    ed.lu.var.list[["Evergreen"]] <- c(evg  , colSums(ncvar_get(ncMT, "Fcomp")[6:8 ,]))
    ed.lu.var.list[["Deciduous"]] <- c(decid, colSums(ncvar_get(ncMT, "Fcomp")[9:11,]))
    ed.lu.var.list[["Grass"    ]] <- c(grass, colSums(ncvar_get(ncMT, "Fcomp")[c(1,5,12:16),]))
    # ----------------------

    nc_close(ncMT)      
  }
  names(ed.lu.var.list) <- c(ed.lu.var, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(ed.lu.var.list)){
    if(s == 1){
      ed.lu[[v]] <- data.frame(ed.lu.var.list[[v]]) 
    } else {
      ed.lu[[v]][,s] <- ed.lu.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(ed.lu) <- c(ed.lu.var, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(ed.lu)){
  names(ed.lu[[i]]) <- site.list
} 
# -----------------------------------


# -----------------------------------
# CLM 4.5 BGC 
# -----------------------------------
clm.bgc <- list() 
for(s in 1:length(site.list)){
  dir.clm.bgc <- file.path(model.dir, "CLM45BGC.v5.1", paste0(site.list[s], ".CLM45BGC"))
  # dir.clm.bgc <- file.path(model.dir, "CLM45.v3", site.list[s])
  files.clm.bgc <- dir(dir.clm.bgc)
  clm.bgc.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.clm.bgc)){
    ncMT <- nc_open(file.path(dir.clm.bgc, files.clm.bgc[i]))
    for(v in 1:length(clm.bgc.var)){
      if(i == 1) temp <- vector() else temp <- clm.bgc.var.list[[v]]
      if(clm.bgc.var[v] %in% var.diversity[1:3]){
      	temp <- c(temp, rowSums(ncvar_get(ncMT, clm.bgc.var[[v]])))
      } else if(clm.bgc.var[v] %in% soil.var[2:3]){
        soil.temp <- (ncvar_get(ncMT, clm.bgc.var[v]))[,soil.clm.bgc.5]
        for(q in 1:ncol(soil.temp)){
        	soil.temp[,q] <- soil.temp[,q]* vol.clm.bgc[soil.clm.bgc.5[q]]/sum(vol.clm.bgc[soil.clm.bgc.5])
        }
        temp <- c(temp, rowSums(soil.temp))
      } else {      
        temp <- c(temp, ncvar_get(ncMT, clm.bgc.var[v]))  }
      clm.bgc.var.list[[v]] <- temp
    }    
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# ----------------------
    if(i == 1){ evg <- decid <- grass <- vector() 
    } else { 
	  evg   <- clm.bgc.var.list[["Evergreen"]]
	  decid <- clm.bgc.var.list[["Deciduous"]]
	  grass <- clm.bgc.var.list[["Grass"]]
    }

    clm.bgc.var.list[["Evergreen"]] <- c(evg  , rowSums(ncvar_get(ncMT, "Fcomp")[,c(2:3, 5:6)]))
    clm.bgc.var.list[["Deciduous"]] <- c(decid, rowSums(ncvar_get(ncMT, "Fcomp")[,c(4,7:9)]))
    clm.bgc.var.list[["Grass"    ]] <- c(grass, rowSums(ncvar_get(ncMT, "Fcomp")[,c(13:17)]))
    # ----------------------
    nc_close(ncMT)      
  }
  names(clm.bgc.var.list) <- c(clm.bgc.var, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(clm.bgc.var.list)){
    if(s == 1){
      clm.bgc[[v]] <- data.frame(clm.bgc.var.list[[v]]) 
    } else {
      clm.bgc[[v]][,s] <- clm.bgc.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(clm.bgc) <- c(clm.bgc.var, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(clm.bgc)){
  names(clm.bgc[[i]]) <- site.list
}
# -----------------------------------


# -----------------------------------
# CLM 4.5 CN 
# -----------------------------------
clm.cn <- list() 
for(s in 1:length(site.list)){
  dir.clm.cn <- file.path(model.dir, "CLM45CN.v3.1", paste0(site.list[s], ".CLM45CN"))
  files.clm.cn <- dir(dir.clm.cn)
  clm.cn.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.clm.cn)){
    ncMT <- nc_open(file.path(dir.clm.cn, files.clm.cn[i]))
    for(v in 1:length(clm.cn.var)){
      if(i == 1) temp <- vector() else temp <- clm.cn.var.list[[v]]
      if(clm.cn.var[v] %in% var.diversity[1:3]){
      	temp <- c(temp, rowSums(ncvar_get(ncMT, clm.cn.var[[v]])))
      } else if(clm.cn.var[v] %in% soil.var[2:3]){
        soil.temp <- (ncvar_get(ncMT, clm.cn.var[v]))[,soil.clm.cn.5]
        for(q in 1:ncol(soil.temp)){
        	soil.temp[,q] <- soil.temp[,q]* vol.clm.cn[soil.clm.cn.5[q]]/sum(vol.clm.cn[soil.clm.cn.5])
        }
        temp <- c(temp, rowSums(soil.temp))
      } else {      
        temp <- c(temp, ncvar_get(ncMT, clm.cn.var[v]))  }
      clm.cn.var.list[[v]] <- temp
    }    
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# ----------------------
    if(i == 1){ evg <- decid <- grass <- vector() 
    } else { 
	  evg   <- clm.cn.var.list[["Evergreen"]]
	  decid <- clm.cn.var.list[["Deciduous"]]
	  grass <- clm.cn.var.list[["Grass"]]
    }

    clm.cn.var.list[["Evergreen"]] <- c(evg  , rowSums(ncvar_get(ncMT, "Fcomp")[,c(2:3, 5:6)]))
    clm.cn.var.list[["Deciduous"]] <- c(decid, rowSums(ncvar_get(ncMT, "Fcomp")[,c(4,7:9)]))
    clm.cn.var.list[["Grass"    ]] <- c(grass, rowSums(ncvar_get(ncMT, "Fcomp")[,c(13:17)]))
    # ----------------------
    nc_close(ncMT)      
  }
  names(clm.cn.var.list) <- c(clm.cn.var, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(clm.cn.var.list)){
    if(s == 1){
      clm.cn[[v]] <- data.frame(clm.cn.var.list[[v]]) 
    } else {
      clm.cn[[v]][,s] <- clm.cn.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(clm.cn) <- c(clm.cn.var, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(clm.cn)){
  names(clm.cn[[i]]) <- site.list
}
# -----------------------------------


# -----------------------------------
# LPJ-GUESS
# -----------------------------------
lpj.g <- list()
lpj.pft <- c(which(lpj.g.var.y=="AGB"), which(lpj.g.var.y=="TotLivBiom"))
for(s in 1:length(site.list)){
  dir.lpj.g <- file.path(model.dir, "LPJ-GUESS.v6", paste(site.list[s], "LPJ-GUESS", sep="_"))
  files.lpj.g <- dir(dir.lpj.g)
  
  index <- gregexpr("month",files.lpj.g[2])[[1]][1] # LPJ-GUESS has separate annual and monthly files & we just want the monthly
  files.lpj.g.m <- files.lpj.g[substr(files.lpj.g, index, index+4)=="month"]
  files.lpj.g.y <- files.lpj.g[substr(files.lpj.g, index, index+5)=="annual"]
  
  lpj.g.var.list <- list()
  #-----------------------------------
  # Monthly Variables
  for(i in 1:length(files.lpj.g.m)){
    ncMT <- nc_open(file.path(dir.lpj.g, files.lpj.g.m[i]))
    for(v in 1:length(lpj.g.var.m)){
      if(i == 1) temp <- vector() else temp <- lpj.g.var.list[[v]]
      if(lpj.g.var.m[v] %in% var.diversity[1:3]){
      	temp <- c(temp, colSums(ncvar_get(ncMT, lpj.g.var.m[[v]])))
      } else if(lpj.g.var[v] %in% soil.var[2]){
        temp <- c(temp, ncvar_get(ncMT, lpj.g.var.m[v])[soil.lpj.g.5,])
      } else {      
        temp <- c(temp, ncvar_get(ncMT, lpj.g.var.m[v]))  }
      lpj.g.var.list[[v]] <- temp
    }    
    nc_close(ncMT)      
  }
  # Yearly Variables
  for(i in 1:length(files.lpj.g.y)){
    ncMT <- nc_open(file.path(dir.lpj.g, files.lpj.g.y[i]))
    for(v in 1:length(lpj.g.var.y)){
      if(i == 1) temp <- vector() else temp <- lpj.g.var.list[[v+length(lpj.g.var.m)]] # this tells us to go past monthly vars
      if(v %in% lpj.pft | lpj.g.var.y[v] %in% var.diversity[1:3]){
        temp <- c(temp, ncvar_get(ncMT, lpj.g.var.y[v])[13,])
      } else {
      temp <- c(temp, ncvar_get(ncMT, lpj.g.var.y[v]))  }
      lpj.g.var.list[[v+length(lpj.g.var.m)]] <- temp
    }    
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# ----------------------
    if(i == 1){ evg <- decid <- grass <- vector() 
    } else { 
	  evg   <- lpj.g.var.list[["Evergreen"]]
	  decid <- lpj.g.var.list[["Deciduous"]]
	  grass <- lpj.g.var.list[["Grass"]]
    }

    lpj.g.var.list[["Evergreen"]] <- c(evg  , colSums(ncvar_get(ncMT, "Fcomp")[c(1:2, 7:9),]))
    lpj.g.var.list[["Deciduous"]] <- c(decid, colSums(ncvar_get(ncMT, "Fcomp")[c(3:6, 10),]))
    lpj.g.var.list[["Grass"    ]] <- c(grass, colSums(ncvar_get(ncMT, "Fcomp")[c(11:12),]))
    # ----------------------
    nc_close(ncMT)      
  }
  names(lpj.g.var.list) <- c(lpj.g.var, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(lpj.g.var.list)){
    if(s == 1){
      lpj.g[[v]] <- data.frame(lpj.g.var.list[[v]]) 
    } else {
      lpj.g[[v]][,s] <- lpj.g.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(lpj.g) <- c(lpj.g.var, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(lpj.g)){
  names(lpj.g[[i]]) <- site.list
}
# -----------------------------------


# -----------------------------------
# LPJ-WSL
# -----------------------------------
lpj.w <- list()
var.pft.lpj.w <- c("LAI", "NPP")
for(s in 1:length(site.list)){
  dir.lpj.w <- file.path(model.dir, "LPJ-WSL.v5")
  files.lpj.w <- dir(dir.lpj.w)
  lpj.w.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  ncMT <- nc_open(file.path(dir.lpj.w, paste(site.list[s], "LPJ-wsl.850.nc", sep=".")))
  for(v in 1:length(lpj.w.var)){
	if(lpj.w.var[v] == "NPP" | lpj.w.var[v] %in% var.diversity[1:3]){
	lpj.w.var.list[[v]] <- rowSums(ncvar_get(ncMT, lpj.w.var[v]))
	} else {
    lpj.w.var.list[[v]] <- ncvar_get(ncMT, lpj.w.var[v])
    }    
  }
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# ----------------------
    w.fcomp <- rowSums(ncvar_get(ncMT, "Fcomp")[c(1,3:4,6),])
    lpj.w.var.list[["Evergreen"]] <- rowSums(ncvar_get(ncMT, "Fcomp")[,c(1,3:4,6)])/w.fcomp
    lpj.w.var.list[["Deciduous"]] <- rowSums(ncvar_get(ncMT, "Fcomp")[,c(2,5,7)])/w.fcomp
    lpj.w.var.list[["Grass"    ]] <- rowSums(ncvar_get(ncMT, "Fcomp")[,c(8:9)])/w.fcomp
    # ----------------------
  nc_close(ncMT)      

  names(lpj.w.var.list) <- c(lpj.w.var, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(lpj.w.var.list)){
    if(s == 1){
      lpj.w[[v]] <- data.frame(lpj.w.var.list[[v]]) 
    } else {
      lpj.w[[v]][,s] <- lpj.w.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(lpj.w) <- c(lpj.w.var, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(lpj.w)){
  names(lpj.w[[i]]) <- site.list
}
# -----------------------------------


# -----------------------------------
# JULES_STATIC
# Notes for fixing the loop: 
#	- currently no soil depth
#	- NPP broken down by PFT
#	- BIOMASS = 0; see note in email
# -----------------------------------
pft.vars <- c("NPP_PFT", "LAI", "Qh", "Qle", "SnowDepth")
jules.s <- list() 
for(s in 1:length(site.list)){
  dir.jules.s <- file.path(model.dir, "JULES.v2", paste(site.list[s], "JULES_STATIC", sep="_"))
  files.jules.s <- dir(dir.jules.s)
  jules.s.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.jules.s)){
    ncMT <- nc_open(file.path(dir.jules.s, files.jules.s[i]))
    for(v in 1:length(jules.s.var)){
      if(i == 1) temp <- vector() else temp <- jules.s.var.list[[v]]
      if(jules.s.var[v] %in% pft.vars){ 
      	temp <- c(temp, rowSums(t(ncvar_get(ncMT, jules.s.var[v]))))}
      else {
      if(jules.s.var[v] %in% soil.var[2:3]){
        soil.temp <- t(ncvar_get(ncMT, jules.s.var[v]))[,soil.jules.s.5]
        for(q in 1:ncol(soil.temp)){
        	soil.temp[,q] <- soil.temp[,q]*vol.jules.s[soil.jules.s.5[q]]/sum(vol.jules.s[soil.jules.s.5])
        }
        temp <- c(temp, rowSums(soil.temp))
      } else {      
        temp <- c(temp, ncvar_get(ncMT, jules.s.var[v]))  } }
      jules.s.var.list[[v]] <- temp
     }   
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# NOTE: Because Jules isn't giving us Biomass, we're doing it off of LAI
	# NOTE: "GRASS" includes both grass & shrubs for Jules
	# ----------------------
    if(i == 1){ evg <- decid <- grass <- vector() 
    } else { 
	  evg   <- jules.s.var.list[["Evergreen"]]
	  decid <- jules.s.var.list[["Deciduous"]]
	  grass <- jules.s.var.list[["Grass"]]
    }
	lai.df  <- data.frame(t(ncvar_get(ncMT, "LAI"))) 
	lai.tot <- rowSums(lai.df)
	
    jules.s.var.list[["Evergreen"]] <- c(evg  , lai.df[,2]/lai.tot)
    jules.s.var.list[["Deciduous"]] <- c(decid, lai.df[,1]/lai.tot)
    jules.s.var.list[["Grass"    ]] <- c(grass, rowSums(lai.df[,3:5])/lai.tot)
    # ----------------------
    nc_close(ncMT)      
  }
  names(jules.s.var.list) <- c(jules.s.var2, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(jules.s.var.list)){
    if(s == 1){
      jules.s[[v]] <- data.frame(jules.s.var.list[[v]]) 
    } else {
      jules.s[[v]][,s] <- jules.s.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(jules.s) <- c(jules.s.var2, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(jules.s)){
  names(jules.s[[i]]) <- site.list
}
# -----------------------------------

# -----------------------------------
# JULES_TRIFFED (Dynamic Veg)
# Notes for fixing the loop: 
#	- currently no soil depth
#	- NPP broken down by PFT
# -----------------------------------
pft.vars <- c("NPP_PFT", "Fcomp", "TotLivBio_PFT", "Height", "LAI", "Qh", "Qle", "SnowDepth")
jules.triff <- list() 
for(s in 1:length(site.list)){
  dir.jules.triff <- file.path(model.dir, "JULES_TRIFFID.v1", paste(site.list[s], "JULES_TRIFFID", sep="_"))
  files.jules.triff <- dir(dir.jules.triff)
  jules.triff.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.jules.triff)){
    ncMT <- nc_open(file.path(dir.jules.triff, files.jules.triff[i]))
    for(v in 1:length(jules.triff.var)){
      if(i == 1) temp <- vector() else temp <- jules.triff.var.list[[v]]
      if(jules.triff.var[v] %in% pft.vars){ 
      	temp <- c(temp, rowSums(t(ncvar_get(ncMT, jules.triff.var[v]))))}
      else {
      if(jules.triff.var[v] %in% soil.var[2:3]){
        soil.temp <- t(ncvar_get(ncMT, jules.triff.var[v]))[,soil.jules.triff.5]
        for(q in 1:ncol(soil.temp)){
        	soil.temp[,q] <- soil.temp[,q]*vol.jules.triff[soil.jules.triff.5[q]]/sum(vol.jules.triff[soil.jules.triff.5])
        }
        temp <- c(temp, rowSums(soil.temp))
      } else {      
        temp <- c(temp, ncvar_get(ncMT, jules.triff.var[v]))  } }
      jules.triff.var.list[[v]] <- temp
     }   
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# NOTE: "GRASS" includes both grass & shrubs for Jules
	# ----------------------
    if(i == 1){ evg <- decid <- grass <- vector() 
    } else { 
	  evg   <- jules.triff.var.list[["Evergreen"]]
	  decid <- jules.triff.var.list[["Deciduous"]]
	  grass <- jules.triff.var.list[["Grass"]]
    }
	lai.tot <- colSums(ncvar_get(ncMT, "LAI"))

	lai.df  <- data.frame(t(ncvar_get(ncMT, "LAI"))) 
	lai.tot <- rowSums(lai.df)
	
    jules.triff.var.list[["Evergreen"]] <- c(evg  , lai.df[,2]/lai.tot)
    jules.triff.var.list[["Deciduous"]] <- c(decid, lai.df[,1]/lai.tot)
    jules.triff.var.list[["Grass"    ]] <- c(grass, rowSums(lai.df[,3:5])/lai.tot)
    # ----------------------
    nc_close(ncMT)      
  }
  names(jules.triff.var.list) <- c(jules.triff.var2, "Evergreen", "Deciduous", "Grass")
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(jules.triff.var.list)){
    if(s == 1){
      jules.triff[[v]] <- data.frame(jules.triff.var.list[[v]]) 
    } else {
      jules.triff[[v]][,s] <- jules.triff.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(jules.triff) <- c(jules.triff.var2, "Evergreen", "Deciduous", "Grass")
for(i in 1:length(jules.triff)){
  names(jules.triff[[i]]) <- site.list
}
# -----------------------------------


# -----------------------------------
# Linkages
# -----------------------------------
linkages <- list()
link.vars <- c("AGB", "TotLivBiomass", "TotSoilCarb", "GWBI", "HeteroResp", "NPP", "NEE", "Evap")

for(i in 1:length(link.vars)){
	linkages[[link.vars[i]]] <- data.frame(Year=850:2010)
}
	linkages[["Evergreen"]] <- data.frame(Year=850:2010)
	linkages[["Deciduous"]] <- data.frame(Year=850:2010)

for(s in 1:length(site.list)){
  dir.linkages <- file.path(model.dir, "LINKAGES.v1.3", paste0(site.list[s], "_LINKAGES"))
  files.linkages <- dir(dir.linkages)
  files.linkages <- files.linkages[which(as.numeric(substr(files.linkages,1,4))>=850)]
  #-----------------------------------
  # File loop extracting time series by variable group
  #-----------------------------------
  for(i in 1:length(files.linkages)){
    ncMT <- nc_open(file.path(dir.linkages, files.linkages[i]))
    for(v in 1:length(link.vars)){
    	linkages[[v]][linkages[[v]]$Year==as.numeric(substr(files.linkages[i], 1, 4)) | linkages[[v]]$Year==as.numeric(substr(files.linkages[i], 1, 4))+1,site.list[s]] <- ncvar_get(ncMT, names(linkages)[v])
    }
	# ----------------------
	# Adding in Fraction Evergreen Tree, Deciduous Tree, Grass
	# ----------------------
    # if(i == 1){ evg <- decid <- grass <- vector() 
    # } else { 
	  # evg   <- linkages[["Evergreen"]]
	  # decid <- linkages[["Deciduous"]]
	  # grass <- linkages[["Grass"]]
    # }
    	linkages[["Evergreen"]][linkages[["Evergreen"]]$Year==as.numeric(substr(files.linkages[i], 1, 4)) | linkages[["Evergreen"]]$Year==as.numeric(substr(files.linkages[i], 1, 4))+1,site.list[s]] <- sum(ncvar_get(ncMT, "Fcomp")[c(3,6:7)])
    	linkages[["Deciduous"]][linkages[["Deciduous"]]$Year==as.numeric(substr(files.linkages[i], 1, 4)) | linkages[["Deciduous"]]$Year==as.numeric(substr(files.linkages[i], 1, 4))+1,site.list[s]] <- sum(ncvar_get(ncMT, "Fcomp")[c(1:2, 4:5, 8:9)])

    # linkages[["Evergreen"]] <- c(evg  , sum(ncvar_get(ncMT, "Fcomp")[c(3,6:7),]))
    # linkages[["Deciduous"]] <- c(decid, sum(ncvar_get(ncMT, "Fcomp")[c(1:2, 4:5, 8:9),]))
    # linkages[["Grass"    ]] <- c(grass, sum(ncvar_get(ncMT, "Fcomp")[3:4,]))
    # ----------------------
    nc_close(ncMT)
  }
} # Close the model loop
# -----------------------------------

# -----------------------------------
# SiBCASA
# -----------------------------------
sib <- list()
sib.diversity <- list()
for(s in 1:length(site.list)){
  dir.sib <- file.path(model.dir, "SiBCASA.v1", paste(site.list[s], "SiBCASA", sep="_"))
  files.sib <- dir(dir.sib)
  
  #  nee.temp <- npp.temp <- rh.temp <- ah.temp <- gpp.temp <- vector()
  sib.var.list <- list()
  div.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  for(i in 1:length(files.sib)){
    ncMT <- nc_open(file.path(dir.sib, files.sib[i]))
    for(v in 1:length(sib.var)){
      if(i == 1) temp <- vector() else temp <- sib.var.list[[v]]
      if(sib.var[v] %in% var.diversity[1:3]){
      	temp <- c(temp, colSums(ncvar_get(ncMT, sib.var[[v]])))
      } else if(sib.var[v] %in% soil.var[2:3]){
        soil.temp <- t(ncvar_get(ncMT, sib.var[v]))[,soil.sib.5]
        for(q in 1:ncol(soil.temp)){
        	soil.temp[,q] <- soil.temp[,q]* vol.sib[soil.sib.5[q]]/sum(vol.sib[soil.sib.5])
        }
        temp <- c(temp, rowSums(soil.temp))
      } else {      
      temp <- c(temp, ncvar_get(ncMT, sib.var[v])) }
      sib.var.list[[v]] <- temp
    }
    nc_close(ncMT)      
  }
  names(sib.var.list) <- sib.var
  #-----------------------------------
  # Adding variable groups to master model list
  for(v in 1:length(sib.var)){
    if(s == 1){
      sib[[v]] <- data.frame(sib.var.list[[v]]) 
    } else {
      sib[[v]][,s] <- sib.var.list[[v]]
    }
  }
} # Close the model loop
# Adding site label to each variable
names(sib) <- c(sib.var)
for(i in 1:length(sib.var)){
  names(sib[[i]]) <- site.list
}
# -----------------------------------




# ------------------------------------------------------------------------
# ORGANIZING MODEL OUTPUTS BY VARIABLE
# ------------------------------------------------------------------------
names(clm.bgc) 
names(clm.cn)
names(ed) 
names(ed.lu) 
names(lpj.g) 
names(lpj.w)
names(jules.s)
names(jules.triff)
names(linkages)
names(sib)

# # ------------------------------------------------------------------------
# # ------------------------------------------------------------------------
# # Monthly Variables for all Models
# # ------------------------------------------------------------------------
# # ------------------------------------------------------------------------
# NEE <- NPP <- HeteroResp <- AutoResp <- GPP <- LAI.m  <- Qs <- Qs <- SoilMoist <- SoilTemp <- SWE <- SnowDepth <- Transp <- list()
# var1 <- c("NEE", "NPP", "HeteroResp", "AutoResp", "GPP", "LAI.m", "Evap", "Qs", "SoilMoist", "SoilTemp", "SWE", "SnowDepth", "Transp")
# for(i in 1:length(var1)){
  # assign(var1[i], list())
# }

# #dim(jules.s[["HeteroResp"]])
# #dim(ed[["HeteroResp"]])

# for(s in 1:length(site.list)){
  # # -------------------------------------
  # # Variables for which we have output from all models
  # HeteroResp[[s]] <- data.frame(cbind(clm.bcg[["HeteroResp"]][,s], ed[["HeteroResp"]][,s], lpj.g[["HeteroResp"]][,s], lpj.w[["HeteroResp"]][,s], c(jules.s[["HeteroResp"]][,s]*exp(-1.3), rep(NA, 12))))  
  # AutoResp[[s]] <- data.frame(cbind(clm[["AutoResp"]][,s], ed[["AutoResp"]][,s], lpj.g[["AutoResp"]][,s], lpj.w[["AutoResp"]][,s]), c(jules.s[["AutoResp"]][,s], rep(NA, 12)))  
  # GPP[[s]] <- data.frame(cbind(clm[["GPP"]][,s], ed[["GPP"]][,s], lpj.g[["GPP"]][,s], lpj.w[["GPP"]][,s]), c(jules.s[["GPP"]][,s], rep(NA, 12)))  
  # SoilMoist[[s]] <- data.frame(cbind(clm[["SoilMoist"]][,s], ed[["SoilMoist"]][,s], lpj.g[["SoilMoist"]][,s], lpj.w[["SoilMoist"]][,s], c(jules.s[["SoilMoist"]][,s], rep(NA, 12))))  
  # Qs[[s]] <- data.frame(cbind(clm[["Qs"]][,s], ed[["Qs"]][,s], lpj.g[["Qs"]][,s], lpj.w[["Qs"]][,s], c(jules.s[["Qs"]][,s], rep(NA, 12))))  
    
  # names(HeteroResp[[s]]) <- names(AutoResp[[s]]) <- names(GPP[[s]]) <- names(SoilMoist[[s]]) <- names(Qs[[s]]) <- c("clm45", "ed2", "lpj.guess", "lpj.wsl", "jules.stat")

  # # -------------------------------------
  # # Missing JULES_STATIC  
  # NEE[[s]] <- data.frame(cbind(clm[["NEE"]][,s], ed[["NEE"]][,s], lpj.g[["NEE"]][,s], lpj.w[["NEE"]][,s]))  
  # Evap[[s]] <- data.frame(cbind(clm[["Evap"]][,s], ed[["Evap"]][,s], lpj.g[["Evap"]][,s], lpj.w[["Evap"]][,s]))  
  # Transp[[s]] <- data.frame(cbind(clm[["Tranp"]][,s], ed[["Transp"]][,s], lpj.g[["Transp"]][,s], lpj.w[["Tranp"]][,s]))  
   
  # names(NEE[[s]]) <- names(Evap[[s]]) <- names(Transp[[s]]) <- c("clm45", "ed2", "lpj.guess", "lpj.wsl")
  
  # # -------------------------------------
  # # Missing LPJ-WSL
  # NPP[[s]] <- data.frame(cbind(clm[["NPP"]][,s], ed[["NPP"]][,s], lpj.g[["NPP"]][,s], c(jules.s[["NPP"]][,s], rep(NA, 12))))   
  # LAI.m[[s]] <- data.frame(cbind(clm[["LAI"]][,s], ed[["LAI"]][,s], lpj.g[["LAI"]][,s], c(jules.s[["LAI"]][,s], rep(NA, 12))))  
  # SoilTemp[[s]] <- data.frame(cbind(clm[["SoilTemp"]][,s], ed[["SoilTemp"]][,s], lpj.g[["SoilTemp"]][,s], c(jules.s[["SoilTemp"]][,s], rep(NA, 12))))  
  # SWE[[s]] <- data.frame(cbind(clm[["SWE"]][,s], ed[["SWE"]][,s], lpj.g[["SWE"]][,s], c(jules.s[["SWE"]][,s], rep(NA, 12))))  


  # names(NPP[[s]]) <- names(LAI.m[[s]]) <- names(SoilTemp[[s]]) <- names(SWE[[s]]) <- c("clm45", "ed2", "lpj.guess", "jules.stat")
  
  # # -------------------------------------
  # # Missing LPJ-WSL and LPJ-GUESS
  # SnowDepth[[s]] <- data.frame(cbind(clm[["SnowDepth"]][,s], ed[["SnowDepth"]][,s], c(jules.s[["SnowDepth"]][,s], rep(NA, 12))))  
  
  # names(SnowDepth[[s]]) <- c("clm45", "ed2", "jules.stat")
  
# }
# names(NEE) <- names(NPP) <- names(HeteroResp) <- names(AutoResp) <- names(GPP) <- names(LAI.m) <- 
  # names(Evap) <- names(Qs) <- names(SoilMoist) <- names(SoilTemp) <- names(SWE) <- names(SnowDepth) <- names(Transp) <- site.list
# summary(NEE)
# summary(NEE[[1]]); 
# summary(NPP[[1]]); 
# summary(HeteroResp[[1]]); 
# summary(AutoResp[[1]]); 
# summary(GPP[[1]]); 
# summary(Evap[[1]]); 
# summary(Transp[[1]])
# summary(Qs[[1]]); 
# summary(SoilMoist[[1]]); 
# summary(LAI.m[[1]])
# summary(SoilTemp[[1]]); 
# summary(SWE[[1]]); 
# summary(SnowDepth[[1]])


# summary(NEE[[1]])
# plot(NEE[[1]][,"ed2"], type="l", ylab="NEE KgC/m2/s", xlab="months since 850-01-01", main=paste(site.list[1], "NEE", sep=": "))
# #plot(NEE[[1]][,"clm45"], type="l", ylim=c(-1,1), ylab="NEE KgC/m2", xlab="months since 850-01-01", main=site.list[1])
# #lines(NEE[[1]][,"ed2"], col="green3")
# lines(NEE[[1]][,"lpj.guess"], col="lightblue", lwd=.1)
# lines(NEE[[1]][,"lpj.wsl"], col="red", lwd=.5)
# lines(NEE[[1]][,"clm45"], col="green3", lwd=.1)
# legend("bottomleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45"), col=c("black", "lightblue", "red", "green3"), lwd=2, bg="white")

# summary(NPP[[1]])
# plot(NPP[[1]][,"ed2"], type="l", ylab="NPP", xlab="months since 850-01-01", main=paste(site.list[1], "NPP", sep=": "))
# lines(NPP[[1]][,"clm45"], col="green3", lwd=0.1)
# lines(NPP[[1]][,"lpj.guess"], col="lightblue", lwd=0.1)
# lines(NPP[[1]][,"jules.stat"], col="orange3", lwd=0.1)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "CLM45", "JULES_STATIC"), col=c("black", "lightblue", "green3", "orange3"), lwd=2, bg="white")

# summary(HeteroResp[[1]])
# plot(HeteroResp[[1]][,"ed2"], type="l", ylab="Heterotrophic Respiration", xlab="months since 850-01-01", main=paste(site.list[1], "HeteroResp", sep=": "))
# lines(HeteroResp[[1]][,"jules.stat"], col="orange3", lwd=0.1)
# lines(HeteroResp[[1]][,"ed2"], col="black", lwd=0.5)
# lines(HeteroResp[[1]][,"clm45"], col="green3", lwd=0.1)
# lines(HeteroResp[[1]][,"lpj.wsl"], col="red", lwd=0.5)
# lines(HeteroResp[[1]][,"lpj.guess"], col="lightblue", lwd=0.1)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45", "JULES_STATIC"), col=c("black", "lightblue", "red", "green3", "orange3"), lwd=2, bg="white")

# summary(AutoResp[[1]])
# plot(AutoResp[[1]][,"ed2"], type="l", ylab="Autotrophic Respiration", xlab="months since 850-01-01", main=paste(site.list[1], "AutoResp", sep=": "))
# lines(AutoResp[[1]][,"lpj.wsl"], col="red", lwd=0.3)
# #lines(AutoResp[[1]][,"ed2"], col="black", lwd=0.5)
# lines(AutoResp[[1]][,"lpj.guess"], col="lightblue", lwd=0.1)
# lines(AutoResp[[1]][,"clm45"], col="green3", lwd=0.1)
# lines(AutoResp[[1]][,"jules.stat"], col="orange3", lwd=0.1)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45", "JULES_STATIC"), col=c("black", "lightblue", "red", "green3", "orange3"), lwd=2, bg="white")

# summary(Evap[[1]])
# plot(Evap[[1]][,"ed2"], type="l", ylab="Evap KgC/m2/s", xlab="months since 850-01-01", main=paste(site.list[1], "Evaporation", sep=": "), ylim=range(Evap[[1]], na.rm=T))
# #plot(Evap[[1]][,"ed2"], type="l", ylab="Evap KgC/m2/s", xlab="months since 850-01-01", main=paste(site.list[1], "Evaporation", sep=": "), ylim=c(-1e-6, 2e-4))
# lines(-Evap[[1]][,"lpj.wsl"], col="red", lwd=.5)
# lines(Evap[[1]][,"clm45"], col="green3", lwd=.1)
# lines(Evap[[1]][,"lpj.guess"], col="lightblue", lwd=.1)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45"), col=c("black", "lightblue", "red", "green3"), lwd=2, bg="white")

# summary(Transp[[1]])
# plot(Transp[[1]][,"ed2"], type="l", ylab="Transp KgC/m2/s", xlab="months since 850-01-01", main=paste(site.list[1], "Transp", sep=": "), ylim=range(Transp[[1]]))
# lines(Transp[[1]][,"clm45"], col="green3", lwd=.1)
# lines(Transp[[1]][,"lpj.wsl"], col="red", lwd=.5)
# lines(Transp[[1]][,"lpj.guess"], col="lightblue", lwd=.1)
# lines(Transp[[1]][,"ed2"], col="black", lwd=.5)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45"), col=c("black", "lightblue", "red", "green3"), lwd=2, bg="white")

# summary(Qs[[1]]) # NOTE: ED is small because I had pretty much turned it off!
# plot(Qs[[1]][,"ed2"], type="l", ylab="Qs KgC/m2/s", xlab="months since 850-01-01", ylim=range(Qs, na.rm=T), col="black", lwd=.5, main=paste(site.list[1], "Runoff", sep=": "))
# lines(Qs[[1]][,"lpj.wsl"], col="red", lwd=.5)
# lines(Qs[[1]][,"lpj.guess"], col="lightblue", lwd=.5)
# lines(Qs[[1]][,"clm45"], col="green3", lwd=.3)
# lines(Qs[[1]][,"jules.stat"], col="orange3", lwd=0.5)
# lines(Qs[[1]][,"ed2"], col="black", lwd=0.5)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45", "JULES_STATIC"), col=c("black", "lightblue", "red", "green3", "orange3"), lwd=2, bg="white")

# summary(SoilMoist[[1]]) # NOTE: ED is wrong currently
# plot(SoilMoist[[1]][,"ed2"], type="l", ylim=range(SoilMoist[[1]], na.rm=T), ylab="SoilMoist Kg/m3", xlab="months since 850-01-01", main=site.list[1])
# # plot(SoilMoist[[1]][,"clm45"], type="l", ylim=range(SoilMoist[[1]], na.rm=T), ylab="Soil Moist KgC/m2", xlab="months since 850-01-01", col="green3", lwd=.5, main=paste(site.list[1], "SoilMoist", sep=": "))
# lines(SoilMoist[[1]][,"jules.stat"], col="orange3", lwd=0.3)
# lines(SoilMoist[[1]][,"clm45"], col="green3", lwd=.3)
# lines(SoilMoist[[1]][,"lpj.guess"], col="lightblue", lwd=.5)
# lines(SoilMoist[[1]][,"ed2"], col="black", lwd=.3)
# lines(SoilMoist[[1]][,"lpj.wsl"], col="red", lwd=.3)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45", "JULES_STATIC"), col=c("black", "lightblue", "red", "green3", "orange3"), lwd=2, bg="white")


# summary(LAI.m[[1]])
# plot(LAI.m[[1]][,"ed2"], type="l", ylab="Leaf Area Index", xlab="months since 850-01-01", main=paste(site.list[1], "Monthly LAI", sep=": "))
# lines(LAI.m[[1]][,"clm45"], col="green3", lwd=0.5)
# lines(LAI.m[[1]][,"ed2"], col="black", lwd=0.5)
# lines(LAI.m[[1]][,"jules.stat"], col="orange3", lwd=0.4)
# lines(LAI.m[[1]][,"lpj.guess"], col="lightblue", lwd=0.5)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "CLM45", "JULES_STATIC"), col=c("black", "lightblue", "green3", "orange3"), lwd=2, bg="white")


# summary(SoilTemp[[1]])
# plot(SoilTemp[[1]][,"lpj.guess"], col="lightblue", type="l", ylab="SoilTemp K", xlab="months since 850-01-01", main=paste(site.list[1], "SoilTemp", sep=": "))
# lines(SoilTemp[[1]][,"jules.stat"], col="orange3", lwd=.4)
# lines(SoilTemp[[1]][,"ed2"], col="black", lwd=.35)
# lines(SoilTemp[[1]][,"clm45"], col="green3", lwd=.3)
# legend("topleft", legend=c("ED2", "LPJ-GUESS", "CLM45"), col=c("black", "lightblue", "green3"), lwd=2, bg="white")


# summary(SWE[[1]])
# plot(SWE[[1]][,"ed2"], col="black", type="l", ylab="SWE", xlab="months since 850-01-01", ylim=range(SWE[[1]][,c("ed2", "jules.stat")]*1.5, na.rm=T), main=paste(site.list[1], "SWE", sep=": "))
# lines(SWE[[1]][,"clm45"], col="green3", lwd=.3)
# lines(SWE[[1]][,"ed2"], col="black", lwd=.4)
# lines(SWE[[1]][,"jules.stat"], col="orange3", lwd=.5)
# legend("topleft", legend=c("ED2", "CLM45", "JULES_STATIC"), col=c("black", "green3", "orange3"), lwd=2, bg="white")

# summary(SnowDepth[[1]])
# plot(SnowDepth[[1]][,"ed2"], col="black", type="l", ylab="SnowDepth", xlab="months since 850-01-01", ylim=range(SnowDepth[[1]][,c("ed2", "jules.stat")]*1.5, na.rm=T), main=paste(site.list[1], "SnowDepth", sep=": "))
# lines(SnowDepth[[1]][,"clm45"], col="green3", lwd=.3)
# #lines(SnowDepth[[1]][,"ed2"], col="black", lwd=.4)
# lines(SnowDepth[[1]][,"jules.stat"], col="orange3", lwd=.5)
# legend("topleft", legend=c("ED2", "CLM45", "JULES_STATIC"), col=c("black", "green3", "orange3"), lwd=2, bg="white")

# ##------------------------------------------------------
# summary(ed)
# pdf("PrelimGraphs/MetDrivers_Tair_Monthly.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2))
# for(i in 1:ncol(ed[["tair"]])){
	# plot(ed[["tair"]][,1], type="l", ylim=range(ed[["tair"]], na.rm=T), lwd=0.25, main=names(ed[["tair"]])[i], ylab="Monthly Temp")
# }
# dev.off()

# pdf("PrelimGraphs/MetDrivers_precipf_Monthly.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2))
# for(i in 1:ncol(ed[["precipf"]])){
	# plot(ed[["precipf"]][,1], type="l", ylim=range(ed[["precipf"]], na.rm=T), lwd=0.25, main=names(ed[["precipf"]])[i], ylab="Monthly Precip")
# }
# dev.off()

# pdf("PrelimGraphs/MetDrivers_psurf_Monthly.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2))
# for(i in 1:ncol(ed[["psurf"]])){
	# plot(ed[["psurf"]][,1], type="l", ylim=range(ed[["psurf"]], na.rm=T), lwd=0.25, main=names(ed[["psurf"]])[i], ylab="Monthly Psurf")
# }
# dev.off()

# pdf("PrelimGraphs/MetDrivers_qair_Monthly.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2))
# for(i in 1:ncol(ed[["qair"]])){
	# plot(ed[["qair"]][,1], type="l", ylim=range(ed[["qair"]], na.rm=T), lwd=0.25, main=names(ed[["qair"]])[i], ylab="Monthly Qair")
# }
# dev.off()

# pdf("PrelimGraphs/MetDrivers_lwdown_Monthly.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2))
# for(i in 1:ncol(ed[["lwdown"]])){
	# plot(ed[["lwdown"]][,1], type="l", ylim=range(ed[["lwdown"]], na.rm=T), lwd=0.25, main=names(ed[["lwdown"]])[i], ylab="Monthly LWdown")
# }
# dev.off()

# pdf("PrelimGraphs/MetDrivers_swdown_Monthly.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2))
# for(i in 1:ncol(ed[["swdown"]])){
	# plot(ed[["swdown"]][,1], type="l", ylim=range(ed[["swdown"]], na.rm=T), lwd=0.25, main=names(ed[["swdown"]])[i], ylab="Monthly swdown")
# }
# dev.off()
