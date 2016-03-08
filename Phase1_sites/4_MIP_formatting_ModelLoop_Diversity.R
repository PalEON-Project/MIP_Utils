
# Doing some EDA to make sure the ED runs are at least somewhat on par with others
library(ncdf4)
library(car)
#setwd("~/Desktop/PalEON_CR/PalEON_MIP_Site/phase1a_model_output")

# ------------------------------------------------
# Setting up to compare the inital data from the models
# ------------------------------------------------
model.dir <- "~/Dropbox/PalEON_CR/PalEON_MIP_Site/phase1a_model_output"

# Models for which we have data
model.list <- dir(model.dir)

# Sites
site.list <- c("PHA", "PHO", "PUN", "PBL", "PDL", "PMB") # Removing PHA for now until I get it working with ED

# useful numbers
yr2sec <- 1/(365*24*60*60)
mo2sec <- 1/(12*24*60*60)

# ------------------------------------------------------------------------
# Extracting Variables names to make life easier
# ------------------------------------------------------------------------
# Setting up directories to pull an example file
dir.ed <- file.path(model.dir, "ED2.v7", site.list[1])
files.ed <- dir(dir.ed)

dir.ed.lu <- file.path(model.dir, "ED2-LU.v5", site.list[1])
files.ed.lu <- dir(dir.ed.lu)

dir.clm.bgc <- file.path(model.dir, "CLM-BGC.v4", site.list[1])
dir.clm.cn <- file.path(model.dir, "CLM-CN.v2", site.list[1])
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

ed.pft <- ncvar_get(ed, "PFT")
lpj.g.pft <- ncvar_get(lpj.g.y, "PFT")
lpj.w.pft <- c("TrBE", "TrBR", "TeNE", "TeBE", "BES", "BNE", "BBS", "C3", "C4")
clm.pft <- c("bare", "TeNE", "BNE", "BNS", "TrBE", "TeBE", "TrBS", "TeBS", "BBS", "SBD", "TeSBS", "BSBS", "C3_arctic", "C3", "C4", "corn", "wheat")
# jules.var <- 


# ------------------------------------------------------------------------
# EXTRACTING MODEL OUTPUTS
# ------------------------------------------------------------------------

ed.fcomp <- ed.lu.fcomp <- lpj.g.fcomp <- lpj.w.fcomp <- clm.fcomp <- list()
for(s in 1:length(site.list)){
  #-----------------------------------  
  # ED
  #-----------------------------------
  dir.ed <- file.path(model.dir, "ED2.v7", site.list[s])   
  files.ed <- dir(dir.ed)
  #ed.var.list <- list()
  # File loop extracting time series by variable group
  for(i in 1:length(files.ed)){
    ncMT <- nc_open(file.path(dir.ed, files.ed[i]))
    npft <- length(ncvar_get(ncMT, "PFT"))
    if(i == 1) ed.fcomp[[s]] <- as.data.frame(t(ncvar_get(ncMT, "Fcomp"))) 
    else ed.fcomp[[s]] <- rbind(ed.fcomp[[s]], t(ncvar_get(ncMT, "Fcomp")))
    nc_close(ncMT)      
  }
  
  #-----------------------------------  
  # ED - Land Use
  #-----------------------------------
  dir.ed.lu <- file.path(model.dir, "ED2-LU.v5", site.list[s])   
  files.ed.lu <- dir(dir.ed.lu)
  #ed.var.list <- list()
  # File loop extracting time series by variable group
  for(i in 1:length(files.ed.lu)){
    ncMT <- nc_open(file.path(dir.ed.lu, files.ed.lu[i]))
    npft <- length(ncvar_get(ncMT, "PFT"))
    if(i == 1) ed.lu.fcomp[[s]] <- as.data.frame(t(ncvar_get(ncMT, "Fcomp"))) 
    else ed.lu.fcomp[[s]] <- rbind(ed.lu.fcomp[[s]], t(ncvar_get(ncMT, "Fcomp")))
    nc_close(ncMT)      
  }
  
   #-----------------------------------
  # CLM45
  #-----------------------------------
  dir.clm <- file.path(model.dir, "CLM-BGC.v4", paste(site.list[s]))
  # dir.clm <- file.path(model.dir, "CLM45.v3", paste(site.list[s]))
  files.clm <- dir(dir.clm)  
#  clm.var.list <- list()
  for(i in 1:length(files.clm)){
    ncMT <- nc_open(file.path(dir.clm, files.clm[i]))
    npft <- nrow(ncvar_get(ncMT, "Fcomp"))
    if(i == 1) clm.fcomp[[s]] <- as.data.frame(ncvar_get(ncMT, "Fcomp"))
    else clm.fcomp[[s]] <- rbind(clm.fcomp[[s]], ncvar_get(ncMT, "Fcomp"))
    nc_close(ncMT)      
  }
 
  #-----------------------------------
  # LPJ-Guess
  #-----------------------------------
  dir.lpj.g <- file.path(model.dir, "LPJ-GUESS.v6", paste(site.list[s], "LPJ-GUESS", sep="_"))   
  # dir.lpj.g <- file.path(model.dir, "LPJ-GUESS", paste(site.list[s], "LPJ-GUESS", sep="_"))   
  files.lpj.g <- dir(dir.lpj.g)
  index <- gregexpr("month",files.lpj.g[2])[[1]][1] # LPJ-GUESS has separate annual and monthly files & we just want the monthly
  #files.lpj.g.m <- files.lpj.g[substr(files.lpj.g, index, index+4)=="month"]
  files.lpj.g.y <- files.lpj.g[substr(files.lpj.g, index, index+5)=="annual"]
  
  #lpj.g.var.list <- list()
  for(i in 1:length(files.lpj.g.y)){
    ncMT <- nc_open(file.path(dir.lpj.g, files.lpj.g.y[i]))
    npft <- nrow(ncvar_get(ncMT, "Fcomp"))
    if(i == 1) lpj.g.fcomp[[s]] <- as.data.frame(t(ncvar_get(ncMT, "Fcomp")))
    else lpj.g.fcomp[[s]] <- rbind(lpj.g.fcomp[[s]], t(ncvar_get(ncMT, "Fcomp")))
    nc_close(ncMT)      
  }
  
  #-----------------------------------
  # LPJ-WSL
  #-----------------------------------
  dir.lpj.w <- file.path(model.dir, "LPJ-WSL.v5")
  files.lpj.w <- dir(dir.lpj.w)
#  lpj.w.var.list <- list()
  #-----------------------------------
  # File loop extracting time series by variable group
  ncMT <- nc_open(file.path(dir.lpj.w, paste(site.list[s], "LPJ-wsl.850.nc", sep=".")))
  lpj.w.fcomp[[s]] <- data.frame(ncvar_get(ncMT, "Fcomp"))
  nc_close(ncMT)
  
  
  } # Close the model loop

# Adding site label to each variable
names(ed.fcomp) <- names(ed.lu.fcomp) <- names(lpj.g.fcomp) <- names(lpj.w.fcomp) <- names(clm.fcomp) <- site.list
for(i in 1:length(site.list)){
  names(ed.fcomp[[i]]) <- ed.pft
  names(ed.lu.fcomp[[i]]) <- ed.pft
  names(clm.fcomp[[i]]) <- clm.pft
  names(lpj.g.fcomp[[i]]) <- lpj.g.pft
  names(lpj.w.fcomp[[i]]) <- lpj.w.pft
}
summary(ed.fcomp[[1]])
summary(ed.lu.fcomp[[1]])
summary(clm.fcomp[[1]])
summary(lpj.g.fcomp[[1]])
summary(lpj.w.fcomp[[1]])

summary(rowSums(lpj.w.fcomp[[1]]))
summary(rowSums(lpj.g.fcomp[[1]][,1:(ncol(lpj.g.fcomp[[1]])-1)]))
summary(rowSums(clm.fcomp[[1]]))
summary(rowSums(ed.fcomp[[1]]))
summary(rowSums(ed.lu.fcomp[[1]]))


# -----------------------------------------------------
# Exploratory Graphing
# -----------------------------------------------------
#-----------------------------------  
# ED
#-----------------------------------
pft.colors.ed <- c("green3", "darkgreen", "darkgoldenrod3", "darkorange3", "red3")
pfts.ed <- c("6-North Pine", "8-Late Conifer", "9-Early Hardwood", "10-Mid Hardwood", "11-Late Hardwood")


par(mfrow=c(1,1))
plot(ed.fcomp[[1]][,6], type="l", ylim=c(0,1), col=pft.colors.ed[1], lwd=3, main=site.list[1])
lines(ed.fcomp[[1]][,c(8)], col=pft.colors.ed[2], lwd=3)
lines(ed.fcomp[[1]][,c(9)], col=pft.colors.ed[3], lwd=3)
lines(ed.fcomp[[1]][,c(10)], col=pft.colors.ed[4], lwd=3)
lines(ed.fcomp[[1]][,c(11)], col=pft.colors.ed[5], lwd=3)
legend(x=8000, y=0.6, legend=pfts.ed, col=pft.colors.ed, bty="n", bg="white", lwd=3)

pdf(width=11, height=8.5, file="../PrelimGraphs/ED_Paleon_Prelim_Fcomp.pdf")
par(mfrow=c(3,2), mar=c(3,5,1,1)+.1)
for(i in 1:length(site.list)){
plot(ed.fcomp[[i]][,6], type="l", ylim=c(0,1), col=pft.colors.ed[1], lwd=2, xlab="", ylab="Fraction of AGB")
lines(ed.fcomp[[i]][,c(8)], col=pft.colors.ed[2], lwd=1.5)
lines(ed.fcomp[[i]][,c(9)], col=pft.colors.ed[3], lwd=1.5)
lines(ed.fcomp[[i]][,c(10)], col=pft.colors.ed[4], lwd=1.5)
lines(ed.fcomp[[i]][,c(11)], col=pft.colors.ed[5], lwd=1.5)
lines(ed.fcomp[[i]][,c(5)], col="navajowhite3", lwd=2)
text(x=150, y=1, site.list[i], cex=1.5, font=2)
#legend("topleft", legend=pfts.ed, col=pft.colors.ed, bty="n", lwd=3)
}
# plot(-5, type="l", ylim=c(0,1), col=pft.colors.ed[1], lwd=0.0001, xlab="", ylab="Fraction of AGB")
legend("topright", legend=c("5-C3 Grass", pfts.ed), col=c("navajowhite3", pft.colors.ed), bty="n", lwd=5, cex=1, ncol=3)
dev.off()


summary(ed.fcomp[[1]])
yr.rows <- seq(1, nrow(ed.fcomp[[1]]), by=12)

ed.fcomp.yr <- data.frame(Year=850:2010)
for(j in 1:ncol(ed.fcomp[[1]])){
	for(i in 1:length(yr.rows)){
		if(i==1) ed.temp <- vector()
		ed.temp <- c(ed.temp, mean(ed.fcomp[[1]][yr.rows[i]:(yr.rows[i]+11),j]))
	}
	ed.fcomp.yr[,names(ed.fcomp[[1]])[j]] <- ed.temp
}
summary(ed.fcomp.yr)

names(ed.fcomp.yr)[2:ncol(ed.fcomp.yr)] <- 1:(ncol(ed.fcomp.yr)-1)

ed.pha.stack <- stack(ed.fcomp.yr[,2:ncol(ed.fcomp.yr)])[,c(2,1)]
names(ed.pha.stack) <- c("PFT", "Fcomp")
ed.pha.stack$Year <- ed.fcomp.yr$Year
# ed.pha.stack <- ed.pha.stack[as.numeric(ed.pha.stack$PFT) %in% c(6, 8:11),]
# ed.pha.stack$PFT <- as.ordered(ed.pha.stack$PFT)
ed.pha.stack$PFT <- as.factor(ifelse(nchar(as.character(ed.pha.stack$PFT))==1, paste0("0", ed.pha.stack$PFT), paste0(ed.pha.stack$PFT)))
levels(ed.pha.stack$PFT) <- names(ed.fcomp[[1]])
summary(ed.pha.stack)

pfts.use <- levels(ed.pha.stack$PFT)[c(6,8:11)]
summary(ed.pha.stack[ed.pha.stack$PFT %in% pfts.use,])
pft.colors.ed <- c("green3", "darkgreen", "darkgoldenrod3", "darkorange3", "red3")

library(ggplot2)
library(grid)
ggplot(data=ed.pha.stack[ed.pha.stack$PFT %in% pfts.use,]) +
	geom_line(aes(x=Year, y=Fcomp, color=PFT), size=2) +
	scale_y_continuous(limits=c(0,1)) +
	scale_color_manual(values=pft.colors.ed) +
	theme_bw() +
	theme(plot.title=element_text(face="bold", size=rel(3))) + theme(legend.position=c(0.6,0.9), legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(2))) + labs(color="PFT", title="ED Fractional PFT Composition",y="Fraction of Biomass") +
	guides(color=guide_legend(ncol=3)) +
	theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=rel(2)), axis.text.y=element_text(color="black", size=rel(2)), axis.title.x=element_text(face="bold", size=rel(2), vjust=-0.5),  axis.title.y=element_text(face="bold", size=rel(2), vjust=1), plot.margin=unit(c(0.1,0.5,0.5,0.1), "lines"))



#-----------------------------------  
# ED - Land Use
#-----------------------------------
pft.colors.ed <- c("green3", "darkgreen", "darkgoldenrod3", "darkorange3", "red3")
pfts.ed <- c("6-North Pine", "8-Late Conifer", "9-Early Hardwood", "10-Mid Hardwood", "11-Late Hardwood")

par(mfrow=c(1,1))
plot(ed.lu.fcomp[[1]][,6], type="l", ylim=c(0,1), col=pft.colors.ed[1], lwd=3, main=site.list[1])
lines(ed.lu.fcomp[[1]][,c(8)], col=pft.colors.ed[2], lwd=3)
lines(ed.lu.fcomp[[1]][,c(9)], col=pft.colors.ed[3], lwd=3)
lines(ed.lu.fcomp[[1]][,c(10)], col=pft.colors.ed[4], lwd=3)
lines(ed.lu.fcomp[[1]][,c(11)], col=pft.colors.ed[5], lwd=3)
legend(x=8000, y=0.6, legend=pfts.ed, col=pft.colors.ed, bty="n", bg="white", lwd=3)

pdf(width=11, height=8.5, file="PrelimGraphs/ED-LU_Paleon_Prelim_Fcomp.pdf")
par(mfrow=c(3,2), mar=c(3,5,1,1)+.1)
for(i in 1:length(site.list)){
plot(ed.lu.fcomp[[i]][,6], type="l", ylim=c(0,1), col=pft.colors.ed[1], lwd=2, xlab="", ylab="Fraction of AGB")
lines(ed.lu.fcomp[[i]][,c(8)], col=pft.colors.ed[2], lwd=1.5)
lines(ed.lu.fcomp[[i]][,c(9)], col=pft.colors.ed[3], lwd=1.5)
lines(ed.lu.fcomp[[i]][,c(10)], col=pft.colors.ed[4], lwd=1.5)
lines(ed.lu.fcomp[[i]][,c(11)], col=pft.colors.ed[5], lwd=1.5)
lines(ed.lu.fcomp[[i]][,c(5)], col="navajowhite3", lwd=2)
text(x=150, y=1, site.list[i], cex=1.5, font=2)
#legend("topleft", legend=pfts.ed, col=pft.colors.ed, bty="n", lwd=3)
}
# plot(-5, type="l", ylim=c(0,1), col=pft.colors.ed[1], lwd=0.0001, xlab="", ylab="Fraction of AGB")
legend("topright", legend=c("5-C3 Grass", pfts.ed), col=c("navajowhite3", pft.colors.ed), bty="n", lwd=5, cex=1, ncol=3)
dev.off()

#-----------------------------------  
# LPJ-GUESS
#-----------------------------------
summary(lpj.g.fcomp[[1]])
summary(lpj.g.fcomp[[2]])
summary(lpj.g.fcomp[[3]])
summary(lpj.g.fcomp[[4]])
summary(lpj.g.fcomp[[5]])
summary(lpj.g.fcomp[[6]])

pft.colors.lpj.g <- c("darkgreen", "green3", "black", "darkorange3", "darkgoldenrod3", "red3")
pfts.lpj.g <- c("BNE", "BINE", "BIBS", "TeBS", "TeIBS", "C3G")

par(mfrow=c(1,1))
plot(lpj.g.fcomp[[1]][,"BNE"], type="l", ylim=c(0,1), col=pft.colors.lpj.g[1], lwd=3, main=site.list[1])
lines(lpj.g.fcomp[[1]][,"BINE"], col=pft.colors.lpj.g[2], lwd=3)
lines(lpj.g.fcomp[[1]][,"BIBS"], col=pft.colors.lpj.g[3], lwd=3)
lines(lpj.g.fcomp[[1]][,"TeBS"], col=pft.colors.lpj.g[4], lwd=3)
lines(lpj.g.fcomp[[1]][,"TeIBS"], col=pft.colors.lpj.g[5], lwd=3)
lines(lpj.g.fcomp[[1]][,"C3G"], col=pft.colors.lpj.g[6], lwd=3)
legend(x=50, y=0.8, legend=pfts.lpj.g, col=pft.colors.lpj.g, bty="n", lwd=3)

pdf(width=11, height=8.5, file="PrelimGraphs/LPJ-GUESS_Paleon_Prelim_Fcomp.pdf")
par(mfrow=c(3,2), mar=c(3,5,1,1)+.1)
for(i in 1:length(site.list)){
  plot(lpj.g.fcomp[[i]][,"BNE"], type="l", ylim=c(0,1), col=pft.colors.lpj.g[1], lwd=3, ylab="Fraction of AGB")
  lines(lpj.g.fcomp[[i]][,"BINE"], col=pft.colors.lpj.g[2], lwd=3)
  lines(lpj.g.fcomp[[i]][,"BIBS"], col=pft.colors.lpj.g[3], lwd=3)
  lines(lpj.g.fcomp[[i]][,"TeBS"], col=pft.colors.lpj.g[4], lwd=3)
  lines(lpj.g.fcomp[[i]][,"TeIBS"], col=pft.colors.lpj.g[5], lwd=3)
  lines(lpj.g.fcomp[[i]][,"C3G"], col=pft.colors.lpj.g[6], lwd=3)
  text(x=100, y=.4, site.list[i], cex=1.5, font=2)
  #legend("topleft", legend=pfts.lpj.g, col=pft.colors.lpj.g, bty="n", lwd=3)
}
# plot(-5, type="l", ylim=c(0,1), col=pft.colors.lpj.g[1], lwd=0.0001, xlab="", ylab="Fraction of AGB")
legend("topleft", legend=pfts.lpj.g, col=pft.colors.lpj.g, bty="n", lwd=3, cex=1, horiz=T)
dev.off()



#-----------------------------------  
# LPJ-WSL
#-----------------------------------
summary(lpj.w.fcomp[[1]])
summary(lpj.w.fcomp[[2]])
summary(lpj.w.fcomp[[3]])
summary(lpj.w.fcomp[[4]])
summary(lpj.w.fcomp[[5]])
#summary(lpj.w.fcomp[[6]])

pft.colors.lpj.w <- c("green3", "darkgreen", "black","red3")
pfts.lpj.w <- c("BES", "BNE", "BBS", "C3")

par(mfrow=c(1,1))
plot(lpj.w.fcomp[[1]][,"BES"], type="l", ylim=c(0,1), col=pft.colors.lpj.w[1], lwd=3, main=site.list[1])
lines(lpj.w.fcomp[[1]][,"BNE"], col=pft.colors.lpj.w[2], lwd=3)
lines(lpj.w.fcomp[[1]][,"BBS"], col=pft.colors.lpj.w[3], lwd=3)
lines(lpj.w.fcomp[[1]][,"C3"], col=pft.colors.lpj.w[4], lwd=3)
legend("topleft", legend=pfts.lpj.w, col=pft.colors.lpj.w, bty="n", lwd=3)

pdf(width=11, height=8.5, file="PrelimGraphs/LPJ-GUESS_Paleon_Prelim_Fcomp.pdf")
par(mfrow=c(3,2), mar=c(3,5,1,1)+.1)
for(i in 1:length(site.list)){
  plot(lpj.w.fcomp[[i]][,"BES"], type="l", ylim=c(0,1), col=pft.colors.lpj.w[1], lwd=3, ylab="Fraction of AGB")
  lines(lpj.w.fcomp[[i]][,"BNE"], col=pft.colors.lpj.w[2], lwd=3)
  lines(lpj.w.fcomp[[i]][,"BBS"], col=pft.colors.lpj.w[3], lwd=3)
  lines(lpj.w.fcomp[[i]][,"C3"], col=pft.colors.lpj.w[4], lwd=3)
  text(x=50, y=.95, site.list[i], cex=1.5, font=2)
  #legend("topleft", legend=pfts.lpj.w, col=pft.colors.lpj.w, bty="n", lwd=3)
  if(i==1) legend(x=0, y=0.85, legend=pfts.lpj.w, col=pft.colors.lpj.w, bty="n", lwd=3, cex=1, horiz=T)

}
# plot(-5, type="l", ylim=c(0,1), col=pft.colors.lpj.w[1], lwd=0.0001, xlab="", ylab="Fraction of AGB")
dev.off()


#-----------------------------------  
# CLM
#-----------------------------------
summary(clm.fcomp[[1]]); summary(rowSums(clm.fcomp[[1]]))
summary(clm.fcomp[[2]]); summary(rowSums(clm.fcomp[[2]]))
summary(clm.fcomp[[3]]); summary(rowSums(clm.fcomp[[3]]))
summary(clm.fcomp[[4]]); summary(rowSums(clm.fcomp[[4]]))
summary(clm.fcomp[[5]]); summary(rowSums(clm.fcomp[[5]]))
summary(clm.fcomp[[6]]); summary(rowSums(clm.fcomp[[6]]))

pft.colors.clm <- c("darkgreen", "green3")
pfts.clm <- c("TeNE", "TeBS")

par(mfrow=c(1,1))
plot(clm.fcomp[[1]][,"TeNE"], type="l", ylim=c(0,1), col=pft.colors.clm[1], lwd=3, main=site.list[1])
lines(clm.fcomp[[1]][,"TeBS"], col=pft.colors.clm[2], lwd=3)
legend("topleft", legend=pfts.clm, col=pft.colors.clm, bty="n", lwd=3)

pdf(width=11, height=8.5, file="PrelimGraphs/CLM45_Paleon_Prelim_Fcomp.pdf")
par(mfrow=c(3,2), mar=c(3,5,1,1)+.1)
for(i in 1:length(site.list)){
  plot(clm.fcomp[[i]][,"TeNE"], type="l", ylim=c(0,1), col=pft.colors.clm[1], lwd=3, ylab="Fraction of AGB")
  lines(clm.fcomp[[i]][,"TeBS"], col=pft.colors.clm[2], lwd=3)
  text(x=700, y=.9, site.list[i], cex=1.5, font=2)
  #legend("topleft", legend=pfts.clm, col=pft.colors.clm, bty="n", lwd=3)
}
#plot(-5, type="l", ylim=c(0,1), col=pft.colors.clm[1], lwd=0.0001, xlab="", ylab="Fraction of AGB")
legend(x=9000, y=.90, legend=pfts.clm, col=pft.colors.clm, bty="n", lwd=5, cex=2)
dev.off()
