# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Yearly Output
# This script starts where 1_ left off, aggregates things to annual & packages them together
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

sec2yr <- 1*60*60*24*365

# ---------------------------------------------- #
# Variables that exist for all models, but need to be aggregated to year for comparision with LPJ

AGB.y <- TotLivBiom <- TotSoilCarb <- LAI.y <- GPP.y <- NEE.y <- NPP.y <- HeteroResp.y <- AutoResp.y <- LAI.y  <- Qs.y <- SoilMoist.y <- SoilTemp.y <- Evap.y <- Transp.y <- Fcomp.y <- Dens.y <- BA.y <- list()

tair.y <- precipf.y <- wind.y <- lwdown.y <- swdown.y <- qair.y <- psurf.y <- list()

var1 <- c("NEE", "NPP", "HeteroResp", "AutoResp", "GPP", "LAI.m", "Evap", "Qs", "SoilMoist", "SoilTemp", "SWE", "Transp", "Fcomp", "Dens", "BA")
for(i in 1:length(var1)){
  assign(var1[i], list())
}

lpj.g.soilC <- data.frame(array(NA, dim=dim(lpj.w[["TotSoilCarb"]]))); names(lpj.g.soilC) <- site.list


for(s in 1:length(site.list)){
  yr.rows <- seq(1, length(ed[["AGB"]][,s]), by=12)

	#---------------------------------------------------------------------
	# DIVERSITY / STRUCTURE
	#---------------------------------------------------------------------
    #-----------------------------------
	# Fcomp (check)
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["Fcomp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["Fcomp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["Fcomp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["Fcomp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  jules.temp <- c(jules.temp, NA)
	  lpj.w.temp <- c(lpj.w.temp, NA)
	  lpj.g.temp <- c(lpj.g.temp, NA)
	  linkages.temp <- c(linkages.temp, NA)
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["Fcomp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  sib.temp <- c(sib.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["Fcomp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["Fcomp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    Fcomp.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["Fcomp"]][,s], lpj.g[["Fcomp"]][,s]), jules.temp, jules.triff.temp, linkages.temp, sib.temp)  
    names(Fcomp.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


    #-----------------------------------
	# Density (total)
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["Dens"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["Dens"]][yr.rows[i]:(yr.rows[i]+11),s]))

	  clm.bgc.temp <- c(clm.bgc.temp, NA)
	  clm.cn.temp <- c(clm.cn.temp, NA)
	  jules.temp <- c(jules.temp, NA)
	  jules.triff.temp <- c(jules.triff.temp, NA)
	  lpj.w.temp <- c(lpj.w.temp, NA)
	  lpj.g.temp <- c(lpj.g.temp, NA)
	  linkages.temp <- c(linkages.temp, NA)
	  sib.temp <- c(sib.temp, NA)
      # clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["Dens"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # jules.temp <- c(jules.temp, mean(jules.s[["Dens"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["Dens"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["Dens"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    Dens.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["Dens"]][,s], lpj.g[["Dens"]][,s]), jules.temp, jules.triff.temp, linkages.temp, sib.temp)  
    names(Dens.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")



    #-----------------------------------
	# Basal Area (total)
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["BA"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["BA"]][yr.rows[i]:(yr.rows[i]+11),s]))

	  clm.bgc.temp <- c(clm.bgc.temp, NA)
	  clm.cn.temp <- c(clm.cn.temp, NA)
	  jules.temp <- c(jules.temp, NA)
	  jules.triff.temp <- c(jules.triff.temp, NA)
	  lpj.w.temp <- c(lpj.w.temp, NA)
	  lpj.g.temp <- c(lpj.g.temp, NA)
	  linkages.temp <- c(linkages.temp, NA)
	  sib.temp <- c(sib.temp, NA)
      # clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["BA"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # jules.temp <- c(jules.temp, mean(jules.s[["BA"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["BA"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["BA"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    BA.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["BA"]][,s], lpj.g[["BA"]][,s]), jules.temp, jules.triff.temp, linkages.temp, sib.temp)  
    names(BA.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#---------------------------------------------------------------------
	# CARBON FLUXES
	#---------------------------------------------------------------------

    #-----------------------------------
	# GPP
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)
	  sib.temp <- c(sib.temp, mean(sib[["GPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    GPP.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(GPP.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


    #-----------------------------------
	# NPP
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  lpj.w.temp <- c(lpj.w.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)
	  sib.temp <- c(sib.temp, mean(sib[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    NPP.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["NPP"]][,s], lpj.g.temp, jules.temp, jules.triff.temp, linkages$NPP[,site.list[s]], sib.temp))  
    names(NPP.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


    #-----------------------------------
	# NEE
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["NEE"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["NEE"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["NEE"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["NEE"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s])-mean(jules.s[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]*0.237))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["NPP"]][yr.rows[i]:(yr.rows[i]+11),s])-mean(jules.triff[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]*0.237))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["NEE"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["NEE"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)      
	  sib.temp <- c(sib.temp, mean(sib[["NEE"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    NEE.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages$NEE[,site.list[s]], sib.temp))  
    names(NEE.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


    #-----------------------------------
	# AutoResp
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)
	  sib.temp <- c(sib.temp, mean(sib[["AutoResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    AutoResp.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(AutoResp.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


    #-----------------------------------
	# HeteroResp
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]*0.237))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]*0.237))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)      
	  sib.temp <- c(sib.temp, mean(sib[["HeteroResp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  }
    HeteroResp.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages$HeteroResp[,site.list[s]], sib.temp))  
    names(HeteroResp.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


	#---------------------------------------------------------------------
	# CARBON POOLS
	#---------------------------------------------------------------------
	
	#-----------------------------------
	# AGB
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))

	  lpj.w.temp <- c(lpj.w.temp, NA)
	  lpj.g.temp <- c(lpj.g.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)      
	  sib.temp <- c(sib.temp, mean(sib[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))

      }
    AGB.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["TotLivBiom"]][,s], lpj.g[["AGB"]][,s], jules.temp, jules.triff.temp, linkages$AGB[,site.list[s]], sib.temp))  
    names(AGB.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


	#-----------------------------------
	# TotLivBiom
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))

	  lpj.w.temp <- c(lpj.w.temp, NA)
	  lpj.g.temp <- c(lpj.g.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)      
	  sib.temp <- c(sib.temp, mean(sib[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))

      }
    TotLivBiom[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["TotLivBiom"]][,s], lpj.g[["TotLivBiom"]][,s], jules.temp, jules.triff.temp, linkages$TotLivBiomass[,site.list[s]], sib.temp))  
    names(TotLivBiom[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


	#-----------------------------------
	# TotSoilCarb
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["TotSoilCarb"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["TotSoilCarb"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["TotSoilCarb"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["TotSoilCarb"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["TotSoilCarb"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["TotSoilCarb"]][yr.rows[i]:(yr.rows[i]+11),s]))

	  lpj.w.temp <- c(lpj.w.temp, NA)
	  lpj.g.temp <- c(lpj.g.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["TotLivBiom"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["AGB"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)      
	  sib.temp <- c(sib.temp, mean(sib[["TotSoilCarb"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  }

	# Extracting Soil from LPJ.G
	lpj.g.soilC <- vector()
    dir.lpj.g <- file.path(model.dir, "LPJ-GUESS.v6", paste(site.list[1], "LPJ-GUESS", sep="_"))
    # dir.lpj.g <- file.path(model.dir, "LPJ-GUESS.v2", paste(site.list[1], "LPJ-GUESS", sep="_"))
	files.lpj.g <- dir(dir.lpj.g)
  
	index <- gregexpr("month",files.lpj.g[2])[[1]][1] # LPJ-GUESS has separate annual and monthly files & we just want the monthly
	files.lpj.g.y <- files.lpj.g[substr(files.lpj.g, index, index+5)=="annual"]
  
	for(i in 1:length(files.lpj.g.y)){
   		ncMT <- nc_open(file.path(dir.lpj.g, files.lpj.g.y[i]))
        lpj.g.soilC <- c(lpj.g.soilC, colSums(ncvar_get(ncMT, "CarbPools")[3:4,]))
      	}    
    nc_close(ncMT)

    TotSoilCarb[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["TotSoilCarb"]][,s], lpj.g.soilC, jules.temp, jules.triff.temp, linkages$TotSoilCarb[,site.list[s]], sib.temp))  
    names(TotSoilCarb[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


	#---------------------------------------------------------------------
	# OTHER
	#---------------------------------------------------------------------
    
	#-----------------------------------
	# LAI
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  lpj.w.temp <- c(lpj.w.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)      
	  sib.temp <- c(sib.temp, mean(sib[["LAI"]][yr.rows[i]:(yr.rows[i]+11),s]))

      }
    LAI.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w[["LAI"]][,s], lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(LAI.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

        
	#-----------------------------------
	# Evap
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)      
	  sib.temp <- c(sib.temp, mean(sib[["Evap"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    Evap.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages$Evap[,site.list[s]], sib.temp))  
    names(Evap.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")
        
	#-----------------------------------
	# Transp
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["Transp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["Transp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["Tranp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["Tranp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  jules.temp <- c(jules.temp, "NA")
	  jules.triff.temp <- c(jules.triff.temp, "NA")
      # jules.temp <- c(jules.temp, mean(jules.s[["Transp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["Tranp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["Transp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["Tranp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  }
    Transp.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(Transp.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#-----------------------------------
	# SoilMoist
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["SoilMoist"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    SoilMoist.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(SoilMoist.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


	#-----------------------------------
	# SoilTemp
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["SoilTemp"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    SoilTemp.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(SoilTemp.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


	#-----------------------------------
	# Qs
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["Qs"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    Qs.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(Qs.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


	#---------------------------------------------------------------------
	# MET CHECKING!
	#---------------------------------------------------------------------

	#-----------------------------------
	# tair
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  lpj.w.temp <- c(lpj.w.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["tair"]][yr.rows[i]:(yr.rows[i]+11),s]))

      }
    tair.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(tair.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#-----------------------------------
	# precipf
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  lpj.w.temp <- c(lpj.w.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["precipf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    precipf.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(precipf.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#-----------------------------------
	# wind
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, NA)
      lpj.g.temp <- c(lpj.g.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["wind"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    wind.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(wind.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#-----------------------------------
	# lwdown
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, NA)
      lpj.g.temp <- c(lpj.g.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["lwdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    lwdown.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(lwdown.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#-----------------------------------
	# swdown
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["swdown"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    swdown.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp,  linkages.temp, sib.temp))  
    names(swdown.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#-----------------------------------
	# qair
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, NA)
      lpj.g.temp <- c(lpj.g.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["qair"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    qair.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp,  linkages.temp, sib.temp))  
    names(qair.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")

	#-----------------------------------
	# psurf
    #-----------------------------------
    for(i in 1:length(yr.rows)){
     if(i==1) ed.temp <- ed.lu.temp <- clm.bgc.temp <- clm.cn.temp <- jules.temp <- jules.triff.temp <- lpj.w.temp <- lpj.g.temp <- linkages.temp <- sib.temp <- vector()
      ed.temp <- c(ed.temp, mean(ed[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      ed.lu.temp <- c(ed.lu.temp, mean(ed.lu[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.bgc.temp <- c(clm.bgc.temp, mean(clm.bgc[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      clm.cn.temp <- c(clm.cn.temp, mean(clm.cn[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.temp <- c(jules.temp, mean(jules.s[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      jules.triff.temp <- c(jules.triff.temp, mean(jules.triff[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      lpj.w.temp <- c(lpj.w.temp, NA)
      lpj.g.temp <- c(lpj.g.temp, NA)
      # lpj.w.temp <- c(lpj.w.temp, mean(lpj.w[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      # lpj.g.temp <- c(lpj.g.temp, mean(lpj.g[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
	  linkages.temp <- c(linkages.temp, NA)            
	  sib.temp <- c(sib.temp, mean(sib[["psurf"]][yr.rows[i]:(yr.rows[i]+11),s]))
      }
    psurf.y[[s]] <- data.frame(cbind(ed.temp, ed.lu.temp, clm.bgc.temp, clm.cn.temp, lpj.w.temp, lpj.g.temp, jules.temp, jules.triff.temp, linkages.temp, sib.temp))  
    names(psurf.y[[s]]) <- c("ed2", "ed2.lu", "clm.bgc", "clm.cn", "lpj.wsl", "lpj.guess", "jules.stat", "jules.triffid", "linkages", "sibcasa")


}


names(AGB.y) <- names(TotLivBiom) <- names(TotSoilCarb) <- names(LAI.y) <- names(GPP.y) <- names(NEE.y) <- names(NPP.y) <- names(HeteroResp.y) <- names(AutoResp.y) <- names(LAI.y)  <- names(Qs.y) <- names(SoilMoist.y) <- names(SoilTemp.y) <- names(Evap.y) <- names(Transp.y) <- names(Fcomp.y) <- names(Dens.y) <- names(BA.y) <- site.list




# ----------------------------------------------------------------------------
# Storing key variables as .nc files
dim.years <- ncdim_def(name="Year", units="Years", vals=1850:(1850+nrow(GPP.y[[1]])-1))
dim.allmods <- ncdim_def(name="Models", units="", vals=1:ncol(GPP.y[[1]]))
# dim.3mods <- ncdim_def(name="Models", units="", vals=1:3)
# dim.4mods <- ncdim_def(name="Models", units="", vals=1:4)

dim.string <- ncdim_def("names", "", 1:24, create_dimvar=FALSE)
dim.allmods2 <- ncdim_def(name="Model Names", units="", vals=1:ncol(GPP.y[[1]]), create_dimvar=FALSE)
# dim.3mods2 <- ncdim_def(name="Model Names", units="", vals=1:3, create_dimvar=FALSE)
# dim.4mods2 <- ncdim_def(name="Model Names", units="", vals=1:4, create_dimvar=FALSE)

units.fluxes     <- "kg m-2 s-1"
units.pools      <- "kgC m-2"
units.lai        <- "m2 m-2"
units.temp       <- "K"
units.smoist <- "kg m-2"
units.radiation  <- "W m-2"
units.psurf      <- "Pa"
units.qair       <- "kg kg-1"
units.wind       <- "m s-1"
models.all <- names(GPP.y[[1]])
# models.4 <- names(tair.y[[1]])

GPP.vars <- AGB.vars <- LAI.vars <-  NEE.vars <- NPP.vars <- AutoResp.vars <- HeteroResp.vars <- TotSoilCarb.vars <- Evap.vars <- Transp.vars <- SoilMoist.vars <- list()

tair.vars <- precipf.vars <- swdown.vars <- lwdown.vars <- wind.vars <- psurf.vars <- qair.vars <- list()

for(i in 1:length(site.list)){
	GPP.vars[[i]] <- ncvar_def(site.list[i], units=units.fluxes, dim=list(dim.allmods, dim.years))
	AGB.vars[[i]] <- ncvar_def(site.list[i], units= units.pools, dim=list(dim.allmods, dim.years))
	LAI.vars[[i]] <- ncvar_def(site.list[i], units=units.lai,    dim=list(dim.allmods, dim.years))
	NEE.vars[[i]] <- ncvar_def(site.list[i], units=units.fluxes, dim=list(dim.allmods, dim.years))
	NPP.vars[[i]] <- ncvar_def(site.list[i], units=units.fluxes, dim=list(dim.allmods, dim.years))

	AutoResp.vars[[i]]    <- ncvar_def(site.list[i], units=units.fluxes, dim=list(dim.allmods, dim.years))
	HeteroResp.vars[[i]]  <- ncvar_def(site.list[i], units=units.fluxes, dim=list(dim.allmods, dim.years))
	TotSoilCarb.vars[[i]] <- ncvar_def(site.list[i], units=units.pools,  dim=list(dim.allmods, dim.years))
	Evap.vars[[i]]        <- ncvar_def(site.list[i], units=units.fluxes, dim=list(dim.allmods, dim.years))
	Transp.vars[[i]]      <- ncvar_def(site.list[i], units=units.fluxes, dim=list(dim.allmods, dim.years))
	SoilMoist.vars[[i]]   <- ncvar_def(site.list[i], units=units.smoist, dim=list(dim.allmods, dim.years))

	# met drivers
	tair.vars[[i]]    <- ncvar_def(site.list[i], units=units.temp,      dim=list(dim.allmods, dim.years))
	precipf.vars[[i]] <- ncvar_def(site.list[i], units=units.fluxes,    dim=list(dim.allmods, dim.years))
	swdown.vars[[i]]  <- ncvar_def(site.list[i], units=units.radiation, dim=list(dim.allmods, dim.years))
	lwdown.vars[[i]]  <- ncvar_def(site.list[i], units=units.radiation, dim=list(dim.allmods, dim.years))
	wind.vars[[i]]    <- ncvar_def(site.list[i], units=units.wind,      dim=list(dim.allmods, dim.years))
	psurf.vars[[i]]   <- ncvar_def(site.list[i], units=units.psurf,     dim=list(dim.allmods, dim.years))
	qair.vars[[i]]    <- ncvar_def(site.list[i], units=units.qair,      dim=list(dim.allmods, dim.years))

}	
GPP.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
AGB.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
LAI.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
NPP.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
NEE.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")

AutoResp.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
HeteroResp.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
TotSoilCarb.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
Evap.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
Transp.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
SoilMoist.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")

tair.vars[[length(site.list)+1]]    <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
precipf.vars[[length(site.list)+1]] <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
swdown.vars[[length(site.list)+1]]  <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
lwdown.vars[[length(site.list)+1]]  <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
wind.vars[[length(site.list)+1]]    <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
psurf.vars[[length(site.list)+1]]   <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")
qair.vars[[length(site.list)+1]]    <- ncvar_def("ModelNames", units="", dim=list(dim.string, dim.allmods2), prec="char")


names(GPP.vars) <- names(AGB.vars) <- names(AutoResp.vars) <- names(HeteroResp.vars) <- names(TotSoilCarb.vars) <- names(Evap.vars) <- names(Transp.vars) <- names(SoilMoist.vars) <-  names(tair.vars) <- names(precipf.vars) <- names(swdown.vars) <- names(lwdown.vars) <- names(wind.vars) <- names(psurf.vars) <- names(qair.vars) <- c(site.list, "ModelNames")
summary(GPP.vars)
summary(NPP.vars)

output.location <- "phase1a_output_variables"
gpp <- nc_create(file.path(output.location, "GPP.annual.nc"), GPP.vars)
agb <- nc_create(file.path(output.location, "AGB.annual.nc"), AGB.vars)
lai <- nc_create(file.path(output.location, "LAI.annual.nc"), LAI.vars)
npp <- nc_create(file.path(output.location, "NPP.annual.nc"), NPP.vars)
nee <- nc_create(file.path(output.location, "NEE.annual.nc"), NEE.vars)
auto.resp   <- nc_create(file.path(output.location, "RespirationAuto.annual.nc"  ), AutoResp.vars)
hetero.resp <- nc_create(file.path(output.location, "RespirationHetero.annual.nc"), HeteroResp.vars)
soilcarb    <- nc_create(file.path(output.location, "TotSoilCarb.annual.nc"      ), TotSoilCarb.vars)
evap        <- nc_create(file.path(output.location, "Evap.annual.nc"             ), Evap.vars)
soilmoist   <- nc_create(file.path(output.location, "SoilMoist.annual.nc"        ), SoilMoist.vars)
# transp <- nc_create(file.path(output.location, "Transp.annual.nc"), Transp.vars)

# All Drivers
tair    <- nc_create(file.path(output.location, "tair.annual.nc"   ), tair.vars   )
precipf <- nc_create(file.path(output.location, "precipf.annual.nc"), precipf.vars)
swdown  <- nc_create(file.path(output.location, "swdown.annual.nc" ), swdown.vars )
lwdown  <- nc_create(file.path(output.location, "lwdown.annual.nc" ), lwdown.vars )
wind    <- nc_create(file.path(output.location, "wind.annual.nc"   ), wind.vars   )
psurf   <- nc_create(file.path(output.location, "psurf.annual.nc"  ), psurf.vars  )
qair    <- nc_create(file.path(output.location, "qair.annual.nc"   ), qair.vars   )


# for(i in 1:5){
for(i in 1:length(site.list)){
	ncvar_put(npp, NPP.vars[[i]], t(NPP.y[[i]]))
	ncvar_put(gpp, GPP.vars[[i]], t(GPP.y[[i]]))
	ncvar_put(agb, AGB.vars[[i]], t(AGB.y[[i]]))
	ncvar_put(lai, LAI.vars[[i]], t(LAI.y[[i]]))
	ncvar_put(nee, NEE.vars[[i]], t(NEE.y[[i]]))

	ncvar_put(auto.resp  , AutoResp.vars[[i]]   , t(AutoResp.y[[i]])   )
	ncvar_put(hetero.resp, HeteroResp.vars[[i]] , t(HeteroResp.y[[i]]) )
	ncvar_put(soilcarb   , TotSoilCarb.vars[[i]], t(TotSoilCarb[[i]])  )
	ncvar_put(evap       , Evap.vars[[i]]       , t(Evap.y[[i]])       )
	ncvar_put(soilmoist  , SoilMoist.vars[[i]]  , t(SoilMoist.y[[i]])  )
	# ncvar_put(transp   , Transp.vars[[i]], t(Transp.y[[i]]))

	# Drivers
	ncvar_put(tair   , tair.vars[[i]]   , t(tair.y[[i]])    )
	ncvar_put(precipf, precipf.vars[[i]], t(precipf.y[[i]]) )
	ncvar_put(swdown , swdown.vars[[i]] , t(swdown.y[[i]])  )
	ncvar_put(lwdown , lwdown.vars[[i]] , t(lwdown.y[[i]])  )
	ncvar_put(wind   , wind.vars[[i]]   , t(wind.y[[i]])    )
	ncvar_put(psurf  , psurf.vars[[i]]  , t(psurf.y[[i]])   )
	ncvar_put(qair   , qair.vars[[i]]   , t(qair.y[[i]])    )

}

# NOTE: Some very weird bug is causing PMB to not save right in the loop, 
#   so we need to repeat it manually here
ncvar_put(npp, NPP.vars[[6]], t(NPP.y[[6]]))

ncvar_put(gpp, GPP.vars[[length(site.list)+1]], models.all)
ncvar_put(agb, AGB.vars[[length(site.list)+1]], models.all)
ncvar_put(lai, LAI.vars[[length(site.list)+1]], models.all)
ncvar_put(npp, NPP.vars[[length(site.list)+1]], models.all)
ncvar_put(nee, NEE.vars[[length(site.list)+1]], models.all)

ncvar_put(auto.resp  , AutoResp.vars[[length(site.list)+1]]   , models.all)
ncvar_put(hetero.resp, HeteroResp.vars[[length(site.list)+1]] , models.all)
ncvar_put(soilcarb   , TotSoilCarb.vars[[length(site.list)+1]], models.all)
ncvar_put(evap       , Evap.vars[[length(site.list)+1]]       , models.all)
ncvar_put(soilmoist  , SoilMoist.vars[[length(site.list)+1]]  , models.all)
# ncvar_put(transp     , Transp.vars[[length(site.list)+1]]       , models.all)

ncvar_put(tair   , tair.vars[[length(site.list)+1]]   , models.all)
ncvar_put(precipf, precipf.vars[[length(site.list)+1]], models.all)
ncvar_put(swdown , swdown.vars[[length(site.list)+1]] , models.all)
ncvar_put(lwdown , lwdown.vars[[length(site.list)+1]] , models.all)
ncvar_put(wind   , wind.vars[[length(site.list)+1]]   , models.all)
ncvar_put(psurf  , psurf.vars[[length(site.list)+1]]  , models.all)
ncvar_put(qair   , qair.vars[[length(site.list)+1]]   , models.all)


nc_close(gpp); nc_close(agb); nc_close(auto.resp); nc_close(hetero.resp); nc_close(soilcarb); nc_close(evap); nc_close(soilmoist); #nc_close(transp)
nc_close(tair); nc_close(precipf); nc_close(swdown); nc_close(lwdown); nc_close(wind); nc_close(psurf); nc_close(qair);

# nc <- nc_open(file.path(file.path(output.location, "GPP.annual.nc")))
# summary(nc$var)
# ncvar_get(nc, "ModelNames")
# ----------------------------------------------------------------------------


# --------------------------------------------------------

# GPP
summary(GPP.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/GPP_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(GPP.y[[s]][,"ed2"]*sec2yr, ylim=range(GPP.y, na.rm=T)*sec2yr, col="black", type="l", lwd=2, ylab="GPP KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "GPP", sep=": "))
	lines(GPP.y[[s]][,"ed2.lu"]*sec2yr, col="gray50", lwd=1.5)
	lines(GPP.y[[s]][,"jules.stat"]*sec2yr, col="orange4", lwd=2.5)
	lines(GPP.y[[s]][,"jules.triffid"]*sec2yr, col="goldenrod3", lwd=1.5)
	lines(GPP.y[[s]][,"lpj.guess"]*sec2yr, col="blue", lwd=2)
	lines(GPP.y[[s]][,"lpj.wsl"]*sec2yr, col="red", lwd=2)
	lines(GPP.y[[s]][,"clm.bgc"]*sec2yr, col="green3", lwd=2)
	lines(GPP.y[[s]][,"clm.cn"]*sec2yr, col="darkolivegreen", lwd=2)
	lines(GPP.y[[s]][,"sibcasa"]*sec2yr, col="salmon", lwd=2)
	if(site.list[s]=="PHA") {
		# From Urbanski et al 2007
		points(x=2000-850, y=14*1e3*1e-4, cex=1.5, col="gray50", pch=19)
		# arrows(x0=2000-850, y0=(14-1.641137817)*1e-6*yr2sec*1e4, x1=2000-850, y1=(14+1.641137817)*1e-6*yr2sec*1e4, length=0, lwd=2, col="gray50") 
		}
	if(s==4) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", ncol=3, cex=0.75)
}
dev.off()


# NPP
summary(NPP.y[[1]]) 
pdf(width=8.5, height=11, file="PrelimGraphs/NPP_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(NPP.y[[s]][,"ed2"], ylim=range(NPP.y, na.rm=T), col="black", type="l", lwd=2, ylab="NPP KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "NPP", sep=": "))
	lines(NPP.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(NPP.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(NPP.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(NPP.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(NPP.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(NPP.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(NPP.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(NPP.y[[s]][,"linkages"], col="purple2", lwd=3)
	lines(NPP.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==6) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()


pdf(width=8.5, height=11, file="PrelimGraphs/NPP_Annual_AllSites_Truncated.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(NPP.y[[s]][,"ed2"], ylim=c(0, 8e-8), col="black", type="l", lwd=2, ylab="NPP KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "NPP", sep=": "))
	lines(NPP.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(NPP.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(NPP.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(NPP.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(NPP.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(NPP.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(NPP.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(NPP.y[[s]][,"linkages"], col="purple2", lwd=3)
	lines(NPP.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==5) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.75)
}
dev.off()

# NEE
summary(NEE.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/NEE_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(-NEE.y[[s]][,"ed2"], type="l", ylim=range(NEE.y, na.rm=T), lwd=2, ylab="NEE KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "NEE", sep=": "))
	lines(NEE.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(NEE.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(NEE.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(NEE.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(NEE.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(NEE.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(NEE.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(NEE.y[[s]][,"linkages"], col="purple2", lwd=2)
	lines(NEE.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(site.list[s]=="PHA") {
		# From Urbanski et al 2007
		points(x=2000-850, y=2.453846154*1e3*yr2sec*1e-4, cex=2, col="gray70", pch=19)
		# arrows(x0=2000-850, y0=(2.064461538)*1e-6*yr2sec*1e4, x1=2000-850, y1=(2.725384615)*1e-6*yr2sec*1e4, length=0, lwd=2, col="gray50") 
		}
	if(s==1) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SibCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.75)
}
dev.off()

pdf(width=8.5, height=11, file="PrelimGraphs/NEE_Annual_AllSites_Truncated.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(-NEE.y[[s]][,"ed2"], type="l", ylim=c(-4e-8, 5e-8), lwd=2, ylab="NEE KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "NEE", sep=": "))
	lines(NEE.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(NEE.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(NEE.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(NEE.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(NEE.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(NEE.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(NEE.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(NEE.y[[s]][,"linkages"], col="purple2", lwd=2)
	lines(NEE.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(site.list[s]=="PHA") {
		# From Urbanski et al 2007
		points(x=2000-850, y=2.453846154*1e3*yr2sec*1e-4, cex=2, col="gray70", pch=19)
		# arrows(x0=2000-850, y0=(2.064461538)*1e-6*yr2sec*1e4, x1=2000-850, y1=(2.725384615)*1e-6*yr2sec*1e4, length=0, lwd=2, col="gray50") 
		}
	if(s==1) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SibCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()

# AutoResp
summary(AutoResp.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/AutoResp_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(AutoResp.y[[s]][,"ed2"], ylim=range(AutoResp.y, na.rm=T), type="l", lwd=2, ylab="AutoResp KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "AutoResp", sep=": "))
	lines(AutoResp.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(AutoResp.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(AutoResp.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(AutoResp.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	# lines(AutoResp.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(AutoResp.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(AutoResp.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(AutoResp.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(AutoResp.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==5) legend("topleft", legend=c("ED2","ED-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "SibCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", cex=0.5, ncol=2)
}
dev.off()

# HeteroResp
summary(HeteroResp.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/HeteroResp_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(HeteroResp.y[[s]][,"ed2"], ylim=range(HeteroResp.y, na.rm=T), type="l", lwd=2, ylab="HeteroResp KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "HeteroResp", sep=": "))
	lines(HeteroResp.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(HeteroResp.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(HeteroResp.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(HeteroResp.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(HeteroResp.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(HeteroResp.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(HeteroResp.y[[s]][,"jules.triffid"], col="goldenrod", lwd=2)
	lines(HeteroResp.y[[s]][,"linkages"], col="purple2", lwd=3)
	lines(HeteroResp.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==6) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3)
}
dev.off()

pdf(width=8.5, height=11, file="PrelimGraphs/HeteroResp_Annual_AllSites_Truncated.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(HeteroResp.y[[s]][,"ed2"], ylim=c(0, 8e-8), type="l", lwd=2, ylab="HeteroResp KgC/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "HeteroResp", sep=": "))
	lines(HeteroResp.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(HeteroResp.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(HeteroResp.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(HeteroResp.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(HeteroResp.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(HeteroResp.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(HeteroResp.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(HeteroResp.y[[s]][,"linkages"], col="purple2", lwd=3)
	lines(HeteroResp.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==5) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="o", bg="white", ncol=3, cex=0.5)
}
dev.off()

# AGB
summary(AGB.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/AGB_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(AGB.y[[s]][,"ed2"], ylim=range(AGB.y, na.rm=T), type="l", lwd=2, ylab="AGB KgC/m2", xlab="years since 850-01-01", main=paste(site.list[s], "AGB", sep=": "))
	lines(AGB.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(AGB.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(AGB.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(AGB.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(AGB.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(AGB.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(AGB.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(AGB.y[[s]][,"linkages"], col="purple2", lwd=3)
	lines(AGB.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(site.list[s]=="PHA") {
		# From Urbanski et al 2007
		points(x=2000-850, y=107.7942857*1e3*1e-4, cex=2, col="gray70", pch=19)
		}
	if(s==5) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()

# TotLivBiom
summary(TotLivBiom[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/TotLivBiom_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(TotLivBiom[[s]][,"ed2"], ylim=range(TotLivBiom, na.rm=T), type="l", lwd=2, ylab="TotLivBiom KgC/m2", xlab="years since 850-01-01", main=paste(site.list[s], "TotLivBiom", sep=": "))
	lines(TotLivBiom[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(TotLivBiom[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(TotLivBiom[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(TotLivBiom[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(TotLivBiom[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(TotLivBiom[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(TotLivBiom[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(TotLivBiom[[s]][,"linkages"], col="purple2", lwd=2)
	lines(TotLivBiom[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==6) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()

# TotSoilCarb
summary(TotSoilCarb[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/TotSoilCarb_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(TotSoilCarb[[s]][,"ed2"], ylim=range(TotSoilCarb, na.rm=T), type="l", lwd=2, ylab="TotSoilCarb KgC/m2", xlab="years since 850-01-01", main=paste(site.list[s], "TotSoilCarb", sep=": "))
	lines(TotSoilCarb[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(TotSoilCarb[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(TotSoilCarb[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(TotSoilCarb[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(TotSoilCarb[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(TotSoilCarb[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(TotSoilCarb[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(TotSoilCarb[[s]][,"linkages"], col="purple2", lwd=2)
	lines(TotSoilCarb[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==2) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=2, cex=0.75)
}
dev.off()

pdf(width=8.5, height=11, file="PrelimGraphs/TotSoilCarb_Annual_AllSites_Truncated.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(TotSoilCarb[[s]][,"ed2"], ylim=c(0,100), type="l", lwd=2, ylab="TotSoilCarb KgC/m2", xlab="years since 850-01-01", main=paste(site.list[s], "TotSoilCarb", sep=": "))
	lines(TotSoilCarb[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(TotSoilCarb[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(TotSoilCarb[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(TotSoilCarb[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(TotSoilCarb[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(TotSoilCarb[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(TotSoilCarb[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(TotSoilCarb[[s]][,"linkages"], col="purple2", lwd=2)
	lines(TotSoilCarb[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==3) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()

# LAI
summary(LAI.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/LAI_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(LAI.y[[s]][,"ed2"], ylim=range(LAI.y, na.rm=T), col="black", type="l", lwd=2, ylab="LAI", xlab="years since 850-01-01", main=paste(site.list[s], "LAI", sep=": "))
	lines(LAI.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(LAI.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(LAI.y[[s]][,"clm.bgc"], col="green3", lwd=0.5)
	lines(LAI.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(LAI.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(LAI.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(LAI.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(LAI.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(LAI.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(site.list[s]=="PHA") {
		# From Urbanski et al 2007
		points(x=2000-850, y=5.225, cex=2, col="gray50", pch=19)
		}
	if(s==4) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", bg="white", ncol=3, cex=0.5)
}
dev.off()

# Evap
summary(Evap.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/Evap_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
#	plot(Evap.y[[s]][,"ed2"], ylim=range(Evap.y, na.rm=T), type="l", lwd=2, ylab="Evap kg/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "Evaporation", sep=": "))
	plot(Evap.y[[s]][,"ed2"], ylim=c(min(Evap.y[[s]], na.rm=T), max(Evap.y[[s]][,1:(ncol(Evap.y$PHA)-1)], na.rm=T)), type="l", lwd=2, ylab="Evap kg/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "Evaporation", sep=": "))
	lines(Evap.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(Evap.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(Evap.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(Evap.y[[s]][,"clm.bgc"], col="green3", lwd=0.5)
	lines(Evap.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(Evap.y[[s]][,"linkages"], col="purple2", lwd=2)
	lines(Evap.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(Evap.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(Evap.y[[s]][,"sibcasa"]/sec2yr, col="salmon", lwd=2)
	if(s==5) legend("topleft", legend=c("ED2", "ED-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "LINKAGES", "SibCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "purple2", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()


# Transp
# summary(Transp.y[[2]])
# pdf(width=8.5, height=11, file="PrelimGraphs/Transp_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(Transp.y[[s]][,"ed2"], ylim=range(Transp.y, na.rm=T), type="l", lwd=2, ylab="Transp kg/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "Transpiration", sep=": "))
	# lines(Transp.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	# lines(Transp.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	# lines(Transp.y[[s]][,"clm.bgc"], col="green3", lwd=0.5)
	# # lines(Transp.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	# if(s==5) legend("topleft", legend=c("ED2", "LPJ-GUESS", "LPJ-WSL", "CLM45"), col=c("black", "blue", "red", "green3"), lwd=2, bty="n")
# }
# dev.off()

# SoilMoist
summary(SoilMoist.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/SoilMoist_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(SoilMoist.y[[s]][,"ed2"], ylim=range(SoilMoist.y, na.rm=T), type="l", lwd=2, ylab="SoilMoist", xlab="years since 850-01-01", main=paste(site.list[s], "Soil Moisture", sep=": "))
	lines(SoilMoist.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(SoilMoist.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(SoilMoist.y[[s]][,"clm.bgc"], col="green3", lwd=0.5)
	lines(SoilMoist.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(SoilMoist.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(SoilMoist.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(SoilMoist.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(SoilMoist.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==1) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()

# SoilTemp
summary(SoilTemp.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/SoilTemp_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(SoilTemp.y[[s]][,"ed2"], ylim=c(270, 285), type="l", lwd=2, ylab="SoilTemp K", xlab="years since 850-01-01", main=paste(site.list[s], "Soil Temperature", sep=": "))
	lines(SoilTemp.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(SoilTemp.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	# lines(SoilTemp.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(SoilTemp.y[[s]][,"clm.bgc"], col="green3", lwd=0.5)
	lines(SoilTemp.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(SoilTemp.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(SoilTemp.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(SoilTemp.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==1) legend("bottomleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "gray50", "blue", "green3", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()

# Qs
summary(Qs.y[[1]])
pdf(width=8.5, height=11, file="PrelimGraphs/Runoff_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(Qs.y[[s]][,"ed2"], ylim=range(Qs.y, na.rm=T), type="l", lwd=2, ylab="Runoff kg/m2/s", xlab="years since 850-01-01", main=paste(site.list[s], "Runoff", sep=": "))
	lines(Qs.y[[s]][,"ed2.lu"], col="gray50", lwd=1.5)
	lines(Qs.y[[s]][,"lpj.guess"], col="blue", lwd=2)
	lines(Qs.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(Qs.y[[s]][,"clm.bgc"], col="green3", lwd=0.5)
	lines(Qs.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(Qs.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(Qs.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(Qs.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	if(s==5) legend("topleft", legend=c("ED2", "ED2-LU", "LPJ-GUESS", "LPJ-WSL", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "gray50", "blue", "red", "green3", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", ncol=3, cex=0.5)
}
dev.off()


# Fcomp
summary(Fcomp.y[[1]])

# Dens
# BA


# --------------------------------------------------------------------------------
# Met Checking
# --------------------------------------------------------------------------------
# Tair

#pdf(width=8.5, height=11, file="PrelimGraphs/Tair_Annual_AllSites.pdf")
summary(tair.y[[1]])
pdf("PrelimGraphs/MetDrivers_Tair_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(tair.y[[s]][,"ed2"], ylim=c(270, max(tair.y[[s]], na.rm=T)), type="l", lwd=2, ylab="Air Temp (k)", xlab="years since 850-01-01", main=paste(site.list[s], "Air Temp", sep=": "))
	lines(tair.y[[s]][,"lpj.guess"], col="blue", lwd=2.5)
#	lines(tair.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(tair.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(tair.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(tair.y[[s]][,"jules.stat"], col="orange4", lwd=1.5)
	lines(tair.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=1.5)
	# lines(tair.y[[s]][,"lpj.guess"], col="blue", lwd=1)
	# lines(tair.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=1)
	lines(tair.y[[s]][,"ed2"], col="black", lwd=1)
	lines(tair.y[[s]][,"sibcasa"], col="salmon", lwd=1)
	lines(tair.y[[s]][,"ed2.lu"], col="gray50", lwd=0.8)
	if(s==6) legend("bottomleft", legend=c("ED2", "ED2-LU","CLM-BGC", "CLM-CN", "LPJ-GUESS", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black","gray50", "green3", "blue", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", ncol=3, cex=0.75)
}
dev.off()

summary(precipf.y[[1]])
pdf("PrelimGraphs/MetDrivers_precipf_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(precipf.y[[s]][,"ed2"], ylim=range(precipf.y, na.rm=T), type="l", lwd=2, ylab="Precip Rate", xlab="years since 850-01-01", main=paste(site.list[s], "Precip Rate", sep=": "))
	lines(precipf.y[[s]][,"ed2.lu"], col="gray50", lwd=2)
	lines(precipf.y[[s]][,"lpj.guess"], col="blue", lwd=2.5)
#	lines(precipf.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(precipf.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(precipf.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(precipf.y[[s]][,"jules.stat"], col="orange4", lwd=1)
	lines(precipf.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=1)
	lines(precipf.y[[s]][,"sibcasa"], col="salmon", lwd=1)
	lines(precipf.y[[s]][,"ed2.lu"], col="gray50", lwd=0.8)
	if(s==5) legend("topleft", legend=c("ED2", "ED2-LU", "CLM-BGC", "CLM-CN", "LPJ-GUESS", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "gray50", "green3", "darkolivegreen", "blue", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", ncol=2)
}
dev.off()

summary(swdown.y[[1]])
pdf("PrelimGraphs/MetDrivers_swdown_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(swdown.y[[s]][,"ed2"], ylim=range(swdown.y[[s]], na.rm=T), type="l", lwd=2, ylab="swdown (W/m2/s)", xlab="years since 850-01-01", main=paste(site.list[s], "Downwelling Shortwave Rad", sep=": "))
	lines(swdown.y[[s]][,"ed2.lu"], col="gray50", lwd=2)
	lines(swdown.y[[s]][,"lpj.guess"], col="blue", lwd=2)
#	lines(tair.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	lines(swdown.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	lines(swdown.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(swdown.y[[s]][,"jules.stat"], col="orange4", lwd=2)
	lines(swdown.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2)
	lines(swdown.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	# if(s==1) legend("bottomleft", legend=c("ED2", "CLM-BGC", "CLM-CN", "LPJ-GUESS", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "green3", "darkolivegreen", "blue", "orange4", "goldenrod3", "salmon"), lwd=2, bty="c", bg="white", ncol=3, cex=0.7)
}
dev.off()

pdf("PrelimGraphs/MetDrivers_lwdown_Annual_AllSites.pdf")
par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
for(s in 1:length(site.list)){
	plot(lwdown.y[[s]][,"ed2"], ylim=range(lwdown.y[[s]], na.rm=T), type="l", lwd=3, ylab="swdown (W/m2/s)", xlab="years since 850-01-01", main=paste(site.list[s], "Downwelling Longwave Rad", sep=": "))
	# lines(lwdown.y[[s]][,"lpj.guess"], col="blue", lwd=2.5)
#	lines(lwdown.y[[s]][,"lpj.wsl"], col="red", lwd=2)
	# lines(lwdown.y[[s]][,"clm.bgc"], col="green3", lwd=2)
	# lines(lwdown.y[[s]][,"clm.cn"], col="darkolivegreen", lwd=2)
	lines(lwdown.y[[s]][,"jules.stat"], col="orange4", lwd=2.75)
	lines(lwdown.y[[s]][,"jules.triffid"], col="goldenrod3", lwd=2.5)
	lines(lwdown.y[[s]][,"sibcasa"], col="salmon", lwd=2)
	lines(lwdown.y[[s]][,"ed2.lu"], col="gray50", lwd=1)
	if(s==1) legend("bottomleft", legend=c("ED2", "ED2-LU", "CLM-BGC", "CLM-CN", "JULES-STATIC", "JULES-TRIFFID", "SiBCASA"), col=c("black", "gray50", "green3", "darkolivegreen", "orange4", "goldenrod3", "salmon"), lwd=2, bty="n", ncol=2, cex=0.5)
}
dev.off()

#------------------------------------------
# # tair
# pdf("PrelimGraphs/MetDrivers_Tair_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(tair.y[[s]][,"ed2"], ylim=c(270, 283), type="l", lwd=1, ylab="Air Temp", xlab="years since 850-01-01", main=paste(site.list[s]))
# }
# dev.off()

# # precipf
# pdf("PrelimGraphs/MetDrivers_precipf_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(precipf.y[[s]][,"ed2"], ylim=range(precipf.y[[s]][["ed2"]], na.rm=T), type="l", lwd=1, ylab="Precip", xlab="years since 850-01-01", main=paste(site.list[s]))
# }
# dev.off()

# # psurf
# pdf("PrelimGraphs/MetDrivers_psurf_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(psurf.y[[s]][,"ed2"], ylim=range(psurf.y[[s]][["ed2"]], na.rm=T), type="l", lwd=1, ylab="Psurf", xlab="years since 850-01-01", main=paste(site.list[s]))
# }
# dev.off()

# # qair
# pdf("PrelimGraphs/MetDrivers_qair_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(qair.y[[s]][,"ed2"], ylim=range(qair.y[[s]][["ed2"]], na.rm=T), type="l", lwd=1, ylab="Qair", xlab="years since 850-01-01", main=paste(site.list[s]))
# }
# dev.off()

# # wind
# pdf("PrelimGraphs/MetDrivers_wind_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(wind.y[[s]][,"ed2"], ylim=range(wind.y[[s]][["ed2"]], na.rm=T), type="l", lwd=1, ylab="Wind", xlab="years since 850-01-01", main=paste(site.list[s]))
# }
# dev.off()

# # co2
# co2.ann.dir <- "~/Dropbox/PalEON CR/phase1a_env_drivers/phase1a_env_drivers_v3/paleon_co2/paleon_annual_co2.nc"
# co2.ann.nc <- nc_open(co2.ann.dir)
# co2.ann <- ncvar_get(co2.ann.nc, "co2")
# nc_close(co2.ann.nc)
# summary(co2.ann)
# pdf("PrelimGraphs/MetDrivers_co2_Annual_AllSites.pdf")
 # par(new=F, mfrow=c(1,1))
 # plot(co2.ann, type="l", lwd=2, main="CO2", xlab="years since 850-01-01")
# dev.off()
# # pdf("PrelimGraphs/MetDrivers_co2_Annual_AllSites.pdf")
# # par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# # for(s in 1:length(site.list)){
	# # plot(co2.y[[s]][,"ed2"], ylim=range(co2.y[[s]][["ed2"]], na.rm=T), type="l", lwd=1, ylab="CO2", xlab="years since 850-01-01", main=paste(site.list[s]))
# # }
# # dev.off()

# # lwdown
# pdf("PrelimGraphs/MetDrivers_lwdown_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(lwdown.y[[s]][,"ed2"], ylim=range(lwdown.y[[s]][["ed2"]], na.rm=T), type="l", lwd=1, ylab="lwdown", xlab="years since 850-01-01", main=paste(site.list[s]))
# }
# dev.off()

# # swdown
# pdf("PrelimGraphs/MetDrivers_swdown_Annual_AllSites.pdf")
# par(mfrow=c(round((length(site.list)+.5)/2, 0), 2), mar=c(4,5,4,1)+0.1)
# for(s in 1:length(site.list)){
	# plot(swdown.y[[s]][,"ed2"], ylim=range(swdown.y[[s]][["ed2"]], na.rm=T), type="l", lwd=1, ylab="swdown", xlab="years since 850-01-01", main=paste(site.list[s]))
# }
# dev.off()
