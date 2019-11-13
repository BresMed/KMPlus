source("UDFs.R")

# Functions ---------------------------------------------------------------

# as a function
GenKMPlusModels <- function(dat,t_var,e_var,S_var = NULL,covs= NULL,facs= NULL,namCovs= NULL,namFacs= NULL,Tcutoff) {
  
  t <- dat[,t_var]
  e <- dat[,e_var]
  # stratification parameters if needed
  if(!is.null(S_var)) {
    S <- dat[,S_var]
    S_uniq <- unique(S)
    S_len <- length(S_uniq)
  } else {
    S <- NULL
    S_uniq <- "Patients"
    S_len <- NULL
  }
  # iteration index for stratified or unstratified
  if (is.null(S)) {
    StratLoop <- 1
  } else {
    StratLoop <- 1:S_len
  }
  # survival distributions
  dists     <- c("exponential","weibull","llogis",      "lnorm",     "gompertz","gengamma")
  
  # generate rebaselined data
  dat_kmplus <- dat[(which(dat[,t_var] >= Tcutoff)),]
  dat_kmplus[,t_var] <- dat_kmplus[,t_var] - Tcutoff
  
  # fit as survival data depending on stratification variable
  if(is.null(S)) {
    sfit <- survfit(Surv(dat_kmplus[,t_var],dat_kmplus[,e_var]) ~ 1, data = dat_kmplus)
    sforms <- list(
      SurvregFormulaGen(
        t = dat_kmplus[, t_var],
        e = dat_kmplus[, e_var],
        covs = covs,
        factors = facs,
        nam_t = t_var,
        nam_e = e_var,
        nam_covs = namCovs,
        nam_factors = namFacs,
        DEBUG = FALSE
      )
    )
  } else {
    sfit <- survfit(Surv(dat_kmplus[,t_var],dat_kmplus[,e_var]) ~ dat_kmplus[,S_var], data = dat_kmplus)
    sforms <- lapply(1:S_len, function(STRATA) {
      SurvregFormulaGen(
        t = dat_kmplus[(which(dat_kmplus[,S_var] == unique(S)[STRATA])),t_var],
        e = dat_kmplus[(which(dat_kmplus[,S_var] == unique(S)[STRATA])),e_var],
        covs = covs,
        factors = facs,
        nam_t = t_var,
        nam_e = e_var,
        nam_covs = namCovs,
        nam_factors = namFacs,
        DEBUG = FALSE
      )
    })
  }
  
  # run the regressions
  MODELS_KMplus <- lapply(StratLoop, function(STRATA) {
    Temp <- lapply(1:6, function(MODEL) {
      mod <- flexsurvreg(sforms[[STRATA]], data = dat_kmplus, dist = dists[MODEL])
    })
    Temp <- setNames(Temp,dists)
  })
  MODELS_KMplus <- setNames(MODELS_KMplus,S_uniq)
  return(MODELS_KMplus)
}


GenKMPlusExtraps <- function(dat,MODELS_KMplus,time_steps,t_var,e_var,S_var=NULL) {
  
  t <- dat[,t_var]
  e <- dat[,e_var]
  # stratification parameters if needed
  if(!is.null(S_var)) {
    S <- dat[,S_var]
    S_uniq <- unique(S)
    S_len <- length(S_uniq)
  } else {
    S <- NULL
    S_uniq <- "Patients"
    S_len <- NULL
  }
  # iteration index for stratified or unstratified
  if (is.null(S)) {
    StratLoop <- 1
  } else {
    StratLoop <- 1:S_len
  }
  # survival distributions
  dists     <- c("exponential","weibull","llogis",      "lnorm",     "gompertz","gengamma")
  
  
  KM <- lapply(StratLoop, function(STRATA) {
    
    # stratify dataset if necessary 
    if (is.null(S)) {
      DF <- dat
    } else {
      DF <- dat[which(dat[,S_var] == S_uniq[STRATA]),]
    }
    
    # make the KM plot data both full and cut off
    KM <- MakeKMPlotData(
      DATA     = DF,
      timevar  = t_var,
      eventvar = e_var
    )
    KM_cutoff <- KM %>% filter(t <= Tcutoff)
    
    OUT <- list(
      KM_full = KM,
      KM_cutoff = KM_cutoff
    )
    
    # return a list of KM plot data and KM plot data that's been cut off
    return(OUT)
    
  })
  KM <- setNames(KM,S_uniq)
  
  extraps_KMplus <- lapply(StratLoop, function(STRATA) {
    data.frame(
      time         = time_steps,
      exponential  = summary(MODELS_KMplus[[STRATA]][[dists[1]]],t=time_steps,ci=FALSE)[[1]]$est,
      weibull      = summary(MODELS_KMplus[[STRATA]][[dists[2]]],t=time_steps,ci=FALSE)[[1]]$est,
      log_logistic = summary(MODELS_KMplus[[STRATA]][[dists[3]]],t=time_steps,ci=FALSE)[[1]]$est,
      log_normal   = summary(MODELS_KMplus[[STRATA]][[dists[4]]],t=time_steps,ci=FALSE)[[1]]$est,
      gompertz     = summary(MODELS_KMplus[[STRATA]][[dists[5]]],t=time_steps,ci=FALSE)[[1]]$est,
      gen_gamma    = summary(MODELS_KMplus[[STRATA]][[dists[6]]],t=time_steps,ci=FALSE)[[1]]$est
    )
  })
  extraps_KMplus <- setNames(extraps_KMplus,S_uniq)
  
  
  extraps_KMplusMod <- lapply(StratLoop, function(STRATA){
    KM_cutoff <- KM[[STRATA]]$KM_cutoff
    rbtime <- max(KM_cutoff$t)
    rbsurv <- min(KM_cutoff$s_t)
    
    # make a data.frame with KM+extrapolation for each extrapolation
    data.frame(
      time = c(KM_cutoff$t, extraps_KMplus[[STRATA]]$time + rbtime),
      exponential = c(KM_cutoff$s_t, extraps_KMplus[[STRATA]]$exponential * rbsurv),
      weibull = c(KM_cutoff$s_t, extraps_KMplus[[STRATA]]$weibull * rbsurv),
      log_logistic = c(KM_cutoff$s_t, extraps_KMplus[[STRATA]]$log_logistic * rbsurv),
      log_normal = c(KM_cutoff$s_t, extraps_KMplus[[STRATA]]$log_normal * rbsurv),
      gompertz = c(KM_cutoff$s_t, extraps_KMplus[[STRATA]]$gompertz * rbsurv),
      gen_gamma = c(KM_cutoff$s_t, extraps_KMplus[[STRATA]]$gen_gamma * rbsurv)
    )
  })
  
  return(extraps_KMplusMod) 
}


GenKMPlusPlots <- function(extraps_KMplusMod,S_var=NULL,Tcutoff,time_steps) {
  
  if(!is.null(S_var)) {
    S <- dat[,S_var]
    S_uniq <- unique(S)
    S_len <- length(S_uniq)
  } else {
    S <- NULL
    S_uniq <- "Patients"
    S_len <- NULL
  }
  
  if (is.null(S)) {
    StratLoop <- 1
  } else {
    StratLoop <- 1:S_len
  }
  TitleString <- lapply(StratLoop,function(STRATA){
    paste0("Strata: ",S_uniq[[STRATA]]," KM with extrapolations starting at t=",Tcutoff)
  })
  
  plots_KMPlus <- lapply(StratLoop, function(STRATA) {
    p <- ggplot(extraps_KMplusMod[[STRATA]],aes(x = time, y = exponential)) + 
      geom_line() +
      geom_line(data = extraps_KMplusMod[[STRATA]], colour = 2,linetype = 2,aes(x = time, y = weibull)) +
      geom_line(data = extraps_KMplusMod[[STRATA]], colour = 3,linetype = 3,aes(x = time, y = log_logistic)) +
      geom_line(data = extraps_KMplusMod[[STRATA]], colour = 4,linetype = 4,aes(x = time, y = log_normal)) +
      geom_line(data = extraps_KMplusMod[[STRATA]], colour = 5,linetype = 5,aes(x = time, y = gompertz)) +
      geom_line(data = extraps_KMplusMod[[STRATA]], colour = 6,linetype = 6,aes(x = time, y = gen_gamma)) + 
      theme_bw() +
      scale_x_continuous(limits = c(0, max(time_steps)), expand = expand_scale(mult = c(0, .05))) +
      scale_y_continuous(limits = c(0, 1), expand = expand_scale(mult = c(0, .05)),labels = scales::percent) +
      theme(
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.text = element_text(size = 10),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        strip.background = element_rect(colour = "black", fill = "#d4d4d4")
      ) +
      guides(colour = guide_legend(nrow = 2)) +
      scale_linetype_manual(
        c(
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash",
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash",
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash",
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash",
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash",
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash",
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash",
          "solid", "dashed", "dotted", "dotdash", "longdash","twodash"
        )
      ) + 
      scale_colour_manual(
        values = 
          c(
            "black","#DB3B93","#45B9D1","#002D5C","#F9161C","#B97E22",
            "#75AC3F","blue","orange","green","brown","yellow",
            "red","brown","grey50","purple","lightblue",
            "black","#DB3B93","#45B9D1","#002D5C","#F9161C","#B97E22",
            "#75AC3F","blue","orange","green","brown","yellow",
            "red","brown","grey50","purple","lightblue",
            "black","#DB3B93","#45B9D1","#002D5C","#F9161C","#B97E22",
            "#75AC3F","blue","orange","green","brown","yellow",
            "red","brown","grey50","purple","lightblue"
          )
      ) +
      scale_size_manual(values = c(1,0.6,0.6,0.6,0.6,0.6,0.6,
                                   1,0.6,0.6,0.6,0.6,0.6,0.6,
                                   1,0.6,0.6,0.6,0.6,0.6,0.6,
                                   1,0.6,0.6,0.6)) +
      labs(x = "time", y = "Survival (%)",title = TitleString[[STRATA]])
    
    return(p)
  })
  
  
  return(do.call("grid.arrange", c(plots_KMPlus,ncol=S_len,nrow=1)))
}


# Testing -----------------------------------------------------------------



MODELS_strat    <- GenKMPlusModels( dat = dat,t_var = t_var,e_var = e_var,S_var = S_var, Tcutoff = Tcutoff)
MODELS_unstrat  <- GenKMPlusModels( dat = dat,t_var = t_var,e_var = e_var,               Tcutoff = Tcutoff)
EXTRAPS_strat   <- GenKMPlusExtraps(dat = dat,t_var = t_var,e_var = e_var,S_var = S_var,time_steps = time_steps,MODELS_KMplus = MODELS_strat)
EXTRAPS_unstrat <- GenKMPlusExtraps(dat = dat,t_var = t_var,e_var = e_var,              time_steps = time_steps,MODELS_KMplus = MODELS_unstrat)

GenKMPlusPlots(extraps_KMplusMod = EXTRAPS_strat  ,S_var = "sex",Tcutoff = Tcutoff, time_steps = time_steps)
GenKMPlusPlots(extraps_KMplusMod = EXTRAPS_unstrat,              Tcutoff = Tcutoff, time_steps = time_steps)



# 
# dat <- lung
# t_var <- "time"
# e_var <- "status"
# S_var <- "sex"
# covs <- NULL
# facs <- NULL
# namCovs <- NULL
# namFacs <- NULL
# Tcutoff <- 500
# time_steps <- 0:1000
# 
# 
# # dat <- ovarian
# # t_var <- "futime"
# # e_var <- "fustat"
# # S_var <- "ecog.ps"
# 
# t <- dat[,t_var]
# e <- dat[,e_var]
# S <- dat[,S_var]
# S_uniq <- unique(dat[, S_var])
# S_len <- length(S_uniq)
# 
# time_steps <- seq(0,2000,1)
# Tcutoff <- 500
# 
