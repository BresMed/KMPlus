# script to generate KM+ estimate based on example data

source("UDFs.R")

dbloader(c(
  "flexsurv",
  "tidyverse",
  "gridExtra"
))

# parameters ----------------------------------------------------------------
dat <- lung
t_var <- "time"
e_var <- "status"
S_var <- "sex"

t <- dat[,t_var]
e <- dat[,e_var]
S <- dat[,S_var]
S_uniq <- unique(dat[, S_var])
S_len <- length(S_uniq)

time_steps <- seq(0,2000,1)
Tcutoff <- 500

dists     <- c("exponential","weibull","llogis",      "lnorm",     "gompertz","gengamma")

# ~ Full fit --------------------------------------------------------------

# from the beginning
sfit <- survfit(Surv(t,e) ~ S, data = dat)

sforms <- lapply(1:S_len, function(STRATA) {
  SurvregFormulaGen(
    t = dat[(which(dat[,S_var] == unique(S)[STRATA])),t_var],
    e = dat[(which(dat[,S_var] == unique(S)[STRATA])),e_var],
    covs = NULL,
    factors = NULL,
    nam_t = "time",
    nam_e = "status",
    nam_covs = NULL,
    nam_factors = NULL,
    DEBUG = FALSE
  )
})


MODELS_orig <- lapply(1:S_len, function(STRATA) {
  Temp <- lapply(1:6, function(MODEL) {
    mod <- flexsurvreg(sforms[[1]], data = dat, dist = dists[MODEL])
  })
  Temp <- setNames(Temp,dists)
})
MODELS_orig <- setNames(MODELS,S_uniq)

extraps_orig <- lapply(1:S_len, function(STRATA) {
  data.frame(
    time         = time_steps,
    exponential  = summary(MODELS_orig[[STRATA]][[dists[1]]],t=time_steps,ci=FALSE)[[1]]$est,
    weibull      = summary(MODELS_orig[[STRATA]][[dists[2]]],t=time_steps,ci=FALSE)[[1]]$est,
    log_logistic = summary(MODELS_orig[[STRATA]][[dists[3]]],t=time_steps,ci=FALSE)[[1]]$est,
    log_normal   = summary(MODELS_orig[[STRATA]][[dists[4]]],t=time_steps,ci=FALSE)[[1]]$est,
    gompertz     = summary(MODELS_orig[[STRATA]][[dists[5]]],t=time_steps,ci=FALSE)[[1]]$est,
    gen_gamma    = summary(MODELS_orig[[STRATA]][[dists[6]]],t=time_steps,ci=FALSE)[[1]]$est
  )
})
extraps_orig <- setNames(extraps,S_uniq)


plots_orig <- lapply(1:S_len, function(STRATA) {
  p <- ggplot(extraps_orig[[STRATA]],aes(x = time, y = exponential)) + 
    geom_line() +
    geom_line(data = extraps_orig[[STRATA]], colour = 2,linetype = 2,aes(x = time, y = weibull)) +
    geom_line(data = extraps_orig[[STRATA]], colour = 3,linetype = 3,aes(x = time, y = log_logistic)) +
    geom_line(data = extraps_orig[[STRATA]], colour = 4,linetype = 4,aes(x = time, y = log_normal)) +
    geom_line(data = extraps_orig[[STRATA]], colour = 5,linetype = 5,aes(x = time, y = gompertz)) +
    geom_line(data = extraps_orig[[STRATA]], colour = 6,linetype = 6,aes(x = time, y = gen_gamma)) + 
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
    labs(x = "time", y = "Survival (%)",title = paste0("Strata ",S_uniq[STRATA]," full dataset"))
  
  return(p)
  
})

# how to programatically plot n plots
# p_orig <- do.call("grid.arrange", c(plots_orig,ncol=S_len))

# ~ KM+ ---------------------------------------------------------------

# ~~ extrapolation based on re-baselined data -----------------------------

# re-baseline the remaining data at that cut-off
dat_kmplus <- dat[(which(dat[,t_var] >= Tcutoff)),]
dat_kmplus[,t_var] <- dat_kmplus[,t_var] - Tcutoff

# perform survival analysis on that re-baselined data

# from the beginning
sfit <- survfit(Surv(dat_kmplus[,t_var],dat_kmplus[,e_var]) ~ dat_kmplus[,S_var], data = dat_kmplus)

sforms <- lapply(1:S_len, function(STRATA) {
  SurvregFormulaGen(
    t = dat_kmplus[(which(dat_kmplus[,S_var] == unique(S)[STRATA])),t_var],
    e = dat_kmplus[(which(dat_kmplus[,S_var] == unique(S)[STRATA])),e_var],
    covs = NULL,
    factors = NULL,
    nam_t = "time",
    nam_e = "status",
    nam_covs = NULL,
    nam_factors = NULL,
    DEBUG = FALSE
  )
})

MODELS_KMplus <- lapply(1:S_len, function(STRATA) {
  Temp <- lapply(1:6, function(MODEL) {
    mod <- flexsurvreg(sforms[[1]], data = dat_kmplus, dist = dists[MODEL])
  })
  Temp <- setNames(Temp,dists)
})
MODELS_KMplus <- setNames(MODELS,S_uniq)

extraps_KMplus <- lapply(1:S_len, function(STRATA) {
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
extraps_KMplus <- setNames(extraps,S_uniq)


plots_rebaselined <- lapply(1:S_len, function(STRATA) {
  p <- ggplot(extraps_KMplus[[STRATA]],aes(x = time, y = exponential)) + 
    geom_line() +
    geom_line(data = extraps_KMplus[[STRATA]], colour = 2,linetype = 2,aes(x = time, y = weibull)) +
    geom_line(data = extraps_KMplus[[STRATA]], colour = 3,linetype = 3,aes(x = time, y = log_logistic)) +
    geom_line(data = extraps_KMplus[[STRATA]], colour = 4,linetype = 4,aes(x = time, y = log_normal)) +
    geom_line(data = extraps_KMplus[[STRATA]], colour = 5,linetype = 5,aes(x = time, y = gompertz)) +
    geom_line(data = extraps_KMplus[[STRATA]], colour = 6,linetype = 6,aes(x = time, y = gen_gamma)) + 
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
    labs(x = "time", y = "Survival (%)",title = paste0("Strata ",S_uniq[STRATA]," re-baselined at t=",Tcutoff))
  
  return(p)
  
})

# how to programatically plot n plots
# p_rebaselined <- do.call("grid.arrange", c(plots_rebaselined,ncol=S_len))


# ~~ Kaplan-Meier ---------------------------------------------------------

# generate KM by strata

KM <- lapply(1:S_len, function(STRATA) {
  
  DF <- dat %>% dplyr::filter(sex == S_uniq[STRATA])
  
  KM <- MakeKMPlotData(
    DATA     = DF,
    timevar  = "time",
    eventvar = "status"
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

# the last row in the KM_cuttoff objects are the beginning of the curves
#   simply multiply suvival by that number and add the time on

extraps_KMplusMod <- lapply(1:S_len, function(STRATA){
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


plots_KMPlus <- lapply(1:S_len, function(STRATA) {
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
    labs(x = "time", y = "Survival (%)",title = paste0("Strata ",S_uniq[STRATA]," KM with extrapolations starting at t=",Tcutoff))
  
  return(p)
  
})


plotlist <- c(plots_orig, plots_rebaselined,plots_KMPlus)
do.call("grid.arrange", c(plotlist,ncol=S_len,nrow=3))

