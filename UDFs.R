# UDFs
message("Server:  Loading user-defined functions (UDFs)")
# Package loading machine -------------------------------------------------
dbloader <- function(packs) {
  system.time(
    if (length(setdiff(packs, rownames(installed.packages()))) > 0) {
      message("you need the following packages to run this script")
      message(setdiff(packs, rownames(installed.packages())))
      message("Installing them now...")
      install.packages(setdiff(packs, rownames(installed.packages())))
      message("Now loading libraries...")
      sapply(packs, require, character.only = TRUE)
    } else {
      message("All packages are installed already")
      message("Loading the specified libraries...")
      sapply(packs, require, character.only = TRUE)
    }
  )
}


# Convenience functions ---------------------------------------------------


# ~ Time unit conversions -------------------------------------------------
# nice convenience function to help converting numbers consistently
UnitsInOneYear <- function(Time_unit){
  switch (Time_unit,
          Days = 365.25,
          Weeks = 52,
          `4-week cycles` = 13,
          Months = 12,
          Years = 1
  )
}

DefaultStepsize <- function(Time_unit) {
  switch (Time_unit,
          Days = 7,
          Weeks = 1,
          `4-week cycles` = 1,
          Months = 1,
          Years = 1/52
  )
}

NatRiskDefaultStepSize <- function(Time_unit) {
  switch (Time_unit,
          Days = 180,
          Weeks = 24,
          `4-week cycles` = 26,
          Months = 6,
          Years = 1/2
  )
}


# ~ Round to nearest X ----------------------------------------------------
# can round to any arbitrary number, which is useful for numbers of cycles generated in this model
round_any <- function(x, accuracy, f=round){f(x/ accuracy) * accuracy}


# ~ Age adjustment --------------------------------------------------------

AgeAdjustExtrapolations <- function(extraps,life_table) {
  # N.B. at this stage, the extrapolations are always a list object with at least 1 element
  # Adjustment for general population survival
  #   1. lapply separating out each model individually and adding column for gpop mortality p(death)
  #   2. use ifelse() to pick the probability to use
  #   3. re-calculate survival in each cycle using foreach() %do% {}
  #   4. stick back together into n(strata) data.frames containing the original time(units) and strata columns
  
  extrapnames <- colnames(extraps[[1]])
  
  extraps <- lapply(1:length(extraps), function(STRATA){
    
    # for each strata, perform a lapply and stick back together columnwise
    DF <- data.frame(extraps[[STRATA]])
    
    # Make the middle columns (the ones that the adjustment applies to)
    # shave off the left and right most columns
    TEMP <- DF[,2:(ncol(DF)-1)]
    NAMS <- colnames(TEMP)
    
    # LOGIC: recalc max pr(death) for each row, for each column vs gpop
    TEMP2 <- do.call(
      'cbind',
      lapply(1:ncol(TEMP),function(COLUMN){
        
        # LOGIC: Take column, make probs, add in gpop vlookup, logical check, replace if true
        # Reduce can be used to quickly calc time-dependency: mutate(s = Reduce(function(a,b) a*p[b], 1:(n()-1), 1, acc=T))
        ThisCol <- data.frame(OrigSurv = TEMP[,COLUMN]) %>%
          mutate(
            DeathProb     = ifelse(row_number() == n(),0,1-(dplyr::lead(OrigSurv)/OrigSurv)),
            DeathProbGpop = ifelse(row_number() == n(),0,1-(dplyr::lead(life_table$SurvCuveWA)/life_table$SurvCuveWA))
          ) %>%
          mutate(
            FinalSurvProb = ifelse(DeathProb < DeathProbGpop,1-DeathProbGpop,1-DeathProb)
          ) %>%
          mutate(
            FinalSurv = Reduce(function(a,b) a*FinalSurvProb[b], 1:(n()-1), 1, acc=T)
          )
        
        #rename the column with its original name again
        ThisCol <- ThisCol %>% select(FinalSurv)
        colnames(ThisCol) <- NAMS[COLUMN]
        return(ThisCol)
      })
    )
    
    #make a sandwich with the results
    OUT <- cbind(DF$time,TEMP2,DF$Strata)
    return(OUT)
  })
  
  
  for (STRATA in 1:length(extraps)) {
    colnames(extraps[[STRATA]]) <- extrapnames
  }
  
  return(extraps)
  
}

# x <- readRDS("~/../Downloads/intRact_SaveState_Oct.15.2019_08.24.19_2.rds")
# 
# extraps <- lapply(1:length(unique(x$DEBUG$patient_flow$Strata)), function(STRATA){
#   x$DEBUG$patient_flow %>% dplyr::filter(Strata == unique(x$DEBUG$patient_flow$Strata)[STRATA])
# })
#   
# 
# AgeAdjustExtrapolations(extraps = extraps,life_table = life_table)


#  ~ string manipulation --------------------------------------------------

#  ~~ row-wise concatenation ----------------------------------------------

concat_rowwise <- function(DF,cols_as_vect) {
  do.call(paste, c(DF[cols_as_vect], sep = ""))
}
# x <- read.csv("~/../Downloads/URGHHHHHH.csv")
# concat_rowwise(x,1:3)
# do.call('paste', c(x[1:3], sep = ""))



# MakeKM - outputs the KM as a dataframe## ----------------------------------------------------
MakeKM <- function(Time, Event, DataSource) {
  x <-
    survfit(Surv(Time,Event) ~ 1,
            data = DataSource,
            type = "kaplan-meier")
  KMSum <- (summary(x))
  y <-
    data.frame(
      "Time" = c(
        0,
        min(DataSource[Event == 1, ]$Event),
        KMSum$time,
        max(Time)
      ),
      "Survival" = c(1, 1, KMSum$surv, min(KMSum$surv))
    )
  y <- filter(y, Time != "Inf", Time != "-Inf")
  y <- mutate(y, timeyrs = Time/365.25)
  y <- cbind(y, "data" = "KM")
  colnames(y) <- c("time", "est", "timeyrs","data")
  return(y)
}



# Make Km plot data -------------------------------------------------------
MakeKMPlotData <- function(DATA, timevar, eventvar) {
  #make KM  testing ground
  #load example data and do survfit
  x <- DATA
  y <- Surv(x[,timevar], x[,eventvar])
  z <- survfit(y ~ 1)
  
  #make adjustments to data required for graph (i.e. double both datasets)
  t   <- rep(z$time,2)
  s_t <- rep(z$surv,2)
  
  #sort parameters as required (by time ascending, by survival descending), adding extra data points
  t   <- append(0,t[order(t)])
  s_t <- append(1,append(1,s_t[order(s_t,decreasing = TRUE)]))[1:length(t)]
  
  #put into a dataframe and return as output of function
  df <- data.frame(t = t, s_t = s_t)
  return(df)
}


# Get AIC and BIC table from flexsurv output ------------------------------
#Make a list, extract names of models, extract AIC, extract BIC and compile into a nice DF
get_aic_bic <- function(models) {
  output <- list()
  output$names  <- c("Exponential", "Weibull",
                     "Log-logistic", "Log-normal",
                     "Gompertz", "Generalised gamma")
  output$AICs   <- lapply(models, FUN = function(x) AIC(x))
  output$BICs   <- lapply(models, FUN = function(x) BIC(x))
  output <- data.frame(Model = output$names,
                       AIC   = unlist(output$AICs),
                       BIC   = unlist(output$BICs))
  output
}


# Generate flexsurv formula -----------------------------------------------

#'function: SurvregFormulaGen
#'Purpose: generate a formula object 
#'Author: Darren Burns
#'Inst.: BresMed
#'Date: 28/08/2019
#'Ver.: 0.0.1
SurvregFormulaGen <- function(t,e,covs,factors,nam_t,nam_e,nam_covs,nam_factors,treat_factor = FALSE,nam_Treat,DEBUG=FALSE) {
  # make sure time and event data is numeric
  t <- t %>% as.numeric()
  e <- e %>% as.numeric()
  
  # simple formula if there are no covariates
  if (all(is.null(nam_covs), is.null(nam_factors))) {
    survreg_formula <- "Surv(t, e) ~ 1"
  } else {
    
    # more complicated formula if there are covariates/factors at play
    #start from the basic:
    survreg_formula <- "Surv(t, e) ~ "
    
    # add in the rest using a logical framework for the model specification
    
    if (all(is.null(nam_covs), !is.null(nam_factors))) {
      # no covariates but categorical (factor) variables included
      # Add in these categoricals one at a time
      for (i in 1:length(nam_factors)) {
        if (i == 1) {
          #The first one should take the formula and add in the first covariate 
          #    WITHOUT a + in front of it
          survreg_formula <- paste0(survreg_formula, "factor(",nam_factors[1],")")
        } else {
          #all the additional covariates will have a + separating it from the others
          survreg_formula <- paste0(survreg_formula, " + factor(", nam_factors[i],")")
        }
      }
    } else if(all(!is.null(nam_covs), is.null(nam_factors))) {
      #no stratifications but continuous variables included
      # Add in the covariates one at a time
      for (i in 1:length(nam_covs)) {
        if (i == 1) {
          #The first one should take the formula and add in the first covariate 
          #    WITHOUT a + in front of it
          survreg_formula <- paste0(survreg_formula, nam_covs[1])
        } else {
          #all the additional covariates will have a + separating it from the others
          survreg_formula <- paste0(survreg_formula, " + ", nam_covs[i])
        }
      }
    } else {
      #both stratifications and covariates
      for (i in 1:length(nam_covs)) {
        if (i == 1) {
          #The first one should take the formula and add in the first covariate 
          #    WITHOUT a + in front of it
          survreg_formula <- paste0(survreg_formula, nam_covs[1])
        } else {
          #all the additional covariates will have a + separating it from the others
          survreg_formula <- paste0(survreg_formula, " + ", nam_covs[i])
        }
        
      }
      # now add in the stratification factors one at a time
      for (i in 1:length(nam_factors)) {
        survreg_formula <- paste0(survreg_formula, " + factor(", nam_factors[i],")")
      }
    }
  }
  
  # Final step: If ARM included as a covariate, include it here as a factor.
  if (treat_factor) {
    #if treatment is covariate, then add it to specification as a catagorical variable
    survreg_formula <- paste0(survreg_formula, " + factor(",nam_Treat,")")
  }
  
  # Debug output: print formula to console
  if (DEBUG) {
    print(survreg_formula)
  }
  
  # Now define the final string as a formula object to be used in flexsurv estimation
  survreg_formula <- as.formula(survreg_formula)
  
  # This formula is then passed into the regression analysis
  return(survreg_formula)
  
}

# # Testing example:
# t <- ovarian[, "futime"] %>% as.numeric()
# e <- ovarian[, "fustat"] %>% as.numeric()
# covs <- ovarian[, "age"] %>% as.numeric()
# Factors <- ovarian[, "rx"] %>% as.factor()
# 
# nam_t <- "futime"
# nam_e <- "fustat"
# nam_covs <- c("age")
# nam_factors <- c("rx","resid.ds")
# SurvregFormulaGen(
#   t = t,
#   e = e,
#   covs = covs,
#   factors = factors,
#   nam_t = nam_t,
#   nam_e = nam_e,
#   nam_covs = nam_covs,
#   nam_factors = nam_factors,
#   DEBUG = TRUE
# )


# Estimate MIXTURE CURE model (MCM) ---------------------------------------------

# The function flexsurvcure() is used directly in the server



# Extrapolate flexsurv fits -----------------------------------------------
#MUCH faster implementation is to simply not estimate the CIs!
extrapolate_models <- function(models,time_steps,Category=NULL,Category_level=NULL){
  # If the context is a stratification of the data (or just the whole dataset)  
  if (is.null(Category)) {
    extrapolations <- data.frame(
      time         = time_steps,
      exponential  = summary(models$exponential,t=time_steps,ci=FALSE)[[1]]$est,
      weibull      = summary(models$weibull,t=time_steps,ci=FALSE)[[1]]$est,
      log_logistic = summary(models$log_logistic,t=time_steps,ci=FALSE)[[1]]$est,
      log_normal   = summary(models$log_normal,t=time_steps,ci=FALSE)[[1]]$est,
      gompertz     = summary(models$gompertz,t=time_steps,ci=FALSE)[[1]]$est,
      gen_gamma    = summary(models$gen_gamma,t=time_steps,ci=FALSE)[[1]]$est
    )
  } else {
    # If the context is a factor variable and we want to see curves by levels of
    #   this variable, then we have to generate a set of extrapolations with
    #   that factor variable set to different values (but based on the same
    #   regression analysis). here we use the Category character value to generate 
    #   a set of fits specifically for that value of the Category_level variable
    # an example with ovarian dataset:
    # x <- ovarian
    # fit <- flexsurvreg(Surv(futime, fustat) ~ resid.ds + sigma(resid.ds), data=x, dist="gengamma")
    # t <- seq(1, 1227, by=1)
    # summary(fit, t=t, newdata=data.frame(resid.ds=1),ci=FALSE)[[1]][,"est"]
    #
    # Applying that here:
    #
    #   Generate newdat and rename the column to the same name as the data used
    newdat <- data.frame(Category_level)
    names(newdat)[1] <- Category
    # perform all extrapolations for that level of the category
    extrapolations <- data.frame(
      time = time_steps,
      exponential = summary(object = models$exponential, t = time_steps,newdata = newdat,ci = FALSE)[[1]]$est,
      weibull = summary(object = models$weibull, t = time_steps, newdata = newdat, ci = FALSE)[[1]]$est,
      log_logistic = summary(object = models$log_logistic, t = time_steps, newdata = newdat, ci = FALSE)[[1]]$est,
      log_normal = summary(object = models$log_normal, t = time_steps, newdata = newdat, ci = FALSE)[[1]]$est,
      gompertz = summary(object = models$gompertz, t = time_steps, newdata = newdat, ci = FALSE)[[1]]$est,
      gen_gamma = summary(object = models$gen_gamma, t = time_steps, newdata = newdat, ci = FALSE)[[1]]$est
    )
    
  }
  
  return(extrapolations)  
}


# median survival --------------------------------------------

# ~ from PLD --------------------------------------------------------------
get_median_survival <- function(survfit_object) {
  read.table(textConnection(capture.output(survfit_object)),
             skip = 2,
             header = TRUE)$median
}

# ~ from fitted models --------------------------------------------------------------
getSurv_Med <- function(models,Category=NULL,Category_level=NULL){
  # If the context is a stratification of the data (or just the whole dataset)  
  if (is.null(Category)) {
    medians <- data.frame(
      exponential  = summary(models$exponential , type = "median", ci=FALSE),
      weibull      = summary(models$weibull     , type = "median", ci=FALSE),
      log_logistic = summary(models$log_logistic, type = "median", ci=FALSE),
      log_normal   = summary(models$log_normal  , type = "median", ci=FALSE),
      gompertz     = summary(models$gompertz    , type = "median", ci=FALSE),
      gen_gamma    = summary(models$gen_gamma   , type = "median", ci=FALSE)
    )
    
  } else {
    # Generate newdat and rename the column to the same name as the data used
    #   See above function (extrapolate_models) for more detailed notes
    newdat <- data.frame(Category_level)
    names(newdat)[1] <- Category
    # perform all extrapolations for that level of the category
    medians <- data.frame(
      exponential  = summary(object = models$exponential , type = "median", newdata = newdat, ci = FALSE),
      weibull      = summary(object = models$weibull     , type = "median", newdata = newdat, ci = FALSE),
      log_logistic = summary(object = models$log_logistic, type = "median", newdata = newdat, ci = FALSE),
      log_normal   = summary(object = models$log_normal  , type = "median", newdata = newdat, ci = FALSE),
      gompertz     = summary(object = models$gompertz    , type = "median", newdata = newdat, ci = FALSE),
      gen_gamma    = summary(object = models$gen_gamma   , type = "median", newdata = newdat, ci = FALSE)
    )
    
  }
  
  # Note the hard-coding of names here. QC point to check (that these are ordered correctly, & that it works with MCM).
  medians <- medians %>%
    rename(
      `Exponential` = est,
      `Weibull` = est.1,
      `Log-logistic` = est.2,
      `Log-normal` = est.3,
      `Gompertz` = est.4,
      `Generalised gamma` = est.5,
    ) %>% 
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Model") %>%
    rename(`Median Survival` = V1)
  
  
  return(medians)  
}



# functions for graphics --------------------------------------------------


# ~ plot functions --------------------------------------------------------

#  ~~ extrapolation of survival models ------------------------------------
# expects a data.frame() with 4 columns: line, time, value, Strata
UI_ExtrapPlot <- function(Plot_data,xlim,TimeLab,nrow=2) {
  ggplot(Plot_data, aes(x = time, y = value, colour = line, size = line)) +
    geom_line() +
    theme_bw() +
    scale_x_continuous(limits = c(0, xlim), expand = c(0, 0)) +
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
    guides(colour = guide_legend(nrow = nrow)) +
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
    labs(x = TimeLab, y = "Survival (%)")
}


# ~~ LCHP -----------------------------------------------------------------
# Expects the same data.frame as the above function for extrapolations!
UI_LCHP <- function(Plot_data,ylim_max_mult,ylim_min_mult,TimeLab,nrow=2) {
  
  # convert the time and survivla values
  Plot_data$time <- log(Plot_data$time)
  Plot_data$value <- log(-log(Plot_data$value))
  
  # draw graph
  # max and min are set by the data, but you can multiply them to see beyond the data!!
  ymax <-
    max(Plot_data[which(Plot_data$line == "Kaplan-Meier" & is.finite(Plot_data$value)), "value"]) +
    abs(max(Plot_data[which(Plot_data$line == "Kaplan-Meier" & is.finite(Plot_data$value)), "value"]) * (ylim_max_mult - 1))
  ymin <-
    min(Plot_data[which(Plot_data$line == "Kaplan-Meier" & is.finite(Plot_data$value)), "value"]) -
    abs(min(Plot_data[which(Plot_data$line == "Kaplan-Meier" & is.finite(Plot_data$value)), "value"]) * (ylim_min_mult - 1))
  
  # draw the plot unstratified, and use facet grid if the data is stratified
  #    or a treatment covariate is being used to spearate the data out
  p <- ggplot(Plot_data, aes(x = time, y = value, colour = line, size = line)) +
    geom_line() +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.text = element_text(size = 10),
      axis.text=element_text(size=10),
      axis.title=element_text(size=12),
      strip.background = element_rect(colour = "black", fill = "#d4d4d4")
    ) +
    guides(colour = guide_legend(nrow = nrow)) +
    scale_y_continuous(limits = c(ymin,ymax)) +
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
    labs(x = paste0("ln(",TimeLab,")"), y = paste0("ln(-ln(S(",TimeLab,")))"))
}


# ~~ Comparative efficacy plot(s) --------------------------------------------
# ~~~ Absolute --------------------------------------------

# expects a data.frame() with 4 columns: line, time, value, Strata
UI_CompAbsPlot <- function(Plot_data,xlim,ylim,TimeLab) {
  
  ggplot(Plot_data, aes(x = time, y = value, colour = Strata)) +
    geom_line() +
    theme_bw() +
    scale_x_continuous(limits = c(0, xlim), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, ylim), expand = expand_scale(mult = c(0, .05))) +
    theme(
      legend.title = element_blank(),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.text = element_text(size = 10),
      axis.text=element_text(size=10),
      axis.title=element_text(size=12),
      strip.background = element_rect(colour = "black", fill = "#d4d4d4")
    ) +
    guides(colour = guide_legend(nrow = 1)) +
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
    labs(x = TimeLab, y = paste0("pr(E), per ",TimeLab))
}

# ~~~ Relative --------------------------------------------
# expects a data.frame() with 4 columns: line, time, value, Strata
UI_CompRelPlot <- function(Plot_data,xlim,ylim,TimeLab) {
  
  
  ggplot(Plot_data, aes(x = time, y = value)) +
    geom_line() +
    theme_bw() +
    scale_x_continuous(limits = c(0, xlim), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, ylim), expand = expand_scale(mult = c(0, .05))) +
    theme(
      legend.title = element_blank(),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.text = element_text(size = 10),
      axis.text=element_text(size=10),
      axis.title=element_text(size=12),
      strip.background = element_rect(colour = "black", fill = "#d4d4d4")
    ) +
    guides(colour = guide_legend(nrow = 1)) +
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
    labs(x = TimeLab, y = paste0("RR(E), per ",TimeLab))
}


# ~ interactive UI elements -------------------------------------------------

# ~~ hover panels for values ----------------------------------------------

# ~~~ Extrapolation plots -------------------------------------------------

# ~~~~ Extrapolations -----------------------------------------------------

# graphical float overlay showing curve values (closest line). wrap in renderUI({}).
UI_HoverValWindow_Extraps <- function(Plot_data,hover_location) {
  hover <- hover_location
  point <- nearPoints(
    Plot_data,
    hover,
    threshold = 5,
    maxpoints = 1,
    addDist = TRUE
  )
  if (nrow(point) == 0) return(NULL)
  
  # calculate point position INSIDE the image as percent of total dimensions
  # from left (horizontal) and from top (vertical)
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  
  # calculate distance from left and bottom side of the picture in pixels
  left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  
  #rgba is used here as it has a 4th dimension which is controlling transaprency
  b_colour <- "rgba(71, 186, 209, 0.81)"
  colour   <- "#ffffff"
  
  # create style property fot tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ",
                  "left:", left_px + 5, "px; top:", top_px + 5, "px;")
  
  # actual tooltip created as wellPanel
  wellPanel(
    style = style,
    p(HTML(
      paste0('<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#9ABAD9;}
.tg td{font-family:Arial, sans-serif;font-size:12px;padding:5px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#444;background-color:#EBF5FF;}
.tg th{font-family:Arial, sans-serif;font-size:12px;font-weight:normal;padding:5px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#fff;background-color:#409cff;}
.tg .tg-1wig{font-weight:bold;text-align:left;vertical-align:top}
.tg .tg-hmp3{background-color:#D2E4FC;text-align:left;vertical-align:top}
</style>
<table class="tg" style="undefined;table-layout: fixed; width: 183px">
<colgroup>
<col style="width: 73px">
<col style="width: 110px">
</colgroup>
  <tr>
    <th class="tg-1wig">Curve:</th>
    <th class="tg-1wig">',point$line,'</th>
  </tr>
  <tr>
    <td class="tg-1wig">Time:</td>
    <td class="tg-hmp3">',round(point$time,3),'</td>
  </tr>
  <tr>
    <td class="tg-1wig">Survival:</td>
    <td class="tg-hmp3">',round(point$value,4)*100,'%</td>
  </tr>
</table>')
    ))
  )
  
}

# ~~~~ LCHP -----------------------------------------------------
UI_HoverValWindow_LCHP <- function(Plot_data,hover_location) {
  graph_data <- Plot_data
  graph_data$time <- log(graph_data$time)
  graph_data$value <- log(-log(graph_data$value))
  
  hover <- hover_location
  point <- nearPoints(
    graph_data,
    hover,
    threshold = 5,
    maxpoints = 1,
    addDist = TRUE
  )
  if (nrow(point) == 0) return(NULL)
  
  # calculate point position INSIDE the image as percent of total dimensions
  # from left (horizontal) and from top (vertical)
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  
  
  # calculate distance from left and bottom side of the picture in pixels
  left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  
  #rgba is used here as it has a 4th dimension which is controlling transaprency
  b_colour <- "rgba(71, 186, 209, 0.81)"
  colour   <- "#ffffff"
  
  # create style property fot tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ",
                  "left:", left_px + 5, "px; top:", top_px + 5, "px;")
  
  wellPanel(
    style = style,
    p(HTML(
      paste0('<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#9ABAD9;}
.tg td{font-family:Arial, sans-serif;font-size:12px;padding:5px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#444;background-color:#EBF5FF;}
.tg th{font-family:Arial, sans-serif;font-size:12px;font-weight:normal;padding:5px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#fff;background-color:#409cff;}
.tg .tg-1wig{font-weight:bold;text-align:left;vertical-align:top}
.tg .tg-hmp3{background-color:#D2E4FC;text-align:left;vertical-align:top}
</style>
<table class="tg" style="undefined;table-layout: fixed; width: 183px">
<colgroup>
<col style="width: 73px">
<col style="width: 110px">
</colgroup>
  <tr>
    <th class="tg-1wig">Curve:</th>
    <th class="tg-1wig">',point$line,'</th>
  </tr>
  <tr>
    <td class="tg-1wig">ln(t):</td>
    <td class="tg-hmp3">',round(point$time,3),'</td>
  </tr>
  <tr>
    <td class="tg-1wig">ln(-ln(S(t))):</td>
    <td class="tg-hmp3">',round(point$value,3),'</td>
  </tr>
</table>')
    ))
  )
}


# ~~~~ digitizer image ----------------------------------------------------
UI_HoverValWindow_Digi <- function(hover_location,Max_X,Max_Y,info_text,direction_text) {
  hover <- hover_location
  
  if (is.null(hover)) return(NULL)
  
  # calculate point position INSIDE the image as percent of total dimensions
  # from left (horizontal) and from top (vertical)
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  
  # scaled
  left_scaled <- left_pct * Max_X
  top_scaled <- (1-top_pct) * Max_Y
  
  # calculate distance from left and bottom side of the picture in pixels
  left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  right_px <- hover$domain$right - hover$x
  top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  bot_px <- hover$domain$bottom - hover$y
  
  pr_from_left <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  pr_from_right <- (hover$domain$right - hover$x) / (hover$domain$right - hover$domain$left)
  pr_from_top <- (-hover$y) / (hover$domain$top - hover$domain$bottom)
  pr_from_bot <- 1-((hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom))
  
  nudge_left  <- (hover$domain$right - hover$domain$left) * -0.355
  nudge_right <- (hover$domain$right - hover$domain$left) * 0.02
  nudge_up    <- (hover$domain$top - hover$domain$bottom) * 0.35
  nudge_down  <- (hover$domain$top - hover$domain$bottom) * -0.01
  
  #rgba is used here as it has a 4th dimension which is controlling transaprency
  b_colour <- "rgba(71, 186, 209, 0.81)"
  colour   <- "#ffffff"
  
  # create style property fot tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  
  if (all(pr_from_left < 0.2,pr_from_top >= 0.2,pr_from_right >= 0.2, pr_from_bot >= 0.2)) {
    # print("too close to left, rendering window to bottom right")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_right, "px; top:", top_px + nudge_down, "px;")
  } else if (all(pr_from_left >= 0.2,pr_from_top < 0.2,pr_from_right >= 0.2, pr_from_bot >= 0.2)) {
    # print("too close to top, rendering window to bottom right")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_right, "px; top:", top_px + nudge_down, "px;")
  } else if (all(pr_from_left >= 0.2,pr_from_top >= 0.2,pr_from_right < 0.2, pr_from_bot >= 0.2)) {
    # print("too close to right, rendering window to bottom left")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_left, "px; top:", top_px + nudge_down, "px;")
  } else if (all(pr_from_left >= 0.2,pr_from_top >= 0.2,pr_from_right >= 0.2, pr_from_bot < 0.2)) {
    # print("too close to bottom, rendering window to top right")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_right, "px; top:", top_px + nudge_up, "px;")
  } else if (all(pr_from_left < 0.2,pr_from_top < 0.2,pr_from_right >= 0.2, pr_from_bot >= 0.2)) {
    # print("too close to top left, rendering window to bottom right")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_right, "px; top:", top_px + nudge_down, "px;")
  } else if (all(pr_from_left >= 0.2,pr_from_top < 0.2,pr_from_right < 0.2, pr_from_bot >= 0.2)) {
    # print("too close to top right, rendering window to bottom left")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_left, "px; top:", top_px + nudge_down, "px;")
  } else if (all(pr_from_left >= 0.2,pr_from_top >= 0.2,pr_from_right < 0.2, pr_from_bot < 0.2)) {
    # print("too close to bottom right, rendering window to top left")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_left, "px; top:", top_px + nudge_up, "px;")
  } else if (all(pr_from_left < 0.2,pr_from_top >= 0.2,pr_from_right >= 0.2, pr_from_bot < 0.2)) {
    # print("too close to bottom left, rendering window to top right")
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_right, "px; top:", top_px + nudge_up, "px;")
  } else {
    # not too close to anything!
    style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ","left:", left_px + nudge_right, "px; top:", top_px + nudge_down, "px;")
  }
  
  
  wellPanel(
    style = style,
    p(HTML(paste0(
      '<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#9ABAD9;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#444;background-color:#EBF5FF;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#fff;background-color:#409cff;}
.tg .tg-mim8{font-size:12px;border-color:#000000;text-align:left;vertical-align:top}
.tg .tg-5unb{font-weight:bold;font-size:16px;border-color:#000000;text-align:center;vertical-align:top}
.tg .tg-0ttd{font-weight:bold;font-size:12px;border-color:#000000;text-align:left;vertical-align:top}
.tg .tg-2kao{background-color:#D2E4FC;font-size:12px;border-color:#000000;text-align:left;vertical-align:top}
.tg .tg-73oq{border-color:#000000;text-align:left;vertical-align:top}
</style>
<table class="tg" style="undefined;table-layout: fixed; width: 387px">
<colgroup>
<col style="width: 124px">
<col style="width: 66px">
<col style="width: 136px">
<col style="width: 61px">
</colgroup>
  <tr>
    <th class="tg-l93j" colspan="4">Digitizer info:<br></th>
  </tr>
  <tr>
    <td class="tg-77x5">Px from left:</td>
    <td class="tg-w66w">',round(hover$x,0),'</td>
    <td class="tg-77x5">Px from top:</td>
    <td class="tg-w66w">',round(hover$y,0),'</td>
  </tr>
  <tr>
    <td class="tg-77x5">Px from right:<br></td>
    <td class="tg-w66w">',round(right_px,0),'</td>
    <td class="tg-77x5">Px from bottom:</td>
    <td class="tg-w66w">',round(bot_px,0),'</td>
  </tr>
  <tr>
    <td class="tg-xsvg">Time:</td>
    <td class="tg-7nlv">',round(left_scaled,5),'</td>
    <td class="tg-xsvg">Survival:</td>
    <td class="tg-7nlv">',round(top_scaled,5),'</td>
  </tr>
  <tr>
    <td class="tg-0lax" colspan="4">',info_text,'</td>
  </tr>
  <tr>
    <td class="tg-73a0" colspan="4">',direction_text,'</td>
  </tr>
</table>'
    )))
  )
}

#  ~~~~ OTHER -------------------------------------------------------------
# graphical float overlay showing curve values (closest line). wrap in renderUI({}).
UI_HoverValWindow_Comparisons <- function(Plot_data,hover_location,ylab = "Value",ByLabel = "none") {
  hover <- hover_location
  point <- nearPoints(
    Plot_data,
    hover,
    threshold = 5,
    maxpoints = 1,
    addDist = TRUE
  )
  if (nrow(point) == 0) return(NULL)
  
  # calculate point position INSIDE the image as percent of total dimensions
  # from left (horizontal) and from top (vertical)
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  
  # calculate distance from left and bottom side of the picture in pixels
  left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  
  #rgba is used here as it has a 4th dimension which is controlling transaprency
  b_colour <- "rgba(71, 186, 209, 0.81)"
  colour   <- "#ffffff"
  
  # create style property fot tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ",
                  "left:", left_px + 5, "px; top:", top_px + 5, "px;")
  
  if (ByLabel == "Strata") {
    TextString <- paste0(
      "<b> Strata: </b>" , point["Strata"], "<br/>",
      "<b> Time: </b>" , round(point$time, 3), "<br/>",
      "<b>", " ",ylab,": ","</b>", round(point$value, 3)
    )
  } else if (ByLabel == "line") {
    TextString <- paste0(
      "<b> Curve: </b>", point$line, "<br/>",
      "<b> Time: </b>" , round(point$time, 3), "<br/>",
      "<b>", " ",ylab,": ","</b>", round(point$value, 3)
    )
  } else {
    TextString <- paste0(
      "<b> Time: </b>" , round(point$time, 3), "<br/>",
      "<b>", " ",ylab,": ","</b>", round(point$value, 3)
    )
  }
  
  
  # actual tooltip created as wellPanel
  wellPanel(style = style,p(HTML(TextString)))
}


# graphical float overlay showing curve values (closest line). wrap in renderUI({}).
UI_HoverValWindow_ExpoFit <- function(hover_location) {
  hover <- hover_location
  
  # calculate point position INSIDE the image as percent of total dimensions
  # from left (horizontal) and from top (vertical)
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  s <- hover$y
  t <- hover$x
  
  
  # calculate distance from left and bottom side of the picture in pixels
  left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  
  #rgba is used here as it has a 4th dimension which is controlling transaprency
  b_colour <- "rgba(71, 186, 209, 0.81)"
  colour   <- "#ffffff"
  
  # create style property fot tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  style <- paste0("position:absolute; z-index:100; background-color: ",b_colour,"; color: ",colour,"; ",
                  "left:", left_px + 5, "px; top:", top_px + 5, "px;")
  
  # actual tooltip created as wellPanel
  wellPanel(
    style = style,
    p(HTML(
      paste0('<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#9ABAD9;}
.tg td{font-family:Arial, sans-serif;font-size:12px;padding:5px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#444;background-color:#EBF5FF;}
.tg th{font-family:Arial, sans-serif;font-size:12px;font-weight:normal;padding:5px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#9ABAD9;color:#fff;background-color:#409cff;}
.tg .tg-1wig{font-weight:bold;text-align:left;vertical-align:top}
.tg .tg-hmp3{background-color:#D2E4FC;text-align:left;vertical-align:top}
</style>
<table class="tg" style="undefined;table-layout: fixed; width: 183px">
<colgroup>
<col style="width: 73px">
<col style="width: 110px">
</colgroup>
  <tr>
    <td class="tg-1wig">Time:</td>
    <td class="tg-hmp3">',round(t,3),'</td>
  </tr>
  <tr>
    <td class="tg-1wig">Survival:</td>
    <td class="tg-hmp3">',round(s,4)*100,'%</td>
  </tr>
</table>')
    ))
  )
  
}



# ~~ Helper functions for hover panels ------------------------------------

#for clicking and hovering. need to enter as arguments the definitions you want for x and Y in the server, or it will default to x and y!
xy_str <- function(e,xdef="x",ydef="y") {
  if(is.null(e)) return("NULL\n")
  paste0(xdef,"=", round(e$x, 3),"; ", ydef," =", round(e$y, 4), "\n")
}
xy_range_str <- function(e,xdef="x",ydef="y") {
  if(is.null(e)) return("NULL\n")
  paste0(xdef,"min=", round(e$xmin, 3),"; ",xdef,"max=", round(e$xmax, 1), 
         "; ",ydef,"min=", round(e$ymin, 4),"; ",ydef,"max=", round(e$ymax, 4))
}

# A function to store the correct information from the digitization plot
Digitizer_xy_df <- function(click_info,scaling_factor_x=1,scaling_factor_y=1,StrataLabel = "Historical Kaplan-Meier") {
  if(is.null(click_info)) return()
  # % left to right * Max value of x-axis = value on X-axis
  # % bottom to top * Max value of Y-axis = value on Y axis
  data.frame(
    x = (1-((click_info$domain$right  - click_info$x) / click_info$domain$right))*scaling_factor_x,
    y = ((click_info$domain$bottom - click_info$y) / click_info$domain$bottom)*scaling_factor_y,
    Strata = StrataLabel
  )
}

# A function to store the correct information from the digitization plot
ExpoFit_xy_df <- function(click_info) {
  if(is.null(click_info)) return()
  data.frame(s = clickinfo$y,t = clickinfo$x)
}


#Credit to https://www.r-bloggers.com/an-enhanced-kaplan-meier-plot/
#’ Create a Kaplan-Meier plot using ggplot2
#’
#’ @param sfit a \code{\link[survival]{survfit}} object
#’ @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#’ @param returns logical: if \code{TRUE}, return an arrangeGrob object
#’ @param xlabs x-axis label
#’ @param ylabs y-axis label
#’ @param ystratalabs The strata labels. \code{Default = levels(summary(sfit)$strata)}
#’ @param ystrataname The legend name. Default = “Strata”
#’ @param timeby numeric: control the granularity along the time-axis
#’ @param main plot title
#’ @param pval logical: add the pvalue to the plot?
#’ @return a ggplot is made. if return=TRUE, then an arrangeGlob object
#’ is returned
#’ @author Abhijit Dasgupta with contributions by Gil Tomas
#’ \url{http://statbandit.wordpress.com/2011/03/08/an-enhanced-kaplan-meier-plot/}
#’ @export
#’ @examples
#’ \dontrun{
#’ data(colon)
#’  fit <- survfit(Surv(time,status)~rx, data=colon)
#'  ggkm(fit, timeby=500)
#' }
#’ Create a Kaplan-Meier plot using ggplot2
#’
#’ @param sfit a \code{\link[survival]{survfit}} object
#’ @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#’ @param returns logical: if \code{TRUE}, return an arrangeGrob object
#’ @param xlabs x-axis label
#’ @param ylabs y-axis label
#’ @param ystratalabs The strata labels. \code{Default = levels(summary(sfit)$strata)}
#’ @param ystrataname The legend name. Default = “Strata”
#’ @param timeby numeric: control the granularity along the time-axis
#’ @param main plot title
#’ @param pval logical: add the pvalue to the plot?
#’ @return a ggplot is made. if return=TRUE, then an arrangeGlob object
#’ is returned
#’ @author Abhijit Dasgupta with contributions by Gil Tomas
#’ \url{http://statbandit.wordpress.com/2011/03/08/an-enhanced-kaplan-meier-plot/}
#’ @export
#’ @examples
#’ \dontrun{
#’ data(colon)
#’  fit <- survfit(Surv(time,status)~rx, data=colon)
#'  ggkm(fit, timeby=500)
#' }
ggkm <- function(sfit, table = TRUE, returns = FALSE,
                 xlabs = "Time", ylabs = "survival probability",
                 ystratalabs = NULL, ystrataname = NULL,
                 timeby = 100, main = "Kaplan-Meier Plot",
                 pval = TRUE, ...) {
  require(ggplot2)
  require(survival)
  require(gridExtra)
  if(is.null(ystratalabs)) {
    ystratalabs <- as.character(levels(summary(sfit)$strata))
  }
  
  m <- max(nchar(ystratalabs))
  
  if(is.null(ystrataname)) ystrataname <- "Strata"
  
  times <- seq(0, max(sfit$time), by = timeby)
  
  .df <- data.frame(time = sfit$time, n.risk = sfit$n.risk,
                    n.event = sfit$n.event, surv = sfit$surv, strata = summary(sfit, censored = T)$strata,
                    upper = sfit$upper, lower = sfit$lower)
  
  levels(.df$strata) <- ystratalabs
  
  zeros <- data.frame(time = 0, surv = 1, strata = factor(ystratalabs, levels=levels(.df$strata)),
                      upper = 1, lower = 1)
  
  .df <- plyr::rbind.fill(zeros, .df)
  
  d <- length(levels(.df$strata))
  
  p <- ggplot(.df, aes(time, surv, group = strata)) +
    geom_step(aes(linetype = strata), size = 0.7) +
    theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    scale_x_continuous(xlabs, breaks = times, limits = c(0, max(sfit$time))) +
    scale_y_continuous(ylabs, limits = c(0, 1)) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position = c(ifelse(m < 10, .28, .35), ifelse(d < 4, .25, .35))) +
    theme(legend.key = element_rect(colour = NA)) +
    labs(linetype = ystrataname) +
    theme(plot.margin = unit(c(0, 1, .5, ifelse(m < 10, 1.5, 2.5)), "lines"))
  # theme(title = main)
  
  ## Create a blank plot for place-holding
  ## .df <- data.frame()
  blank.pic <- ggplot(.df, aes(time, surv)) +
    geom_blank() +
    theme_bw() +
    theme(axis.text.x = element_text(), axis.text.y = element_text(),
          axis.title.x = element_text(), axis.title.y = element_text(),
          axis.ticks = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank())
  
  if(pval) {
    sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <- pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
    pvaltxt <- ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
    p <- p + annotate("text", x = 0.6 * max(sfit$time), y = 0.1, label = pvaltxt)
  }
  
  if(table) {
    ## Create table graphic to include at-risk numbers
    risk.data <- data.frame(strata = summary(sfit, times = times, extend = TRUE)$strata,
                            time = summary(sfit, times = times, extend = TRUE)$time,
                            n.risk = summary(sfit, times = times, extend = TRUE)$n.risk)
    data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      #, color = strata)) +
      geom_text(size = 3.5) +
      theme_bw() +
      scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = ystratalabs) +
      # scale_y_discrete(#format1ter = abbreviate,
      # breaks = 1:3,
      # labels = ystratalabs) +
      scale_x_continuous("Numbers at risk", limits = c(0, max(sfit$time))) +
      theme(axis.title.x = element_text(size = 10, vjust = 1), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.border = element_blank(),
            axis.text.x = element_text(), axis.ticks = element_blank(),
            axis.text.y = element_text(face = "bold", hjust = 1))
    data.table <- data.table + theme(legend.position = "none") +
      xlab(NULL) + ylab(NULL)
    data.table <- data.table +
      theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5) - 0.28 * m), "lines"))
    ## Plotting the graphs
    ## p <- ggplotGrob(p)
    ## p <- addGrob(p, textGrob(x = unit(.8, "npc"), y = unit(.25, "npc"), label = pvaltxt,
    ## gp = gpar(fontsize = 12)))
    grid.arrange(p, blank.pic, data.table,
                 clip = FALSE, nrow = 3, ncol = 1,
                 heights = unit(c(2, .1, .25),c("null", "null", "null")))
    if(returns) {
      a <- arrangeGrob(p, blank.pic, data.table, clip = FALSE,
                       nrow = 3, ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null")))
      return(a)
    }
  } else {
    ## p <- ggplotGrob(p)
    ## p <- addGrob(p, textGrob(x = unit(0.5, "npc"), y = unit(0.23, "npc"),
    ## label = pvaltxt, gp = gpar(fontsize = 12)))
    print(p)
    if(returns) return(p)
  }
}


# Plot of extrapolations --------------------------------------------------
ExtrapPlotDat_FromCEMod <- function(EXTRAPS,CE_model_loc,DistNameNamedRangeName) {
  # pass the data to more manageable objects
  t <- EXTRAPS$time
  S <- EXTRAPS$Strata
  exDF <- data.frame(EXTRAPS)
  exDF_ExOnly <- exDF %>% dplyr::select(-time,-Strata)
  colnam_orig <- as.character(readWorkbook(xlsxFile = CE_model_loc,namedRegion = DistNameNamedRangeName,colNames = FALSE,rowNames = FALSE))
  
  # process the data
  graph_data <- do.call(
    "rbind",
    purrr::map(
      .x = 1:(ncol(exDF_ExOnly)),
      .f = function(COL) {
        
        # extract data
        Thiscol <- data.frame(time = t, value = exDF_ExOnly[,COL], Strata = S)
        nam <- rep(colnam_orig[COL],nrow(exDF))
        
        # build data.frame to compile together
        Thiscol$line <- nam
        
        Thiscol <- data.frame(
          line = Thiscol$line,
          time = Thiscol$time,
          value = Thiscol$value,
          Strata = Thiscol$Strata
        )
        
        return(Thiscol)
      }
    )
  )
  return(graph_data)
}


# ~ digitzation plot ------------------------------------------------------
Digi_CropImage <- function(usr_image,up_by,down_by,left_by,right_by,Apply_grid = FALSE) {
  # get size of the image AFTER zooming it when loading it up
  width  <- usr_image %>% image_info() %>% select(width) %>% (function(x) x * 1) %>% as.numeric()
  height <- usr_image %>% image_info() %>% select(height) %>% (function(x) x * 1) %>% as.numeric()
  
  # create a string to use in the crop command that responds to user inputs
  CropString <- paste0(
    width - sum(left_by,right_by), 
    "x",
    height - sum(up_by, down_by),
    "-",
    left_by,
    "+",
    up_by
  )
  
  
  # load the ZOOMED image (i.e. load then resize) and then crop it using the user inputs
  OUT <- usr_image %>%
    image_resize(paste0(width,"x",height)) %>%
    image_crop(CropString) %>%
    image_ggplot() +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.1)) +
    labs(x = "Pixels from top left", y = "Pixels from top left")
  
  if (Apply_grid) {
    OUT <- OUT + 
      geom_vline(xintercept=seq(1, width, by=5),size = 0.01,alpha = 0.1)+
      geom_hline(yintercept=seq(1, height, by=5),size = 0.01,alpha = 0.1) 
  }
  
  return(OUT)
  
}

#Digi_CropImage(image_read("~/../Downloads/ovarianecog.png"),0,0,0,0,TRUE)



# Functions for table layouts ----------------------------------------------------

BasicTableWithCopy <- function(DF) {
  datatable(
    data = DF,
    extensions = 'Buttons',
    rownames = FALSE,
    options = list(
      scrollX = TRUE,
      lengthChange = FALSE,
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      ordering = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf')
    )
  )
}
BasicTableWithpaging <- function(DF) {
  maxLen <- nrow(DF)
  datatable(
    data = DF,
    rownames = FALSE,
    options = list(
      scrollX = TRUE,
      lengthChange = TRUE,
      lengthMenu = c(10,50,100,maxLen),
      paging = TRUE,
      searching = FALSE,
      info = FALSE
    )
  )
}

BasicTableWithXscroll <- function(DF) {
  datatable(
    data = DF,
    rownames = FALSE,
    options = list(
      scrollX = TRUE,
      lengthChange = FALSE,
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      ordering = FALSE
    )
  )
}


# Excel survival output model ---------------------------------------------
# A function to take either stratified or unstratified models and generate
#   the necessary outputs within an excel book
GenerateExcelOutput <- function(KM_table,Extrap_table,AIC_BIC_table,Params_table,vcov_table,DEBUG=FALSE) {
  require(openxlsx)
  
  # print data coming in for debugging
  if (DEBUG) {
    print(KM_table)
    print(Extrap_table)
    print(AIC_BIC_table)
    print(Params_table)
  }  
  
  # data manipulation
  
  ## Extrapolation table
  # The N list elements have been rbound according to strata. we need to split them up again
  #   It is set up to be a list even if the data isn't stratified (i.e. list of length 1)
  #   That way, the code here doesn't need to be sensitive to the dataset treatment
  #   if it's unstratified no action is needed
  if (!is.na(length(unique(Extrap_table$Strata)))) {
    Extrap_table_Excel <-
      lapply(1:length(unique(Extrap_table$Strata)), function(STRATA) {
        Extrap_table %>% filter(Strata == unique(Extrap_table$Strata)[STRATA])
      })
  } else {
    Extrap_table_Excel <- Extrap_table
  }
  
  ## LCHP KM: log-cumulative hazard of reported KM graph data
  KM_LCHP_table <- lapply(1:length(KM_table), function(STRATA){
    data.frame(
      ln_t = log(KM_table[[STRATA]]$t),
      ln_negln_s_t = log(-log(KM_table[[STRATA]]$s_t)),
      Strata = KM_table[[STRATA]]$Strata
    )
  })
  
  
  ## LCHP extrapolations table
  if (!is.na(length(unique(Extrap_table$Strata)))) {
    LCHP_table <-
      lapply(1:length(unique(Extrap_table$Strata)), function(STRATA) {
        DF <- Extrap_table_Excel[[STRATA]]
        # Bind it back together into a data.frame() after the inner lapply is done
        DF <- do.call(
          'cbind',
          # don't want to try to log the strata column!
          lapply(1:7, function(COL){
            if (COL == 1) {
              # time column
              log(DF[,COL])
            } else {
              # Other columns (extrapolated survival)
              log(-log(DF[,COL]))
            }
          })
        )
      })
  } else {
    LCHP_table <-
      list(
        do.call(
          'cbind',
          # don't want to try to log the strata column!
          lapply(1:7, function(COL){
            if (COL == 1) {
              # time column
              log(Extrap_table_Excel[[1]][,COL])
            } else {
              # Other columns (extrapolated survival)
              log(-log(Extrap_table_Excel[[1]][,COL]))
            }
          })
        )
      )
  }
  
  ### Params_table
  # simply separate into a list so we can programatically write multiple tables from that list
  if (!is.null(Params_table$stratification)) {
    Params_table <- lapply(1:length(unique(Params_table$stratification)), function(STRATA){
      Params_table %>% filter(stratification == unique(Params_table$stratification)[STRATA])
    })
  } else {
    Params_table <- list(Params_table = Params_table)
  }
  
  # Excel workbook generation
  # Make empty excel workbook
  wb <- createWorkbook()
  
  # aesthetics
  options("openxlsx.borderColour" = "#4F80BD")
  options("openxlsx.borderStyle" = "thin")
  modifyBaseFont(wb, fontSize = 10, fontName = "Arial Narrow")
  
  # Make necessary worksheets
  addWorksheet(wb, sheetName = "KaplanMeier", gridLines = FALSE)
  addWorksheet(wb, sheetName = "Plot_extrap_Fit", gridLines = FALSE)
  addWorksheet(wb, sheetName = "LCHP", gridLines = FALSE)
  addWorksheet(wb, sheetName = "Params", gridLines = FALSE)
  addWorksheet(wb, sheetName = "varianceCovariance", gridLines = FALSE)
  
  
  ## sheet: KaplanMeier
  freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) ## freeze first row and column
  
  ### programatically generate N tables of KM plot data where N is number of strata
  for (STRATA in 1:length(KM_table)) {
    writeDataTable(
      wb = wb,
      sheet = "KaplanMeier",
      x = data.frame(KM_table[[STRATA]]),
      startCol = 1+((STRATA-1)*4),
      colNames = TRUE,
      rowNames = FALSE,
      tableStyle = "TableStyleLight9",
      tableName = paste0("KMTable_",STRATA)
    )
  }
  
  ## Sheet: Plot_extrap_Fit
  ### Extrapolation data
  writeDataTable(
    wb,
    sheet = "Plot_extrap_Fit",
    AIC_BIC_table,
    startCol = "B",
    startRow = 3,
    tableStyle = "TableStyleLight9",
    tableName =  "AIC_BIC_Table"
  )
  
  ## Sheet: Plot_extrap_Fit
  ### programatically generate N tables of extrapolation data where N is number of strata
  for (STRATA in 1:length(Extrap_table_Excel)) {
    # Kaplan Meier for easy graphing
    writeDataTable(
      wb,
      sheet = "Plot_extrap_Fit",
      data.frame(KM_table[[STRATA]]),
      startCol = 7+((STRATA-1)*13),
      startRow = 3,
      tableStyle = "TableStyleLight9",
      tableName = paste0("KMExtrapTable_",STRATA)
    )
    # Extrapolations 
    writeDataTable(
      wb,
      sheet = "Plot_extrap_Fit",
      data.frame(Extrap_table_Excel[[STRATA]]),
      startCol = 10+((STRATA-1)*13),
      startRow = 3,
      tableStyle = "TableStyleLight9",
      tableName = paste0("ExtrapTable_",STRATA)
    )
  }
  
  ## Sheet: LCHP
  for (STRATA in 1:length(Extrap_table_Excel)) {
    # add in the LCHP for the KM
    writeDataTable(
      wb = wb,
      sheet = "LCHP",
      x = data.frame(KM_LCHP_table[[STRATA]]),
      startCol = 2+((STRATA-1)*13),
      startRow = 3,
      colNames = TRUE,
      rowNames = FALSE,
      tableStyle = "TableStyleLight9",
      tableName = paste0("KMLCHPTable_",STRATA)
    )
    # add in the LCHP for the extrapolations
    writeDataTable(
      wb,
      sheet = "LCHP",
      data.frame(LCHP_table[[STRATA]]),
      startCol = 5+((STRATA-1)*13),
      startRow = 3,
      tableStyle = "TableStyleLight9",
      tableName = paste0("LCHPTable_",STRATA)
    )
  }
  
  ## Sheet: Params
  for (STRATA in 1:length(Params_table)) {
    writeDataTable(
      wb,
      sheet = "Params",
      data.frame(Params_table[[STRATA]]),
      startCol = 2+((STRATA-1)*7),
      startRow = 3,
      tableStyle = "TableStyleLight9",
      tableName = paste0("ParamsTable_",STRATA)
    )
  }
  
  ## Sheet: variance-covariance matrices
  
  # loop through tables required
  if (class(vcov_table[[1]]) == "list") {
    # the analysis is stratified (as going 1 deep in the list still shows a list of results)
    for (STRATA in 1:length(vcov_table)) {
      # gather some information useful for spacing of all the tables
      # gen gamma is always the biggest table and is last in the list
      LargestVcov <- dim(vcov_table[[STRATA]][[6]])[1]
      
      # now go through the individual models spacing them nicely
      for (MODEL in 1:length(vcov_table[[STRATA]])) {
        ThisVcovSize <- dim(vcov_table[[STRATA]][[MODEL]])[1]
        ThisStrataName <- paste0("Strata: ",STRATA)
        ThisModelName <- c("exponential","weibull","log_logistic","log_normal","gompertz","gen_gamm")[MODEL]
        # identifiers
        # openxlsx::writeData(
        #   wb = wb,
        #   sheet = "varianceCovariance",
        #   x = cbind(ThisStrataName,ThisModelName),
        #   startCol = 2+((STRATA-1)*(LargestVcov + 1)),
        #   startRow = (2+((MODEL-1)*(ThisVcovSize + 5)))-2
        # )
        
        # data table
        writeDataTable(
          wb,
          sheet = "varianceCovariance",
          data.frame(vcov_table[[STRATA]][[MODEL]]),
          startCol = 2+((STRATA-1)*(LargestVcov + 1)),
          startRow = 4+((MODEL-1)*(ThisVcovSize + 5)),
          tableStyle = "TableStyleLight9",
          tableName = paste0("VcovTable_",STRATA,"_",MODEL)
        )
      }
    }
  } else {
    # the analysis is not stratified, list only goes one level deep
    # gather some information useful for spacing of all the tables
    # gen gamma is always the biggest table and is last in the list
    
    # now go through the individual models spacing them nicely
    for (MODEL in 1:length(vcov_table)) {
      # adjust for size of matrix
      ThisVcovSize <- dim(vcov_table[[MODEL]])[1]
      ThisModelName <- c("exponential","weibull","log_logistic","log_normal","gompertz","gen_gamm")[MODEL]
      # identifiers
      # openxlsx::writeData(
      #   wb = wb,
      #   sheet = "varianceCovariance",
      #   x = ThisModelName,
      #   startCol = 2,
      #   startRow = (3+((MODEL-1)*(ThisVcovSize + 2)))-2,
      # )
      # data
      writeDataTable(
        wb,
        sheet = "varianceCovariance",
        data.frame(vcov_table[[MODEL]]),
        startCol = 2,
        startRow = 3+((MODEL-1)*(ThisVcovSize + 3)),
        tableStyle = "TableStyleLight9",
        tableName = paste0("VcovTable_",STRATA,"_",MODEL)
      )
    }
  }
  
  return(wb)
  
}

# #TESTING AREA
#
# KM_table <- 
# AIC_BIC_table <- read.csv("C:/Users/dburns/Downloads/AIC_BIC.csv") %>% dplyr::select(-X)
# Params_table <- read.csv("C:/Users/dburns/Downloads/PARAM.csv") %>% dplyr::select(-X)
# Extrap_table <- read.csv("C:/Users/dburns/Downloads/plotdata.csv") %>% dplyr::select(-X)
# GenerateExcelOutput(Extrap_table = Extrap_table,AIC_BIC_table = AIC_BIC_table,Params_table = Params_table)

# if (!is.na(unique(Extrap_table$Strata))) {
#   KM_data <- lapply(1:length(unique(Extrap_table$Strata)), function(STRATA){
#     Extrap_table %>%
#       filter(line == "Kaplan-Meier") %>% 
#       dplyr::select(-line) %>% 
#       filter(Strata == unique(Extrap_table$Strata)[STRATA])
#   })
# } else {
#   KM_data <- list(Extrap_table %>%
#     filter(line == "Kaplan-Meier") %>% 
#     dplyr::select(-line))
# }
# 
# # These can now be put into excel as separate tables in the output excel book programatically
# 
# # The same is true of the extrapolations as well
# # processing of extrapolation data into a form for excel
# OutPutTable <- Extrap_table %>% filter(line != "Kaplan-Meier")
# # separate back out into a list
# #  Note: as the KM has been taken out above, the times are identical for all models
# if (!is.na(unique(Extrap_table$Strata))) {
#   OutPutTable <- lapply(1:length(unique(OutPutTable$Strata)), function(STRATA){
#     lapply(0:length(unique(OutPutTable$line)), function(MODEL){
#       if (MODEL == 0) {
#         # pull out only the time column for the first one
#         DF <- OutPutTable %>%
#           filter(Strata == unique(OutPutTable$Strata)[STRATA])
#         DF <- data.frame(time = unique(DF$time))
#       } else {
#         # pull out the associated distribution with identifier
#         DF <- OutPutTable %>%
#           filter(Strata == unique(OutPutTable$Strata)[STRATA]) %>%
#           filter(line == unique(OutPutTable$line)[MODEL])
#         DF <- DF %>% dplyr::select(line,value)
#       }
#     })
#   })
# } else {
#   # if the data is not stratified at all, then Strata == NA is true, just create one table
#   OutPutTable <- lapply(0:length(unique(OutPutTable$line)), function(MODEL){
#     if (MODEL == 0) {
#       # pull out only the time column for the first one
#       DF <- OutPutTable
#       DF <- data.frame(time = unique(DF$time))
#     } else {
#       # pull out the associated distribution with identifier
#       DF <- OutPutTable %>%
#         filter(line == unique(OutPutTable$line)[MODEL])
#       DF <- DF %>% dplyr::select(line,value)
#     }
#   })
# }
# 
# # compile list of available distributions
# if (!is.na(unique(Extrap_table$Strata))) {
#   # stratified setting: look in first element of 1 time element + 6 model extrapolation elements
#   Dists <-
#     do.call('rbind', lapply(1:length(OutPutTable[[1]]), function(MODEL) {
#       as.character(OutPutTable[[1]][[MODEL]][1, 1])
#     }))
#   Dists[1,1] <- "time"
# } else {
#   # unstratified setting: depth level of list is reduced by 1
#   Dists <-
#     do.call('rbind', lapply(1:length(OutPutTable), function(MODEL) {
#       as.character(OutPutTable[[MODEL]][1, 1])
#     }))
#   Dists[1,1] <- "time"
# }
#   
# 
# 
# # reformat everything back into separate dataframes, stick them together
# #   and rename properly
# if (!is.na(unique(Extrap_table$Strata))) {
#   for (j in 1:length(OutPutTable)) {
#     for (i in 2:7) {
#       OutPutTable[[j]][[i]] <- OutPutTable[[j]][[i]][,2] %>% data.frame()
#     }
#     OutPutTable[[j]] <- do.call('cbind',OutPutTable[[j]])
#     names(OutPutTable[[j]]) <- as.character(Dists[,1])
#     OutPutTable$Strata <- unique(Extrap_table)[j]
#   }
# } else {
#   for (i in 2:7) {
#     OutPutTable[[i]] <- OutPutTable[[i]][,2] %>% data.frame()
#   }
#   OutPutTable <- do.call('cbind',OutPutTable)
#   names(OutPutTable) <- as.character(Dists[,1])
#   OutPutTable <- list(OutPutTable)
# }
# 
# ### Params_table
# # simply separate into a list so we can programatically write multiple tables from that list
# if (!is.null(Params_table$stratification)) {
#   Params_table <- lapply(1:length(unique(Params_table$stratification)), function(STRATA){
#     Params_table %>% filter(stratification == unique(Params_table$stratification)[STRATA])
#   })
# } else {
#   Params_table <- list(Params_table = Params_table)
# }
# 
# # Excel workbook generation
# # Make empty excel workbook
# wb <- createWorkbook()
# 
# # aesthetics
# options("openxlsx.borderColour" = "#4F80BD")
# options("openxlsx.borderStyle" = "thin")
# modifyBaseFont(wb, fontSize = 10, fontName = "Arial Narrow")
# 
# # Make necessary worksheets
# addWorksheet(wb, sheetName = "KaplanMeier", gridLines = FALSE)
# addWorksheet(wb, sheetName = "Plot_extrap_Fit", gridLines = FALSE)
# addWorksheet(wb, sheetName = "Params", gridLines = FALSE)
# 
# ## sheet: KaplanMeier
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) ## freeze first row and column
# 
# # programatically generate N tables of KM plot data where N is number of strata
# for (STRATA in 1:length(KM_data)) {
#   writeDataTable(
#     wb = wb,
#     sheet = "KaplanMeier",
#     x = data.frame(KM_data[[STRATA]]),
#     startCol = 1+((STRATA-1)*4),
#     colNames = TRUE,
#     rowNames = FALSE,
#     tableStyle = "TableStyleLight9"
#   )
# }
# 
# ## Sheet: Plot_extrap_Fit
# ### Extrapolation data
# writeDataTable(wb, sheet = "Plot_extrap_Fit", AIC_BIC_table, startCol = "B", startRow = 3,tableStyle = "TableStyleLight9")
# 
# ### Extraps
# # programatically generate N tables of extrapolation data where N is number of strata
# for (STRATA in 1:length(KM_data)) {
#   writeDataTable(
#     wb,
#     sheet = "Plot_extrap_Fit",
#     OutPutTable[[STRATA]],
#     startCol = 7+((STRATA-1)*9),
#     startRow = 3,
#     tableStyle = "TableStyleLight9"
#   )
# }
# 
# ## Sheet: Params
# writeDataTable(wb, sheet = "Params", PARAM, startCol = "B", startRow = 3,tableStyle = "TableStyleLight9")
# 
# saveWorkbook(wb = wb, file = "~/TEST.xlsx")




# Generate save state -----------------------------------------------------


# ~ List of empty lists ---------------------------------------------------
# simple: generate a list using a vector of names of the items as strings
GenSaveStateList <- function(Items_as_stringvec) {
  OUT <- list()
  OUT <- lapply(1:length(Items_as_stringvec),function(LIST_ITEM){
    OUT[[Items_as_stringvec[LIST_ITEM]]] <- list()
  })
  names(OUT) <- Items_as_stringvec
  return(OUT)
}




# UI functions ------------------------------------------------------------

# creates a big blue button that takes a slot in the middle
UI_NiceActionButton <- function(inputId,label,icon) {
  fluidRow(style = "padding-bottom:25px;padding-top:25px;",
           column(
             12,
             align = "center",
             actionBttn(
               inputId = inputId,
               label = label,
               style = "material-flat",
               icon = icon,
               color = "primary"
             )
           ))
}





#use: get_median_survival(survfit(Surv(t, e) ~ 1))

# testing (please comment out) --------------------------------------------
# # #get data
# x <- ovarian
# t <- x[, "futime"] %>% as.numeric()
# e <- x[, "fustat"] %>% as.numeric()
# #run models
# models <- list(
#   exponential  = flexsurvreg(Surv(t, e) ~ resid.ds, data = x, dist = "exponential"),
#   weibull      = flexsurvreg(Surv(t, e) ~ resid.ds, data = x, dist = "weibull")    ,
#   log_logistic = flexsurvreg(Surv(t, e) ~ resid.ds, data = x, dist = "llogis")     ,
#   log_normal   = flexsurvreg(Surv(t, e) ~ resid.ds, data = x, dist = "lnorm")      ,
#   gompertz     = flexsurvreg(Surv(t, e) ~ resid.ds, data = x, dist = "gompertz")   ,
#   gen_gamma    = flexsurvreg(Surv(t, e) ~ resid.ds, data = x, dist = "gengamma")
# )
# 
# extrapolate_models(models = models,time_steps = 1:10,Category = "resid.ds",Category_level = 1)



# #extrapolate (these are the two steps we do in the shiny model as the models have already been run by this point!)
# #extraps <- extrapolate_models(models,0:5000)
# 
# 
# results <- lapply(models,function(x) as.data.frame(x$res.t[,c("est","se")]))
# #because exponential is annoying
# if (ncol(results$exponential) == 1) {
#   results$exponential <- t(results$exponential)
#   rownames(results$exponential) <- "rate"
# }
# do.call(rbind,lapply(results,function(x) x))
# 
# models$exponential$res
# models$exponential$coefficients
# 
# lapply(models,function(x) as.data.frame(x$coef[,c("est","se")]))


