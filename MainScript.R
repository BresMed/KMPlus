# script to generate KM+ estimate based on example data

source("UDFs.R")

dbloader(c(
  "flexsurv",
  "tidyverse"
))

dat <- lung

# from the beginning
sfit <- survfit(Surv(dat$time,dat$status) ~ dat$sex, data = dat)


sforms <- lapply(1:length(unique(dat[,"sex"])), function(STRATA) {
  SurvregFormulaGen(
    t = dat %>% dplyr::filter(sex == unique(dat$sex)[1]) %>% dplyr::select(time) %>% (function(x) x$time),
    e = dat %>% dplyr::filter(sex == unique(dat$sex)[STRATA]) %>% select(status) %>% (function(x) x$status),
    covs = NULL,
    factors = NULL,
    nam_t = "time",
    nam_e = "status",
    nam_covs = NULL,
    nam_factors = NULL,
    DEBUG = FALSE
  )
})

dists     <- c("exponential","weibull","llogis",      "lnorm",     "gompertz","gengamma")


MODELS <- lapply(1:length(unique(dat[,"sex"])), function(STRATA) {
  lapply(1:6,function(MODEL){
    flexsurvreg(sforms[[1]], data = dat, dist = dists[MODEL])
  })
})



extraps <- lapply(1:length(unique(dat[,"sex"])), function(STRATA) {
  extrapolate_models(
    models = MODELS[[STRATA]],
    time_steps = seq(
      0,
      2000,
      1
    )
  ) %>%
    mutate(Strata = unique(dat[,"sex"])[STRATA])
})

time_steps <- seq(0,2000,1)

extrapolations1 <- data.frame(
  time         = time_steps,
  exponential  = summary(MODELS[[1]][[1]],t=time_steps,ci=FALSE)[[1]]$est,
  weibull      = summary(MODELS[[1]][[2]],t=time_steps,ci=FALSE)[[1]]$est,
  log_logistic = summary(MODELS[[1]][[3]],t=time_steps,ci=FALSE)[[1]]$est,
  log_normal   = summary(MODELS[[1]][[4]],t=time_steps,ci=FALSE)[[1]]$est,
  gompertz     = summary(MODELS[[1]][[5]],t=time_steps,ci=FALSE)[[1]]$est,
  gen_gamma    = summary(MODELS[[1]][[6]],t=time_steps,ci=FALSE)[[1]]$est
)

extrapolations2 <- data.frame(
  time         = time_steps,
  exponential  = summary(MODELS[[2]][[1]],t=time_steps,ci=FALSE)[[1]]$est,
  weibull      = summary(MODELS[[2]][[2]],t=time_steps,ci=FALSE)[[1]]$est,
  log_logistic = summary(MODELS[[2]][[3]],t=time_steps,ci=FALSE)[[1]]$est,
  log_normal   = summary(MODELS[[2]][[4]],t=time_steps,ci=FALSE)[[1]]$est,
  gompertz     = summary(MODELS[[2]][[5]],t=time_steps,ci=FALSE)[[1]]$est,
  gen_gamma    = summary(MODELS[[2]][[6]],t=time_steps,ci=FALSE)[[1]]$est
)


EXTRAPS <- list(
  extrapolations1,
  extrapolations2
)

ggplot(EXTRAPS[[1]],aes(x = time, y = exponential)) + geom_line()

plot(sfit)


# KM to a certain point and then an exponential line