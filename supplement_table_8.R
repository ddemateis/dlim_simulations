library(ggplot2)
library(dlim)
library(dlnm)
library(mgcv)
library(dplyr)
library(xtable)

### ---------------------- Functions ---------------------- ###

sim_pois_data <- function(x, L=NULL, modifiers, type=2, ncovariates=0, gamma=1){
  
  #create lagged structure
  if(is.vector(x)){
    X <- Lag(x,0:L)[-c(1:L),]
    modifiers <- modifiers[-c(1:L)]
  }else{
    L <- ncol(x)-1
    X <- x
    modifiers <- modifiers
  }
  
  #Create Betas
  betas <- sim_dlf(L,modifiers,type)/10 #scaling so mean count is reasonable
  betas_cumul <- colSums(betas)
  
  y_mean <- colSums(t(X)*betas)
  
  #Create gammas and covariates
  if(ncovariates!=0){
    gammas <- c(gamma,matrix(rnorm(ncovariates),ncol=1))
    Z <- matrix(rnorm(length(modifiers)*(ncovariates)), ncol = ncovariates)
    mod_Z <- cbind(modifiers,Z)
    mu <- exp(y_mean + mod_Z%*%gammas)
    y <- rpois(length(modifiers), lambda = mu)
  }else{
    Z <- NULL
    gammas <- gamma/10 #scaling so mean count is reasonable
    mu <- exp(y_mean + matrix(modifiers,ncol=1)*gammas)
    y <- rpois(length(modifiers), lambda = mu)
  }
  
  result <- list(x=X,
                 L=L,
                 modifiers=modifiers,
                 y=y,
                 betas=betas,
                 betas_cumul=betas_cumul,
                 Z = Z,
                 gammas = gammas)
  
  return(result)
}

df <- c()
n_sim=100
array_idx <- as.numeric(commandArgs(trailingOnly=TRUE)) #read in array values from the shell file, which is just the task number
type <- array_idx
set.seed(013023)
print(type)
plot_df_cumul <- c()
plot_df_pw <- c()

for(i in 1:n_sim){
  print(i)
  dta <- sim_pois_data(x = exposure,
                       modifiers = c(0.25,0.5,0.75,runif(997)),
                       type = type)
  fit <- dlim(y = dta$y,
              x = dta$x,
              modifiers = dta$modifiers,
              df_m = 20,
              df_l = 20,
              penalize = T,
              method = "REML", #arguments for gam below
              family = "poisson")
  pred <- predict(fit)
  
  cb_dlm <- crossbasis(x=dta$x,argvar=list(fun="lin"),arglag = list(fun="ps",df=20))
  penalty <- cbPen(cb_dlm)
  model_dlm <- bam(dta$y~cb_dlm+dta$modifiers,family = "poisson", paraPen = list(cb_dlm = penalty), method="REML")
  dlm_crosspred <- crosspred(cb_dlm,model_dlm,at=rep(1,37),cen = F)
  beta_by_lag <- matrix(rep(dlm_crosspred$matfit, 1000),nrow=37)
  cumul_betas <- dlm_crosspred$allfit
  z <- qnorm(1 - (1 - 0.95)/2)
  cumul_lb <- dlm_crosspred$allfit - z * dlm_crosspred$allse 
  cumul_ub <- dlm_crosspred$allfit + z * dlm_crosspred$allse
  dlm_lb <- matrix(rep(dlm_crosspred$matfit - z * dlm_crosspred$matse, 1000), nrow=37) #for some reason $matlow is gone
  dlm_ub <- matrix(rep(dlm_crosspred$matfit + z * dlm_crosspred$matse, 1000), nrow=37) #for some reason $mathigh is gone
  
  cumul_RMSE <- sqrt(mean((dta$betas_cumul - pred$est_dlim$betas_cumul)^2))
  cumul_RMSE_dlm <- sqrt(mean((dta$betas_cumul - cumul_betas)^2))
  cumul_coverage <- mean(pred$est_dlim$cumul_LB < dta$betas_cumul & dta$betas_cumul < pred$est_dlim$cumul_UB)
  cumul_coverage_dlm <- mean(cumul_lb < dta$betas_cumul & dta$betas_cumul < cumul_ub)
  cumul_width <- mean(pred$est_dlim$cumul_UB - pred$est_dlim$cumul_LB)
  cumul_width_dlm <- mean(cumul_ub - cumul_lb)
  PW_RMSE <- sqrt(mean((dta$betas - t(pred$est_dlim$betas))^2))
  PW_RMSE_dlm <- sqrt(mean((dta$betas - beta_by_lag)^2))
  PW_coverage <- mean(pred$est_dlim$LB < t(dta$betas) & t(dta$betas) < pred$est_dlim$UB)
  PW_coverage_dlm <- mean(dlm_lb < dta$betas & dta$betas < dlm_ub)
  PW_width <- mean(pred$est_dlim$UB - pred$est_dlim$LB)
  PW_width_dlm <- mean(dlm_ub - dlm_lb)
  
  df <- rbind(df, data.frame(Scenario = rep(type,2),
                             Model = c("DLIM", "DLM"),
                             Cumul_RMSE = c(cumul_RMSE, cumul_RMSE_dlm),
                             Cumul_coverage = c(cumul_coverage, cumul_coverage_dlm),
                             Cumul_width = c(cumul_width, cumul_width_dlm),
                             PW_RMSE = c(PW_RMSE, PW_RMSE_dlm),
                             PW_coverage = c(PW_coverage, PW_coverage_dlm),
                             PW_width = c(PW_width, PW_width_dlm)))  
}

save(df, file=paste0("df", array_idx, ".rda"))


summ_table <- df %>% 
  group_by(Scenario, Model) %>%
  dplyr::summarize(Cumul_RMSE = round(mean(Cumul_RMSE),3),
                   Cumul_coverage = round(mean(Cumul_coverage),2),
                   Cumul_width = round(mean(Cumul_width), 3),
                   PW_RMSE = round(mean(PW_RMSE),3),
                   PW_coverage = round(mean(PW_coverage),2),
                   PW_width = round(mean(PW_width),3)) 


xtable(summ_table, align=c("|c","|c","|c","|c","|c","|c","|c|","|c|","|c|"), digits = c(0,0,0,3,2,3,3,2,3))
