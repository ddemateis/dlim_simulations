library(dlim)
library(dplyr)
library(xtable)
library(dlnm)
library(mgcv)

### ------------------ Functions ------------------ ###


fit_models <- function(type, SNR, n, df, model_type, 
                       penalize=T, new_mods=NULL, 
                       pen_fn = "ps", mod_fn = list(fn="runif"), 
                       mod_arg = NULL, lag_arg = NULL,  ...){
  
  library(dlnm)
  library(mgcv)
  source("../Functions/bias.R")
  source("../Functions/cover.R")
  source("../Functions/RMSE.R")
  
  data("exposure")
  
  summary_list <- function(X, truth, func) {
    do.call(func,list(X-truth))
  }
  
  gamma_bias_sq <- c()
  gamma_bias <- c()
  gamma_cover <- c()
  CE_RMSE <- c()
  CE_bias <- c()
  PW_RMSE <- c()
  PW_bias <- c()
  CE_cover <- c()
  CE_width <- c()
  PW_cover <- c()
  PW_width <- c()
  
  for(j in 1:n){
    print(j)
    #initialize metrics for data set j (not actually necessary)
    gamma_ests <- c()
    gamma_SE <- c()
    betas <- list()
    betas_cumul <- list()
    
    #generate data
    modifiers <- runif(1000)
    
    if(is.null(new_mods)){
      new_mods <- modifiers
    }
    data <- sim_data(x=exposure, modifiers=modifiers, type=type,
                     SNR=SNR, ncovariates = 3, gamma = 1)
    true_betas <- sim_dlf(L = 36, modifiers = new_mods, type = type)
    true_betas_cumul <- colSums(true_betas)
    
    # Models 2-4
    if(model_type=="linear" | model_type=="quad" | model_type=="standard"){
      model4 <- sim_dlim(data=data, df_m=df[1], df_l=df[2], penalize=penalize,
                         fit_dlm=F, model_type=model_type, pen_fn = pen_fn, 
                         mod_arg = mod_arg, lag_arg = lag_arg, ...)
      #print(model4$fit$method)
      gamma_idx <- which(names(model4$fit$coefficients)=="Zmodifiers")
      gamma_ests[1] <- model4$fit$coefficients[gamma_idx]
      gamma_SE[1] <- sqrt(vcov(model4$fit)[gamma_idx, gamma_idx])
      pred4 <- predict(model4, newdata = new_mods, type=c("DLF", "CE"))
      betas[[1]] <- pred4$est_dlim$betas
      betas_LB <- pred4$est_dlim$LB
      betas_UB <- pred4$est_dlim$UB
      betas_cumul[[1]] <- pred4$est_dlim$betas_cumul
      betas_cumul_LB <- pred4$est_dlim$cumul_LB
      betas_cumul_UB <-pred4$est_dlim$cumul_UB
    } else if(model_type=="DLM"){
      if(penalize){
        cb1 <- crossbasis(x=data$x,argvar=list(fun="lin"),arglag = list(fun="ps",df=10))
        penalty <- cbPen(cb1)
        model1 <- gam(data$y ~ data$modifiers + cb1 + data$Z, paraPen = list(cb1 = penalty), ...)
        #print(model1$method)
      }else{
        cb1 <- crossbasis(x=data$x,argvar=list(fun="lin"),arglag = list(fun="ns",df=5))
        model1 <- gam(data$y ~ data$modifiers + cb1 + data$Z)
      }
      gamma_idx <- which(names(model1$coefficients)=="data$modifiers")
      gamma_ests[1] <- model1$coefficients[gamma_idx]
      gamma_SE[1] <- sqrt(vcov(model1)[gamma_idx, gamma_idx])
      dlm_crosspred <- crosspred(cb1,model1,at=rep(1,37),cen = F)
      betas[[1]] <- matrix(rep(dlm_crosspred$matfit,length(new_mods)),nrow=length(new_mods),byrow=T)
      betas_LB <- matrix(rep(dlm_crosspred$matlow,length(new_mods)),nrow=length(new_mods),byrow=T)
      betas_UB <- matrix(rep(dlm_crosspred$mathigh,length(new_mods)),nrow=length(new_mods),byrow=T)
      betas_cumul[[1]] <- rep(dlm_crosspred$allfit,length(new_mods))
      betas_cumul_LB <- rep(dlm_crosspred$alllow,length(new_mods))
      betas_cumul_UB <- rep(dlm_crosspred$allhigh,length(new_mods))
    }
    
    #compute gamma squared bias, bias, cover for data set j
    gamma_bias_sq <- rbind(gamma_bias_sq,(gamma_ests - rep(1,1))^2)
    gamma_bias <- rbind(gamma_bias,gamma_ests - rep(1,1))
    gamma_cover <- rbind(gamma_cover, apply(cbind(gamma_ests - qnorm(0.975)*gamma_SE, rep(1,1), gamma_ests + qnorm(0.975)*gamma_SE),1, FUN=cover))
    
    #compute cumulative effect RMSE, bias, cover for data set j
    CE_RMSE <- rbind(CE_RMSE, unlist(lapply(betas_cumul, summary_list, truth=true_betas_cumul, func="RMSE")))
    CE_bias <- rbind(CE_bias, unlist(lapply(betas_cumul, summary_list, truth=true_betas_cumul, func="bias")))
    CE_width <- rbind(CE_width, mean(apply(cbind(betas_cumul_LB, betas_cumul_UB), 1, width)))
    CE_cover <- rbind(CE_cover, mean(apply(cbind(betas_cumul_LB, true_betas_cumul, betas_cumul_UB), 1, cover)))
    
    #compute pointwise effect RMSE, bias, cover for data set j
    PW_RMSE <- rbind(PW_RMSE, unlist(lapply(betas, summary_list, truth=t(true_betas), func="RMSE")))
    PW_bias <- rbind(PW_bias, unlist(lapply(betas, summary_list, truth=t(true_betas), func="bias")))
    PW_width <- rbind(PW_width, mean(betas_UB - betas_LB))
    PW_cover <- rbind(PW_cover, mean(betas_LB < t(true_betas) & t(true_betas) < betas_UB))
  }
  
  if(model_type=="DLM"){
    model_name <- "DLM"
  }else if(model_type=="linear"){
    model_name <- "DLIM-linear"
  }else if(model_type=="quad"){
    model_name <- "DLIM-quad"
  }else if(model_type=="standard"){
    model_name <- paste0("DLIM(",df[1],",",df[2],")")
  }
  
  results <- data.frame(SNR=rep(SNR,n),
                        Type=rep(type,n),
                        Model=rep(model_name,n),
                        Gamma_RMSE=signif(sqrt(gamma_bias_sq),4),
                        Gamma_bias=signif(gamma_bias,4),
                        Gamma_coverage=signif(gamma_cover,4),
                        Cumul_RMSE=signif(CE_RMSE,4),
                        Cumul_bias=signif(CE_bias,4),
                        Cumul_width=signif(CE_width,4),
                        Cumul_cover=signif(CE_cover,4),
                        DLF_RMSE=signif(PW_RMSE,4),
                        DLF_bias=signif(PW_bias,4),
                        DLF_width=signif(PW_width,4),
                        DLF_cover=signif(PW_cover,4))
  
  
  return(results)
  
}

### ------------------ Run PS(2,2) simulation ------------------ ###

## Initialize results table
table_results <- c()

## Run penalized simulation
n=200
method="REML"
penalize = T
pen_fn <- "ps"
lag_arg <- NULL #second order difference matrix for time dim
mod_arg <- NULL #second order difference matrix for modifier dim

for(SNR in c(0.5, 1, 10)){ #0.5, 1, 10
  print(SNR)
  
  for(type in c(1, 2, 3, 4)){ #1, 2, 3, 4
    print(type)
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLIM(20,20)
    table_results <- rbind(table_results, 
                           fit_models(type=type, 
                                      SNR=SNR, 
                                      n=n, 
                                      df=c(20,20), 
                                      model_type="standard", 
                                      method=method, 
                                      penalize=penalize, 
                                      pen_fn = pen_fn,
                                      mod_arg = mod_arg,
                                      lag_arg = lag_arg))
    save(table_results, file=paste0("nested_models_results_ps22",".rda"))
  }
}


### ------------------ Run PS(2,1) simulation ------------------ ###

## Initialize results table
table_results <- c()

## Run penalized simulation
n=200
method="REML"
penalize = T
pen_fn <- "ps"
lag_arg <- NULL #second order difference matrix for time dim
mod_arg <- list(diff = 1) #first order difference matrix for modifier dim

for(SNR in c(0.5, 1, 10)){ #0.5, 1, 10
  print(SNR)
  
  for(type in c(1, 2, 3, 4)){ #1, 2, 3, 4
    print(type)
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLIM(20,20)
    table_results <- rbind(table_results, 
                           fit_models(type=type, 
                                      SNR=SNR, 
                                      n=n, 
                                      df=c(20,20), 
                                      model_type="standard", 
                                      method=method, 
                                      penalize=penalize, 
                                      pen_fn = pen_fn,
                                      mod_arg = mod_arg,
                                      lag_arg = lag_arg))
    save(table_results, file=paste0("nested_models_results_ps21",".rda"))
  }
}

### ------------------ Run PS(1,2) simulation ------------------ ###

## Initialize results table
table_results <- c()

## Run penalized simulation
n=200
method="REML"
penalize = T
pen_fn <- "ps"
lag_arg <- list(diff = 1) #first order difference matrix for time dim
mod_arg <- NULL #second order difference matrix for modifier dim

for(SNR in c(0.5, 1, 10)){ #0.5, 1, 10
  print(SNR)
  
  for(type in c(1, 2, 3, 4)){ #1, 2, 3, 4
    print(type)
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLIM(20,20)
    table_results <- rbind(table_results, 
                           fit_models(type=type, 
                                      SNR=SNR, 
                                      n=n, 
                                      df=c(20,20), 
                                      model_type="standard", 
                                      method=method, 
                                      penalize=penalize, 
                                      pen_fn = pen_fn,
                                      mod_arg = mod_arg,
                                      lag_arg = lag_arg))
    save(table_results, file=paste0("nested_models_results_ps12",".rda"))
  }
}

### ------------------ Run CR simulation ------------------ ###

## Initialize results table
table_results <- c()

## Run penalized simulation
n=200
method="REML"
penalize = T
pen_fn <- "cr"
lag_arg <- NULL 
mod_arg <- NULL 

for(SNR in c(0.5, 1, 10)){ #0.5, 1, 10
  print(SNR)
  
  for(type in c(1, 2, 3, 4)){ #1, 2, 3, 4
    print(type)
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLIM(20,20)
    table_results <- rbind(table_results, 
                           fit_models(type=type, 
                                      SNR=SNR, 
                                      n=n, 
                                      df=c(20,20), 
                                      model_type="standard", 
                                      method=method, 
                                      penalize=penalize, 
                                      pen_fn = pen_fn,
                                      mod_arg = mod_arg,
                                      lag_arg = lag_arg))
    save(table_results, file=paste0("nested_models_results_cr",".rda"))
  }
}

### ----------------- Summarize in table ----------------- ###

#use dplyr to summarize over the 200 simulated data sets for 
#each penalization method




