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




### ------------------ Run simulation ------------------ ###

## Initialize results table
table_results <- c()

## Run penalized simulation
n=200
method="REML"
penalize = T
pen_fn <- "ps"
lag_arg <- NULL 
mod_arg <- NULL

for(SNR in c(0.5, 1, 10)){ #0.5, 1, 10
  print(SNR)
  
  for(type in c(1, 2, 3, 4)){ #1, 2, 3, 4
    print(type)
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLM
    table_results <- rbind(table_results, 
                           fit_models(type=type, 
                                      SNR=SNR, 
                                      n=n, 
                                      df=c(10,10), 
                                      model_type="DLM", 
                                      method=method, 
                                      penalize=penalize,
                                      pen_fn = pen_fn,
                                      mod_arg = mod_arg,
                                      lag_arg = lag_arg))
    save(table_results, file=paste0("nested_models_results",".rda"))
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLIM-linear
    table_results <- rbind(table_results, 
                           fit_models(type=type,
                                      SNR=SNR, 
                                      n=n, 
                                      df=c(10,10), 
                                      model_type="linear", 
                                      method=method, 
                                      penalize=penalize, 
                                      pen_fn = pen_fn,
                                      mod_arg = mod_arg,
                                      lag_arg = lag_arg))
    save(table_results, file=paste0("nested_models_results",".rda"))
    
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
    save(table_results, file=paste0("nested_models_results",".rda"))
  }
}

### ------------------ Process into table ------------------ ###

#load("../Paper_simulations/nested_models_results.rda")

space <- "fixed"#default is fixed
table_results_data <- table_results
colnames(table_results_data)[which(colnames(table_results_data)=="Type")] <- "Scenario"
table_results_data$SNR[table_results_data$SNR==0.5] <- "low"
table_results_data$SNR[table_results_data$SNR==1] <- "medium"
table_results_data$SNR[table_results_data$SNR==10] <- "high"
table_results_data$SNR <- factor(table_results_data$SNR, levels = c("low", "medium", "high"))
table_results_data$Model <- factor(table_results_data$Model, levels = c("DLIM(20,20)","DLIM(10,10)", "DLIM(5,5)", "DLIM-linear", "DLM"))
table_sum_results_main <- table_results_data %>% filter(Model !="DLIM(10,10)") %>%
  group_by(Scenario, SNR, Model) %>%
  dplyr::summarize(Cumul_RMSE = round(mean(Cumul_RMSE),3),
                   Cumul_coverage = round(mean(Cumul_cover),2),
                   DLF_RMSE = round(mean(DLF_RMSE),3),
                   DLF_coverage = round(mean(DLF_cover),2)) 


xtable(table_sum_results_main, align=c("|c","|c","|c","|c","|c","|c","|c","|c|"), digits = c(0,0,0,0,3,2,3,2))

### ------------------ Model comparison ------------------ ###


### run on HPC with array_idx from 0:2399

source("model_comparison.R")
library(dlim)
library(dlnm)
library(mgcv)

array_idx <- as.numeric(commandArgs(trailingOnly=TRUE)) #read in array values from the shell file, which is just the task number

#set up grids
type_grid <- 1:4 #length 4
SNR_grid <- c(0.5, 1, 10) #length 3
sim_sets <- 1:200 #length 200
n_type <- length(type_grid)
n_SNR <- length(SNR_grid)
N <- length(sim_sets)

#choose type, SNR
type <- type_grid[floor(array_idx/(n_SNR*N))+1]
SNR <- SNR_grid[(floor(array_idx/N))%%n_SNR+1]
sim_set <- sim_sets[array_idx%%N+1]

B <- 1000 #number of bootstrap samples
rejects <- c()

#simulate data
#since I have no way of recovering each individual 
#data set from the main simulations study, need to
#simulate until I get the same data set
set.seed(013023)
for(i in 1:sim_set){
  dta <- sim_data(x = exposure, 
                  modifiers = runif(1000),
                  SNR=SNR,
                  ncovariates = 3,
                  type = type)
}

#fit full model
model_full <- dlim(y = dta$y,
                   x = dta$x,
                   modifiers = dta$modifiers,
                   z = dta$Z,
                   df_m = 20,
                   df_l = 20,
                   penalize = T,
                   method = "REML")

#model comparison to DLM
decision <- model_comparison(fit = model_full,
                             x = exposure, 
                             null = "DLM",
                             B = B)

#record whether decision was rejection
rejects[1] <- ifelse(decision == "reject", T, F)

#model comparison to DLIM-linear
decision <- model_comparison(fit = model_full,
                             x = exposure, 
                             null = "linear",
                             B = B)

#record whether decision was rejection
rejects[2] <- ifelse(decision == "reject", T, F)

#fit full DLIM-linear
linear_ful <- dlim(y = dta$y,
                   x = dta$x,
                   modifiers = dta$modifiers,
                   z = dta$Z,
                   df_m = 20,
                   df_l = 20,
                   penalize = T,
                   model_type = "linear",
                   method = "REML")

#model comparison to DLM
decision <- model_comparison(fit = linear_ful,
                             x = exposure, 
                             null = "DLM",
                             B = B)

#record whether decision was rejection
rejects[3] <- ifelse(decision == "reject", T, F)

#create results dataframe
results <- data.frame(Type = rep(type,3),
                      SNR = rep(SNR,3),
                      Model = c("DLIM & DLM", 
                                "DLIM & DLIM-linear",
                                "DLIM-linear & DLM"),
                      Reject = matrix(rejects,ncol=1))

save(results, file=paste0("results/model_comparison_results_",array_idx,".rda"))

sum_table <- c()
for(idx in 0:2399){
  load(paste0("~/AA CSU/Research/DLIM/Post_submission_work/model_comparison_results/model_comparison_results_",idx,".rda"))
  sum_table <- rbind(sum_table, results)
}
View(sum_table %>% 
       group_by(Type, SNR, Model) %>%
       dplyr::summarize(Prob_Reject = round(mean(Reject),2)) )
