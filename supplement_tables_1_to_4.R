library(dlim)
library(knitr)
library(ggplot2)
library(dplyr)
library(xtable)


## ----Initialize results table-------------------------------------------------------------------
table_results <- c()

### ------------------ Run penalized simulation ------------------ ###

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
    
    #DLIM(10,10)
    table_results <- rbind(table_results, 
                           fit_models(type=type, 
                                      SNR=SNR, 
                                      n=n, 
                                      df=c(10,10), 
                                      model_type="standard", 
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


## ----Run non-penalized simualtions-------------------------------------------------------------------
#set.seed(013023)
n=200
method="REML"
penalize = F

for(SNR in c(0.5, 1, 10)){ #0.5, 1, 10
  print(SNR)
  
  for(type in c(1, 2, 3, 4)){ #1, 2, 3, 4
    print(type)
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLM
    table_results <- rbind(table_results, fit_models(type=type, SNR=SNR, n=n, df=c(5,5), model_type="DLM", method=method, penalize=penalize))
    save(table_results, file="nested_models_results_ns.rda")
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLIM-linear
    table_results <- rbind(table_results, fit_models(type=type, SNR=SNR, n=n, df=c(5,5), model_type="linear", method=method, penalize=penalize))
    save(table_results, file="nested_models_results_ns.rda")
    
    #reset seed for every model to use the same 200 data sets
    set.seed(013023)
    
    #DLIM(5,5)
    table_results <- rbind(table_results, fit_models(type=type, SNR=SNR, n=n, df=c(5,5), model_type="standard", method=method, penalize=penalize))
    save(table_results, file="nested_models_results_ns.rda")
  }
}


## -------Prep for table----------------------------------------------------------------------
#load penalized results
load("nested_models_results.rda")
table_results_data_ps <- table_results

#load non-penalized results
load("nested_models_results_ns.rda")
table_results_data_ns <- table_results
table_results_data_ns$Model[table_results_data_ns$Model=="DLIM-linear"] <- "DLIM-linear-ns"
table_results_data_ns$Model[table_results_data_ns$Model=="DLM"] <- "DLM-ns"
table_results_data_ns$Model[table_results_data_ns$Model=="DLIM(5,5)"] <- "DLIM-ns(5,5)"

#combine penalized and non-penalized results
table_results_supp <- rbind(table_results_data_ps,table_results_data_ns)

#proceed as before with cleaning
space <- "fixed"#default is fixed
colnames(table_results_supp)[which(colnames(table_results_supp)=="Type")] <- "Scenario"
table_results_supp$SNR[table_results_supp$SNR==0.5] <- "low"
table_results_supp$SNR[table_results_supp$SNR==1] <- "medium"
table_results_supp$SNR[table_results_supp$SNR==10] <- "high"
table_results_supp$SNR <- factor(table_results_supp$SNR, levels = c("low", "medium", "high"))

table_results_supp$Model <- factor(table_results_supp$Model, levels = c("DLIM(20,20)", "DLIM(10,10)", "DLIM-linear", "DLM", "DLIM-ns(5,5)", "DLIM-linear-ns", "DLM-ns"))



## --------Create table with penalized and non-penalized results---------------------------------------------------------------------
table_sum_results_supp <- table_results_supp %>%
  group_by(Scenario, SNR, Model) %>%
  dplyr::summarize(Cumul_RMSE = round(mean(Cumul_RMSE),3),
                   Cumul_coverage = round(mean(Cumul_cover),2),
                   Cumul_width = round(mean(Cumul_width),3),
                   DLF_RMSE = round(mean(DLF_RMSE),3),
                   DLF_coverage = round(mean(DLF_cover),2),
                   DLF_width = round(mean(DLF_width),3))


xtable(table_sum_results_supp, align=c("|c","|c","|c","|c","|c","|c","|c","|c","|c","|c|"))
