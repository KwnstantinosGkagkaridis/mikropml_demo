install.packages("mikropml")
library(mikropml)


source("code/genus_process.R")

srn_genus_data <- composite %>% 
  select(group, taxonomy, rel_abund, srn) %>%
  pivot_wider(names_from = taxonomy, values_from = rel_abund) %>%
  select(-group) %>%
  mutate(srn = if_else(srn, "srn", "healthy")) %>%
  select(srn, everything())



#practice <- composite %>%
#  select(group, fit_result, site, gender, srn, weight) %>%
#  distinct() %>%
#  mutate(weight = na_if(weight, 0)) %>%
#  select(-group) #%>%
  #mutate(near_zero_variance = 23) %>% #create a column with no variance in order to be removed later
  #mutate(perfect_corr = fit_result) #create a column with perfect correlation with another column in order to be removed later


#preprocessed <- preprocess_data(practice, outcome_colname = "srn")
#summary(preprocessed$dat_transformed)




srn_genus_preprocess <- preprocess_data(srn_genus_data, 
                                       outcome_colname = "srn")$dat_transformed





srn_genus_results_no_pp <- run_ml(srn_genus_data, 
                            method = "glmnet",
                            outcome_colname = "srn", #what it wants to predict
                            kfold = 5,
                            cv_times = 100,
                            training_frac = 0.8,
                            seed = 19760620) 




test_hp <- list(alpha = 0, 
                lambda = c(0.1, 1, 2, 3, 4, 5, 10))

get_srn_genus_results <- function(seed){
  
  srn_genus_results <- run_ml(srn_genus_preprocess, #instead of srn_genus_data
                              method = "glmnet",
                              outcome_colname = "srn", #what it wants to predict
                              kfold = 5,
                              cv_times = 100,
                              training_frac = 0.8,
                              hyperparameters = test_hp,
                              seed = seed) 
}



iterative_run_ml_results <- map(c(1,2,3), get_srn_genus_results)


performance <- iterative_run_ml_results %>%
  map(pluck, "trained_model") %>%
  combine_hp_performance()

plot_hp_performance(performance$dat, lambda, AUC)


# how to know the different hyperparameters options and their default values
get_hyperparams_list(srn_genus_preprocess, "glmnet")
#get_hyperparams_list(srn_genus_preprocess, "rf")
