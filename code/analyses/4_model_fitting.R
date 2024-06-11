##############################################################################################
#######################---------------- Model fitting ----------------########################
##############################################################################################

### In this script: 
# (1) Preprocess data
# (2) Model fitting
# (3) Convergence check
# (4) Model comparison
# (5) Posterior predictive checks

# Set working directory
here::i_am("github/semaglutide-study/code/analyses/4_model_fitting.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/plot_funs.R")
source("github/semaglutide-study/code/functions/model_preprocess_fun.R")
source("github/semaglutide-study/code/functions/model_convergence_check_fun.R")
source("github/semaglutide-study/code/functions/model_comparison_fun.R")
source("github/semaglutide-study/code/functions/parameter_estimates_fun.R")
source("github/semaglutide-study/code/functions/extract_posterior_predictions_fun.R")


# source dataset
main_data <- readRDS("data/processed_data/main_data.RDS")
non_diabetic_data <- readRDS("data/processed_data/non_diabetic_matched.RDS")
non_diabetic_normal_weight_data <- readRDS("data/processed_data/non_diabetic_normal_weight_matched.RDS")

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, 
                 writexl, lubridate, magrittr, pushoverr, nlme, gridExtra)
library(cmdstanr)
# run model fitting? 
model_fitting <- FALSE


### (1) Preprocess data -----------------------------------------------

task_data <- main_data$task_data %>% 
  filter(phase == "game") %>% 
  select(subj_id:group, trial:choice) %>% 
  rename(effort_a = effort, amount_a = reward) %>%
  add_column(effort_b = 0, amount_b = 1, 
             .after = "amount_a")

# Treatment S1
task_data_treat_s1 <- task_data %>% 
  filter(group == "treatment", 
         session == 1)
model_dat_treat_s1 <- model_preprocessing(raw_data = task_data_treat_s1, 
                                 retest = FALSE,
                                 subjs = unique(task_data_treat_s1 %>% .$subj_id), 
                                 n_subj = length(unique(task_data_treat_s1 %>% .$subj_id)), 
                                 t_subjs = aggregate(trial ~ subj_id, FUN = max, data = task_data_treat_s1)[,2], 
                                 t_max = max(aggregate(trial ~ subj_id, FUN = max, data = task_data_treat_s1)[,2]))

# Treatment S2
task_data_treat_s2 <- task_data %>% 
  filter(group == "treatment", 
         session == 2)
model_dat_treat_s2 <- model_preprocessing(raw_data = task_data_treat_s2, 
                                              retest = FALSE,
                                              subjs = unique(task_data_treat_s2 %>% .$subj_id), 
                                              n_subj = length(unique(task_data_treat_s2 %>% .$subj_id)), 
                                              t_subjs = aggregate(trial ~ subj_id, FUN = max, data = task_data_treat_s2)[,2], 
                                              t_max = max(aggregate(trial ~ subj_id, FUN = max, data = task_data_treat_s2)[,2]))

# Control
task_data_control <- task_data %>% 
  filter(group == "control")
model_dat_control <- model_preprocessing(raw_data = task_data_control, 
                                              retest = FALSE,
                                              subjs = unique(task_data_control %>% .$subj_id), 
                                              n_subj = length(unique(task_data_control %>% .$subj_id)), 
                                              t_subjs = aggregate(trial ~ subj_id, FUN = max, data = task_data_control)[,2], 
                                              t_max = max(aggregate(trial ~ subj_id, FUN = max, data = task_data_control)[,2]))
 
# Non-diabetics
task_data_non_diabetics <- non_diabetic_data$task_data %>% 
  filter(phase == "game") %>% 
  select(subj_id, trial:choice) %>% 
  rename(effort_a = offerEffort, amount_a = offerReward, trial_type = trialType) %>%
  add_column(effort_b = 0, amount_b = 1, 
             .after = "amount_a")

model_dat_non_diabetics <- model_preprocessing(raw_data = task_data_non_diabetics, 
                                         retest = FALSE,
                                         subjs = unique(task_data_non_diabetics %>% .$subj_id), 
                                         n_subj = length(unique(task_data_non_diabetics %>% .$subj_id)), 
                                         t_subjs = aggregate(trial ~ subj_id, FUN = max, data = task_data_non_diabetics)[,2], 
                                         t_max = max(aggregate(trial ~ subj_id, FUN = max, data = task_data_non_diabetics)[,2]))

# Non-overweight
task_data_normal_weight <- non_diabetic_normal_weight_data$task_data %>% 
  filter(phase == "game") %>% 
  select(subj_id, trial:choice) %>% 
  rename(effort_a = offerEffort, amount_a = offerReward, trial_type = trialType) %>%
  add_column(effort_b = 0, amount_b = 1, 
             .after = "amount_a")

model_dat_normal_weight <- model_preprocessing(raw_data = task_data_normal_weight, 
                                               retest = FALSE,
                                               subjs = unique(task_data_normal_weight %>% .$subj_id), 
                                               n_subj = length(unique(task_data_normal_weight %>% .$subj_id)), 
                                               t_subjs = aggregate(trial ~ subj_id, FUN = max, data = task_data_normal_weight)[,2], 
                                               t_max = max(aggregate(trial ~ subj_id, FUN = max, data = task_data_normal_weight)[,2]))


### (2) Model fitting -----------------------------------------------

if(model_fitting){
  # Load models
  m1_parabolic_stan_model <- cmdstanr::cmdstan_model("github/semaglutide-study/code/stan/models_parabolic/ed_m1_parabolic.stan")
  m2_parabolic_stan_model <- cmdstanr::cmdstan_model("github/semaglutide-study/code/stan/models_parabolic/ed_m2_parabolic.stan")
  m3_parabolic_stan_model <- cmdstanr::cmdstan_model("github/semaglutide-study/code/stan/models_parabolic/ed_m3_parabolic.stan")
  m1_linear_stan_model <- cmdstanr::cmdstan_model("github/semaglutide-study/code/stan/models_linear/ed_m1_linear.stan")
  m2_linear_stan_model <- cmdstanr::cmdstan_model("github/semaglutide-study/code/stan/models_linear/ed_m2_linear.stan")
  m3_linear_stan_model <- cmdstanr::cmdstan_model("github/semaglutide-study/code/stan/models_linear/ed_m3_linear.stan")
  
  ### Treatment group - Session 1 -----  
  ## Parabolic discounting models
  # Model 1
  m1_para_treat_s1_fit <- m1_parabolic_stan_model$sample(
    data = model_dat_treat_s1, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_para_treat_s1_check <- convergence_check(m1_para_treat_s1_fit, 
                                              params = c("kE", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m1_para_treat_s1_check$trace_plot
  saveRDS(list(m1_para_treat_s1_check$Rhat, m1_para_treat_s1_check$ess), 
          here::here("data/model_fits/treatment_s1/m1_para_treat_s1_check.RDS"))
  # LOO for model comparisons
  m1_para_treat_s1_loo <- m1_para_treat_s1_fit$loo()
  saveRDS(m1_para_treat_s1_loo, here::here("data/model_fits/treatment_s1/m1_para_treat_s1_loo.RDS"))
  # Parameter estimates
  m1_para_treat_s1_params <- get_params(subj_id = unique(task_data_treat_s1$subj_id), 
                                        model_fit = m1_para_treat_s1_fit, 
                                        n_subj = length(unique(task_data_treat_s1$subj_id)), 
                                        n_params = 2, 
                                        param_names = c("kE", "a"))
  saveRDS(m1_para_treat_s1_params, here::here("data/model_fits/treatment_s1/m1_para_treat_s1_params.RDS"))
  
  # Model 2
  m2_para_treat_s1_fit <- m2_parabolic_stan_model$sample(
    data = model_dat_treat_s1, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_para_treat_s1_check <- convergence_check(m2_para_treat_s1_fit, 
                                     params = c("kE", "kR"), 
                                     Rhat = TRUE, ess = TRUE,
                                     trace_plot = TRUE, rank_hist = FALSE)
  m2_para_treat_s1_check$trace_plot
  saveRDS(list(m2_para_treat_s1_check$Rhat, m2_para_treat_s1_check$ess), 
          here::here("data/model_fits/treatment_s1/m2_para_treat_s1_check.RDS"))
  # LOO for model comparisons
  m2_para_treat_s1_loo <- m2_para_treat_s1_fit$loo()
  saveRDS(m2_para_treat_s1_loo, here::here("data/model_fits/treatment_s1/m2_para_treat_s1_loo.RDS"))
  # Parameter estimates
  m2_para_treat_s1_params <- get_params(subj_id = unique(task_data_treat_s1$subj_id), 
                               model_fit = m2_para_treat_s1_fit, 
                               n_subj = length(unique(task_data_treat_s1$subj_id)), 
                               n_params = 2, 
                               param_names = c("kE", "kR"))
  saveRDS(m2_para_treat_s1_params, here::here("data/model_fits/treatment_s1/m2_para_treat_s1_params.RDS"))
  
  # Model 3
  m3_para_treat_s1_fit <- m3_parabolic_stan_model$sample(
    data = model_dat_treat_s1, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_para_treat_s1_check <- convergence_check(m3_para_treat_s1_fit, 
                                              params = c("kE", "kR", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m3_para_treat_s1_check$trace_plot
  saveRDS(list(m3_para_treat_s1_check$Rhat, m3_para_treat_s1_check$ess), 
          here::here("data/model_fits/treatment_s1/m3_para_treat_s1_check.RDS"))
  # LOO for model comparisons
  m3_para_treat_s1_loo <- m3_para_treat_s1_fit$loo()
  saveRDS(m3_para_treat_s1_loo, here::here("data/model_fits/treatment_s1/m3_para_treat_s1_loo.RDS"))
  # Posterior Predictions (for target model only)
  m3_para_treat_s1_PPC_dat <- posterior_predictions(csv_paths = m3_para_treat_s1_fit$output_files(),
                                           n_chains = 4,
                                           n_iter = (6000),
                                           n_subj = 58, 
                                           n_trials = 64,
                                           real_dat = task_data_treat_s1) 
  saveRDS(m3_para_treat_s1_PPC_dat, here::here("data/model_fits/treatment_s1/m3_para_treat_s1_PPC.RDS"))
  # Parameter estimates
  m3_para_treat_s1_params <- get_params(subj_id = unique(task_data_treat_s1$subj_id), 
                                        model_fit = m3_para_treat_s1_fit, 
                                        n_subj = length(unique(task_data_treat_s1$subj_id)), 
                                        n_params = 3, 
                                        param_names = c("kE", "kR", "a"))
  saveRDS(m3_para_treat_s1_params, here::here("data/model_fits/treatment_s1/m3_para_treat_s1_params.RDS"))

  ## Linear discounting models
  # Model 1
  m1_lin_treat_s1_fit <- m1_linear_stan_model$sample(
    data = model_dat_treat_s1, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_lin_treat_s1_check <- convergence_check(m1_lin_treat_s1_fit, 
                                              params = c("kE", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m1_lin_treat_s1_check$trace_plot
  saveRDS(list(m1_lin_treat_s1_check$Rhat, m1_lin_treat_s1_check$ess), 
          here::here("data/model_fits/m1_lin_treat_s1_check.RDS"))
  # LOO for model comparisons
  m1_lin_treat_s1_loo <- m1_lin_treat_s1_fit$loo()
  saveRDS(m1_lin_treat_s1_loo, here::here("data/model_fits/treatment_s1/m1_lin_treat_s1_loo.RDS"))
  # Parameter estimates
  m1_lin_treat_s1_params <- get_params(subj_id = unique(task_data_treat_s1$subj_id), 
                                        model_fit = m1_lin_treat_s1_fit, 
                                        n_subj = length(unique(task_data_treat_s1$subj_id)), 
                                        n_params = 2, 
                                        param_names = c("kE", "a"))
  saveRDS(m1_lin_treat_s1_params, here::here("data/model_fits/treatment_s1/m1_lin_treat_s1_params.RDS"))
  
  # Model 2
  m2_lin_treat_s1_fit <- m2_linear_stan_model$sample(
    data = model_dat_treat_s1, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_lin_treat_s1_check <- convergence_check(m2_lin_treat_s1_fit, 
                                              params = c("kE", "kR"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m2_lin_treat_s1_check$trace_plot
  saveRDS(list(m2_lin_treat_s1_check$Rhat, m2_lin_treat_s1_check$ess), 
          here::here("data/model_fits/treatment_s1/m2_lin_treat_s1_check.RDS"))
  # LOO for model comparisons
  m2_lin_treat_s1_loo <- m2_lin_treat_s1_fit$loo()
  saveRDS(m2_lin_treat_s1_loo, here::here("data/model_fits/treatment_s1/m2_lin_treat_s1_loo.RDS"))
  # Parameter estimates
  m2_lin_treat_s1_params <- get_params(subj_id = unique(task_data_treat_s1$subj_id), 
                                        model_fit = m2_lin_treat_s1_fit, 
                                        n_subj = length(unique(task_data_treat_s1$subj_id)), 
                                        n_params = 2, 
                                        param_names = c("kE", "kR"))
  saveRDS(m2_lin_treat_s1_params, here::here("data/model_fits/treatment_s1/m2_lin_treat_s1_params.RDS"))
  
  # Model 3
  m3_lin_treat_s1_fit <- m3_linear_stan_model$sample(
    data = model_dat_treat_s1, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_lin_treat_s1_check <- convergence_check(m3_lin_treat_s1_fit, 
                                              params = c("kE", "kR", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m3_lin_treat_s1_check$trace_plot
  saveRDS(list(m3_lin_treat_s1_check$Rhat, m3_lin_treat_s1_check$ess), 
          here::here("data/model_fits/treatment_s1/m3_lin_treat_s1_check.RDS"))
  # LOO for model comparisons
  m3_lin_treat_s1_loo <- m3_lin_treat_s1_fit$loo()
  saveRDS(m3_lin_treat_s1_loo, here::here("data/model_fits/treatment_s1/m3_lin_treat_s1_loo.RDS"))
  # Parameter estimates
  m3_lin_treat_s1_params <- get_params(subj_id = unique(task_data_treat_s1$subj_id), 
                                        model_fit = m3_lin_treat_s1_fit, 
                                        n_subj = length(unique(task_data_treat_s1$subj_id)), 
                                        n_params = 3, 
                                        param_names = c("kE", "kR", "a"))
  saveRDS(m3_lin_treat_s1_params, here::here("data/model_fits/treatment_s1/m3_lin_treat_s1_params.RDS"))
  
  ### Treatment group - Session 2 -----  
  ## Parabolic discounting models
  # Model 1
  m1_para_treat_s2_fit <- m1_parabolic_stan_model$sample(
    data = model_dat_treat_s2, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_para_treat_s2_check <- convergence_check(m1_para_treat_s2_fit, 
                                              params = c("kE", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m1_para_treat_s2_check$trace_plot
  saveRDS(list(m1_para_treat_s2_check$Rhat, m1_para_treat_s2_check$ess), 
          here::here("data/model_fits/treatment_s2/m1_para_treat_s2_check.RDS"))
  # LOO for model comparisons
  m1_para_treat_s2_loo <- m1_para_treat_s2_fit$loo()
  saveRDS(m1_para_treat_s2_loo, here::here("data/model_fits/treatment_s2/m1_para_treat_s2_loo.RDS"))
  # Parameter estimates
  m1_para_treat_s2_params <- get_params(subj_id = unique(task_data_treat_s2$subj_id), 
                                        model_fit = m1_para_treat_s2_fit, 
                                        n_subj = length(unique(task_data_treat_s2$subj_id)), 
                                        n_params = 2, 
                                        param_names = c("kE", "a"))
  saveRDS(m1_para_treat_s2_params, here::here("data/model_fits/treatment_s2/m1_para_treat_s2_params.RDS"))
  
  # Model 2
  m2_para_treat_s2_fit <- m2_parabolic_stan_model$sample(
    data = model_dat_treat_s2, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_para_treat_s2_check <- convergence_check(m2_para_treat_s2_fit, 
                                              params = c("kE", "kR"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m2_para_treat_s2_check$trace_plot
  saveRDS(list(m2_para_treat_s2_check$Rhat, m2_para_treat_s2_check$ess), 
          here::here("data/model_fits/treatment_s2/m2_para_treat_s2_check.RDS"))
  # LOO for model comparisons
  m2_para_treat_s2_loo <- m2_para_treat_s2_fit$loo()
  saveRDS(m2_para_treat_s2_loo, here::here("data/model_fits/treatment_s2/m2_para_treat_s2_loo.RDS"))
  # Parameter estimates
  m2_para_treat_s2_params <- get_params(subj_id = unique(task_data_treat_s2$subj_id), 
                                        model_fit = m2_para_treat_s2_fit, 
                                        n_subj = length(unique(task_data_treat_s2$subj_id)), 
                                        n_params = 2, 
                                        param_names = c("kE", "kR"))
  saveRDS(m2_para_treat_s2_params, here::here("data/model_fits/treatment_s2/m2_para_treat_s2_params.RDS"))
  
  # Model 3
  m3_para_treat_s2_fit <- m3_parabolic_stan_model$sample(
    data = model_dat_treat_s2, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_para_treat_s2_check <- convergence_check(m3_para_treat_s2_fit, 
                                              params = c("kE", "kR", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m3_para_treat_s2_check$trace_plot
  saveRDS(list(m3_para_treat_s2_check$Rhat, m3_para_treat_s2_check$ess), 
          here::here("data/model_fits/treatment_s2/m3_para_treat_s2_check.RDS"))
  # LOO for model comparisons
  m3_para_treat_s2_loo <- m3_para_treat_s2_fit$loo()
  saveRDS(m3_para_treat_s2_loo, here::here("data/model_fits/treatment_s2/m3_para_treat_s2_loo.RDS"))
  # Posterior Predictions (for target model only)
  m3_para_treat_s2_PPC_dat <- posterior_predictions(csv_paths = m3_para_treat_s2_fit$output_files(),
                                                    n_chains = 4,
                                                    n_iter = (6000),
                                                    n_subj = 25, 
                                                    n_trials = 64,
                                                    real_dat = task_data_treat_s2) 
  saveRDS(m3_para_treat_s2_PPC_dat, here::here("data/model_fits/treatment_s2/m3_para_treat_s2_PPC.RDS"))
  # Parameter estimates
  m3_para_treat_s2_params <- get_params(subj_id = unique(task_data_treat_s2$subj_id), 
                                        model_fit = m3_para_treat_s2_fit, 
                                        n_subj = length(unique(task_data_treat_s2$subj_id)), 
                                        n_params = 3, 
                                        param_names = c("kE", "kR", "a"))
  saveRDS(m3_para_treat_s2_params, here::here("data/model_fits/treatment_s2/m3_para_treat_s2_params.RDS"))
  
  ## Linear discounting models
  # Model 1
  m1_lin_treat_s2_fit <- m1_linear_stan_model$sample(
    data = model_dat_treat_s2, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_lin_treat_s2_check <- convergence_check(m1_lin_treat_s2_fit, 
                                                  params = c("kE", "a"), 
                                                  Rhat = TRUE, ess = TRUE,
                                                  trace_plot = TRUE, rank_hist = FALSE)
  m1_lin_treat_s2_check$trace_plot
  saveRDS(list(m1_lin_treat_s2_check$Rhat, m1_lin_treat_s2_check$ess), 
          here::here("data/model_fits/treatment_s2/m1_lin_treat_s2_check.RDS"))
  # LOO for model comparisons
  m1_lin_treat_s2_loo <- m1_lin_treat_s2_fit$loo()
  saveRDS(m1_lin_treat_s2_loo, here::here("data/model_fits/treatment_s2/m1_lin_treat_s2_loo.RDS"))
  # Parameter estimates
  m1_lin_treat_s2_params <- get_params(subj_id = unique(task_data_treat_s2$subj_id), 
                                            model_fit = m1_lin_treat_s2_fit, 
                                            n_subj = length(unique(task_data_treat_s2$subj_id)), 
                                            n_params = 2, 
                                            param_names = c("kE", "a"))
  saveRDS(m1_lin_treat_s2_params, here::here("data/model_fits/treatment_s2/m1_lin_treat_s2_params.RDS"))
  
  # Model 2
  m2_lin_treat_s2_fit <- m2_linear_stan_model$sample(
    data = model_dat_treat_s2, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_lin_treat_s2_check <- convergence_check(m2_lin_treat_s2_fit, 
                                                  params = c("kE", "kR"), 
                                                  Rhat = TRUE, ess = TRUE,
                                                  trace_plot = TRUE, rank_hist = FALSE)
  m2_lin_treat_s2_check$trace_plot
  saveRDS(list(m2_lin_treat_s2_check$Rhat, m2_lin_treat_s2_check$ess), 
          here::here("data/model_fits/treatment_s2/m2_lin_treat_s2_check.RDS"))
  # LOO for model comparisons
  m2_lin_treat_s2_loo <- m2_lin_treat_s2_fit$loo()
  saveRDS(m2_lin_treat_s2_loo, here::here("data/model_fits/treatment_s2/m2_lin_treat_s2_loo.RDS"))
  # Parameter estimates
  m2_lin_treat_s2_params <- get_params(subj_id = unique(task_data_treat_s2$subj_id), 
                                            model_fit = m2_lin_treat_s2_fit, 
                                            n_subj = length(unique(task_data_treat_s2$subj_id)), 
                                            n_params = 2, 
                                            param_names = c("kE", "kR"))
  saveRDS(m2_lin_treat_s2_params, here::here("data/model_fits/treatment_s2/m2_lin_treat_s2_params.RDS"))
  
  # Model 3
  m3_lin_treat_s2_fit <- m3_linear_stan_model$sample(
    data = model_dat_treat_s2, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_lin_treat_s2_check <- convergence_check(m3_lin_treat_s2_fit, 
                                                  params = c("kE", "kR", "a"), 
                                                  Rhat = TRUE, ess = TRUE,
                                                  trace_plot = TRUE, rank_hist = FALSE)
  m3_lin_treat_s2_check$trace_plot
  saveRDS(list(m3_lin_treat_s2_check$Rhat, m3_lin_treat_s2_check$ess), 
          here::here("data/model_fits/treatment_s2/m3_lin_treat_s2_check.RDS"))
  # LOO for model comparisons
  m3_lin_treat_s2_loo <- m3_lin_treat_s2_fit$loo()
  saveRDS(m3_lin_treat_s2_loo, here::here("data/model_fits/treatment_s2/m3_lin_treat_s2_loo.RDS"))
  # Parameter estimates
  m3_lin_treat_s2_params <- get_params(subj_id = unique(task_data_treat_s2$subj_id), 
                                            model_fit = m3_lin_treat_s2_fit, 
                                            n_subj = length(unique(task_data_treat_s2$subj_id)), 
                                            n_params = 3, 
                                            param_names = c("kE", "kR", "a"))
  saveRDS(m3_lin_treat_s2_params, here::here("data/model_fits/treatment_s2/m3_lin_treat_s2_params.RDS"))
  
  ### Control group -----  
  ## Parabolic discounting models
  # Model 1
  m1_para_control_fit <- m1_parabolic_stan_model$sample(
    data = model_dat_control, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_para_control_check <- convergence_check(m1_para_control_fit, 
                                              params = c("kE", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m1_para_control_check$trace_plot
  saveRDS(list(m1_para_control_check$Rhat, m1_para_control_check$ess), 
          here::here("data/model_fits/controls/m1_para_control_check.RDS"))
  # LOO for model comparisons
  m1_para_control_loo <- m1_para_control_fit$loo()
  saveRDS(m1_para_control_loo, here::here("data/model_fits/controls/m1_para_control_loo.RDS"))
  # Parameter estimates
  m1_para_control_params <- get_params(subj_id = unique(task_data_control$subj_id), 
                                        model_fit = m1_para_control_fit, 
                                        n_subj = length(unique(task_data_control$subj_id)), 
                                        n_params = 2, 
                                        param_names = c("kE", "a"))
  saveRDS(m1_para_control_params, here::here("data/model_fits/controls/m1_para_control_params.RDS"))
  
  # Model 2
  m2_para_control_fit <- m2_parabolic_stan_model$sample(
    data = model_dat_control, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_para_control_check <- convergence_check(m2_para_control_fit, 
                                              params = c("kE", "kR"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m2_para_control_check$trace_plot
  saveRDS(list(m2_para_control_check$Rhat, m2_para_control_check$ess), 
          here::here("data/model_fits/controls/m2_para_control_check.RDS"))
  # LOO for model comparisons
  m2_para_control_loo <- m2_para_control_fit$loo()
  saveRDS(m2_para_control_loo, here::here("data/model_fits/controls/m2_para_control_loo.RDS"))
  # Parameter estimates
  m2_para_control_params <- get_params(subj_id = unique(task_data_control$subj_id), 
                                        model_fit = m2_para_control_fit, 
                                        n_subj = length(unique(task_data_control$subj_id)), 
                                        n_params = 2, 
                                        param_names = c("kE", "kR"))
  saveRDS(m2_para_control_params, here::here("data/model_fits/controls/m2_para_control_params.RDS"))
  
  # Model 3
  m3_para_control_fit <- m3_parabolic_stan_model$sample(
    data = model_dat_control, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_para_control_check <- convergence_check(m3_para_control_fit, 
                                              params = c("kE", "kR", "a"), 
                                              Rhat = TRUE, ess = TRUE,
                                              trace_plot = TRUE, rank_hist = FALSE)
  m3_para_control_check$trace_plot
  saveRDS(list(m3_para_control_check$Rhat, m3_para_control_check$ess), 
          here::here("data/model_fits/controls/m3_para_control_check.RDS"))
  # LOO for model comparisons
  m3_para_control_loo <- m3_para_control_fit$loo()
  saveRDS(m3_para_control_loo, here::here("data/model_fits/controls/m3_para_control_loo.RDS"))
  # Posterior Predictions (for target model only)
  m3_para_control_PPC_dat <- posterior_predictions(csv_paths = m3_para_control_fit$output_files(),
                                                    n_chains = 4,
                                                    n_iter = (6000),
                                                    n_subj = 54, 
                                                    n_trials = 64,
                                                    real_dat = task_data_control) 
  saveRDS(m3_para_control_PPC_dat, here::here("data/model_fits/controls/m3_para_controls_PPC.RDS"))
  # Parameter estimates
  m3_para_control_params <- get_params(subj_id = unique(task_data_control$subj_id), 
                                        model_fit = m3_para_control_fit, 
                                        n_subj = length(unique(task_data_control$subj_id)), 
                                        n_params = 3, 
                                        param_names = c("kE", "kR", "a"))
  saveRDS(m3_para_control_params, here::here("data/model_fits/controls/m3_para_control_params.RDS"))

  ## Linear discounting models
  # Model 1
  m1_lin_control_fit <- m1_linear_stan_model$sample(
    data = model_dat_control, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_lin_control_check <- convergence_check(m1_lin_control_fit, 
                                             params = c("kE", "a"), 
                                             Rhat = TRUE, ess = TRUE,
                                             trace_plot = TRUE, rank_hist = FALSE)
  m1_lin_control_check$trace_plot
  saveRDS(list(m1_lin_control_check$Rhat, m1_lin_control_check$ess), 
          here::here("data/model_fits/controls/m1_lin_control_check.RDS"))
  # LOO for model comparisons
  m1_lin_control_loo <- m1_lin_control_fit$loo()
  saveRDS(m1_lin_control_loo, here::here("data/model_fits/controls/m1_lin_control_loo.RDS"))
  # Parameter estimates
  m1_lin_control_params <- get_params(subj_id = unique(task_data_control$subj_id), 
                                       model_fit = m1_lin_control_fit, 
                                       n_subj = length(unique(task_data_control$subj_id)), 
                                       n_params = 2, 
                                       param_names = c("kE", "a"))
  saveRDS(m1_lin_control_params, here::here("data/model_fits/controls/m1_lin_control_params.RDS"))
  
  # Model 2
  m2_lin_control_fit <- m2_linear_stan_model$sample(
    data = model_dat_control, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_lin_control_check <- convergence_check(m2_lin_control_fit, 
                                             params = c("kE", "kR"), 
                                             Rhat = TRUE, ess = TRUE,
                                             trace_plot = TRUE, rank_hist = FALSE)
  m2_lin_control_check$trace_plot
  saveRDS(list(m2_lin_control_check$Rhat, m2_lin_control_check$ess), 
          here::here("data/model_fits/controls/m2_lin_control_check.RDS"))
  # LOO for model comparisons
  m2_lin_control_loo <- m2_lin_control_fit$loo()
  saveRDS(m2_lin_control_loo, here::here("data/model_fits/controls/m2_lin_control_loo.RDS"))
  # Parameter estimates
  m2_lin_control_params <- get_params(subj_id = unique(task_data_control$subj_id), 
                                       model_fit = m2_lin_control_fit, 
                                       n_subj = length(unique(task_data_control$subj_id)), 
                                       n_params = 2, 
                                       param_names = c("kE", "kR"))
  saveRDS(m2_lin_control_params, here::here("data/model_fits/controls/m2_lin_control_params.RDS"))
  
  # Model 3
  m3_lin_control_fit <- m3_linear_stan_model$sample(
    data = model_dat_control, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_lin_control_check <- convergence_check(m3_lin_control_fit, 
                                             params = c("kE", "kR", "a"), 
                                             Rhat = TRUE, ess = TRUE,
                                             trace_plot = TRUE, rank_hist = FALSE)
  m3_lin_control_check$trace_plot
  saveRDS(list(m3_lin_control_check$Rhat, m3_lin_control_check$ess), 
          here::here("data/model_fits/controls/m3_lin_control_check.RDS"))
  # LOO for model comparisons
  m3_lin_control_loo <- m3_lin_control_fit$loo()
  saveRDS(m3_lin_control_loo, here::here("data/model_fits/controls/m3_lin_control_loo.RDS"))
  # Parameter estimates
  m3_lin_control_params <- get_params(subj_id = unique(task_data_control$subj_id), 
                                       model_fit = m3_lin_control_fit, 
                                       n_subj = length(unique(task_data_control$subj_id)), 
                                       n_params = 3, 
                                       param_names = c("kE", "kR", "a"))
  saveRDS(m3_lin_control_params, here::here("data/model_fits/controls/m3_lin_control_params.RDS"))
  
  ### Non-diabetic group -----  
  ## Parabolic discounting models
  # Model 1
  m1_para_non_diabetics_fit <- m1_parabolic_stan_model$sample(
    data = model_dat_non_diabetics, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_para_non_diabetics_check <- convergence_check(m1_para_non_diabetics_fit, 
                                             params = c("kE", "a"), 
                                             Rhat = TRUE, ess = TRUE,
                                             trace_plot = TRUE, rank_hist = FALSE)
  m1_para_non_diabetics_check$trace_plot
  saveRDS(list(m1_para_non_diabetics_check$Rhat, m1_para_non_diabetics_check$ess), 
          here::here("data/model_fits/non_diabetics/m1_para_non_diabetics_check.RDS"))
  # LOO for model comparisons
  m1_para_non_diabetics_loo <- m1_para_non_diabetics_fit$loo()
  saveRDS(m1_para_non_diabetics_loo, here::here("data/model_fits/non_diabetics/m1_para_non_diabetics_loo.RDS"))
  # Parameter estimates
  m1_para_non_diabetics_params <- get_params(subj_id = unique(task_data_non_diabetics$subj_id), 
                                       model_fit = m1_para_non_diabetics_fit, 
                                       n_subj = length(unique(task_data_non_diabetics$subj_id)), 
                                       n_params = 2, 
                                       param_names = c("kE", "a"))
  saveRDS(m1_para_non_diabetics_params, here::here("data/model_fits/non_diabetics/m1_para_non_diabetics_params.RDS"))
  
  # Model 2
  m2_para_non_diabetics_fit <- m2_parabolic_stan_model$sample(
    data = model_dat_non_diabetics, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_para_non_diabetics_check <- convergence_check(m2_para_non_diabetics_fit, 
                                             params = c("kE", "kR"), 
                                             Rhat = TRUE, ess = TRUE,
                                             trace_plot = TRUE, rank_hist = FALSE)
  m2_para_non_diabetics_check$trace_plot
  saveRDS(list(m2_para_non_diabetics_check$Rhat, m2_para_non_diabetics_check$ess), 
          here::here("data/model_fits/non_diabetics/m2_para_non_diabetics_check.RDS"))
  # LOO for model comparisons
  m2_para_non_diabetics_loo <- m2_para_non_diabetics_fit$loo()
  saveRDS(m2_para_non_diabetics_loo, here::here("data/model_fits/non_diabetics/m2_para_non_diabetics_loo.RDS"))
  # Parameter estimates
  m2_para_non_diabetics_params <- get_params(subj_id = unique(task_data_non_diabetics$subj_id), 
                                       model_fit = m2_para_non_diabetics_fit, 
                                       n_subj = length(unique(task_data_non_diabetics$subj_id)), 
                                       n_params = 2, 
                                       param_names = c("kE", "kR"))
  saveRDS(m2_para_non_diabetics_params, here::here("data/model_fits/non_diabetics/m2_para_non_diabetics_params.RDS"))
  
  # Model 3
  m3_para_non_diabetics_fit <- m3_parabolic_stan_model$sample(
    data = model_dat_non_diabetics, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_para_non_diabetics_check <- convergence_check(m3_para_non_diabetics_fit, 
                                             params = c("kE", "kR", "a"), 
                                             Rhat = TRUE, ess = TRUE,
                                             trace_plot = TRUE, rank_hist = FALSE)
  m3_para_non_diabetics_check$trace_plot
  saveRDS(list(m3_para_non_diabetics_check$Rhat, m3_para_non_diabetics_check$ess), 
          here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_check.RDS"))
  # LOO for model comparisons
  m3_para_non_diabetics_loo <- m3_para_non_diabetics_fit$loo()
  saveRDS(m3_para_non_diabetics_loo, here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_loo.RDS"))
  # Posterior Predictions (for target model only)
  m3_para_non_diabetics_PPC_dat <- posterior_predictions(csv_paths = m3_para_non_diabetics_fit$output_files(),
                                                   n_chains = 4,
                                                   n_iter = (6000),
                                                   n_subj = 59, 
                                                   n_trials = 64,
                                                   real_dat = task_data_non_diabetics) 
  saveRDS(m3_para_non_diabetics_PPC_dat, here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_PPC.RDS"))
  # Parameter estimates
  m3_para_non_diabetics_params <- get_params(subj_id = unique(task_data_non_diabetics$subj_id), 
                                       model_fit = m3_para_non_diabetics_fit, 
                                       n_subj = length(unique(task_data_non_diabetics$subj_id)), 
                                       n_params = 3, 
                                       param_names = c("kE", "kR", "a"))
  saveRDS(m3_para_non_diabetics_params, here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_params.RDS"))
  
  ## Linear discounting models
  # Model 1
  m1_lin_non_diabetics_fit <- m1_linear_stan_model$sample(
    data = model_dat_non_diabetics, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m1_lin_non_diabetics_check <- convergence_check(m1_lin_non_diabetics_fit, 
                                            params = c("kE", "a"), 
                                            Rhat = TRUE, ess = TRUE,
                                            trace_plot = TRUE, rank_hist = FALSE)
  m1_lin_non_diabetics_check$trace_plot
  saveRDS(list(m1_lin_non_diabetics_check$Rhat, m1_lin_non_diabetics_check$ess), 
          here::here("data/model_fits/non_diabetics/m1_lin_non_diabetics_check.RDS"))
  # LOO for model comparisons
  m1_lin_non_diabetics_loo <- m1_lin_non_diabetics_fit$loo()
  saveRDS(m1_lin_non_diabetics_loo, here::here("data/model_fits/non_diabetics/m1_lin_non_diabetics_loo.RDS"))
  # Parameter estimates
  m1_lin_non_diabetics_params <- get_params(subj_id = unique(task_data_non_diabetics$subj_id), 
                                      model_fit = m1_lin_non_diabetics_fit, 
                                      n_subj = length(unique(task_data_non_diabetics$subj_id)), 
                                      n_params = 2, 
                                      param_names = c("kE", "a"))
  saveRDS(m1_lin_non_diabetics_params, here::here("data/model_fits/non_diabetics/m1_lin_non_diabetics_params.RDS"))
  
  # Model 2
  m2_lin_non_diabetics_fit <- m2_linear_stan_model$sample(
    data = model_dat_non_diabetics, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m2_lin_non_diabetics_check <- convergence_check(m2_lin_non_diabetics_fit, 
                                            params = c("kE", "kR"), 
                                            Rhat = TRUE, ess = TRUE,
                                            trace_plot = TRUE, rank_hist = FALSE)
  m2_lin_non_diabetics_check$trace_plot
  saveRDS(list(m2_lin_non_diabetics_check$Rhat, m2_lin_non_diabetics_check$ess), 
          here::here("data/model_fits/non_diabetics/m2_lin_non_diabetics_check.RDS"))
  # LOO for model comparisons
  m2_lin_non_diabetics_loo <- m2_lin_non_diabetics_fit$loo()
  saveRDS(m2_lin_non_diabetics_loo, here::here("data/model_fits/non_diabetics/m2_lin_non_diabetics_loo.RDS"))
  # Parameter estimates
  m2_lin_non_diabetics_params <- get_params(subj_id = unique(task_data_non_diabetics$subj_id), 
                                      model_fit = m2_lin_non_diabetics_fit, 
                                      n_subj = length(unique(task_data_non_diabetics$subj_id)), 
                                      n_params = 2, 
                                      param_names = c("kE", "kR"))
  saveRDS(m2_lin_non_diabetics_params, here::here("data/model_fits/non_diabetics/m2_lin_non_diabetics_params.RDS"))
  
  # Model 3
  m3_lin_non_diabetics_fit <- m3_linear_stan_model$sample(
    data = model_dat_non_diabetics, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_lin_non_diabetics_check <- convergence_check(m3_lin_non_diabetics_fit, 
                                            params = c("kE", "kR", "a"), 
                                            Rhat = TRUE, ess = TRUE,
                                            trace_plot = TRUE, rank_hist = FALSE)
  m3_lin_non_diabetics_check$trace_plot
  saveRDS(list(m3_lin_non_diabetics_check$Rhat, m3_lin_non_diabetics_check$ess), 
          here::here("data/model_fits/non_diabetics/m3_lin_non_diabetics_check.RDS"))
  # LOO for model comparisons
  m3_lin_non_diabetics_loo <- m3_lin_non_diabetics_fit$loo()
  saveRDS(m3_lin_non_diabetics_loo, here::here("data/model_fits/non_diabetics/m3_lin_non_diabetics_loo.RDS"))
  # Parameter estimates
  m3_lin_non_diabetics_params <- get_params(subj_id = unique(task_data_non_diabetics$subj_id), 
                                      model_fit = m3_lin_non_diabetics_fit, 
                                      n_subj = length(unique(task_data_non_diabetics$subj_id)), 
                                      n_params = 3, 
                                      param_names = c("kE", "kR", "a"))
  saveRDS(m3_lin_non_diabetics_params, here::here("data/model_fits/non_diabetics/m3_lin_non_diabetics_params.RDS"))
  
  ### Non-diabetic normal weight group -----  
  ## Parabolic discounting models
  # Model 3
  m3_para_normal_weight_fit <- m3_parabolic_stan_model$sample(
    data = model_dat_normal_weight, 
    refresh = 0, chains = 4, parallel_chains = 4, 
    iter_warmup = 2000, iter_sampling = 6000, 
    adapt_delta = 0.8, step_size = 1, max_treedepth = 10, save_warmup = TRUE, 
    output_dir = NULL
  )
  # Convergence check
  m3_para_normal_weight_check <- convergence_check(m3_para_normal_weight_fit, 
                                                   params = c("kE", "kR", "a"), 
                                                   Rhat = TRUE, ess = TRUE,
                                                   trace_plot = TRUE, rank_hist = FALSE)
  m3_para_normal_weight_check$trace_plot
  saveRDS(list(m3_para_normal_weight_check$Rhat, m3_para_normal_weight_check$ess), 
          here::here("data/model_fits/low_bmi/m3_para_normal_weight_check.RDS"))
  # LOO for model comparisons
  m3_para_normal_weight_loo <- m3_para_normal_weight_fit$loo()
  saveRDS(m3_para_normal_weight_loo, here::here("data/model_fits/low_bmi/m3_para_normal_weight_loo.RDS"))
  # Posterior Predictions (for target model only)
  m3_para_normal_weight_PPC_dat <- posterior_predictions(csv_paths = m3_para_normal_weight_fit$output_files(),
                                                         n_chains = 4,
                                                         n_iter = (6000),
                                                         n_subj = 59, 
                                                         n_trials = 64,
                                                         real_dat = task_data_normal_weight) 
  saveRDS(m3_para_normal_weight_PPC_dat, here::here("data/model_fits/low_bmi/m3_para_normal_weight_PPC_dat.RDS"))
  # Parameter estimates
  m3_para_normal_weight_params <- get_params(subj_id = unique(task_data_normal_weight$subj_id), 
                                             model_fit = m3_para_normal_weight_fit, 
                                             n_subj = length(unique(task_data_normal_weight$subj_id)), 
                                             n_params = 3, 
                                             param_names = c("kE", "kR", "a"))
  saveRDS(m3_para_normal_weight_params, here::here("data/model_fits/low_bmi/m3_para_normal_weight_params.RDS"))
  

} else {
  ### Read in models if already fitted
  ## Treatment group - Session 1
  # Parabolic
  m1_para_treat_s1_loo <- readRDS(here::here("data/model_fits/treatment_s1/m1_para_treat_s1_loo.RDS"))
  m1_para_treat_s1_check <- readRDS(here::here("data/model_fits/treatment_s1/m1_para_treat_s1_check.RDS"))
  m1_para_treat_s1_params <- readRDS(here::here("data/model_fits/treatment_s1/m1_para_treat_s1_params.RDS"))
  m2_para_treat_s1_loo <- readRDS(here::here("data/model_fits/treatment_s1/m2_para_treat_s1_loo.RDS"))
  m2_para_treat_s1_check <- readRDS(here::here("data/model_fits/treatment_s1/m2_para_treat_s1_check.RDS"))
  m2_para_treat_s1_params <- readRDS(here::here("data/model_fits/treatment_s1/m2_para_treat_s1_params.RDS"))
  m3_para_treat_s1_loo <- readRDS(here::here("data/model_fits/treatment_s1/m3_para_treat_s1_loo.RDS"))
  m3_para_treat_s1_check <- readRDS(here::here("data/model_fits/treatment_s1/m3_para_treat_s1_check.RDS"))
  m3_para_treat_s1_params <- readRDS(here::here("data/model_fits/treatment_s1/m3_para_treat_s1_params.RDS"))
  
  ## Treatment group - Session 2
  # Parabolic
  m1_para_treat_s2_loo <- readRDS(here::here("data/model_fits/treatment_s2/m1_para_treat_s2_loo.RDS"))
  m1_para_treat_s2_check <- readRDS(here::here("data/model_fits/treatment_s2/m1_para_treat_s2_check.RDS"))
  m1_para_treat_s2_params <- readRDS(here::here("data/model_fits/treatment_s2/m1_para_treat_s2_params.RDS"))
  m2_para_treat_s2_loo <- readRDS(here::here("data/model_fits/treatment_s2/m2_para_treat_s2_loo.RDS"))
  m2_para_treat_s2_check <- readRDS(here::here("data/model_fits/treatment_s2/m2_para_treat_s2_check.RDS"))
  m2_para_treat_s2_params <- readRDS(here::here("data/model_fits/treatment_s2/m2_para_treat_s2_params.RDS"))
  m3_para_treat_s2_loo <- readRDS(here::here("data/model_fits/treatment_s2/m3_para_treat_s2_loo.RDS"))
  m3_para_treat_s2_check <- readRDS(here::here("data/model_fits/treatment_s2/m3_para_treat_s2_check.RDS"))
  m3_para_treat_s2_params <- readRDS(here::here("data/model_fits/treatment_s2/m3_para_treat_s2_params.RDS"))
  
  ## Control group 
  # Parabolic
  m1_para_control_loo <- readRDS(here::here("data/model_fits/controls/m1_para_control_loo.RDS"))
  m1_para_control_check <- readRDS(here::here("data/model_fits/controls/m1_para_control_check.RDS"))
  m1_para_control_params <- readRDS(here::here("data/model_fits/controls/m1_para_control_params.RDS"))
  m2_para_control_loo <- readRDS(here::here("data/model_fits/controls/m2_para_control_loo.RDS"))
  m2_para_control_check <- readRDS(here::here("data/model_fits/controls/m2_para_control_check.RDS"))
  m2_para_control_params <- readRDS(here::here("data/model_fits/controls/m2_para_control_params.RDS"))
  m3_para_control_loo <- readRDS(here::here("data/model_fits/controls/m3_para_control_loo.RDS"))
  m3_para_control_check <- readRDS(here::here("data/model_fits/controls/m3_para_control_check.RDS"))
  m3_para_control_params <- readRDS(here::here("data/model_fits/controls/m3_para_control_params.RDS"))
  
  ## Non diabetic group 
  # Parabolic
  m1_para_non_diabetics_loo <- readRDS(here::here("data/model_fits/non_diabetics/m1_para_non_diabetics_loo.RDS"))
  m1_para_non_diabetics_check <- readRDS(here::here("data/model_fits/non_diabetics/m1_para_non_diabetics_check.RDS"))
  m1_para_non_diabetics_params <- readRDS(here::here("data/model_fits/non_diabetics/m1_para_non_diabetics_params.RDS"))
  m2_para_non_diabetics_loo <- readRDS(here::here("data/model_fits/non_diabetics/m2_para_non_diabetics_loo.RDS"))
  m2_para_non_diabetics_check <- readRDS(here::here("data/model_fits/non_diabetics/m2_para_non_diabetics_check.RDS"))
  m2_para_non_diabetics_params <- readRDS(here::here("data/model_fits/non_diabetics/m2_para_non_diabetics_params.RDS"))
  m3_para_non_diabetics_loo <- readRDS(here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_loo.RDS"))
  m3_para_non_diabetics_check <- readRDS(here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_check.RDS"))
  m3_para_non_diabetics_params <- readRDS(here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_params.RDS"))
}

### (3) Convergence check -----------------------------------------------

# Rhats maximum
max(m1_para_treat_s1_check[[1]][6], m2_para_treat_s1_check[[1]][6], m3_para_treat_s1_check[[1]][6],
    m1_para_treat_s2_check[[1]][6], m2_para_treat_s2_check[[1]][6], m3_para_treat_s2_check[[1]][6],
    m1_para_control_check[[1]][6], m2_para_control_check[[1]][6], m3_para_control_check[[1]][6],
    m1_para_non_diabetics_check[[1]][6], m2_para_non_diabetics_check[[1]][6], m3_para_non_diabetics_check[[1]][6])

# ESS minimum
min(m1_para_treat_s1_check[[2]][1], m2_para_treat_s1_check[[2]][1], m3_para_treat_s1_check[[2]][1],
    m1_para_treat_s2_check[[2]][1], m2_para_treat_s2_check[[2]][1], m3_para_treat_s2_check[[2]][1],
    m1_para_control_check[[2]][1], m2_para_control_check[[2]][1], m3_para_control_check[[2]][1],
    m1_para_non_diabetics_check[[2]][1], m2_para_non_diabetics_check[[2]][1], m3_para_non_diabetics_check[[2]][1])

### (4) Model comparison -----------------------------------------------

## Treatment group - Session 1
comparison_treat_s1 <- model_comparison(loo_paths = list.files(path = here::here("data/model_fits/treatment_s1"), 
                                                      pattern = "loo.RDS", full.names = TRUE), 
                               model_names = c("lin 1", "para 1", 
                                               "lin 2", "para 2", 
                                               "lin 3", "para 3"), 
                               LOO_ELPD_title = "")

## Treatment group - Session 2
comparison_treat_s2 <- model_comparison(loo_paths = list.files(path = here::here("data/model_fits/treatment_s2"), 
                                                               pattern = "loo.RDS", full.names = TRUE), 
                                        model_names = c("lin 1", "para 1", 
                                                        "lin 2", "para 2", 
                                                        "lin 3", "para 3"), 
                                        LOO_ELPD_title = "")

## Control group
comparison_controls <- model_comparison(loo_paths = list.files(path = here::here("data/model_fits/controls"), 
                                                               pattern = "loo.RDS", full.names = TRUE), 
                                        model_names = c("lin 1", "para 1", 
                                                        "lin 2", "para 2", 
                                                        "lin 3", "para 3"), 
                                        LOO_ELPD_title = "")

## Non-diabetics group 
comparison_non_diabetics <- model_comparison(loo_paths = list.files(path = here::here("data/model_fits/non_diabetics"), 
                                                               pattern = "loo.RDS", full.names = TRUE), 
                                        model_names = c("lin 1", "para 1", 
                                                        "lin 2", "para 2", 
                                                        "lin 3", "para 3"), 
                                        LOO_ELPD_title = "")

pdf(file = here::here("output/figures/R_plots/model_comparison_main.pdf"),  
    width = 10, # The width of the plot in cm (transformed to inches)
    height = 4) # The height of the plot in cm (transformed to inches)
par(mar=c(0,4,0.5,0.5))


ggarrange(
  comparison_treat_s1$LOO_ELPD_plot +
    scale_x_discrete(labels = 
                       c(expression("Linear" ~beta[E]~alpha), expression("Linear" ~beta[E]~beta[R]), expression("Linear" ~beta[E]~beta[R]~alpha), 
                         expression("Parabolic" ~beta[E]~alpha), expression("Parabolic" ~beta[E]~beta[R]), expression("Parabolic" ~beta[E]~beta[R]~alpha))),
  comparison_controls$LOO_ELPD_plot +
    scale_x_discrete(labels = 
                       c(expression("Linear" ~beta[E]~alpha), expression("Linear" ~beta[E]~beta[R]), expression("Linear" ~beta[E]~beta[R]~alpha), 
                         expression("Parabolic" ~beta[E]~alpha), expression("Parabolic" ~beta[E]~beta[R]), expression("Parabolic" ~beta[E]~beta[R]~alpha))),  
  ncol = 2, nrow = 1, 
  labels = c("Treatment group", "Control group"), 
  hjust = c(-1.25, -1.6))

dev.off()

### (5) Posterior predictive checks -----------------------------------------------

## Treatment group - Session 1
m3_para_treat_s1_PPC_dat <- readRDS(here::here("data/model_fits/treatment_s1/m3_para_treat_s1_PPC.RDS")) 
# Plot
ppc_treat_s1_plot <- ppc_plots(m3_para_treat_s1_PPC_dat, indiv_plot_title = "")

## Treatment group - Session 2
m3_para_treat_s2_PPC_dat <- readRDS(here::here("data/model_fits/treatment_s2/m3_para_treat_s2_PPC.RDS")) 
# Plot
ppc_treat_s2_plot <- ppc_plots(m3_para_treat_s2_PPC_dat, indiv_plot_title = "")

## Control group 
m3_para_controls_PPC_dat <- readRDS(here::here("data/model_fits/controls/m3_para_controls_PPC.RDS")) 
# Plot
ppc_controls_plot <- ppc_plots(m3_para_controls_PPC_dat, indiv_plot_title = "")

## Non-diabetic group
m3_para_non_diabetics_PPC_dat <- readRDS(here::here("data/model_fits/non_diabetics/m3_para_non_diabetics_PPC.RDS")) 
# Plot
ppc_non_diabetics_plot <- ppc_plots(m3_para_non_diabetics_PPC_dat, indiv_plot_title = "")


pdf(file = here::here("output/figures/R_plots/ppc_main.pdf"),  
    width = 10, # The width of the plot in cm (transformed to inches)
    height = 4) # The height of the plot in cm (transformed to inches)
par(mar=c(0,4,0.5,0.5))

ggarrange(ppc_treat_s1_plot$indiv_plot,
          ppc_controls_plot$indiv_plot,
          ncol = 2, nrow = 1, 
          labels = c("Treatment group", "Control group"), 
          hjust = c(-1.25, -1.6))


dev.off()


