#################################################################################################
###################---------------- Within subject comparison----------------####################
#################################################################################################

### In this script: 
# (1) Prepare data
# (2) Task data - model based

# Set working directory
here::i_am("github/semaglutide-study/code/analyses/7_within_subject_comparison.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/plot_funs.R")

# source dataset
main_data <- readRDS("data/processed_data/main_data.RDS")

# source parameter estimates
m3_para_treat_s1_params <- readRDS(here::here("data/model_fits/treatment_s1/m3_para_treat_s1_params.RDS"))

m3_para_treat_s2_params <- readRDS(here::here("data/model_fits/treatment_s2/m3_para_treat_s2_params.RDS"))

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, PupillometryR,
                 writexl, lubridate, magrittr, pushoverr, nlme, gridExtra, cmdstanr, rstanarm, bayestestR)

# Color pallet 
color_pal <- c("#E94D36", "#5B9BD5", "#71AB48", "#FDC219", "#8456B8", "#FF7236", "#1FD5B3", "#F781BE")

### (1) Prepare data  -----------------------------------------------

within_analyses_ids <- main_data$glp_data %>% filter(session == 2) %>% .$subj_id

# Merge datasets for analyses
data <- rbind(
  # Session 1
  main_data$demographic_data %>% 
  filter(group == "treatment") %>% 
  select(subj_id, group, age, gender, bmi) %>%
  left_join(main_data$glp_data %>% 
              filter(session == 1 & group == "treatment") %>% 
              select(subj_id, start_date_glp, side_effects_glp, glp_dose_mg, hours_since_injection), 
            by = "subj_id") %>% 
  left_join(main_data$questionnaire_data %>% 
              filter(session == 1 & group == "treatment") %>% 
              select(subj_id, daq_sumScore, shaps_sumScore, tfeq_cr_sumScore, tfeq_ue_sumScore, tfeq_ee_sumScore, 
                    last_meal_time, last_meal_size, snack, snack_time, hunger_rating, cgl, cgl_measure,cgl_unit),
            by = "subj_id") %>% 
  left_join(rbind(m3_para_treat_s1_params$individual_params %>% 
                    pivot_wider(id_cols = subj_id, names_from = parameter, 
                                values_from = c(estimate, hdi_lower, hdi_upper))), 
            by = "subj_id") %>% 
  filter(subj_id %in% within_analyses_ids)  %>% 
  add_column(session = 1, 
             .after = "subj_id"),
  # Session 2
  main_data$demographic_data %>% 
    filter(group == "treatment" & subj_id %in% within_analyses_ids) %>% 
    select(subj_id, group, age, gender, bmi) %>%
    left_join(main_data$glp_data %>% 
                filter(session == 2 & group == "treatment") %>% 
                select(subj_id, start_date_glp, side_effects_glp, glp_dose_mg, hours_since_injection), 
              by = "subj_id") %>% 
    left_join(main_data$questionnaire_data %>% 
                filter(session == 2 & group == "treatment") %>% 
                select(subj_id, daq_sumScore, shaps_sumScore, tfeq_cr_sumScore, tfeq_ue_sumScore, tfeq_ee_sumScore, 
                       last_meal_time, last_meal_size, snack, snack_time, hunger_rating, cgl, cgl_measure,cgl_unit),
              by = "subj_id") %>% 
    left_join(rbind(m3_para_treat_s2_params$individual_params %>% 
                      pivot_wider(id_cols = subj_id, names_from = parameter, 
                                  values_from = c(estimate, hdi_lower, hdi_upper))), 
              by = "subj_id") %>% 
    add_column(session = 2, 
               .after = "subj_id"))

# Parameters
# Scale parameters to be between 0 and 1
data %<>% 
  ungroup %>% 
  mutate(across(c(estimate_kE, estimate_kR, estimate_a), rescale)) 

# Visualize distributions
ggplot(gather(data %>% select(c(estimate_kE:estimate_a))) %>% na.omit(), 
       aes(value)) + 
  geom_histogram(bins = 15) + 
  facet_wrap(~key, scales = 'free_x')

data %<>% 
  # positively skewed
  mutate_at(c("estimate_a"), sqrt)

# Questionnaires
# Scale parameters to be between 0 and 1
data %<>% 
  ungroup %>% 
  mutate_at(colnames(data)[c(6, 11:15)], rescale)

# Visualize distributions
ggplot(gather(data %>% select(colnames(data)[c(6, 11:15)])) %>% na.omit(), 
       aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

data %<>% 
  # positively skewed
  mutate_at(c("bmi", "daq_sumScore", "shaps_sumScore", "tfeq_cr_sumScore"), sqrt) 

# Add group by time since injection
data %<>% 
  ungroup() %>%
  mutate(injection_group = case_when(hours_since_injection < hours(80) ~ 1, 
                                     hours_since_injection > hours(80) ~ 2), 
         .after = "subj_id")

### (2) Task data - model based -----------------------------------------------

data %<>% 
  mutate(injection_group = factor(injection_group)) 

### Effort sensitivity
t.test(estimate_kE ~ injection_group, data = data, paired = TRUE)

### Reward sensitivity
t.test(estimate_kR ~ injection_group, data = data, paired = TRUE)

### Choice bias
t.test(estimate_a ~ injection_group, data = data, paired = TRUE)




