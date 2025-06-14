#################################################################################################
#######################---------------- Primary analyses ----------------########################
#################################################################################################

### In this script: 
# (1) Prepare data
# (2) Task data - model based
# (3) Questionnaire data

# Set working directory
here::i_am("github/semaglutide-study/code/analyses/5_primary_analyses.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/plot_funs.R")

# source dataset
main_data <- readRDS("data/processed_data/main_data.RDS")

# source parameter estimates
m3_para_treat_s1_params <- readRDS(here::here("github/semaglutide-study/data/model_fits/treatment_s1/m3_para_treat_s1_params.RDS"))
m3_para_control_params <- readRDS(here::here("github/semaglutide-study/data/model_fits/controls/m3_para_control_params.RDS"))

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, PupillometryR,
                 writexl, lubridate, magrittr, pushoverr, nlme, gridExtra, cmdstanr, rstanarm, bayestestR, hms)

# Color pallet 
color_pal <- c("#E94D36", "#5B9BD5", "#71AB48", "#FDC219", "#8456B8", "#FF7236", "#1FD5B3", "#F781BE")

### (1) Prepare data  -----------------------------------------------

# Merge datasets for analyses
data <- main_data$demographic_data %>% 
  select(subj_id, group, age, gender, bmi, antidepressant) %>%
  mutate(antidepressant = ifelse(is.na(antidepressant), 0, antidepressant)) %>% 
  left_join(main_data$glp_data %>% 
              filter(session == 1) %>% 
              select(subj_id, start_date_glp, side_effects_glp, glp_dose_mg, hours_since_injection, testing_day, local_testing_time) %>% 
              mutate(time_on_glp = abs(difftime(as_date(start_date_glp),  testing_day, units="days"))), 
            by = "subj_id") %>% 
  left_join(main_data$questionnaire_data %>% 
              filter(session == 1) %>% 
              select(subj_id, aes_sumScore, bdi_sumScore, findrisc_sumScore, mcq_discounting_rate, mctq_MSF_SC, meq_sumScore, 
                     ocir_sumScore, daq_sumScore, shaps_sumScore, tfeq_cr_sumScore, tfeq_ue_sumScore, tfeq_ee_sumScore, 
                     ipaq_sumScore, last_meal_time, last_meal_size, snack, snack_time, hunger_rating, cgl, cgl_measure,cgl_unit, hunger_rating),
            by = "subj_id") %>% 
  left_join(main_data$task_meta_data %>% 
              filter(session == 1) %>% 
              select(subj_id, start_time)) %>% 
  left_join(rbind(m3_para_treat_s1_params$individual_params %>% 
                    pivot_wider(id_cols = subj_id, names_from = parameter, 
                                values_from = c(estimate, hdi_lower, hdi_upper)), 
                  m3_para_control_params$individual_params %>% 
                    pivot_wider(id_cols = subj_id, names_from = parameter, 
                                values_from = c(estimate, hdi_lower, hdi_upper))), 
            by = "subj_id") %>% 
  # make MCTQ result numeric (minutes since 00:00)
  add_column(mctq_continuous = period_to_seconds(hm(.$mctq_MSF_SC))/60, 
             .before = "mctq_MSF_SC") 

# Parameters
# Scale parameters to be between 0 and 1
data %<>% 
  ungroup %>% 
  mutate(across(c(estimate_kE, estimate_kR, estimate_a), rescale)) 

# Visualize distributions
ggplot(gather(data %>% select(c(estimate_kE:estimate_a))) %>% na.omit(), 
       aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

# Questionnaires
# Scale parameters to be between 0 and 1
data %<>% 
  ungroup %>% 
  mutate_at(colnames(data)[c(14:18, 21:27)], rescale)

# Visualize distributions
ggplot(gather(data %>% select(colnames(data)[c(5, 14:18, 21:27)])) %>% na.omit(), 
       aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

data %<>% 
  # positively skewed
  mutate_at(c("bdi_sumScore", "daq_sumScore", 
              "mcq_discounting_rate", "mctq_continuous",
              "ocir_sumScore", "shaps_sumScore"), sqrt) %>% 
  # negatively skewed
  mutate_at(c("aes_sumScore"), norm_neg_skew)

### (2) Task data - model based -----------------------------------------------

### Effort sensitivity
data %>% 
  group_by(group) %>% 
  summarise(mean_kE = mean(estimate_kE),
            sd_kE = sd(estimate_kE))
# Check if variances are equal
var.test(data$estimate_kE ~ data$group)
# => variances are equal -> t test
t.test(data$estimate_kE ~ data$group, var.equal = TRUE)
kE_glm <- stan_glm(estimate_kE ~ group, data = data, 
                  iter = 10000, seed = 123)
kE_glm$coefficients
hdi(kE_glm)
kE_plot <- raincloud_plot(dat = data, title = "", 
                          xlab = " ", ylab = "Effort sensitivity", 
                          predictor_var = "group", outcome_var = "estimate_kE", 
                          predictor_tick_lab = c("control", "treatment"), col = c(color_pal[1], color_pal[2]), 
                          include_grouping = FALSE, direction = "horizontal", scale_seq = c(-0, 1,0.25)) + 
  theme(legend.position = "none") +
  ggtitle("Effort sensitivity")

### Reward sensitivity
data %>% 
  group_by(group) %>% 
  summarise(mean_kR = mean(estimate_kR),
            sd_kR = sd(estimate_kR))
# Check if variances are equal
var.test(data$estimate_kR ~ data$group)
# => variances are not equal -> welch test
t.test(data$estimate_kR ~ data$group, var.equal = FALSE)
kR_glm <- stan_glm(estimate_kR ~ group, data = data, 
                   iter = 10000, seed = 123)
kR_glm$coefficients
hdi(kR_glm)
kR_plot <- raincloud_plot(dat = data, title = "", 
                         xlab = " ", ylab = "Reward sensitivity", 
                         predictor_var = "group", outcome_var = "estimate_kR", 
                         predictor_tick_lab = c("control", "treatment"), col = c(color_pal[1], color_pal[2]), 
                         include_grouping = FALSE, direction = "horizontal", scale_seq = c(-0, 1,0.25)) + 
  theme(legend.position = "none") +
  ggtitle("Reward sensitivity")

### Choice bias 
data %>% 
  group_by(group) %>% 
  summarise(mean_a = mean(estimate_a),
            sd_a = sd(estimate_a))
# Check if variances are equal
var.test(data$estimate_a ~ data$group)
# => variances are equal -> t test
t.test(data$estimate_a ~ data$group, var.equal = FALSE)
# Bayesian GLM
a_glm <- stan_glm(estimate_a ~ group, data = data, 
                        iter = 100000, seed = 123)
a_glm$coefficients
hdi(a_glm)
a_plot <- raincloud_plot(dat = data, title = "", 
                                xlab = " ", ylab = "Choice bias", 
                                predictor_var = "group", outcome_var = "estimate_a", 
                                predictor_tick_lab = c("control", "treatment"), col = c(color_pal[1], color_pal[2]), 
                                include_grouping = FALSE, direction = "horizontal", scale_seq = c(-0, 1, 0.25)) + 
  theme(legend.position = "none") +
  ggtitle("Choice bias")


### (3) Questionnaire data -----------------------------------------------

### Psychiatric questionnaires ----
# Correlation between questionnaires
psych_cor <- data %>% 
  select(bmi, aes_sumScore, bdi_sumScore, ocir_sumScore, daq_sumScore, shaps_sumScore) %>% 
  rename(aes = aes_sumScore, bdi = bdi_sumScore, ocir = ocir_sumScore, daq = daq_sumScore, shaps = shaps_sumScore) %>% 
  cor()
corrplot::corrplot(psych_cor, method="circle")

# Group comparison
# AES
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_aes = mean(aes_sumScore),
            sd_aes = sd(aes_sumScore))
t.test(aes_sumScore ~ group, data = data)
aes_f_glm <- glm(aes_sumScore ~ group, data = data, family = "gaussian")
summary(aes_f_glm)
# bayesian
aes_b_glm <- stan_glm(aes_sumScore ~ group, data = data, 
                  iter = 10000, seed = 123)
aes_b_glm$coefficients
hdi(aes_b_glm)

# BDI
data %>% 
  group_by(group) %>% 
  summarise(mean_bdi = mean(bdi_sumScore),
            sd_bdi = sd(bdi_sumScore))
t.test(bdi_sumScore ~ group, data = data)
bdi_f_glm <- glm(bdi_sumScore ~ group, data = data, family = "gaussian")
summary(bdi_f_glm)
# bayesian
bdi_b_glm <- stan_glm(bdi_sumScore ~ group, data = data, 
                      iter = 10000, seed = 123)
bdi_b_glm$coefficients
hdi(bdi_b_glm)

# OCIR
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_ocir = mean(ocir_sumScore),
            sd_ocir = sd(ocir_sumScore))
t.test(ocir_sumScore ~ group, data = data)
ocir_f_glm <- glm(ocir_sumScore ~ group, data = data, family = "gaussian")
summary(ocir_f_glm)
# bayesian
ocir_b_glm <- stan_glm(ocir_sumScore ~ group, data = data, 
                      iter = 10000, seed = 123)
ocir_b_glm$coefficients
hdi(ocir_b_glm)

# DAQ
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_daq = mean(daq_sumScore),
            sd_daq = sd(daq_sumScore))
t.test(daq_sumScore ~ group, data = data)
daq_f_glm <- glm(daq_sumScore ~ group, data = data, family = "gaussian")
summary(daq_f_glm)
# bayesian
daq_b_glm <- stan_glm(daq_sumScore ~ group, data = data, 
                       iter = 10000, seed = 123)
daq_b_glm$coefficients
hdi(daq_b_glm)

# SHAPS
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_shaps = mean(shaps_sumScore),
            sd_shaps = sd(shaps_sumScore))
t.test(shaps_sumScore ~ group, data = data)
shaps_f_glm <- glm(shaps_sumScore ~ group, data = data, family = "gaussian")
summary(shaps_f_glm)
# bayesian
shaps_b_glm <- stan_glm(shaps_sumScore ~ group, data = data, 
                      iter = 10000, seed = 123)
shaps_b_glm$coefficients
hdi(shaps_b_glm)

# Monetary discounting
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_mcq = mean(mcq_discounting_rate),
            sd_mcq = sd(mcq_discounting_rate))
t.test(mcq_discounting_rate ~ group, data = data)
mcq_glm <- glm(mcq_discounting_rate ~ group, data = data, family = "gaussian")
summary(mcq_glm)
# bayesian
mcq_glm <- stan_glm(mcq_discounting_rate ~ group, data = data, 
                        iter = 10000, seed = 123)
mcq_glm$coefficients
hdi(mcq_glm)

### Three factor eating questionnaire ----
# Correlation between questionnaires
eat_cor <- data %>% 
  select(tfeq_cr_sumScore, tfeq_ue_sumScore, tfeq_ee_sumScore) %>% 
  rename(cognitive_restraint = tfeq_cr_sumScore, uncontrolled_eating = tfeq_ue_sumScore, emotional_eating = tfeq_ee_sumScore) %>% 
  cor()
corrplot::corrplot(eat_cor, method="circle")

# Cognitive restrained
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_cr = mean(tfeq_cr_sumScore),
            sd_cr = sd(tfeq_cr_sumScore))
t.test(tfeq_cr_sumScore ~ group, data = data)
tfeq_cr_f_glm <- glm(tfeq_cr_sumScore ~ group, data = data, family = "gaussian")
summary(tfeq_cr_f_glm)
# bayesian
tfeq_cr_glm <- stan_glm(tfeq_cr_sumScore ~ group, data = data, 
                      iter = 10000, seed = 123)
tfeq_cr_glm$coefficients
hdi(tfeq_cr_glm)

# Uncontrolled eating
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_ue = mean(tfeq_ue_sumScore),
            sd_ue = sd(tfeq_ue_sumScore))
t.test(tfeq_ue_sumScore ~ group, data = data)
tfeq_ue_f_glm <- glm(tfeq_ue_sumScore ~ group, data = data, family = "gaussian")
summary(tfeq_ue_f_glm)
# bayesian
tfeq_ue_glm <- stan_glm(tfeq_ue_sumScore ~ group, data = data, 
                        iter = 10000, seed = 123)
tfeq_ue_glm$coefficients
hdi(tfeq_ue_glm)

tfeq_ue_plot <- raincloud_plot(dat = data, title = "", 
                         xlab = " ", ylab = "Uncontrolled Eating", 
                         predictor_var = "group", outcome_var = "tfeq_ue_sumScore", 
                         predictor_tick_lab = c("control", "treatment"), col = c(color_pal[1], color_pal[2]), 
                         include_grouping = FALSE, direction = "horizontal", scale_seq = c(-0, 1,0.25)) + 
  theme(legend.position = "none") +
  ggtitle("Uncontrolled Eating")

# Emotional eating
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_ee = mean(tfeq_ee_sumScore),
            sd_ee = sd(tfeq_ee_sumScore))
t.test(tfeq_ee_sumScore ~ group, data = data)
tfeq_ee_f_glm <- glm(tfeq_ee_sumScore ~ group, data = data, family = "gaussian")
summary(tfeq_ee_f_glm)
# bayesian
tfeq_ee_glm <- stan_glm(tfeq_ee_sumScore ~ group, data = data, 
                        iter = 10000, seed = 123)
tfeq_ee_glm$coefficients
hdi(tfeq_ee_glm)

# Controlling for antidepressants
tfeq_cr_f_glm <- glm(tfeq_cr_sumScore ~ group + antidepressant, data = data, family = "gaussian")
summary(tfeq_cr_f_glm)

tfeq_ue_f_glm <- glm(tfeq_ue_sumScore ~ group + antidepressant, data = data, family = "gaussian")
summary(tfeq_ue_f_glm)

tfeq_ee_f_glm <- glm(tfeq_ee_sumScore ~ group + antidepressant, data = data, family = "gaussian")
summary(tfeq_ee_f_glm)

# Effect of time on medication in Ozempic group?
time_on_glp_ue_f_glm <- glm(tfeq_ue_sumScore ~ time_on_glp, 
                            data = data %>% filter(group == "treatment"), 
                            family = "gaussian")
summary(time_on_glp_ue_f_glm)

time_on_glp_ue_glm <- stan_glm(tfeq_ue_sumScore ~ time_on_glp, data = data %>% filter(group == "treatment"), 
                               iter = 10000, seed = 123)
time_on_glp_ue_glm$coefficients
hdi(time_on_glp_ue_glm)

ggplot(data %>% filter(group == "treatment"),
       aes(time_on_glp, tfeq_ue_sumScore)) +            
  geom_point() +                                      
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 


### Circadian questionnaires ----
# Correlation between questionnaires
circ_cor <- data %>% 
  select(mctq_continuous, meq_sumScore) %>% 
  rename(MCTQ = mctq_continuous, MEQ = meq_sumScore) %>% 
  cor()
corrplot::corrplot(circ_cor, method="circle")

# MCTQ
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_mctq = mean(mctq_continuous),
            sd_mctq = sd(mctq_continuous))
t.test(mctq_continuous ~ group, data = data)
mctq_f_glm <- glm(mctq_continuous ~ group, data = data, family = "gaussian")
summary(mctq_f_glm)
# bayesian
mctq_glm <- stan_glm(mctq_continuous ~ group, data = data, 
                        iter = 10000, seed = 123)
mctq_glm$coefficients
hdi(mctq_glm)

# MEQ
# frequentist
data %>% 
  group_by(group) %>% 
  summarise(mean_meq = mean(meq_sumScore),
            sd_meq = sd(meq_sumScore))
t.test(meq_sumScore ~ group, data = data)
meq_f_glm <- glm(meq_sumScore ~ group, data = data, family = "gaussian")
summary(meq_f_glm)
# bayesian
meq_glm <- stan_glm(meq_sumScore ~ group, data = data, 
                        iter = 10000, seed = 123)
meq_glm$coefficients
hdi(meq_glm)

# Make chronotype
data %<>% 
  mutate(chronotype = case_when(meq_sumScore > 58 & hm(mctq_MSF_SC) < hm("02:30") ~ "early",
                                meq_sumScore < 42 & hm(mctq_MSF_SC) > hm("05:30") ~ "late",
                                is.na(mctq_MSF_SC) ~ "NA",
                                .default = "intermediate")) 

data %>% 
  tabyl(group, chronotype) %>% 
  chisq.test()
  
# Test relationship between chronotype and time of testing
data %<>% 
  mutate(local_testing_time = as.POSIXct(paste("01jan2000 ", local_testing_time, ":00", sep = ""), 
                                         format = "%d%b%Y %H:%M:%S")) %>%   
  mutate(local_testing_time = case_when(as.POSIXlt(local_testing_time)$hour < 4 ~ .$local_testing_time + days(1), 
                                       .default = .$local_testing_time ))

chronotype_testing_time <- glm(data$local_testing_time %>% as.numeric() ~ 
                          data$chronotype, family = "gaussian")
summary(chronotype_testing_time)

time_chrono_plot <- raincloud_plot(dat = data, title = "", 
                         xlab = " ", ylab = "Testing Time", 
                         predictor_var = "chronotype", outcome_var = "local_testing_time", 
                         predictor_tick_lab = c("early", "intermediate", "late"), col = c(color_pal[1], color_pal[2], color_pal[3]), 
                         include_grouping = FALSE, direction = "horizontal", scale_seq = NULL) + 
  theme(legend.position = "none") +
  ggtitle("Time of testing by Chronotype")

time_meq_plot <- ggplot(data = data, 
                    aes_string(y = "meq_sumScore", x = "local_testing_time")) + 
  geom_point(aes(color = color_pal[3])) + 
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              aes(color = color_pal[3])) 
             









