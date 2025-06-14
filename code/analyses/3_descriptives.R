####################################################################################################
#######################---------------- SAMPLE DESCRIPTIVES ----------------########################
####################################################################################################

### In this script: 
# (1) Sample description
# (2) Demographics
# (3) Psychiatric comorbidities
# (4) Confirm effort discounting effects

# Set working directory
here::i_am("github/semaglutide-study/code/analyses/3_descriptives.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")

# source datasets
main_data <- readRDS("github/semaglutide-study/data/processed_data/main_data.RDS")
non_diabetic_data_new <- readRDS("github/semaglutide-study/data/processed_data/non_diabetic_matched.RDS")
non_diabetic_normal_weight_data_new <- readRDS("github/semaglutide-study/data/processed_data/non_diabetic_normal_weight_matched.RDS")


# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, 
                 writexl, lubridate, purrr, magrittr, nlme)

### (1) Sample description -----------------------------------------------

main_data$demographic_data %>% 
  tabyl(group)
# 58 treatment participants
# 54 control participants

main_data$glp_data %>% 
  filter(session == 2) %>% 
  select(last_injection_days) %>% 
  tabyl(last_injection_days)
# 25 participants for within subject comparison

non_diabetic_data$demographic_data
# 58 non-diabetic controls

### (2) Demographics -----------------------------------------------

all_demographics <- bind_rows(
  main_data$demographic_data %>% 
    select(subj_id, group, age, gender, ethnicity, ses, bmi, psych_neurdev, psych_neurdev_condition, 
           psych_neurdev_condition_other, antidepressant, antidepressant_type, antidepressant_type_other, 
           chronic_disease, chronic_disease_condition, chronic_disease_condition_other),
  non_diabetic_data$demographic_data %>% 
    select(subj_id, age, gender, ses, bmi, psych_neurdev, psych_neurdev_condition, antidepressant, 
           antidepressant_type) %>% 
    add_column(group = "non_diabetic", 
               .after = "subj_id"), 
  non_diabetic_normal_weight_data$demographic_data %>% 
    select(subj_id, age, gender, ses, bmi, psych_neurdev, psych_neurdev_condition, antidepressant, 
           antidepressant_type) %>% 
    add_column(group = "non_diabetic_normalweight", 
               .after = "subj_id")) %>% 
  left_join(bind_rows(
    main_data$questionnaire_data %>% 
      filter(session == 1) %>% 
      select(subj_id, group, ipaq_sumScore),
    non_diabetic_data$questionnaire_data %>% 
      select(subj_id, ipaq_sumScore) %>% 
      add_column(group = "non_diabetic", 
                 .after = "subj_id"),
    non_diabetic_normal_weight_data$questionnaire_data %>% 
      select(subj_id, ipaq_sumScore) %>% 
      add_column(group = "non_diabetic_normalweight", 
                 .after = "subj_id")))

all_demographics %>% 
  group_by(group) %>%
  summarise(N=n(),
            mean_age = mean(age), sd_age = sd(age), min_age = min(age), max_age = max(age),
            median_ses = median(ses), iqr_upper_ses = quantile(ses, 0.25), iqr_lower_ses = quantile(ses, 0.75), 
            mean_bmi = mean(bmi), sd_bmi = sd(bmi), min_bmi = min(bmi), max_bmi = max(bmi),
            mean_ipaq = mean(ipaq_sumScore), sd_ipaq = sd(ipaq_sumScore), min_ipaq = min(ipaq_sumScore), max_ipaq = max(ipaq_sumScore))

all_demographics %>% 
  tabyl(group, gender)

# Test matching

# Age
summary(aov(age ~ factor(group), data = all_demographics))

# Gender
all_demographics %>% 
  janitor::tabyl(group, gender) %>% 
  chisq.test()

# BMI
summary(aov(bmi ~ factor(group), data = all_demographics %>% filter(group != "non_diabetic_normalweight")))

# IPAQ
summary(aov(ipaq_sumScore ~ factor(group), data = all_demographics))


### (3) Psychiatric comorbidities -----------------------------------------------

# Any psychiatric comorbidity
all_demographics %>% 
  tabyl(group, psych_neurdev) %>% 
  chisq.test()

# Depression
all_demographics %>% 
  group_by(group) %>%
  summarise(N=(str_subset(psych_neurdev_condition, 
                          regex("Major depressive disorder|Depression|MDD", ignore_case = TRUE)) %>% 
                  length()), 
            perc = (N/n())*100)

# Anxiety
all_demographics %>% 
  group_by(group) %>%
  summarise(N=(str_subset(psych_neurdev_condition, 
                          regex("Social anxiety disorder|Generalised anxiety disorder", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

# Addiction
all_demographics %>% 
  group_by(group) %>%
  summarise(N=(str_subset(psych_neurdev_condition, 
                          regex("Substance dependence or addictive disorder", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

# OCD
all_demographics %>% 
  group_by(group) %>%
  summarise(N=(str_subset(psych_neurdev_condition, 
                          regex("Obsessive compulsive disorder", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

# Antidepressant use
# Any
all_demographics %>%
  mutate(antidepressant = ifelse(is.na(antidepressant), 0, antidepressant)) %>% 
  tabyl(group, antidepressant) %>% 
  chisq.test()

# Type
all_demographics %>% 
  group_by(group) %>%
  summarise(N=(str_subset(antidepressant_type, 
                          regex("Fluoxetine", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

all_demographics %>% 
  group_by(group) %>%
  summarise(N=(str_subset(antidepressant_type, 
                          regex("Venlafaxine", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

### (4) Weight loss interventions -----------------------------------------------

main_data$demographic_data %>% 
  tabyl(group, weight_loss)

# Exercise
main_data$demographic_data %>% 
  group_by(group) %>% 
  summarise(N=(str_subset(weight_loss_interventions, 
                          regex("Exercise", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

# Dieting
main_data$demographic_data %>% 
  group_by(group) %>% 
  summarise(N=(str_subset(weight_loss_interventions, 
                          regex("Dieting", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

# Meal replacements
main_data$demographic_data %>% 
  group_by(group) %>% 
  summarise(N=(str_subset(weight_loss_interventions, 
                          regex("Meal replacements", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

# Supplements
main_data$demographic_data %>% 
  group_by(group) %>% 
  summarise(N=(str_subset(weight_loss_interventions, 
                          regex("Supplements", ignore_case = TRUE)) %>% 
                 length()), 
            perc = (N/n())*100)

### (5) Confirm effort discounting effects -----------------------------------------------

### Treatment
# Mixed-effects ANOVA
choice_data <- main_data$task_data %>% 
  filter(phase == "game", 
         group == "treatment") %>%
  group_by(subj_id, effort, reward) %>% 
  summarize(mean_choice = mean(choice))

anova(lme(mean_choice ~ effort * reward, random= ~1|subj_id, 
          data = choice_data))

# Post hoc anovas
# per reward level
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 1)))
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 2)))
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 3)))
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 4)))

# per effort level
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 1)))
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 2)))
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 3)))
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 4)))


### Controls
# Mixed-effects ANOVA
choice_data <- main_data$task_data %>% 
  filter(phase == "game", 
         group == "control") %>%
  group_by(subj_id, effort, reward) %>% 
  summarize(mean_choice = mean(choice))

anova(lme(mean_choice ~ effort * reward, random= ~1|subj_id, 
          data = choice_data))

# Post hoc anovas
# per reward level
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 1)))
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 2)))
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 3)))
anova(lme(mean_choice ~ effort, random= ~1|subj_id, 
          data = choice_data %>% filter(reward == 4)))

# per effort level
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 1)))
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 2)))
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 3)))
anova(lme(mean_choice ~ reward, random= ~1|subj_id, 
          data = choice_data %>% filter(effort == 4)))


### Non-diabetics
# Mixed-effects ANOVA
choice_data <- non_diabetic_data$task_data %>% 
  filter(phase == "game") %>%
  group_by(subj_id, offerEffort, offerReward) %>% 
  summarize(mean_choice = mean(choice))

anova(lme(mean_choice ~ offerEffort * offerReward, random= ~1|subj_id, 
          data = choice_data))

# Post hoc anovas
# per reward level
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 1)))
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 2)))
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 3)))
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 4)))

# per effort level
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 1)))
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 2)))
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 3)))
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 4)))


### Non-diabetics, low BMI
# Mixed-effects ANOVA
choice_data <- non_diabetic_normal_weight_data$task_data %>% 
  filter(phase == "game") %>%
  group_by(subj_id, offerEffort, offerReward) %>% 
  summarize(mean_choice = mean(choice))

anova(lme(mean_choice ~ offerEffort * offerReward, random= ~1|subj_id, 
          data = choice_data))

# Post hoc anovas
# per reward level
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 1)))
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 2)))
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 3)))
anova(lme(mean_choice ~ offerEffort, random= ~1|subj_id, 
          data = choice_data %>% filter(offerReward == 4)))

# per effort level
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 1)))
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 2)))
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 3)))
anova(lme(mean_choice ~ offerReward, random= ~1|subj_id, 
          data = choice_data %>% filter(offerEffort == 4)))

