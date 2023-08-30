########################################################################################################################
#######################---------------- SCREENING PARTICIPANTS FOR MAIN TESTING ----------------########################
########################################################################################################################

### In this script: 
# (1) Parse data
# (2) Exclusion
# (3) Identify treatment and control group
# (4) Adjust data based on Prolific messages
# (5) Check for unrealistic values
# (6) Randomize to testing schedule and identify testing days
# (7) Matching control participants

# Set working directory
here::i_am("github/semaglutide-study/code/screening.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/screener_parsing_fun.R")

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, writexl)


# (1) Parse data ---------------------------------------------------

data_files <- list.files(path = here::here("data/raw_data/screening"), 
                         pattern = ".txt", full.names = TRUE)
meta_files <- list.files(path = here::here("data/raw_data/screening"), 
                         pattern = ".csv", full.names = TRUE)

screening_dat <- screener_parsing(file = data_files[5:6], 
                                  meta_file = meta_files[5:6], 
                                  display_progress = TRUE)

# (2) Exclusion ---------------------------------------------------

# Neurological condition - any severe disorders? 
screening_dat$screening_dat %>%
  filter(neurological == 1) %>% 
  select(prolific_id, neurological_condition, neurological_condition_other) %>% 
  print(n = Inf)

# Exclude subjects with English level < B2 and severe neurological disorders
screening_dat$screening_dat %<>%
  add_column(exclusion = case_when(
    .$english == "A1/A2" | .$english == "B1" ~ 1,
    grepl("Multiple Sclerosis|Traumatic brain injury|
          Stroke|Parkinson's disease", 
          .$neurological_condition, 
          ignore.case = TRUE) ~ 1,
    grepl("brain tumor|Hyperadrenergic POTS|Essential Tremor|
          Petit mal seizures|Cerebral Vasculitis|pseudotumor cerebri|
          seziers|Cataplexy|Conversion disorder|FNES|Fibromyalgia", 
          .$neurological_condition_other, 
          ignore.case = TRUE) ~ 1,
    .default = 0
  ))

# if subj has epilepsy, check what medication they are on 
screening_dat$screening_dat %>% 
  add_column(epilepsy_check = case_when(
    grepl("Epilepsy", 
          .$neurological_condition, 
          ignore.case = TRUE) ~ 1
  )) %>% 
  filter(!is.na(epilepsy_check)) %>%
  select(prolific_id, neurological_condition, medication, other_medication, other_medication_type)

# (3) Identify treatment and control group ---------------------------------------------------

# Add residence (US or UK)
screening_dat$screening_dat %<>% left_join(screening_dat$prolific_dat %>% 
                                             select(prolific_id, residence, screening_day), 
                                           by = "prolific_id")

# Add group (treatment or control)
screening_dat$screening_dat %<>% 
  add_column(group = case_when(
    .$diabetes == "type 2" & 
      .$glp_treatment == "Current" & 
      .$type_glp == "Semaglutide_Ozempic" & 
      .$administration_glp == "Injection" &
      .$schedule_glp == "Weekly" &
      .$exclusion == 0 ~ "treatment",
    .$diabetes == "type 2" & 
      .$glp_treatment == "Current" & 
      .$type_glp == "Semaglutide_Wegovy" & 
      .$administration_glp == "Injection" &
      .$schedule_glp == "Weekly" &
      .$exclusion == 0 ~ "treatment",
    .$diabetes == "type 2" & 
      .$glp_treatment == "No" & 
      is.na(.$type_glp) &
      is.na(.$administration_glp) &
      is.na(.$schedule_glp) &
      .$exclusion == 0 ~ "control", 
    .default = "none"
  ))

screening_dat$screening_dat %>% 
  janitor::tabyl(group)

# Double check everyone in the treatment group have been on Ozempic for long enough

screening_dat$screening_dat %>% 
  filter(group == "treatment") %>% 
  select(prolific_id, start_date_glp, screening_day) %>% 
  mutate(across(start_date_glp, ymd)) %>% 
  mutate(start_date_study = start_date_glp + weeks(4)) 


# (4) Adjust data based on Prolific messages ---------------------------------------------------

## Subject 63f226e602d56e2ab3623b2f
# by chat reported to suffers from anxiety and made a mistake in IPAQ

screening_dat$screening_dat %<>% 
  mutate(psych_neurdev = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                1, psych_neurdev)) %>%
  mutate(psych_neurdev_condition = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                          "[\"Generalised anxiety disorder\",\"Current\"]", psych_neurdev_condition)) %>%
  mutate(ipaq_1 = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                  "None", ipaq_1)) %>% 
  mutate(ipaq_2 = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                  NA, ipaq_2)) %>% 
  mutate(ipaq_3 = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                  "5", ipaq_3)) %>% 
  mutate(ipaq_4 = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                  "0:40", ipaq_4)) %>% 
  mutate(ipaq_5 = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                  "7", ipaq_5)) %>% 
  mutate(ipaq_6 = ifelse(prolific_id == "63f226e602d56e2ab3623b2f", 
                                  "2:0", ipaq_6)) %>%
  
  ## Subject 5e649c198b69732234df7db8
  # by chat reported to have made a mistake in IPAQ
  mutate(ipaq_1 = ifelse(prolific_id == "5e649c198b69732234df7db8", 
                                  "None", ipaq_1)) %>%
  mutate(ipaq_2 = ifelse(prolific_id == "5e649c198b69732234df7db8", 
                                  NA, ipaq_2))
  

# (5) Check for unrealistic values ---------------------------------------------------

# BMI 

screening_dat$screening_dat %>% 
  filter(group == "treatment") %>%
  select(prolific_id, height_cm, weight_kg, bmi, height_cm_raw, weight_kg_raw, group) %>% 
  arrange(bmi) %>% 
  slice_tail(n=20)

# flagged treatment participants: 
# 6095c58f673ef9aa5dafd4e2: height=180, weight=416, bmi=128
# 5dedf76558776a4a4a845c68: height=175, weight=285, bmi=92.7
# 6457e462c16f6f27b5c0519f: height=188, weight=254, bmi=71.9

# => message for clarification when they complete main testing

screening_dat$screening_dat %<>% 
  rowwise() %>%
  mutate(across(height_cm:bmi, ~ case_when(is.infinite(bmi) ~ NA,
                                         is.nan(bmi) ~ NA,
                                         bmi <= 18.5 ~ NA, 
                                         bmi >= 150 ~ NA, 
                                         height_cm > 240 ~ NA,
                                         height_cm < 140 ~ NA,
                                         weight_kg < 30 ~ NA,
                                         weight_kg > 400 ~ NA,
                                        .default = .))) 


# IPAQ
# data cleaning according to manual

suppressMessages({
  suppressWarnings({
    # Make all NA if time spent exercising > 16h per day
    screening_dat$screening_dat %<>% 
      rowwise() %>%
      mutate(across(ipaq_1:ipaq_sumScore, ~ case_when(sum(period_to_seconds(hm(ipaq_2)), 
                                                          period_to_seconds(hm(ipaq_4)), 
                                                          period_to_seconds(hm(ipaq_6)), 
                                                          na.rm = TRUE) > 57600 ~ NA, 
                                                      .default = .))) 
    
    # Round exercise < 10 minutes down to 00:00
    screening_dat$screening_dat %<>% 
      mutate(across(c(ipaq_2, ipaq_4, ipaq_6), ~ case_when(is.na(.) ~ NA,
                                                           period_to_seconds(hm(.)) < 600 ~ "00:00", 
                                                           .default = .)))
    
    # round exercise > 180 minutes down to 03:00
    screening_dat$screening_dat %<>% 
      mutate(across(c(ipaq_2, ipaq_4, ipaq_6), ~ case_when(is.na(.) ~ NA,
                                                           period_to_seconds(hm(.)) > 10800 ~ "03:00", 
                                                           .default = .)))
  })
})

screening_dat$screening_dat %>%
  filter(group == "treatment") %>% 
  select(ipaq_1:ipaq_sumScore) %>% 
  print(n=Inf)

# no treatment group participants have NA IPAQ data

# (6) Randomize to testing schedule and identify testing days  ---------------------------------------------------

screening_dat$screening_dat %>% 
  filter(group == "treatment") %>% 
  add_column(schedule_group = sample(1:2, nrow(.), TRUE)) %>%
  select(prolific_id, screening_day, injection_day_glp, schedule_group, residence, start_date_glp) %>%
  mutate(across(start_date_glp, ymd)) %>% 
  mutate(start_date_study = start_date_glp + weeks(4)) %>% 
  rowwise() %>% 
  mutate(session_1_day = case_when(schedule_group == 1 ~  weekdays(Sys.Date() + 
                                                                    match(injection_day_glp, 
                                                                          weekdays(Sys.Date()+1:7)) - 1),
                                   schedule_group == 2 ~  weekdays(Sys.Date() + 
                                                                    match(injection_day_glp, 
                                                                          weekdays(Sys.Date()+1:7)) + 1),
                                   .default = NA
  )) %>% 
  mutate(session_2_day = case_when(schedule_group == 1 ~  weekdays(Sys.Date() + 
                                                                     match(injection_day_glp, 
                                                                           weekdays(Sys.Date()+1:7)) + 1),
                                   schedule_group == 2 ~  weekdays(Sys.Date() + 
                                                                     match(injection_day_glp, 
                                                                           weekdays(Sys.Date()+1:7)) - 1),
                                   .default = NA
  ))





# (7) Matching control participants ---------------------------------------------------

# To Do: 
# - Identify treatment participants for that have already have a match 
# - Exclude control participants from matching that have any missing data on the matching variables
# - Revisit the BMI issue
# => check if any treatment people have BMI values worth double checking and do so before matching them
# => check I don't match any controls with unrealistic BMI values












