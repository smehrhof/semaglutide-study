##################################################################################################################
#######################---------------- PROCESSING DATA FROM MAIN TESTING ----------------########################
##################################################################################################################

### In this script: 
# (1) Parse data
# (2) Adjust data based on Prolific messages
# (3) Merge with screening data
# (4) Exclusion
# (5) Estimate hours since injection
# (6) Save processed data 
# (7) Find and match non-diabetic controls


# Set working directory
here::i_am("github/semaglutide-study/code/analyses/2_data_processing.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/parsing_fun.R")

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, writexl, tidyr)
# (1) Parse data ---------------------------------------------------

### Load file with extra dose data ------------------------------------------------------------------------------------------

raw_doses_dat <- read.csv(here::here("github/semaglutide-study/data/raw_data/treatment_group/doses.csv"), 
                          sep = ";") %>% as_tibble() 

last_injection_dat <- read.csv(here::here("github/semaglutide-study/data/raw_data/treatment_group/last_injection.csv"), 
                               sep = ";") %>% as_tibble() 

raw_dat <- readRDS("github/semaglutide-study/data/raw_data/treatment_group/session_1/raw_dat.RDS")

raw_meta_dat <- readRDS("github/semaglutide-study/data/raw_data/treatment_group/session_1/raw_meta_dat.RDS")


s1_data <- parsing(raw_dat = raw_dat, 
                   raw_meta_dat = raw_meta_dat, 
                   doses_dat = raw_doses_dat, 
                   last_injection_dat = last_injection_dat,
                   session = "s1", 
                   display_progress = TRUE)

### Treatment - Session 2

raw_doses_dat <- read.csv(here::here("github/semaglutide-study/data/raw_data/treatment_group/doses.csv"), 
                          sep = ";") %>% as_tibble() 

last_injection_dat <- tibble("prolific_id" = NA)

raw_dat <- readRDS("github/semaglutide-study/data/raw_data/treatment_group/session_2/raw_dat.RDS")

raw_meta_dat <- readRDS("github/semaglutide-study/data/raw_data/treatment_group/session_2/raw_meta_dat.RDS")


s2_data <- parsing(raw_dat = raw_dat, 
                   raw_meta_dat = raw_meta_dat, 
                   doses_dat = raw_doses_dat, 
                   last_injection_dat = last_injection_dat,
                   session = "s2", 
                   display_progress = TRUE)



### Controls

raw_doses_dat <- read.csv(here::here("github/semaglutide-study/data/raw_data/control_group/doses.csv"), 
                          sep = ";") %>% as_tibble() 

raw_dat <- readRDS("github/semaglutide-study/data/raw_data/control_group/raw_dat.RDS")

raw_meta_dat <- readRDS("github/semaglutide-study/data/raw_data/control_group/raw_meta_dat.RDS")

controls_data <- parsing(raw_dat = raw_dat, 
                         raw_meta_dat = raw_meta_dat, 
                         doses_dat = raw_doses_dat, 
                         last_injection_dat = NA,
                         session = "control", 
                         display_progress = TRUE)

# Correct "other medication" input bug
controls_data$demographic_dat %>% 
  add_column(
    other_med_text = { ifelse(.$other_med == 0, 
                              NA, 
                              { ifelse(.$prolific_id == "9585dd25f3", 
                                       "vitamin D, tumeric", 
                                       NA)})}, 
    .after = "other_med_type") %>% 
  mutate(contraceptive_med = case_when(prolific_id == "52f6024277" ~ "1",
                                    prolific_id == "9585dd25f3" ~ "1",
                                    .default = other_med_type)) %>%
  mutate(antidepressant_med = case_when(prolific_id == "52f6024277" ~ "1",
                                       prolific_id == "9585dd25f3" ~ "1",
                                       .default = other_med_type)) %>%
  mutate(antidepressant_med_type = case_when(prolific_id == "52f6024277" ~ "bupropion",
                                        prolific_id == "9585dd25f3" ~ " ",
                                        .default = other_med_type)) %>%
  mutate(other_med_type = case_when(prolific_id == "52f6024277" ~ NA,
                                    prolific_id == "9585dd25f3" ~ "Vitamin or supplement",
                                    .default = other_med_type))


# (2) Adjust data based on Prolific messages  ---------------------------------------------------

# subject f9361326d2 having made an error in glp check -> is on Ozempic treatment
s1_data$demographic_dat %<>%
  mutate(across(glp_check, ~ case_when(prolific_id == "f9361326d2" ~ 1, .default = .)),
         across(glp_past, ~ case_when(prolific_id == "f9361326d2" ~ NA, .default = .)),
         across(glp_dose_mg, ~ case_when(prolific_id == "f9361326d2" ~ "1", .default = .)),
         across(last_injection_days, ~ case_when(prolific_id == "f9361326d2" ~ 1, .default = .)),
         across(last_injection_time, ~ case_when(prolific_id == "f9361326d2" ~ "10am - 2pm", .default = .))) 


# (3) Merge with screening data ---------------------------------------------------

screening_dat <- readRDS("github/semaglutide-study/data/processed_data/screening_dat.RDS")

screening_dat$screening_dat %>% 
  tabyl(group)

s1_data %<>% 
  list_modify(screening_dat = screening_dat$screening_dat %>%
                filter(prolific_id %in% s1_data$demographic_dat$prolific_id))
s2_data %<>% 
  list_modify(screening_dat = screening_dat$screening_dat %>%
                filter(prolific_id %in% s2_data$demographic_dat$prolific_id))
controls_data %<>% 
  list_modify(screening_dat = screening_dat$screening_dat %>%
                filter(prolific_id %in% controls_data$demographic_dat$prolific_id))

# (4) Exclusion ---------------------------------------------------

### Treatment - Session 1
s1_exclusion <- c()

# Treatment group conditions: (possibly) not on treatment or not diabetic
s1_exclusion <- append(s1_exclusion, 
                       s1_data$demographic_dat %>% filter(glp_check == 0 |
                                                            # participant reported not to be diabetic (only insulin resistant)
                                                            prolific_id == "8f3602f44a" |
                                                            # participant reported conflicting information about what semaglutide treatment they are on
                                                            prolific_id == "8e16120561") %>% 
                         .$prolific_id)
# 10 participants excluded

# Task based
# no offers accepted
s1_exclusion <- append(s1_exclusion, 
                       s1_data$task_dat %>% 
                         filter(phase == "game") %>% 
                         reframe(acceptance_rate = mean(choice), .by = prolific_id) %>% 
                         filter(acceptance_rate == 0) %>% 
                         .$prolific_id)
# 0 participants excluded

# large difference between minimum and maximum clicking speed during calibration
# -> >3SD based on large sample: 128.9755
s1_exclusion <- append(s1_exclusion, 
                       s1_data$task_dat %>% 
                         filter(phase == "calibration" & 
                                  trial != 1) %>% 
                         reframe(calibration_difference = diff(clicks) %>% abs(), 
                                 .by = prolific_id) %>% 
                         filter(calibration_difference > 128.9755) %>% 
                         .$prolific_id)
# 0 participants excluded

# large change in clicking calibration pre- to post-task
# -> >3SD based on large sample: 133.7341
s1_exclusion <- append(s1_exclusion, 
                       s1_data$task_dat %>% 
                         filter(phase == "calibration") %>% 
                         reframe(pre_max = max(clicks[1:3]),
                                 post_max = clicks[4],
                                 .by = prolific_id) %>% 
                         mutate(calibration_change = abs(pre_max- post_max)) %>%
                         filter(calibration_change > 133.7341) %>% 
                         .$prolific_id)
# 0 participants excluded

# Questionnaire based: catch questions
s1_exclusion <- append(s1_exclusion, 
                       s1_data$questionnaire_dat %>% 
                         filter(catch_questions_pass == 0) %>% 
                         .$prolific_id)
# 2 participants excluded

# Exclude
s1_data <- lapply(s1_data, 
                  function(df){
                    df %>%
                      filter(!prolific_id %in% s1_exclusion)
                  })

### Treatment - Session 2
s2_exclusion <- c()

# Treatment group conditions: (possibly) not on treatment or not diabetic
s2_exclusion <- append(s2_exclusion, 
                       s2_data$demographic_dat %>% filter(glp_check == 0 |
                                                            # participant reported not to be diabetic (only insulin resistant)
                                                            prolific_id == "8f3602f44a") %>% 
                         .$prolific_id)
# 2 participants excluded

# Task based
# no offers accepted
s2_exclusion <- append(s2_exclusion, 
                       s2_data$task_dat %>% 
                         filter(phase == "game") %>% 
                         reframe(acceptance_rate = mean(choice), .by = prolific_id) %>% 
                         filter(acceptance_rate == 0) %>% 
                         .$prolific_id)
# 0 participants excluded

# large difference between minimum and maximum clicking speed during calibration
# -> >3SD based on large sample: 128.9755
s2_exclusion <- append(s2_exclusion, 
                       s2_data$task_dat %>% 
                         filter(phase == "calibration" & 
                                  trial != 1) %>% 
                         reframe(calibration_difference = diff(clicks) %>% abs(), 
                                 .by = prolific_id) %>% 
                         filter(calibration_difference > 128.9755) %>% 
                         .$prolific_id)
# 0 participants excluded

# large change in clicking calibration pre- to post-task
# -> >3SD based on large sample: 133.7341
s2_exclusion <- append(s2_exclusion, 
                       s2_data$task_dat %>% 
                         filter(phase == "calibration") %>% 
                         reframe(pre_max = max(clicks[1:3]),
                                 post_max = clicks[4],
                                 .by = prolific_id) %>% 
                         mutate(calibration_change = abs(pre_max- post_max)) %>%
                         filter(calibration_change > 133.7341) %>% 
                         .$prolific_id)
# 0 participants excluded

# Questionnaire based: catch questions
s2_exclusion <- append(s2_exclusion, 
                       s2_data$questionnaire_dat %>% 
                         filter(catch_questions_pass == 0) %>% 
                         .$prolific_id)
# 0 participants excluded

# Subjects excluded from within subject comparison due to changes in injection day
s2_exclusion <- append(s2_exclusion, 
                       c("6435dd8940", "0da634002b", "f6e6324b81", 
                         "16164003f2", "7fc627dcf2"))

# Exclude
s2_data <- lapply(s2_data, 
                  function(df){
                    df %>%
                      filter(!prolific_id %in% s2_exclusion)
                  })

### Controls

# Do any collected control participants now qualify for treatment condition?
controls_data$demographic_dat %>% 
  filter(glp_treatment == "Current" & 
           type_glp == "Semaglutide_Ozempic" & 
           administration_glp == "Injection" & 
           schedule_glp == "Weekly")
# make sure they aren't excluded for other reasons
controls_data$demographic_dat %>% 
  filter(prolific_id == "9955ce5d61") %>%
  select(start_date_glp)
controls_data$task_dat %>% 
  filter(prolific_id == "9955ce5d61") %>%
  filter(phase == "game") %>% 
  reframe(acceptance_rate = mean(choice), 
          pre_max = max(clicks[1:3]),
          post_max = clicks[4])
controls_data$task_dat %>% 
  filter(prolific_id == "9955ce5d61") %>%
  filter(phase == "calibration" & 
           trial != 1) %>% 
  reframe(calibration_difference = diff(clicks) %>% abs())
controls_data$questionnaire_dat %>% 
  filter(prolific_id == "9955ce5d61") %>%
  select(catch_questions_pass)

# Can be included in treatment group -> add to treatment data
# demographics
s1_data$demographic_dat %<>% 
  bind_rows(controls_data$demographic_dat %>% 
              filter(prolific_id == "9955ce5d61") %>% 
              select(prolific_id, glp_treatment, glp_dose_mg, last_injection_days, last_injection_time, 
                     diabetes_med:time_response) %>% 
              rename(glp_check = glp_treatment) %>% 
              mutate(glp_check = 1)) 
# meta task data
s1_data$task_meta_dat %<>% 
  bind_rows(controls_data$task_meta_dat %>% 
              filter(prolific_id == "9955ce5d61")) 
# task data
s1_data$task_dat %<>% 
  bind_rows(controls_data$task_dat %>% 
              filter(prolific_id == "9955ce5d61")) 
# questionnaire data
s1_data$questionnaire_dat %<>% 
  bind_rows(controls_data$questionnaire_dat %>% 
              filter(prolific_id == "9955ce5d61")) 
# prolific data
s1_data$prolific_dat %<>% 
  bind_rows(controls_data$prolific_dat %>% 
              filter(prolific_id == "9955ce5d61")) 
# screening data - update as their medication status has changes since
s1_data$screening_dat %<>% 
  bind_rows(controls_data$screening_dat %>% 
              filter(prolific_id == "9955ce5d61") %>% 
              mutate(across(.cols = c(glp_treatment:side_effects_glp_symptoms_other)), 
                     controls_data$demographic_dat %>%   
                       filter(prolific_id == "9955ce5d61") %>% 
                       select(glp_treatment:side_effects_glp_symptoms_other)))

# Exclusion for controls
controls_exclusion <- c()

# Treatment group conditions: (possibly) not on treatment or not diabetic
controls_exclusion <- append(controls_exclusion, 
                             controls_data$demographic_dat %>% filter(glp_treatment != "No") %>% 
                               .$prolific_id)
# 9 participants excluded

# Task based
# no offers accepted
controls_exclusion <- append(controls_exclusion, 
                             controls_data$task_dat %>% 
                               filter(phase == "game") %>% 
                               reframe(acceptance_rate = mean(choice), .by = prolific_id) %>% 
                               filter(acceptance_rate == 0) %>% 
                               .$prolific_id)
# 0 participants excluded

# large difference between minimum and maximum clicking speed during calibration
# -> >3SD based on large sample: 128.9755
controls_exclusion <- append(controls_exclusion, 
                             controls_data$task_dat %>% 
                               filter(phase == "calibration" & 
                                        trial != 1) %>% 
                               reframe(calibration_difference = diff(clicks) %>% abs(), 
                                       .by = prolific_id) %>% 
                               filter(calibration_difference > 128.9755) %>% 
                               .$prolific_id)
# 0 participants excluded

# large change in clicking calibration pre- to post-task
# -> >3SD based on large sample: 133.7341
controls_exclusion <- append(controls_exclusion, 
                             controls_data$task_dat %>% 
                               filter(phase == "calibration") %>% 
                               reframe(pre_max = max(clicks[1:3]),
                                       post_max = clicks[4],
                                       .by = prolific_id) %>% 
                               mutate(calibration_change = abs(pre_max- post_max)) %>%
                               filter(calibration_change > 133.7341) %>% 
                               .$prolific_id)
# 0 participants excluded

# Questionnaire based: catch questions
controls_exclusion <- append(controls_exclusion,
                             controls_data$questionnaire_dat %>% 
                               filter(catch_questions_pass == 0) %>% 
                               .$prolific_id)
# 1 participants excluded

# Exclude
controls_data <- lapply(controls_data, 
                  function(df){
                    df %>%
                      filter(!prolific_id %in% controls_exclusion)
                  })

### Across groups: exclude for BMI
# exclude participants who report BMI > median + 2 * sd
bmi_excl <- c(s1_data$screening_dat$bmi, controls_data$screening_dat$bmi) %>% 
  median() + 
  (2 * c(s1_data$screening_dat$bmi, controls_data$screening_dat$bmi) %>% 
     sd())

# cut off: 36.70662 + 2*(13.04223) = 62.79108
controls_exclusion <- append(controls_exclusion, 
                             controls_data$screening_dat %>% 
                               filter(bmi >= bmi_excl) %>% 
                               .$prolific_id)
s1_exclusion <- append(s1_exclusion, 
                       s1_data$screening_dat %>% 
                         filter(bmi >= bmi_excl) %>% 
                         .$prolific_id)

# Exclude
controls_data <- lapply(controls_data, 
                  function(df){
                    df %>%
                      filter(!prolific_id %in% controls_exclusion)
                  })
s1_data <- lapply(s1_data, 
                  function(df){
                    df %>%
                      filter(!prolific_id %in% s1_exclusion)
                  })

# (5) Structure and combine data ---------------------------------------------------

# Arrange tibbles to be in the same order by prolific id
s1_data <- lapply(s1_data, 
                        function(df){
                          df %>%
                            arrange(prolific_id) %>% 
                            rename(subj_id = prolific_id)
                        })
s2_data <- lapply(s2_data, 
                        function(df){
                          df %>%
                            arrange(prolific_id) %>% 
                            rename(subj_id = prolific_id)
                        })
controls_data <- lapply(controls_data, 
                        function(df){
                          df %>%
                            arrange(prolific_id) %>% 
                            rename(subj_id = prolific_id)
                        })

### Demographics
demographic_data <- bind_rows(
  # Treatment
  s1_data$screening_dat %>% 
    select(subj_id:diabetes_other, 
           diabetes_medication:chronic_disease_condition_other, 
           height_cm:weight_loss_interventions_other, 
           residence) %>% 
    left_join(s1_data$prolific_dat %>% 
                select(subj_id, sex)) %>% 
    left_join(s1_data$demographic_dat %>% 
                select(subj_id, games)) %>% 
    # update medication information if indicated
    left_join(s1_data$demographic_dat %>% 
                select(subj_id, diabetes_med:antidepressant_med_type)) %>% 
    # diabetes medication
    mutate(diabetes_medication = ifelse(diabetes_med == 1, diabetes_med, diabetes_medication)) %>% 
    mutate(diabetes_medication_type = ifelse(diabetes_med == 1, diabetes_med_type, diabetes_medication_type)) %>% 
    mutate(diabetes_medication_type_other = ifelse(diabetes_med == 1, diabetes_med_type_other, diabetes_medication_type_other)) %>%
    # other medication
    mutate(medication = ifelse(other_med == 1, other_med, medication)) %>% 
    mutate(contraceptive = ifelse(other_med == 1, contraceptive_med, contraceptive)) %>% 
    mutate(antidepressant = ifelse(other_med == 1, antidepressant_med, antidepressant)) %>% 
    mutate(antidepressant_type = ifelse(other_med == 1, antidepressant_med_type, antidepressant_type)) %>% 
    mutate(other_medication_type = ifelse(other_med == 1, other_med_type, other_medication_type)) %>% 
    add_column(group = "treatment", .after = "subj_id") %>%
    select(-(diabetes_med:antidepressant_med_type)),
  
  # Controls
  controls_data$screening_dat %>% 
    select(subj_id:diabetes_other, 
           diabetes_medication:chronic_disease_condition_other, 
           height_cm:weight_loss_interventions_other, 
           residence) %>% 
    left_join(controls_data$prolific_dat %>% 
                select(subj_id, sex)) %>% 
    left_join(controls_data$demographic_dat %>% 
                select(subj_id, games)) %>% 
    # update medication information if indicated
    left_join(controls_data$demographic_dat %>% 
                select(subj_id, diabetes_med:antidepressant_med_type)) %>% 
    # diabetes medication
    mutate(diabetes_medication = ifelse(diabetes_med == 1, diabetes_med, diabetes_medication)) %>% 
    mutate(diabetes_medication_type = ifelse(diabetes_med == 1, diabetes_med_type, diabetes_medication_type)) %>% 
    mutate(diabetes_medication_type_other = ifelse(diabetes_med == 1, diabetes_med_type_other, diabetes_medication_type_other)) %>%
    # other medication
    mutate(medication = ifelse(other_med == 1, other_med, medication)) %>% 
    mutate(contraceptive = ifelse(other_med == 1, contraceptive_med, contraceptive)) %>% 
    mutate(antidepressant = ifelse(other_med == 1, antidepressant_med, antidepressant)) %>% 
    mutate(antidepressant_type = ifelse(other_med == 1, antidepressant_med_type, antidepressant_type)) %>%
    mutate(other_medication_type = ifelse(other_med == 1, other_med_type, other_medication_type)) %>% 
    add_column(group = "control", .after = "subj_id") %>% 
    select(-(diabetes_med:antidepressant_med_type)),
)

### GLP-1 treatment
glp_data <- bind_rows(
  # Treatment
  # session 1
  s1_data$screening_dat %>% 
    select(subj_id, glp_treatment, type_glp, start_date_glp, administration_glp, schedule_glp, 
           injection_day_glp:side_effects_glp_symptoms_other) %>%
    left_join(s1_data$demographic_dat %>% 
                select(subj_id, glp_check, glp_dose_mg, last_injection_days, last_injection_time, time_response) %>% 
                rename(local_testing_time = time_response)) %>% 
    left_join(s1_data$prolific_dat %>% 
                select(subj_id, study_start) %>% 
                rename(testing_day = study_start)) %>% 
      add_column(group = "treatment", 
                 session = 1, 
                 .after = "subj_id"),
  # session 2 
  s2_data$screening_dat %>%
    select(subj_id, glp_treatment, type_glp, start_date_glp, administration_glp, schedule_glp, 
           injection_day_glp:side_effects_glp_symptoms_other) %>%
    left_join(s2_data$demographic_dat %>% 
                select(subj_id, glp_check, glp_dose_mg, last_injection_days, last_injection_time, time_response) %>% 
                rename(local_testing_time = time_response)) %>% 
    left_join(s2_data$prolific_dat %>% 
                select(subj_id, study_start) %>% 
                rename(testing_day = study_start)) %>% 
    add_column(group = "treatment", 
               session = 2, 
               .after = "subj_id"), 
  # controls
  controls_data$screening_dat %>% 
    select(subj_id, glp_treatment, type_glp, start_date_glp, administration_glp, schedule_glp, 
           injection_day_glp:side_effects_glp_symptoms_other) %>%
    # add session 1 glp data for all participants
    left_join(controls_data$demographic_dat %>% 
                select(subj_id, glp_treatment, glp_dose_mg, last_injection_days, last_injection_time, time_response) %>% 
                rename(glp_check = glp_treatment, local_testing_time = time_response)) %>% 
    left_join(controls_data$prolific_dat %>% 
                select(subj_id, study_start) %>% 
                rename(testing_day = study_start)) %>% 
    mutate(glp_check = case_when(glp_check == "No" ~ 0, 
                                    glp_check == "Yes" ~ 1, 
                                    .default = NA)) %>% 
    add_column(group = "control", 
               session = 1, 
               .after = "subj_id")) %>% 
  mutate(local_testing_time_d = ymd_hm(paste("2000-01-01", local_testing_time, sep = " ")), 
         last_injection_time_d = case_when(last_injection_time == "6am - 10am" ~ ymd_hm("2000-01-01 08:00"),
                                         last_injection_time == "10am - 2pm" ~ ymd_hm("2000-01-01 12:00"), 
                                         last_injection_time == "2pm - 6pm" ~ ymd_hm("2000-01-01 16:00"), 
                                         last_injection_time == "6pm - 10pm" ~ ymd_hm("2000-01-01 20:00"), 
                                         last_injection_time == "after 10pm" ~ ymd_hm("2000-01-01 23:59"))) %>% 
  mutate(hours_since_injection = last_injection_days * 24 + difftime(local_testing_time_d, last_injection_time_d, units="hours")) %>%
  mutate(local_testing_time_d = NULL, last_injection_time_d = NULL)
  
### Task - meta data
task_meta_data <- bind_rows(
  s1_data$task_meta_dat %>% 
    add_column(session = 1,
               group = "treatment", 
               .after = "subj_id"), 
  s2_data$task_meta_dat %>% 
    add_column(session = 2,
               group = "treatment", 
               .after = "subj_id"), 
  controls_data$task_meta_dat %>% 
    add_column(session = 1,
               group = "control", 
               .after = "subj_id")) 

# Task - trial data
task_data <- bind_rows(
  s1_data$task_dat %>% 
    add_column(session = 1,
               group = "treatment", 
               .after = "subj_id"), 
  s2_data$task_dat %>% 
    add_column(session = 2,
               group = "treatment", 
               .after = "subj_id"), 
  controls_data$task_dat %>% 
    add_column(session = 1,
               group = "control", 
               .after = "subj_id")) 

# Questionnaires
questionnaire_data <- bind_rows(
  s1_data$questionnaire_dat %>% 
    add_column(session = 1,
               group = "treatment", 
               .after = "subj_id") %>%
    select(subj_id:tfeq_ee_sumScore) %>% 
    left_join(s1_data$screening_dat %>% 
                select(subj_id, ipaq_1:ipaq_sumScore)) %>%
    left_join(s1_data$questionnaire_dat %>% 
                select(subj_id, catch_question_1:other_comments)), 
  
  s2_data$questionnaire_dat %>% 
    add_column(session = 2,
               group = "treatment", 
               .after = "subj_id"), 
  
  controls_data$questionnaire_dat %>% 
    add_column(session = 1,
               group = "control", 
               .after = "subj_id") %>%
    select(subj_id:tfeq_ee_sumScore) %>% 
    left_join(controls_data$screening_dat %>% 
                select(subj_id, ipaq_1:ipaq_sumScore)) %>%
    left_join(controls_data$questionnaire_dat %>% 
                select(subj_id, catch_question_1:other_comments)))

data <- list(demographic_data = demographic_data, glp_data = glp_data, 
             task_meta_data = task_meta_data, task_data = task_data, 
             questionnaire_data = questionnaire_data)

# (6) Save processed data  ---------------------------------------------------
setwd(here::here())
saveRDS(data, "github/semaglutide-study/data/processed_data/main_data.RDS")

# (7) Load groups of non-diabetic controls  ---------------------------------------------------

# Data from https://github.com/smehrhof/2024_effort_study

non_diabetic_data <- readRDS("github/semaglutide-study/data/raw_data/non_diabetic_groups/non_diabetic_matched.RDS")
non_diabetic_normal_weight_data <- readRDS("github/semaglutide-study/data/raw_data/non_diabetic_groups/non_diabetic_normal_weight_matched.RDS")
