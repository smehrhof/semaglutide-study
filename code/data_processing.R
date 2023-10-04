##################################################################################################################
#######################---------------- PROCESSING DATA FROM MAIN TESTING ----------------########################
##################################################################################################################

# DATA LAST UPDATED: 30.08.23
# s1 = 31
# s2 = 21

### In this script: 
# (1) Parse data
# (2) Exclusion
# (3) Adjust data based on Prolific messages
# (4) Merge with screening data
# (5) Estimate hours since injection


# Set working directory
here::i_am("github/semaglutide-study/code/data_processing.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/parsing_fun.R")

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, writexl)


# (1) Parse data ---------------------------------------------------

### Session 1

data_files <- list.files(path = here::here("data/raw_data/main/session_1"), 
                         pattern = ".txt", full.names = TRUE)

doses_file <- list.files(path = here::here("data/raw_data/main"), 
                         pattern = "doses.csv", full.names = TRUE)

meta_files <- list.files(path = here::here("data/raw_data/main/session_1"), 
                         pattern = ".csv", full.names = TRUE)

s1_data <- parsing(file = data_files, 
                   meta_file = meta_files, 
                   doses_file = doses_file, 
                   session = "s1", 
                   display_progress = TRUE)

### Session 2

data_files <- list.files(path = here::here("data/raw_data/main/session_2"), 
                         pattern = ".txt", full.names = TRUE)

doses_file <- list.files(path = here::here("data/raw_data/main"), 
                         pattern = "doses.csv", full.names = TRUE)

meta_files <- list.files(path = here::here("data/raw_data/main/session_2"), 
                         pattern = ".csv", full.names = TRUE)

s2_data <- parsing(file = data_files, 
                   meta_file = meta_files, 
                   doses_file = doses_file, 
                   session = "s2", 
                   display_progress = TRUE)


# (2) Exclusion ---------------------------------------------------

### Session 1
s1_exclusion <- c()

# Treatment group conditions: not on treatment or not diabetic
s1_exclusion <- append(s1_exclusion, 
                       s1_data$demographic_dat %>% filter(glp_check == 0 |
                                                            prolific_id == "60ff243a844cae10907dfd18") %>% 
                         .$prolific_id)

# Task based
# no offers accepted
s1_exclusion <- append(s1_exclusion, 
                       s1_data$task_dat$task_data %>% 
                         filter(phase == "game") %>% 
                         reframe(acceptance_rate = mean(choice), .by = prolific_id) %>% 
                         filter(acceptance_rate == 0) %>% 
                         .$prolific_id)

# large difference between minimum and maximum clicking speed during calibration
# -> >3SD based on large sample: 128.9755
s1_exclusion <- append(s1_exclusion, 
                       s1_data$task_dat$task_data %>% 
                         filter(phase == "calibration" & 
                                  trial != 1) %>% 
                         reframe(calibration_difference = diff(clicks) %>% abs(), 
                                 .by = prolific_id) %>% 
                         filter(calibration_difference > 128.9755) %>% 
                         .$prolific_id)

# large change in clicking calibration pre- to post-task
# -> >3SD based on large sample: 133.7341
s1_exclusion <- append(s1_exclusion, 
                       s1_data$task_dat$task_data %>% 
                         filter(phase == "calibration") %>% 
                         reframe(pre_max = max(clicks[1:3]),
                                 post_max = clicks[4],
                                 .by = prolific_id) %>% 
                         mutate(calibration_change = abs(pre_max- post_max)) %>%
                         filter(calibration_change > 133.7341) %>% 
                         .$prolific_id)

# Questionnaire based: catch questions
s1_exclusion <- append(s1_exclusion, 
                       s1_data$questionnaire_dat %>% 
                         filter(catch_questions_pass == 0) %>% 
                         .$prolific_id)

### Session 2
s2_exclusion <- c()

# Treatment group conditions: not on treatment or not diabetic
s2_exclusion <- append(s2_exclusion, 
                       s2_data$demographic_dat %>% filter(glp_check == 0 |
                                                            prolific_id == "60ff243a844cae10907dfd18") %>% 
                         .$prolific_id)

# Task based
# no offers accepted
s2_exclusion <- append(s2_exclusion, 
                       s2_data$task_dat$task_data %>% 
                         filter(phase == "game") %>% 
                         reframe(acceptance_rate = mean(choice), .by = prolific_id) %>% 
                         filter(acceptance_rate == 0) %>% 
                         .$prolific_id)

# large difference between minimum and maximum clicking speed during calibration
# -> >3SD based on large sample: 128.9755
s2_exclusion <- append(s2_exclusion, 
                       s2_data$task_dat$task_data %>% 
                         filter(phase == "calibration" & 
                                  trial != 1) %>% 
                         reframe(calibration_difference = diff(clicks) %>% abs(), 
                                 .by = prolific_id) %>% 
                         filter(calibration_difference > 128.9755) %>% 
                         .$prolific_id)

# large change in clicking calibration pre- to post-task
# -> >3SD based on large sample: 133.7341
s2_exclusion <- append(s2_exclusion, 
                       s2_data$task_dat$task_data %>% 
                         filter(phase == "calibration") %>% 
                         reframe(pre_max = max(clicks[1:3]),
                                 post_max = clicks[4],
                                 .by = prolific_id) %>% 
                         mutate(calibration_change = abs(pre_max- post_max)) %>%
                         filter(calibration_change > 133.7341) %>% 
                         .$prolific_id)

# Questionnaire based: catch questions
s2_exclusion <- append(s2_exclusion, 
                       s2_data$questionnaire_dat %>% 
                         filter(catch_questions_pass == 0) %>% 
                         .$prolific_id)

# (3) Adjust data based on Prolific messages  ---------------------------------------------------

# subjects 60ff243a844cae10907dfd18, 63f7e098e99a7db2d59ac6ce, and 603e56588fd0a9ea9f3c710a 
# corrected "days since last injection"
s1_data$demographic_dat %<>% 
  mutate(across(last_injection_days, ~ case_when(prolific_id == "60ff243a844cae10907dfd18" ~ 6, 
                                                 prolific_id == "63f7e098e99a7db2d59ac6ce" ~ 1, 
                                                 prolific_id == "603e56588fd0a9ea9f3c710a" ~ 1, 
                                                 .default = .)))



# estimate hours since injection

s1_data$demographic_dat %>% 
  filter(!prolific_id %in% s1_exclusion) %>% 
  janitor::tabyl(last_injection_days)

s1_data$demographic_dat %>% 
  filter(!prolific_id %in% s1_exclusion) %>% 
  select(prolific_id, last_injection_days, last_injection_time, time_response) %>% 
  print(n = Inf)






