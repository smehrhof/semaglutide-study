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
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, writexl, readxl)


# (1) Parse data ---------------------------------------------------

data_files <- list.files(path = here::here("data/raw_data/screening"), 
                         pattern = ".txt", full.names = TRUE)
meta_files <- list.files(path = here::here("data/raw_data/screening"), 
                         pattern = ".csv", full.names = TRUE)

screening_dat <- screener_parsing(file = data_files, 
                                  meta_file = meta_files, 
                                  display_progress = TRUE)

# (2) Exclusion ---------------------------------------------------
 
# Neurological condition - any severe disorders? 
screening_dat$screening_dat %>%
  filter(neurological == 1) %>% 
  select(prolific_id, neurological_condition, neurological_condition_other) %>% 
  print(n = Inf, width = Inf)

# If subject has epilepsy, check if they are medicated
screening_dat$screening_dat %>% 
  add_column(epilepsy_check = case_when(
    grepl("Epilepsy", 
          .$neurological_condition, 
          ignore.case = TRUE) ~ 1
  )) %>% 
  filter(!is.na(epilepsy_check)) %>%
  select(prolific_id, neurological_condition, medication, other_medication, other_medication_type) %>% 
  print(n = Inf, width = Inf)

# All excluded except for 5db8a4f15f7e1a000a20d09c and 614c6a6334c83a831c18dd06

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
    grepl("Epilepsy", 
          .$neurological_condition, 
          ignore.case = TRUE) & 
     !.$prolific_id %in% c("5db8a4f15f7e1a000a20d09c", "614c6a6334c83a831c18dd06") ~ 1,
    .default = 0
  ))

# (3) Identify treatment and control group ---------------------------------------------------

# Add residence (US or UK) and screening date
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
  mutate(start_date_study = start_date_glp + weeks(4)) %>% 
  filter(start_date_study > now())


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
# Make completely unrealistic values NA
screening_dat$screening_dat %<>% 
  rowwise() %>%
  mutate(across(height_cm:bmi, ~ case_when(is.infinite(bmi) ~ NA,
                                           is.nan(bmi) ~ NA,
                                           bmi <= 12 ~ NA, 
                                           bmi >= 200 ~ NA, 
                                           height_cm > 250 ~ NA,
                                           height_cm < 140 ~ NA,
                                           weight_kg < 20 ~ NA,
                                           weight_kg > 600 ~ NA,
                                           .default = .))) %>%
  # to end durn rowws_df back into normal tibble
  ungroup()

# Flag participants that have are underweight or have a BMI > 2 SD from sample mean
bmi_flag <- screening_dat$screening_dat %>% 
  select(prolific_id, height_cm, weight_kg, bmi, height_cm_raw, weight_kg_raw, group) %>% 
  na.omit() %>%
  filter(bmi < 18.5 | bmi > (mean(.$bmi)+2*sd(.$bmi))) %>% 
  .$prolific_id

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

# Check if any treatment group participants have NA IPAQ data
screening_dat$screening_dat %>%
  filter(group == "treatment") %>% 
  select(ipaq_1:ipaq_sumScore) %>% 
  print(n=Inf)


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
  )) %>% arrange(screening_day) %>% print(n = Inf)

# Are any treatment participants flagged for BMI? 
screening_dat$screening_dat %>% 
  filter(group == "treatment") %>% 
  filter(prolific_id %in% bmi_flag)

# (7) Matching control participants ---------------------------------------------------

### US group
# Already matched participants
control_match_US <- read_excel(here::here("data_collection/main/controls_B1.xlsx"), sheet = "US controls")
control_match_US %<>%
  filter(treatment_id != "TOTAL") 

US_match_dat <- screening_dat$screening_dat %>% 
  # US only
  filter(residence == "United States") %>%
  # Treatment and control groups
  filter(group != "none") %>% 
  # Only include treatment participants that have already been collected
  filter((group == "treatment" & prolific_id %in% control_match_US$treatment_id[control_match_US$control_completed == 0]) |
           
  # Exclude control participants that have already been matched and collected
           (group == "control" & !prolific_id %in% na.omit(control_match_US$control_id))) %>%
  # Variables to screen by
  select(prolific_id, age, gender, bmi, ipaq_sumScore, group)

US_match_dat %<>% 
  # Make group variable binary
  mutate(group = { ifelse(group == "treatment", 1, 0) }) %>% 
  # Use sex for non-binary gender, due to low numbers
  left_join(screening_dat$prolific_dat %>% select(prolific_id, sex),
            join_by(prolific_id == prolific_id)) %>% 
  mutate(gender = { ifelse(gender == "Non-binary", 
                           sex, 
                           gender) }) %>%
  # Make gender variable factor
  mutate(gender = as.factor(gender)) %>% 
  # Make prolific_if rownames
  column_to_rownames(var = "prolific_id") %>% 
  # Exclude participants that have NAs on any matching variables
  na.omit()

# matching
US_matched <- matchit(group ~ age + gender + bmi + ipaq_sumScore, data = US_match_dat,
                   method = "nearest", distance = "glm", link = "probit")
summary(US_matched, un = FALSE)

# Matches
US_matched$match.matrix %<>% 
  .[,1] %>% 
  enframe(name = "treatment_id", value = "control_id") 

# Add new matches to tibble
control_match_US %<>% 
  mutate(control_id = ifelse(control_completed == 1, control_id, NA)) %>% 
  left_join(US_matched$match.matrix, by = "treatment_id") %>% 
  unite("control_id", control_id.x, control_id.y, na.rm = TRUE, remove = FALSE) %>% 
  select(-c(control_id.x, control_id.y))

### UK group
# Already matched participants
control_match_UK <- read_excel(here::here("data_collection/main/controls_B1.xlsx"), sheet = "UK controls")
control_match_UK %<>%
  filter(treatment_id != "TOTAL") 

UK_match_dat <- screening_dat$screening_dat %>% 
  # UK only
  filter(residence == "United Kingdom") %>%
  # Treatment and control groups
  filter(group != "none") %>% 
  # Only include treatment participants that have already been collected
  filter((group == "treatment" & prolific_id %in% control_match_UK$treatment_id[control_match_UK$control_completed == 0]) |
           
           # Exclude control participants that have already been matched and collected
           (group == "control" & !prolific_id %in% na.omit(control_match_UK$control_id))) %>%
  # Variables to screen by
  select(prolific_id, age, gender, bmi, ipaq_sumScore, group)

UK_match_dat %<>% 
  # Make group variable binary
  mutate(group = { ifelse(group == "treatment", 1, 0) }) %>% 
  # Use sex for non-binary gender, due to low numbers
  left_join(screening_dat$prolific_dat %>% select(prolific_id, sex),
            join_by(prolific_id == prolific_id)) %>% 
  mutate(gender = { ifelse(gender == "Non-binary", 
                           sex, 
                           gender) }) %>%
  # Make gender variable factor
  mutate(gender = as.factor(gender)) %>% 
  # Make prolific_if rownames
  column_to_rownames(var = "prolific_id") %>% 
  # Exclude participants that have NAs on any matching variables
  na.omit()

# matching
UK_matched <- matchit(group ~ age + gender + bmi + ipaq_sumScore, data = UK_match_dat,
                      method = "nearest", distance = "glm", link = "probit")
summary(UK_matched, un = FALSE)

# Matches
UK_matched$match.matrix %<>% 
  .[,1] %>% 
  enframe(name = "treatment_id", value = "control_id") 

# Add new matches to tibble
control_match_UK %<>% 
  mutate(control_id = ifelse(control_completed == 1, control_id, NA)) %>% 
  left_join(UK_matched$match.matrix, by = "treatment_id") %>% 
  unite("control_id", control_id.x, control_id.y, na.rm = TRUE, remove = FALSE) %>% 
  select(-c(control_id.x, control_id.y))

### Write new excel file
write_xlsx(setNames(list(control_match_US, control_match_UK), 
                    c("US controls", "UK controls")), 
           path=here::here("data_collection/main/controls_B2.xlsx"))

### Write text file to copy into Prolific
write.table(c(US_matched$match.matrix$control_id, 
              UK_matched$match.matrix$control_id),
            file = here::here("data_collection/main/controls_B2.txt"), 
            sep = ";", row.names = FALSE, col.names = FALSE, quote = FALSE)


