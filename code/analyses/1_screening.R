########################################################################################################################
#######################---------------- SCREENING PARTICIPANTS FOR MAIN TESTING ----------------########################
########################################################################################################################

### In this script: 
# (1) Parse data
# (2) Exclusion
# (3) Identify treatment and control group
# (4) Adjust data based on Prolific messages
# (5) Check for unrealistic values
# (6) Save processed data
# (7) Randomize to testing schedule and identify testing days
# (8) Matching control participants

# Set working directory
here::i_am("github/semaglutide-study/code/analyses/1_screening.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/screener_parsing_fun.R")

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, writexl, readxl)

# (1) Parse data ---------------------------------------------------

raw_dat <- readRDS("github/semaglutide-study/data/raw_data/screening/raw_dat.RDS")

raw_meta_dat <- readRDS("github/semaglutide-study/data/raw_data/screening/raw_meta_dat.RDS")

screening_dat <- screener_parsing(raw_dat = raw_dat, 
                                  raw_meta_dat = raw_meta_dat, 
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

# All excluded except for 58f5dabf41 and 3c661644a3

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
     !.$prolific_id %in% c("58f5dabf41", "3c661644a3") ~ 1,
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

## Subject 02e632f266
# by chat reported to suffers from anxiety and made a mistake in IPAQ

screening_dat$screening_dat %<>% 
  mutate(psych_neurdev = ifelse(prolific_id == "02e632f266", 
                                1, psych_neurdev)) %>%
  mutate(psych_neurdev_condition = ifelse(prolific_id == "02e632f266", 
                                          "[\"Generalised anxiety disorder\",\"Current\"]", psych_neurdev_condition)) %>%
  mutate(ipaq_1 = ifelse(prolific_id == "02e632f266", 
                                  "None", ipaq_1)) %>% 
  mutate(ipaq_2 = ifelse(prolific_id == "02e632f266", 
                                  NA, ipaq_2)) %>% 
  mutate(ipaq_3 = ifelse(prolific_id == "02e632f266", 
                                  "5", ipaq_3)) %>% 
  mutate(ipaq_4 = ifelse(prolific_id == "02e632f266", 
                                  "0:40", ipaq_4)) %>% 
  mutate(ipaq_5 = ifelse(prolific_id == "02e632f266", 
                                  "7", ipaq_5)) %>% 
  mutate(ipaq_6 = ifelse(prolific_id == "02e632f266", 
                                  "2:0", ipaq_6)) %>%
  
  ## Subject 8415e96bc9
  # by chat reported to have made a mistake in IPAQ
  mutate(ipaq_1 = ifelse(prolific_id == "8415e96bc9", 
                                  "None", ipaq_1)) %>%
  mutate(ipaq_2 = ifelse(prolific_id == "8415e96bc9", 
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

# (6) Save processed data  ---------------------------------------------------
setwd(here::here())
saveRDS(screening_dat, "github/semaglutide-study/data/processed_data/screening_dat.RDS")

# (7) Randomize to testing schedule and identify testing days  ---------------------------------------------------

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

# (8) Matching control participants ---------------------------------------------------

# participants excluded from being controls: 
                    # past GLP: 
exclude_controls <- c("c5684e199b", "9585dd25f3", "f975555db0", 
                    # current GLP (but not Ozempic)
                      "52f6024277", "a196346d89", "4a06350ad5",
                    # returned
                    "e8e61a6fc1", "6d8624e15d", "11b6348241", 
                    "ba35fd1637", "99e64946fe", "2d260ef7e2",
                    # timed out
                    "11063adcaa", 
                    # failed catch questions
                    "5d463c442f", "5f85700310", "5bf6185c90",
                    
                    "93564808db", "95b5b01656", "42460ccf0b", 
                    "365621c1d2", "6d1631e969", "b7c60fe1b3"
                    )  

### US group
# Already matched participants
# B1
control_match_US_B1 <- read_excel(here::here("data_collection/main/controls_B1.xlsx"), sheet = "US controls")
control_match_US_B1 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id)) 

exclude_controls %<>%
  append(control_match_US_B1 %>%
           .$control_id)

# B2 
control_match_US_B2 <- read_excel(here::here("data_collection/main/controls_B2.xlsx"), sheet = "US controls")
control_match_US_B2 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id)) 

exclude_controls %<>%
  append(control_match_US_B2 %>% 
           .$control_id)

# B3 
control_match_US_B3 <- read_excel(here::here("data_collection/main/controls_B3.xlsx"), sheet = "US controls")
control_match_US_B3 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id)) 

exclude_controls %<>%
  append(control_match_US_B3 %>% 
           .$control_id)

# B4 
control_match_US_4 <- read_excel(here::here("data_collection/main/controls_B4.xlsx"), sheet = "US controls")
control_match_US_4 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id)) 

exclude_controls %<>%
  append(control_match_US_4 %>% 
           .$control_id)

# B5 
control_match_US_5 <- read_excel(here::here("data_collection/main/controls_B5.xlsx"), sheet = "US controls")
control_match_US_5 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id)) 

exclude_controls %<>%
  append(control_match_US_5 %>% 
           .$control_id)

# B6 
control_match_US_6 <- read_excel(here::here("data_collection/main/controls_B6.xlsx"), sheet = "US controls")
control_match_US_6 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id)) 

exclude_controls %<>%
  append(control_match_US_6 %>% 
           .$control_id)

# B7 
control_match_US_7 <- read_excel(here::here("data_collection/main/controls_B7.xlsx"), sheet = "US controls")
control_match_US_7 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_US_7 %>% 
           .$control_id)

# B8 - current
control_match_US <- read_excel(here::here("data_collection/main/controls_B8.xlsx"), sheet = "US controls")
control_match_US %<>%
  filter(treatment_id != "TOTAL") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id)) 

exclude_controls %<>%
  append(control_match_US %>% 
           .$control_id)


# participant 9955ce5d61 was collected as control participant but started
# Ozempic since screening => qualified as treatment group subject
#screening_dat$screening_dat %<>% 
#  mutate(group = case_when(prolific_id == "9955ce5d61" ~ "treatment", 
#                           .default = group))

US_match_dat <- screening_dat$screening_dat %>% 
  # US only
  filter(residence == "United States") %>%
  # Treatment and control groups
  filter(group != "none") %>% 
  # Only include treatment participants that have already been collected and that don't have a control yet
  filter((group == "treatment" & prolific_id %in% control_match_US$treatment_id[control_match_US$control_completed == 0]) |
           
  # Exclude control participants that have already been matched and collected
  # or are blocked because they previousely didn't respond to study or responded with non-eligibility
           (group == "control" & 
              !prolific_id %in% c(na.omit(control_match_US$control_id), exclude_controls) & 
              # no participants that have missing values
              !is.na(bmi) & !is.na(ipaq_sumScore))) %>%
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
# B1
control_match_UK_B1 <- read_excel(here::here("data_collection/main/controls_B1.xlsx"), sheet = "UK controls")
control_match_UK_B1 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK_B1 %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

# B2 
control_match_UK_B2 <- read_excel(here::here("data_collection/main/controls_B2.xlsx"), sheet = "UK controls")
control_match_UK_B2 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK_B2 %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

# B3 
control_match_UK_B3 <- read_excel(here::here("data_collection/main/controls_B3.xlsx"), sheet = "UK controls")
control_match_UK_B3 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK_B3 %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

# B4 
control_match_UK_4 <- read_excel(here::here("data_collection/main/controls_B4.xlsx"), sheet = "UK controls")
control_match_UK_4 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK_4 %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

# B5
control_match_UK_5 <- read_excel(here::here("data_collection/main/controls_B5.xlsx"), sheet = "UK controls")
control_match_UK_5 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK_5 %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

# B6 
control_match_UK_6 <- read_excel(here::here("data_collection/main/controls_B6.xlsx"), sheet = "UK controls")
control_match_UK_6 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK_6 %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

# B7 
control_match_UK_7 <- read_excel(here::here("data_collection/main/controls_B7.xlsx"), sheet = "UK controls")
control_match_UK_7 %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK_7 %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

# B8 - current
control_match_UK <- read_excel(here::here("data_collection/main/controls_B8.xlsx"), sheet = "UK controls")
control_match_UK %<>%
  filter(treatment_id != "TOTAL" & control_completed == 0 & control_id != "NA") %>%
  rowwise %>%
  mutate(treatment_id = id_shuffle(treatment_id)) %>%
  mutate(control_id = id_shuffle(control_id))

exclude_controls %<>%
  append(control_match_UK %>% 
           filter(control_completed == 0 & control_id != "NA") %>% 
           .$control_id)

UK_match_dat <- screening_dat$screening_dat %>% 
  # UK only
  filter(residence == "United Kingdom") %>%
  # Treatment and control groups
  filter(group != "none") %>% 
  # Only include treatment participants that have already been collected
  filter((group == "treatment" & prolific_id %in% control_match_UK$treatment_id[control_match_UK$control_completed == 0]) |
           
           # Exclude control participants that have already been matched and collected
           # or are blocked because they previousely didn't respond to study or responded with non-eligibility
           (group == "control" & 
              !prolific_id %in% c(na.omit(control_match_UK$control_id), exclude_controls) & 
              # no participants that have missing values
              !is.na(bmi) & !is.na(ipaq_sumScore))) %>%
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
           path=here::here("data_collection/main/controls_B8.xlsx"))

### Write text file to copy into Prolific
write.table(c(US_matched$match.matrix$control_id, 
              UK_matched$match.matrix$control_id),
            file = here::here("data_collection/main/controls_B8.txt"), 
            sep = ";", row.names = FALSE, col.names = FALSE, quote = FALSE)




