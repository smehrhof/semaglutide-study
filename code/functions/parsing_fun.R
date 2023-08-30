########################################################################################################################
############################---------------- PARSING SCREENER DATA FUNCTION ----------------############################
########################################################################################################################

# Set working directory
here::i_am("github/semaglutide-study/code/functions/parsing_fun.R")
setwd(here::here())

file <- list.files(path = here::here("data/raw_data/main/session_2"), 
                         pattern = ".txt", full.names = TRUE)

doses_file <- list.files(path = here::here("data/raw_data/main"), 
                   pattern = "doses.csv", full.names = TRUE)

meta_file <- list.files(path = here::here("data/raw_data/main/session_2"), 
                        pattern = ".csv", full.names = TRUE)

s2_data <- parsing(file = file, 
                   meta_file = meta_file, 
                   doses_file = doses_file, 
                   session = "s2", 
                   display_progress = TRUE)

s1_data$demographic_dat

# TO DO 
# exclude 60ff243a844cae10907dfd18
# exclude 64415cccab771017a37d3d93
  
## Parsing data from screening -----------------------------------------
# @ file: vector of relative path(s) to .txt data file(s) or existing R object
# @ meta_file: vector of relative path(s) to .csv data file(s) with Prolific meta data
  
parsing <- function(file,
                    meta_file, 
                    doses_file, 
                    session = "s1", 
                    display_progress = FALSE){
  
  suppressMessages({
    
    # Required packages
    require(tidyverse)
    require(jsonlite)
    require(MatchIt)
    require(magrittr)
    require(lubridate)
    
    ### Load raw JSON data -----------------------------------------------------------------------------------------------
    if(length(file) == 1){
      raw_dat <- jsonlite::fromJSON(sprintf('[%s]', paste(readLines(file, encoding="UTF-8", warn=F), collapse = ',')))
    } else if(length(file) > 1){
      raw_dat <- list()
      for(i in 1:length(file)){
        raw_dat <- append(raw_dat,  jsonlite::fromJSON(sprintf('[%s]', paste(readLines(file[i], encoding="UTF-8", warn=F), collapse = ','))))
      }
    }
    
    ### Load Prolific meta data ------------------------------------------------------------------------------------------
    if(length(meta_file) == 1){
      raw_meta_dat <- read.csv(meta_file)
    } else if(length(meta_file) > 1){
      raw_meta_dat <- matrix(ncol = length(read.csv(meta_file[1])), nrow = 0)
      colnames(raw_meta_dat) <-  colnames(read.csv(meta_file[1]))
      
      for(i in 1:length(meta_file)){
        raw_meta_dat <- rbind(raw_meta_dat, read.csv(meta_file[i]))
      }
    }
    
    ### Load file with extra dose data ------------------------------------------------------------------------------------------
    raw_doses_dat <- read.csv(doses_file, sep = ";") %>% as_tibble()
    
    
    ### Identify component indices ---------------------------------------------------------------------------------------
    comp_1_id <- raw_dat %>%
      lmap(\(x) list(any(grepl("Welcome to the experiment", x, fixed=T)))) %>%
      unlist(.) %>% which(. == TRUE)
    comp_2_id <- raw_dat %>%
      lmap(\(x) list(any(grepl("This is the start of the game", x, fixed=T)))) %>%
      unlist(.) %>% which(. == TRUE)
    comp_3_id <- raw_dat %>%
      lmap(\(x) list(any(grepl("Well done for completing the game!", x, fixed=T)))) %>%
      unlist(.) %>% which(. == TRUE)
    comp_4_id <- raw_dat %>%
      lmap(\(x) list(any(grepl("Please think back to the last", x, fixed=T)))) %>%
      unlist(.) %>% which(. == TRUE)
    
    
    ### Prolific data -------------------------------------------------------------------------------------
    # Contains Prolific meta data
    
    # Only approved participants
    meta_data <- raw_meta_dat %>%
      filter(Status == "APPROVED") %>%
      select(Participant.id, Age, Sex, Ethnicity.simplified, Country.of.residence, Started.at, Completed.at) %>%
      as_tibble() %>% 
      rename(prolific_id = Participant.id, age = Age, sex = Sex, 
             ethnicity = Ethnicity.simplified, residence = Country.of.residence, 
             study_start = Started.at, study_end = Completed.at) %>%
      mutate(study_start = as_datetime(study_start)) %>%
      mutate(study_end = as_datetime(study_end)) 
 
    # List of prolific IDs for bug fixing
    if(session == "s1"){
      bug_fix_ids <- meta_data %>% 
        filter(study_start < as_datetime("2023-07-28 00:00:00")) %>% 
        select(prolific_id) %>% as_vector() 
    }
    
    ### COMPONENT 1: Medical data -------------------------------------------------------------------------------------------------
    # Contains data on medication data (Semaglutide related and change in diabetes and other daily medication)
    # as well as local time and question of computer games
    
    comp_1_dat <- tibble()
    
    # Progress bar
    if(display_progress == TRUE){
      print("Parsing component 1/3:")
      progress_bar = txtProgressBar(min=0, max=length(comp_1_id), style = 1, char="=")
    }
    
    for(subj in comp_1_id){
      
      # Unconditional questions
      comp_1_dat_subj <- raw_dat[[subj]] %>%
        { if(session == "s1"){
          select(., c(prolific_id, glp_check, diabetes_med, other_med, games, time_response)) 
        } else if(session == "s2"){
          select(., c(prolific_id, glp_check, time_response)) 
        } } %>% 
        summarise(across(prolific_id:time_response, .fns=~na.omit(unique(.)))) %>%
        as_tibble()

      
        # Local time big fix
        if(comp_1_dat_subj$prolific_id %in% bug_fix_ids){
          comp_1_dat_subj %<>% 
            mutate(time_response = {if(raw_dat[[subj]] %>% 
                                       select(responses) %>% 
                                       na.omit() %>% slice_tail(n = 1) %>% 
                                       str_detect(., '"time_pm"')){
              # convert to pm on 24:00
              format((strptime(paste(time_response, "pm"), "%I:%M %p")), format = "%H:%M")
            }
            })
        }

      # GLP follow up questions
      comp_1_dat_subj %<>% 
        add_column(
          glp_dose_mg = { ifelse(.$glp_check == 1,
                                 { ifelse(.$prolific_id %in% raw_doses_dat$prolific_id, 
                                          raw_doses_dat %>% 
                                            filter(prolific_id == comp_1_dat_subj$prolific_id) %>% 
                                            select(all_of(session)) %>% 
                                            as.character(), 
                                          na.omit(raw_dat[[subj]]["glp_dose"])[1,1])},
                              NA) }, 
          last_injection_days = { ifelse(.$glp_check == 1, 
                                             na.omit(raw_dat[[subj]]["last_injection"])[1,1], 
                                             NA) }, 
          last_injection_time = { ifelse(.$glp_check == 1, 
                                         case_when(na.omit(raw_dat[[subj]]["last_injection_time"]) == 0 ~ "before 6am", 
                                                   na.omit(raw_dat[[subj]]["last_injection_time"]) == 1 ~ "6am - 10am", 
                                                   na.omit(raw_dat[[subj]]["last_injection_time"]) == 2 ~ "10am - 2pm", 
                                                   na.omit(raw_dat[[subj]]["last_injection_time"]) == 3 ~ "2pm - 6pm", 
                                                   na.omit(raw_dat[[subj]]["last_injection_time"]) == 4 ~ "6pm - 10pm", 
                                                   na.omit(raw_dat[[subj]]["last_injection_time"]) == 5 ~ "after 10pm"),
                                         NA) }, 
                   .after = "glp_check") %>% 
        add_column(
          glp_past = { ifelse(.$glp_check == 0, 
                              # fix bug in s2 js (glp_past response under glp_check)
                              { ifelse(session == "s2", 
                                       na.omit(raw_dat[[subj]]["glp_check"])[2,1],
                                       na.omit(raw_dat[[subj]]["glp_past"])[1,1]) },
                                 NA) }, 
          .after = "glp_check") %>% 
        add_column(
          glp_last_date = { ifelse(.$glp_past == 1, 
                              na.omit(raw_dat[[subj]]["glp_last_date"])[1,1], 
                              NA) }, 
          .after = "glp_past")
      
      # Medication change follow up questions
      if(session == "s1"){
        comp_1_dat_subj %<>% 
          # diabetes medication
          mutate(diabetes_med = ifelse(.$diabetes_med == 1, 
                                       ifelse(any(str_detect(na.omit(raw_dat[[subj]]["diabetes_med"])$diabetes_med, 
                                                             "None")), 
                                              0, 1),
                                       0)) %>%
          add_column(
            diabetes_med_type = { ifelse(.$diabetes_med == 1, 
                                         paste(
                                           unlist(
                                             jsonlite::parse_json(
                                               na.omit(raw_dat[[subj]]["diabetes_med_type"])[1,1])), 
                                           collapse=', '),
                                         NA) }, 
            diabetes_med_type_other = { ifelse(.$diabetes_med == 1, 
                                               ifelse(na.omit(raw_dat[[subj]]["diabetes_med_type_other"]), 
                                                      na.omit(raw_dat[[subj]]$diabetes_med_type_other_text), 
                                                      NA), 
                                               NA) },
            .after = "diabetes_med") %>%
          # other daily medication
          add_column(
            other_med_type = { ifelse(.$other_med == 1, 
                                      paste(
                                        unlist(
                                          jsonlite::parse_json(
                                            # due to js error of missing comma separation
                                            na.omit(raw_dat[[subj]]["other_med_type"])[1,1] %>% 
                                              str_replace_all(., '\"\"', '\",\"'))), 
                                        collapse=', '),
                                      NA) }, 
            other_med_text = { ifelse(.$other_med == 1, 
                                      paste(
                                        unlist(
                                          jsonlite::parse_json(
                                            # due to js error of missing comma separation
                                            na.omit(raw_dat[[subj]]["other_med_text"])[1,1] %>% 
                                              str_replace_all(., '\"\"', '\",\"'))), 
                                        collapse=', '),
                                      NA) }, 
            .after = "other_med")  
      }
      
      # Combine
      comp_1_dat %<>%
        bind_rows(., comp_1_dat_subj)
      rm(comp_1_dat_subj)
      
      # Progress bar
      if(display_progress == TRUE){
        setTxtProgressBar(progress_bar, value = subj)
      }
    }
    
    # Progress bar
    if(display_progress == TRUE){
      close(progress_bar)
    }
    
    # remove duplicates
    comp_1_dat <- comp_1_dat[!duplicated(comp_1_dat$prolific_id, fromLast = TRUE), ]
    
    
    ### COMPONENT 2: Task data -------------------------------------------------------------------------------------
    # Contains meta and main task data
    
    comp_2_dat_meta <- comp_2_dat_task <- tibble()

    if(display_progress == TRUE){
      print("Parsing component 2/3:")
      progress_bar = txtProgressBar(min=0, max=length(comp_2_id), style = 1, char="=")
    }
    
    for(subj in comp_2_id){
      
      # meta data
      comp_2_dat_meta_subj <- tibble(prolific_id = raw_dat[[subj-1]] %>%
                                       select(prolific_id) %>%
                                       unique(), 
                                     game_id = raw_dat[[subj]]$subjID,
                                     start_time = raw_dat[[subj]]$date[1] %>% as_datetime(),
                                     end_time = raw_dat[[subj]]$date[2] %>% as_datetime()) %>%
        mutate(completion_time = difftime(start_time, end_time) %>% abs())

      
      # task data
      
      # avoid error due to double-registration of click
      if(raw_dat[[subj]]$rt %>% length != raw_dat[[subj]]$phase %>% length){
        raw_dat[[subj]][c("choice", "rt")] %<>% 
          as_tibble() %>% 
          filter(rt == 999 |
                   abs(rt - lag(rt)) > 5 |
                   is.na(abs(rt - lag(rt)))) %>% 
          as.list() 
      } 
        
      # adjust variables if subject had to re-do calibration
      if(raw_dat[[subj]]$clicks %>% length > raw_dat[[subj]]$rt %>% length){
        
        for(i in 1:(abs(raw_dat[[subj]]$clicks %>% length - raw_dat[[subj]]$rt %>% length) / 2)){
          
          if(mean(raw_dat[[subj]]$clicks[(i+1):(i+2)]) < 7){
            
            raw_dat[[subj]][4:13]$phase %<>% 
              append(c("calibration", "calibration"), after = i+2)
            raw_dat[[subj]][4:13]$trial %<>% 
              append(2:3, after = i+2)
            raw_dat[[subj]][4:13]$points %<>% 
              append(c(0, 0), after = i+2)
            
            
            raw_dat[[subj]][4:13][c("trialType", "offerEffort", "offerReward", "choice", "rt", "goalClicks")] %<>% 
              map(\(x) append(x, c(999, 999), after = i+2))
            
          } else {
            stop(paste("Check task data for subject with game ID", raw_dat[[subj]]$subjID, sep = " "))
          }
        }
      }
      
      comp_2_dat_task_subj <- raw_dat[[subj]][4:13] %>% 
        as_tibble() %>%
        rename(trial_type = trialType, effort = offerEffort, reward = offerReward, goal_clicks = goalClicks) %>%
        mutate(across(trial_type:goal_clicks, ~na_if(., 999))) %>% 
        add_column(
          prolific_id = raw_dat[[subj-1]] %>%
            select(prolific_id) %>%
            unique(),
          .before = "phase"
        )
    

      # Combine
      comp_2_dat_meta %<>%
        bind_rows(., comp_2_dat_meta_subj)
      rm(comp_2_dat_meta_subj)
      
      comp_2_dat_task %<>%
        bind_rows(., comp_2_dat_task_subj)
      rm(comp_2_dat_task_subj)
      
      if(display_progress == TRUE){
        setTxtProgressBar(progress_bar, value = subj)
      }
    }
    
    comp_2_dat <- list("task_meta_data" = comp_2_dat_meta, 
                       "task_data" = comp_2_dat_task)
    
    if(display_progress == TRUE){
      close(progress_bar)
    }
    
    
    ### COMPONENT 3: Questionnaire data -------------------------------------------------------------------------------------
    # Contains data on all questionnaires
    
    comp_3_dat <- tibble()
    
    if(display_progress == TRUE){
      print("Parsing component 3/3:")
      progress_bar = txtProgressBar(min = 0, max = length(comp_3_id), style = 1, char = "=")
    }
    
    if(session == "s1"){
      # To estimate delay discounting paramerer k for MCQ
      reward_ll <- c(55, 75, 25, 85, 25, 50, 35, 60, 80, 55, 
                     30, 75, 35, 50, 85, 60, 85, 35, 80, 30,
                     50, 30, 75, 60, 80, 25, 55)
      reward_ss <- c(54, 55, 19, 31, 14, 47, 15, 25, 78, 40,
                     11, 67, 34, 27, 69, 49, 80, 24, 33, 28,
                     34, 25, 41, 54, 54, 22, 20) 
    }
    
    for(subj in comp_3_id){
      
      comp_3_dat_subj <- raw_dat[[subj]] %>%
        select(prolific_id) %>% 
        summarise(across(prolific_id, .fns=~na.omit(unique(.)))) %>%
        as_tibble() 
    
      if(session == "s1"){
        
        comp_3_dat_subj <- bind_cols(
          
        comp_3_dat_subj, 
        
        # AES
        raw_dat[[subj]] %>% 
          filter(questionnaire == "AES") %>% 
          select(aes_response) %>% 
          data.frame(variables = paste("aes", 1:18, sep = "_"), .) %>% 
          pivot_wider(names_from = variables, values_from = aes_response) %>%
          mutate(aes_sumScore = rowSums(.)),
        
        # BDI 
        raw_dat[[subj]] %>% 
          filter(questionnaire == "BDI") %>% 
          select(bdi_response) %>% 
          mutate(bdi_response = as.integer(bdi_response)) %>% 
          data.frame(variables = paste("bdi", 1:21, sep = "_"), .) %>% 
          pivot_wider(names_from = variables, values_from = bdi_response) %>%
          mutate(bdi_sumScore = rowSums(.)),
        
        # FINDRISC
        raw_dat[[subj]] %>% 
          filter(questionnaire == "FINDRISC") %>% 
          select(findrisc_response) %>% 
          .$findrisc_response %>% unlist() %>% as.integer() %>% 
          # fix js bug 
          { if(length(.) == 8) {
            data.frame(variables = paste("findrisc", 1:8, sep = "_"), findrisc_response = .)
          } else {
            data.frame(variables = paste("findrisc", 1:8, sep = "_"), findrisc_response = .[c(1, 4:10)])
          }} %>% 
          pivot_wider(names_from = variables, values_from = findrisc_response) %>%
          mutate(findrisc_sumScore = rowSums(.)), 
        
        # MCQ
        raw_dat[[subj]] %>% 
          filter(questionnaire == "MCQ") %>% 
          select(mcq_response) %>% 
          mutate(mcq_response = as.integer(mcq_response)) %>% 
          data.frame(variables = paste("mcq", 1:27, sep = "_"), .) %>% 
          pivot_wider(names_from = variables, values_from = mcq_response), 
        # fit GLM to compute dd parameter
        suppressWarnings({
          raw_dat[[subj]] %>% 
            filter(questionnaire == "MCQ") %>% 
            select(mcq_response) %>% 
            mutate(mcq_response = as.integer(mcq_response)) %>% 
            add_column(ratio = 1 - (1 / (reward_ss / reward_ll))) %>%
            add_column(delay = c(117, 61, 53, 7, 19, 160, 13, 14, 162, 62,
                                 7, 119, 186, 21, 91, 89, 157, 29, 14, 179, 
                                 30, 80, 20, 111, 30, 136, 7)) %>% 
            glm(mcq_response ~ 0 + ratio + delay,
                data = ., 
                family = "binomial") %>% 
            tidy() %>% as_tibble() %>% 
            mutate(mcq_discounting_rate = estimate[2] / estimate[1]) %>% 
            filter(term == "ratio") %>% 
            select(mcq_discounting_rate)
        }),
        
        # MCTQ
        raw_dat[[subj]] %>% 
          filter(questionnaire == "MCTQ") %>% 
          select(question_no, mctq_response) %>% 
          # remove corrected answers
          filter(!duplicated(.[["question_no"]], fromLast = TRUE)) %>%
          select(!question_no) %>% 
          data.frame(variables = paste("mctq", 1:16, sep = "_"), .) %>% 
          pivot_wider(names_from = variables, values_from = mctq_response) %>% 
          
          # Sleep onset: "get ready to fall asleep" + "time to fall asleep"
          mutate(mctq_SO_w = format(as.POSIXct(parse_date_time(mctq_3, "HM") + minutes(mctq_4)), 
                                    format = "%H:%M"), 
                 mctq_SO_f = format(as.POSIXct(parse_date_time(mctq_10, "HM") + minutes(mctq_11)), 
                                    format = "%H:%M"), 
                 
                 # Getting out of bed: "waking up" + sleep inertia
                 mctq_GU_w = format(as.POSIXct(parse_date_time(mctq_5, "HM") + minutes(mctq_6)), 
                                    format = "%H:%M"),
                 mctq_GU_f = format(as.POSIXct(parse_date_time(mctq_12, "HM") + minutes(mctq_13)), 
                                    format = "%H:%M")
          ) %>% 
          
          # Sleep duration: time difference between "waking up" and sleep onset
          mutate(mctq_SD_w = { if((parse_date_time(mctq_5, "HM") - parse_date_time(mctq_SO_w, "HM")) < 0) {
            24 + (parse_date_time(mctq_5, "HM") - parse_date_time(mctq_SO_w, "HM"))
          } else {
            (parse_date_time(mctq_5, "HM") - parse_date_time(mctq_SO_w, "HM"))
          }}) %>% 
          mutate(mctq_SD_f = { if((parse_date_time(mctq_12, "HM") - parse_date_time(mctq_SO_f, "HM")) < 0) {
            24 + (parse_date_time(mctq_12, "HM") - parse_date_time(mctq_SO_f, "HM"))
          } else {
            (parse_date_time(mctq_12, "HM") - parse_date_time(mctq_SO_f, "HM"))
          }}) %>% 
          
          # Total time in bed: time difference between "getting up" and "going to bed" 
          mutate(mctq_TBT_w = { if((parse_date_time(mctq_GU_w, "HM") - parse_date_time(mctq_2, "HM")) < 0) {
            24 + (parse_date_time(mctq_GU_w, "HM") - parse_date_time(mctq_2, "HM"))
          } else {
            (parse_date_time(mctq_GU_w, "HM") - parse_date_time(mctq_2, "HM"))
          }}) %>% 
          mutate(mctq_TBT_f = { if((parse_date_time(mctq_GU_f, "HM") - parse_date_time(mctq_9, "HM")) < 0) {
            24 + (parse_date_time(mctq_GU_f, "HM") - parse_date_time(mctq_9, "HM"))
          } else {
            (parse_date_time(mctq_GU_f, "HM") - parse_date_time(mctq_9, "HM"))
          }}) %>% 
          
          # Mid point of sleep: sleep onset + half of sleep duration
          mutate(mctq_MS_w = format(as.POSIXct((parse_date_time(mctq_SO_w, "HM") + (mctq_SD_w/2))),
                                    format = "%H:%M"), 
                 mctq_MS_f = format(as.POSIXct((parse_date_time(mctq_SO_f, "HM") + (mctq_SD_f/2))), 
                                    format = "%H:%M")
          ) %>% 
          
          # Weekly average sleep duration: work / free day sleep duration weighted by the number of work / free days
          mutate(mctq_SD_week = ((mctq_SD_w * as.numeric(mctq_1)) + (mctq_SD_f * (7-as.numeric(mctq_1))) / 7)
          ) %>% 
          
          # Chronotype measure: 
          # Equals mid point of sleep on free days if free day sleep duration shorter than on work days
          mutate(mctq_MSF_SC = { if(mctq_SD_f <= mctq_SD_w) {
            mctq_MS_f
          } else {
            # Otherwise is corrected by difference in sleep duration between free and work days
            format(as.POSIXct(parse_date_time(mctq_MS_f, "HM") - ((mctq_SD_f - mctq_SD_w)/2)), 
                   format = "%H:%M")
          }}),
        
        
        # MEQ
        raw_dat[[subj]] %>% 
          filter(questionnaire == "MEQ") %>% 
          select(meq_response) %>% 
          data.frame(variables = paste("meq", 1:19, sep = "_"), .) %>% 
          pivot_wider(names_from = variables, values_from = meq_response) %>%
          mutate(meq_sumScore = rowSums(.)), 
        
        # OCIR
        raw_dat[[subj]] %>% 
          filter(questionnaire == "OCIR") %>% 
          select(ocir_response) %>% 
          data.frame(variables = paste("ocir", 1:18, sep = "_"), .) %>% 
          pivot_wider(names_from = variables, values_from = ocir_response) %>%
          mutate(ocir_sumScore = rowSums(.))
        )
      }
      
      comp_3_dat_subj <- bind_cols(
        
        comp_3_dat_subj, 
      
      # DAQ
      raw_dat[[subj]] %>% 
        filter(questionnaire == "DAQ") %>% 
        # fix js bug 
        { if(comp_3_dat_subj$prolific_id %in% bug_fix_ids) {
          select(., aes_response) %>%
          rename(., daq_response = aes_response)
          } else {
            select(., daq_response)
            }
          } %>% 
        data.frame(variables = paste("daq", 1:6, sep = "_"), .) %>% 
        # because DAQ has range 1-7, not 0-6
        mutate(daq_response = daq_response + 1) %>% 
        pivot_wider(names_from = variables, values_from = daq_response) %>%
        mutate(daq_sumScore = rowSums(.)), 
      
      # SHAPS
      raw_dat[[subj]] %>% 
        filter(questionnaire == "SHAPS") %>% 
        select(shaps_response) %>% 
        data.frame(variables = paste("shaps", 1:14, sep = "_"), .) %>% 
        pivot_wider(names_from = variables, values_from = shaps_response) %>%
        mutate(shaps_sumScore = rowSums(.)), 
      
      # TFEQ
      raw_dat[[subj]] %>% 
        filter(questionnaire == "TFEQ") %>% 
        select(tfeq_response) %>% 
        data.frame(variables = paste("tfeq", 1:18, sep = "_"), .) %>% 
        # because TFEQ has range 1-4, not 0-3
        mutate(tfeq_response = tfeq_response + 1) %>% 
        pivot_wider(names_from = variables, values_from = tfeq_response) %>%
        mutate(tfeq_cr_sumScore = sum(c(tfeq_2, tfeq_11, tfeq_12, tfeq_15, tfeq_16, tfeq_18)), 
               tfeq_ue_sumScore = sum(c(tfeq_1, tfeq_4, tfeq_5, tfeq_7, tfeq_8, tfeq_9, tfeq_13, tfeq_14, tfeq_17)),
               tfeq_ee_sumScore = sum(c(tfeq_3, tfeq_6, tfeq_10))),
      
      # Catch questions
      raw_dat[[subj]] %>% 
        filter(questionnaire == "CATCH_QUESTIONS") %>% 
        select(question_no, responses) %>% 
        data.frame(variables = paste("catch_question", .$question_no, sep = "_"), .) %>% 
        arrange(., question_no) %>%
        select(!question_no) %>% 
        mutate(responses = as.numeric(gsub("\\D", "", responses))) %>% 
        pivot_wider(names_from = variables, values_from = responses) %>% 
        { if(session == "s1"){
          add_column(., 
                     catch_hard_1_pass = ifelse(.$catch_question_1 == 2 | .$catch_question_1 == 3, 1, 0),
                     catch_hard_2_pass = ifelse(.$catch_question_2 == 2 | .$catch_question_2 == 3, 1, 0),
                     catch_easy_3_pass = ifelse(.$catch_question_3 == 0, 1, 0),
                     catch_easy_4_pass = ifelse(.$catch_question_4 == 6, 1, 0))
        } else if(session == "s2"){
          add_column(., 
                     catch_hard_1_pass = ifelse(.$catch_question_1 == 2 | .$catch_question_1 == 3, 1, 0),
                     catch_hard_2_pass = ifelse(.$catch_question_2 == 2 | .$catch_question_2 == 3, 1, 0),
                     catch_easy_3_pass = ifelse(.$catch_question_3 == 2, 1, 0),
                     catch_easy_4_pass = ifelse(.$catch_question_4 == 6, 1, 0))
        }
          } %>% 
        add_column(catch_questions_pass = case_when(((.$catch_hard_1_pass == 0) & (.$catch_hard_2_pass == 0)) ~ 0,
                                                    .$catch_easy_3_pass == 0 ~ 0, 
                                                    .$catch_easy_4_pass == 0 ~ 0,
                                                    .default = 1))

      )
      
      # Combine
      comp_3_dat %<>%
        bind_rows(., comp_3_dat_subj)
      rm(comp_3_dat_subj)
      
      if(display_progress == TRUE){
        setTxtProgressBar(progress_bar, value = subj)
      }
    }
    
    if(display_progress == TRUE){
      close(progress_bar)
    }
    
    
    ### COMPONENT 4: Other questions data -------------------------------------------------------------------------------------
    # Contains data on hunger, CGM, and last meal questions, as well as the final survey
    
    comp_4_dat <- tibble()
    
    if(display_progress == TRUE){
      print("Parsing component 4/4:")
      progress_bar = txtProgressBar(min = 0, max = length(comp_3_id), style = 1, char = "=")
    }

    
    for(subj in comp_4_id){
      
      comp_4_dat_subj <- raw_dat[[subj]] %>%
        select(prolific_id) %>% 
        summarise(across(prolific_id, .fns=~na.omit(unique(.)))) %>%
        as_tibble() 
      
      comp_4_dat_subj <- bind_cols(
        
        comp_4_dat_subj, 
       
        # Last meal
        raw_dat[[subj]] %>% 
          filter(questionnaire == "MEAL", question_no %in% 1:3) %>% 
          select(meal_response) %>% 
          data.frame(variables = c("last_meal_time", "last_meal_size", "snack"),  .) %>% 
          pivot_wider(names_from = variables, values_from = meal_response) %>% 
          mutate(snack_time = { ifelse(snack == 0, NA, 
                                       raw_dat[[subj]] %>% 
                                         filter(questionnaire == "MEAL") %>%
                                         filter(question_no == 4) %>% 
                                         { if(comp_4_dat_subj$prolific_id %in% bug_fix_ids) { 
                                           select(., snack_response) %>% as.character()
                                         } else {
                                           select(., meal_response) %>% as.character()
                                         }
                                         }) }),
        
        
        # Hunger
        raw_dat[[subj]] %>% 
          filter(questionnaire == "HUNGER") %>% 
          # fix js bug mistake
          { if(comp_4_dat_subj$prolific_id %in% bug_fix_ids) {
            select(., tfeq_response) %>%
              rename(., hunger_response = tfeq_response)
          } else {
            select(., hunger_response)
          }
          } %>% 
          data.frame(variables = "hunger_rating",  .) %>% 
          pivot_wider(names_from = variables, values_from = hunger_response),
       
       
        # CGL
        raw_dat[[subj]] %>% 
          filter(questionnaire == "CGL") %>% 
          filter(question_no == 1) %>% 
          select(cgl_response) %>% 
          data.frame(variables = "cgl",  .) %>% 
          pivot_wider(names_from = variables, values_from = cgl_response) %>% 
          mutate(cgl = as.numeric(cgl)) %>% 
          mutate(cgl_measure = case_when(cgl == 0 ~ NA, 
                                         cgl == 1 ~ as.numeric(raw_dat[[subj]] %>% 
                                           filter(questionnaire == "CGL") %>% 
                                           filter(question_no == 2) %>% 
                                           select(cgl_response)))) %>% 
          mutate(cgl_unit = case_when(cgl == 0 ~ NA,
                                      cgl == 1 ~ as.character(raw_dat[[subj]] %>% 
                                                                filter(questionnaire == "CGL") %>% 
                                                                  select(., responses) %>% 
                                                                    slice_tail() %>% 
                                                                    { ifelse(str_detect(., "mgdL"), "mgd_L", "mmol_L") }
                                                                ))),
          
      
        
        # Final survey
        raw_dat[[subj]] %>% 
          select(enjoyment:other_comments) %>% 
          filter_all(any_vars(!is.na(.)))
        
      )
      
      # Combine
      comp_4_dat %<>%
        bind_rows(., comp_4_dat_subj)
      rm(comp_4_dat_subj)
      
      if(display_progress == TRUE){
        setTxtProgressBar(progress_bar, value = subj)
      }
    }
    
    if(display_progress == TRUE){
      close(progress_bar)
    }
    
    ### OUTPUT -------------------------------------------------------------------------------------------------

    return(list(
      # (1): Demographics and medication data
      "demographic_dat" = comp_1_dat, 
      # (2): Task data
      "task_dat" = comp_2_dat, 
      # (3): Questionnaire data
      "questionnaire_dat" = left_join(comp_3_dat, comp_4_dat, by = "prolific_id"), 
      # (4): Prolific data
      "prolific_dat" = meta_data %>%
        filter(prolific_id %in% comp_1_dat$prolific_id)))
    
    
  })
}







