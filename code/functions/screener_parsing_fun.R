########################################################################################################################
############################---------------- PARSING SCREENER DATA FUNCTION ----------------############################
########################################################################################################################

## Parsing data from screening -----------------------------------------
# @ file: vector of relative path(s) to .txt data file(s) or existing R object
# @ meta_file: vector of relative path(s) to .csv data file(s) with Prolific meta data

screener_parsing <- function(file,
                             meta_file, 
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
  
  raw_meta_dat %<>%
    rowwise %>%
    mutate("Submission.id" = id_shuffle(Submission.id)) %>%
    mutate("Participant.id" = id_shuffle(Participant.id))

    
  ### Identify component indices ---------------------------------------------------------------------------------------
  comp_1_id <- raw_dat %>%
    lmap(\(x) list(any(grepl("Welcome to the study", x, fixed=T)))) %>%
    unlist(.) %>% which(. == TRUE)
  comp_2_id <- raw_dat %>%
    lmap(\(x) list(any(grepl("more specific questions regarding your health", x, fixed=T)))) %>%
    unlist(.) %>% which(. == TRUE)
  comp_3_id <- raw_dat %>%
    lmap(\(x) list(any(grepl("The next questions will ask you about the time you spent being physically active", x, fixed=T)))) %>%
    unlist(.) %>% which(. == TRUE)
  
  ### TIBBLE 1: Main data -------------------------------------------------------------------------------------------------
  # Contains demographic data, brief medical data, medication data, current weight loss interventions, BMI, and IPAQ
  
  ## Component 1: Basic demographic questions -------
  comp_1_dat <- tibble()
  
  if(display_progress == TRUE){
    print("Parsing component 1/3:")
    progress_bar = txtProgressBar(min=0, max=length(comp_1_id), style = 1, char="=")
  }
  
  for(subj in comp_1_id){
    comp_1_dat_subj <- raw_dat[[subj]] %>%
      mutate(prolific_id = id_shuffle(prolific_id)) %>%
      select(c(prolific_id, age:diabetes_other)) %>% 
      summarise(across(c(prolific_id, age:diabetes_other), .fns=~na.omit(unique(.)))) %>%
      as_tibble() %>%
      mutate(diabetes = ifelse(!"diabetes" %in% names(.), "0", 
                               ifelse(diabetes == 0, "0", diabetes))) %>% 
      mutate(diabetes_other = ifelse(!diabetes_other, NA, na.omit(raw_dat[[subj]]$diabetes_other_text)))
    
    # Combine
    comp_1_dat %<>%
      bind_rows(., comp_1_dat_subj)
    rm(comp_1_dat_subj)
    
    if(display_progress == TRUE){
    setTxtProgressBar(progress_bar, value = subj)
    }
  }
  
  if(display_progress == TRUE){
  close(progress_bar)
  }
  
  # remove duplicates
  comp_1_dat <- comp_1_dat[!duplicated(comp_1_dat$prolific_id, fromLast = TRUE), ]
  
  ## Component 2: Treatment and medical data -------
  # Contains treatment data (GLP, diabetes medication, and other daily medication), and data on medical conditions
  
  comp_2_dat <- tibble()
  
  if(display_progress == TRUE){
    print("Parsing component 2/3:")
    progress_bar = txtProgressBar(min=0, max=length(comp_2_id), style = 1, char="=")
  }
   
  for(subj in comp_2_id){
    # Unconditional questions
    comp_2_dat_subj <- raw_dat[[subj]] %>%
      mutate(prolific_id = id_shuffle(prolific_id)) %>%
      select(c(prolific_id, glp_treatment, medication, neurological, psych_neurdev, chronic_disease)) %>% 
      summarise(across(.fns=~max(.x, na.rm=TRUE))) %>%
      as_tibble()
    
    # GLP follow up questions
    comp_2_dat_subj %<>% 
      mutate(glp_treatment = ifelse(.$glp_treatment == "Current", 
                                    ifelse(str_detect(na.omit(raw_dat[[subj]]["glp_type"])[1,1], 
                                                      "None"), 
                                           "No", "Current"),
                                    .$glp_treatment)) %>% 
      add_column(
        type_glp = { ifelse(.$glp_treatment == "Current", 
                            na.omit(raw_dat[[subj]]["glp_type"])[1,1], 
                            NA) }, 
        type_glp_other = { ifelse(.$glp_treatment == "Current", 
                                  ifelse(na.omit(raw_dat[[subj]]["glp_type_other"]), 
                                         na.omit(raw_dat[[subj]]$glp_type_other_text), 
                                         NA),
                                  NA) }, 
        start_date_glp = { ifelse(.$glp_treatment == "Current", 
                                  # correct impossible date to prevent lubridate error
                                  if(na.omit(raw_dat[[subj]]["glp_start_date"])[1,1] == "30/2/2019"){
                                    as.character(dmy("01/03/2019"))
                                  } else {
                                    as.character(dmy(na.omit(raw_dat[[subj]]["glp_start_date"])))
                                  }, 
                                  NA) },
        administration_glp = { ifelse(.$glp_treatment == "Current", 
                                      na.omit(raw_dat[[subj]]["glp_administration"])[1,1], 
                                      NA) },
        administration_glp_other = { ifelse(.$glp_treatment == "Current", 
                                            ifelse(na.omit(raw_dat[[subj]]["glp_administration_other"]), 
                                                   na.omit(raw_dat[[subj]]$glp_administration_other_text), 
                                                   NA),
                                            NA) }, 
        schedule_glp = { ifelse(.$glp_treatment == "Current", 
                                { ifelse(.$prolific_id == "93563ddbec", 
                                         "Other", 
                                         na.omit(raw_dat[[subj]]["glp_schedule"])[1,1]) },
                                NA) },
        .after = "glp_treatment") %>% 
      add_column(
        schedule_glp_other = { ifelse(.$schedule_glp == "Other", 
                                      na.omit(raw_dat[[subj]]["responses"])[5,] %>% 
                                        str_replace_all('glp_schedule', '') %>% 
                                        str_replace_all('other', '') %>% 
                                        str_replace_all('_text', '') %>% 
                                        str_replace_all('[[:punct:]]', ''),
                                     NA) },
        .after = "schedule_glp") %>% 
      add_column(
        injection_day_glp = { ifelse(.$schedule_glp == "Weekly", 
                                     na.omit(raw_dat[[subj]]["glp_injection_day"])[1,1], 
                                     NA) },
        side_effects_glp = { ifelse(.$glp_treatment == "Current", 
                                    ifelse(na.omit(raw_dat[[subj]]["glp_side_effects"])[1,1] == "Yes" |
                                             na.omit(raw_dat[[subj]]["glp_side_effects"])[1,1] == "Maybe", 
                                           ifelse(str_detect(na.omit(raw_dat[[subj]]["side_effects_symptoms"])[1,1], 
                                                             "None"), 
                                                  "No", na.omit(raw_dat[[subj]]["glp_side_effects"])[1,1]), "No"),
                                    NA) },
        .after = "schedule_glp_other") %>% 
      add_column(
        side_effects_glp_symptoms = { ifelse(.$side_effects_glp != "No", 
                                             paste(
                                               unlist(
                                                 jsonlite::parse_json(
                                                   na.omit(raw_dat[[subj]]["side_effects_symptoms"])[1,1])), 
                                               collapse=', '),
                                             NA) }, 
        side_effects_glp_symptoms_other = { ifelse(.$side_effects_glp != "No", 
                                                   ifelse(na.omit(raw_dat[[subj]]["side_effects_symptoms_other"]), 
                                                          na.omit(raw_dat[[subj]]$side_effects_symptoms_other_text), 
                                                          NA), 
                                                   NA) },
        .after = "side_effects_glp")
    
    # correct Ozempic / Rybelsus bug in batch 1: 
    
    comp_2_dat_subj %<>%
      mutate(across(type_glp, ~ case_when((type_glp == "Semaglutide" & 
                                             administration_glp == "Injection" & 
                                             schedule_glp == "Weekly") ~ "Semaglutide_Ozempic", 
                                          (type_glp == "Semaglutide" & 
                                             administration_glp == "Pill" & 
                                             schedule_glp == "Daily") ~ "Semaglutide_Rybelsus", 
                                          .default = .)))
    
    
    # Medication follow up questions
    comp_2_dat_subj %<>% 
      add_column(
        diabetes_medication =  { ifelse(comp_1_dat %>%
                                          filter(prolific_id == comp_2_dat_subj$prolific_id) %>%
                                          select(diabetes) %>%
                                          as.character() == "type 2", 
                                        ifelse(na.omit(raw_dat[[subj]]["diabetes_medication"])[1,1] == 1, 
                                               ifelse(str_detect(na.omit(raw_dat[[subj]]["diabetes_medication_type"])[1,1], 
                                                                 "None"), 
                                                      0, na.omit(raw_dat[[subj]]["diabetes_medication"])[1,1]), 
                                               NA), 
                                        NA) },
        .before = "medication") %>% 
      add_column(
        diabetes_medication_type = { ifelse(.$diabetes_medication == 1, 
                                            paste(
                                              unlist(
                                                jsonlite::parse_json(
                                                  na.omit(raw_dat[[subj]]["diabetes_medication_type"])[1,1])), 
                                              collapse=', '),
                                            NA) }, 
        diabetes_medication_type_other = { ifelse(.$diabetes_medication == 1, 
                                                  ifelse(na.omit(raw_dat[[subj]]["diabetes_medication_type_other"]), 
                                                         na.omit(raw_dat[[subj]]$diabetes_medication_type_other_text), 
                                                         NA), 
                                                  NA) },
        .before = "medication") %>% 
      add_column(
        contraceptive = { ifelse(.$medication == 1, 
                                 max(raw_dat[[subj]]["contraceptive"], na.rm = TRUE),
                                 NA) }, 
        antidepressant = { ifelse(.$medication == 1, 
                                  max(raw_dat[[subj]]["antidepressant"], na.rm = TRUE),
                                  NA) },
        .after = "medication") %>% 
      add_column(
        antidepressant_type = { ifelse(.$antidepressant == 1, 
                                       paste(
                                         unlist(
                                           jsonlite::parse_json(
                                             na.omit(raw_dat[[subj]]["antidepressant_type"])[1,1])), 
                                         collapse=', '),
                                       NA) }, 
        antidepressant_type_other = { ifelse(.$antidepressant == 1, 
                                             ifelse(na.omit(raw_dat[[subj]]["antidepressant_type_other"]), 
                                                    na.omit(raw_dat[[subj]]$antidepressant_type_other_text), 
                                                    NA), 
                                             NA) },
        other_medication = { ifelse(.$medication == 1, 
                                    max(raw_dat[[subj]]["other_medication"], na.rm = TRUE),
                                    NA) },
        .after = "antidepressant") %>% 
      add_column(
        other_medication_type = { ifelse(.$other_medication == 1, 
                                         paste(
                                           unlist(
                                             jsonlite::parse_json(
                                               na.omit(raw_dat[[subj]]["other_medication_type"])[1,1])), 
                                           collapse=', '),
                                         NA) }, 
        .after = "other_medication")
    
    # Medical conditions follow up questions
    comp_2_dat_subj %<>% 
      mutate(neurological = ifelse(.$neurological == 1, 
                                   ifelse(str_detect(na.omit(raw_dat[[subj]]["neurological_condition"])[1,1], 
                                                     "None"), 
                                          0, 1),
                                   0)) %>% 
      add_column(
        neurological_condition =  { ifelse(.$neurological == 1, 
                                           paste(
                                             unlist(
                                               jsonlite::parse_json(
                                                 na.omit(raw_dat[[subj]]["neurological_condition"])[1,1])), 
                                             collapse=', '),
                                           NA) }, 
        .after = "neurological") %>% 
      add_column(
        neurological_condition_other = { ifelse(.$neurological == 1, 
                                                ifelse(na.omit(raw_dat[[subj]]["neurological_condition_other"]), 
                                                       na.omit(raw_dat[[subj]]$neurological_condition_other_text), 
                                                       NA), 
                                                NA) }, 
        .after = "neurological_condition") %>% 
      mutate(psych_neurdev = ifelse(.$psych_neurdev == 1, 
                                    ifelse(str_detect(na.omit(raw_dat[[subj]]["psych_neurdev_condition"])[1,1], 
                                                      "None"), 
                                           0, 1),
                                    0)) %>% 
      add_column(
        psych_neurdev_condition =  { ifelse(.$psych_neurdev == 1, 
                                            paste(
                                              unlist(
                                                jsonlite::parse_json(
                                                  na.omit(raw_dat[[subj]]["psych_neurdev_condition"])[1,1])), 
                                              collapse=', '),
                                            NA) }, 
        .after = "psych_neurdev") %>% 
      add_column(
        psych_neurdev_condition_other = { ifelse(.$psych_neurdev == 1, 
                                                 ifelse(na.omit(raw_dat[[subj]]["psych_neurdev_condition_other"]), 
                                                        na.omit(raw_dat[[subj]]$psych_neurdev_condition_other_text), 
                                                        NA), 
                                                 NA) }, 
        .after = "psych_neurdev_condition")  %>% 
      mutate(chronic_disease = ifelse(.$chronic_disease == 1, 
                                      ifelse(str_detect(na.omit(raw_dat[[subj]]["chronic_disease_condition"])[1,1], 
                                                        "None"), 
                                             0, 1),
                                      NA)) %>% 
      add_column(
        chronic_disease_condition =  { ifelse(.$chronic_disease == 1, 
                                              paste(
                                                unlist(
                                                  jsonlite::parse_json(
                                                    na.omit(raw_dat[[subj]]["chronic_disease_condition"])[1,1])), 
                                                collapse=', '),
                                              NA) }, 
        .after = "chronic_disease") %>% 
      add_column(
        chronic_disease_condition_other = { ifelse(.$chronic_disease == 1, 
                                                   ifelse(na.omit(raw_dat[[subj]]["chronic_disease_condition_other"]), 
                                                          na.omit(raw_dat[[subj]]$chronic_disease_condition_other_text), 
                                                          NA), 
                                                   NA) }, 
        .after = "chronic_disease_condition") 
    
    # Combine
    comp_2_dat %<>%
      bind_rows(., comp_2_dat_subj)
    rm(comp_2_dat_subj)
    
    if(display_progress == TRUE){
      setTxtProgressBar(progress_bar, value = subj)
    }
  }
  
  if(display_progress == TRUE){
    close(progress_bar)
  }
  
  # remove duplicates
  comp_2_dat <- comp_2_dat[!duplicated(comp_2_dat$prolific_id, fromLast = TRUE), ]
  
  ## Component 3: Metabolic questions data -------
  # Contains data on current weight loss interventions, BMI, and IPAQ scores
  comp_3_dat <- tibble()
  
  if(display_progress == TRUE){
    print("Parsing component 3/3:")
    progress_bar = txtProgressBar(min=0, max=length(comp_3_id), style = 1, char="=")
  }
  
  for(subj in comp_3_id){
    # IPAQ
    comp_3_dat_subj <- raw_dat[[subj]] %>% 
      mutate(prolific_id = id_shuffle(prolific_id)) %>%
      select(ipaq_response) %>% 
      na.omit() %>%
      data.frame(variables = paste("ipaq", 1:7, sep = "_"), .) %>%
      pivot_wider(names_from = variables, values_from = ipaq_response) %>%
      mutate(across(c(ipaq_1, ipaq_3, ipaq_5), 
                    ~ ifelse(.x == "None", 0, as.double(.)))) %>% 
      mutate(across(c(ipaq_2, ipaq_4, ipaq_6), ~na_if(., "NA"))) %>% 
      add_column(
        prolific_id = { raw_dat[[subj]] %>%
            select(prolific_id) %>% 
            summarise(across(prolific_id, .fns=~na.omit(unique(.x)))) %>% 
            .[1,1]
        },
        .before = "ipaq_1"
      ) %>% 
      add_column(
        ipaq_vigorous_MET = { ifelse(.$ipaq_1 == 0, 
                                     0, 
                                     8 * period_to_seconds(hm(.$ipaq_2))/60 * as.numeric(.$ipaq_1)) }, 
        ipaq_moderate_MET = { ifelse(.$ipaq_3 == 0, 
                                     0, 
                                     4 * period_to_seconds(hm(.$ipaq_4))/60 * as.numeric(.$ipaq_3)) }, 
        ipaq_walking_MET = { ifelse(.$ipaq_5 == 0, 
                                    0, 
                                    3.3 * period_to_seconds(hm(.$ipaq_6))/60 * as.numeric(.$ipaq_5)) }) %>%
      add_column(
        ipaq_sumScore = { sum(.$ipaq_walking_MET, .$ipaq_moderate_MET, .$ipaq_vigorous_MET) })
    
    
    # BMI
    comp_3_dat_subj %<>% 
      add_column(
        height_cm = { 
          raw_dat[[subj]] %>%
            select(bmi_response) %>%
            reframe(across(.fns=~na.omit(.x))) %>% .[1,] %>%
            as.numeric()}, 
        weight_kg = { 
          raw_dat[[subj]] %>%
            select(bmi_response) %>%
            reframe(across(.fns=~na.omit(.x))) %>% .[2,] %>%
            as.numeric()}, 
        height_cm_raw = { 
          raw_dat[[subj]] %>%
            select(bmi_response_raw) %>%
            reframe(across(.fns=~na.omit(.x))) %>% .[1,]}, 
        weight_kg_raw = { 
          raw_dat[[subj]] %>%
            select(bmi_response_raw) %>%
            reframe(across(.fns=~na.omit(.x))) %>% .[2,]}
      ) %>% 
      add_column(
        bmi = { 
          .$weight_kg / ((.$height_cm/100)^2)
        }
      ) 
    
    # Weight loss interventions
    comp_3_dat_subj %<>% 
      add_column(
        weight_loss = { ifelse(na.omit(raw_dat[[subj]]["weight_loss"])[1,1] == "Yes" |
                                 na.omit(raw_dat[[subj]]["weight_loss"])[1,1] == "Maybe", 
                               ifelse(str_detect(na.omit(raw_dat[[subj]]["weight_loss_interventions"])[1,1], 
                                                 "None"), 
                                      "No", na.omit(raw_dat[[subj]]["weight_loss"])[1,1]), "No") }) %>% 
      add_column(
        weight_loss_interventions = { ifelse(.$weight_loss != "No", 
                                             paste(
                                               unlist(
                                                 jsonlite::parse_json(
                                                   na.omit(raw_dat[[subj]]["weight_loss_interventions"])[1,1])), 
                                               collapse=', '),
                                             NA) }, 
        weight_loss_interventions_other = { ifelse(.$weight_loss != "No", 
                                                   ifelse(na.omit(raw_dat[[subj]]["weight_loss_interventions_other"]), 
                                                          na.omit(raw_dat[[subj]]$weight_loss_interventions_other_text), 
                                                          NA), 
                                                   NA) })
    
    
    ## Combine -------
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
  
  # remove duplicates
  comp_3_dat <- comp_3_dat[!duplicated(comp_3_dat$prolific_id, fromLast = TRUE), ]
  
  suppressMessages(
    main_data <- bind_cols(comp_1_dat, comp_2_dat, comp_3_dat) %>% 
      select_if(!duplicated(sub("\\.\\.\\..*", "", names(.)))) %>%
      rename(prolific_id = "prolific_id...1") 
  )


  
  ### TIBBLE 2: Prolific meta data -------------------------------------------------------------------------------------
  # Contains Prolific meta data
  
  # Only approved participants
  meta_data <- raw_meta_dat %>%
    filter(Status == "APPROVED") %>%
    select(Participant.id, Age, Sex, Ethnicity.simplified, Country.of.residence, Started.at) %>%
    as_tibble() %>% 
    rename(prolific_id = Participant.id, age = Age, sex = Sex, 
           ethnicity = Ethnicity.simplified, residence = Country.of.residence, 
           screening_day = Started.at) %>%
    mutate(screening_day = as_datetime(screening_day)) %>%
    filter(prolific_id %in% main_data$prolific_id)
  
  ### OUTPUT -------------------------------------------------------------------------------------------------
  
  return(list("screening_dat" = main_data, 
              "prolific_dat" = meta_data))
  
  })
}
