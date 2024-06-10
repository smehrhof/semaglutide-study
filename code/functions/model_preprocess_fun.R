################################################################################
###############------ PREPROCESS FUNCTION FOR MODEL FITTING -----###############
################################################################################

## Proprocessing task data for cmd stan modelling ------------------------------
# @ raw_data: raw task data to be modelled
# @ subjs: character vector of subj IDs
# @ n_subj: number of subjects
# @ t_subjs: numeric vector with number of trails per subject
# @ t_max: maximum number of trials

model_preprocessing <- function(raw_data, 
                                retest = FALSE,
                                subjs, 
                                n_subj, 
                                t_subjs, 
                                t_max) {
  # Currently class(raw_data) == "data.table"
  
  if(retest == FALSE){
    if(length(raw_data) == 2){
      warning("Two datasets provided but retest = FALSE selected; Only the first dataframe will be used.")
      raw_data <- raw_data[[1]]
    }
    
    # Initialize (model-specific) data arrays
    effort_a <- array( 0, c(n_subj, t_max))
    amount_a <- array( 0, c(n_subj, t_max))
    effort_b <- array( 0, c(n_subj, t_max))
    amount_b <- array( 0, c(n_subj, t_max))
    choice <- array(-1, c(n_subj, t_max))
    
    # Write from raw_data to the data arrays
    for (i in 1:n_subj) {
      subj <- subjs[i]
      t <- t_subjs[i]
      DT_subj <- raw_data[raw_data$subj_id == subj,]
      
      effort_a[i, 1:t]   <- DT_subj$effort_a
      amount_a[i, 1:t]  <- DT_subj$amount_a
      effort_b[i, 1:t]  <- DT_subj$effort_b
      amount_b[i, 1:t] <- DT_subj$amount_b
      choice[i, 1:t]        <- DT_subj$choice
    }
    
    # Wrap into a list for Stan
    data_list <- list(
      N             = n_subj,
      T             = t_max,
      Tsubj         = t_subjs,
      effort_a   = effort_a,
      amount_a  = amount_a,
      effort_b  = effort_b,
      amount_b = amount_b,
      choice        = choice
    )
    
    
  } else if(retest == TRUE){
    if(is.data.frame(raw_data)){
      stop("Please provide dataframes for at least two sessions as a list.")
    }
    n_time <- length(raw_data)
    t_subjs <- array(t_max, c(n_subj, n_time))
    
    # Initialize (model-specific) data arrays
    effort_a <- array(0, c(n_subj, n_time, t_max))
    amount_a <- array( 0, c(n_subj, n_time, t_max))
    effort_b <- array( 0, c(n_subj, n_time, t_max))
    amount_b <- array( 0, c(n_subj, n_time, t_max))
    choice <- array(-1, c(n_subj, n_time, t_max))
    
    # Write from raw_data to the data arrays
    for(time in 1:n_time){
      for (i in 1:n_subj) {
        subj <- subjs[i]
        t <- t_subjs[i, time]
        
        DT_subj <- raw_data[[time]][raw_data[[time]]$subj_id == subj,]
        
        effort_a[i, time, 1:t]   <- DT_subj$effort_a
        amount_a[i, time, 1:t]  <- DT_subj$amount_a
        effort_b[i, time, 1:t]  <- DT_subj$effort_b
        amount_b[i, time, 1:t] <- DT_subj$amount_b
        choice[i, time, 1:t]        <- DT_subj$choice
      }
    }
    
    # Wrap into a list for Stan
    data_list <- list(
      N = n_subj,
      T = t_max,
      N_time = n_time,
      T_subj = t_subjs,
      effort_a = effort_a,
      amount_a = amount_a,
      effort_b = effort_b,
      amount_b = amount_b,
      choice = choice
    )
  }
  
  # Returned data_list will directly be passed to Stan
  return(data_list)
}





