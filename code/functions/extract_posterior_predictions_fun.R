################################################################################
###############------ POSTERIOR PREDICTIVE CHECK FUNCTION -------###############
################################################################################

## Extracts generated quantities in a memory saving manner ---------------------
# @ csv_paths: character vector with paths to csv  files
# @ n_chains: number of chains
# @ n_iter: number of iterations
# @ n_subj: number of subjects
# @ n_trials: number of trials
# @ real_dat: real data 


posterior_predictions <- function(csv_paths,
                                  n_chains,
                                  n_iter,
                                  n_subj, 
                                  n_trials,
                                  real_dat) {
  
  ppc_results <- list()
  
  pred_var_names <- as.vector(
    sapply(1:n_subj, function(i)
      sapply(1:n_trials, function(j) paste0("y_pred", "[", i, ",", j, "]"))
    )
  )
  
  indiv_vars <- split(
    pred_var_names, ceiling(seq_along(pred_var_names)/n_trials)
  )
  
  rm(pred_var_names)
  
  posterior_predictions_trial_type <- data.frame()
  posterior_predictions_effort <- data.frame()
  posterior_predictions_reward <- data.frame()
  
  for(id in 1:n_subj){
    
    real_dat_indiv <- real_dat[real_dat$subj_id == unique(real_dat$subj_id)[id],]
    
    all_indiv_draws <- list()
    
    all_indiv_draws_df <- data.frame(matrix("", ncol = 0, nrow = n_iter))  
    for(o in seq_along(csv_paths)){
      all_indiv_draws[[o]] <- 
        cmdstanr::read_cmdstan_csv(csv_paths[o], variables = indiv_vars[[id]],
                                   format = "draws_list")$post_warmup_draws
      
      all_indiv_draws_df <- cbind(all_indiv_draws_df, rbind.data.frame(all_indiv_draws[[o]]))
    }  
    
    colnames(all_indiv_draws_df) <- paste0(rep(1:n_chains, each = n_trials), "_y_pred_", 1:n_trials)
    rm(all_indiv_draws)
    
    # get mean prediction and 95% hdi for each effort reward combination
    
    indiv_ppc <- data.frame()
    for(eff in seq_along(sort(unique(real_dat_indiv$effort_a)))){
      for(rew in seq_along(sort(unique(real_dat_indiv$amount_a)))){
        trial_type <- which(real_dat_indiv$effort_a == sort(unique(real_dat_indiv$effort_a))[eff] 
                            & real_dat_indiv$amount_a == sort(unique(real_dat_indiv$amount_a))[rew])
        trial_type_name <- paste0(rep(1:n_chains, each = length(trial_type)), "_y_pred_", trial_type)
        trial_type_preds <- all_indiv_draws_df[,trial_type_name]
        
        summed_trial_type_preds <- rowSums(trial_type_preds)/dim(trial_type_preds)[2]
        
        trial_type_hdi <- hBayesDM::HDIofMCMC(summed_trial_type_preds, credMass = 0.95)
        
        indiv_trial_type_ppc <- data.frame("subj_id" = unique(real_dat_indiv$subj_id),
                                           "effort_a" = sort(unique(real_dat_indiv$effort_a))[eff],
                                           "amount_a" = sort(unique(real_dat_indiv$amount_a))[rew],
                                           "observation" = mean(real_dat_indiv[trial_type,]$choice),
                                           "prediction_mean" = mean(summed_trial_type_preds),
                                           "prediction_hdi_lower" = as.numeric(trial_type_hdi[1]),
                                           "prediction_hdi_higher" = as.numeric(trial_type_hdi[2]))
        indiv_ppc <- rbind(indiv_ppc, indiv_trial_type_ppc)
      }
    }
    posterior_predictions_trial_type <- rbind(posterior_predictions_trial_type, indiv_ppc)
    ppc_results$posterior_predictions_trial_type <- posterior_predictions_trial_type
    
    indiv_ppc <- data.frame()
    for(eff in seq_along(sort(unique(real_dat_indiv$effort_a)))){
      
      trial_type <- which(real_dat_indiv$effort_a == sort(unique(real_dat_indiv$effort_a))[eff])
      trial_type_name <- paste0(rep(1:n_chains, each = length(trial_type)), "_y_pred_", trial_type)
      trial_type_preds <- all_indiv_draws_df[,trial_type_name]
      
      summed_trial_type_preds <- rowSums(trial_type_preds)/dim(trial_type_preds)[2]
      
      
      trial_type_hdi <- hBayesDM::HDIofMCMC(summed_trial_type_preds, credMass = 0.95)
      
      indiv_trial_type_ppc <- data.frame("subj_id" = unique(real_dat_indiv$subj_id),
                                         "effort_a" = sort(unique(real_dat_indiv$effort_a))[eff],
                                         "observation" = mean(real_dat_indiv[trial_type,]$choice),
                                         "prediction_mean" = mean(summed_trial_type_preds),
                                         "prediction_hdi_lower" = as.numeric(trial_type_hdi[1]),
                                         "prediction_hdi_higher" = as.numeric(trial_type_hdi[2]))
      indiv_ppc <- rbind(indiv_ppc, indiv_trial_type_ppc)
      
    }
    posterior_predictions_effort <- rbind(posterior_predictions_effort, indiv_ppc)
    ppc_results$posterior_predictions_effort <- posterior_predictions_effort
    
    indiv_ppc <- data.frame()
    for(rew in seq_along(sort(unique(real_dat_indiv$amount_a)))){
      
      trial_type <- which(real_dat_indiv$amount_a == sort(unique(real_dat_indiv$amount_a))[rew])
      trial_type_name <- paste0(rep(1:n_chains, each = length(trial_type)), "_y_pred_", trial_type)
      trial_type_preds <- all_indiv_draws_df[,trial_type_name]
      
      summed_trial_type_preds <- rowSums(trial_type_preds)/dim(trial_type_preds)[2]
      
      
      trial_type_hdi <- hBayesDM::HDIofMCMC(summed_trial_type_preds, credMass = 0.95)
      
      indiv_trial_type_ppc <- data.frame("subj_id" = unique(real_dat_indiv$subj_id),
                                         "amount_a" = sort(unique(real_dat_indiv$amount_a))[rew],
                                         "observation" = mean(real_dat_indiv[trial_type,]$choice),
                                         "prediction_mean" = mean(summed_trial_type_preds),
                                         "prediction_hdi_lower" = as.numeric(trial_type_hdi[1]),
                                         "prediction_hdi_higher" = as.numeric(trial_type_hdi[2]))
      indiv_ppc <- rbind(indiv_ppc, indiv_trial_type_ppc)
      
    }
    posterior_predictions_reward <- rbind(posterior_predictions_reward, indiv_ppc)
    ppc_results$posterior_predictions_reward <- posterior_predictions_reward
    
  }
  
  return(ppc_results)
}






