################################################################################
###############----------- CONVERGENCE CHECK FUNCTION -----------###############
################################################################################

## Computes multiple model convergence checks ----------------------------------
# @ model_fit: Model fit object
# @ params: Parameters of the model
# @ Rhat: Include Rhat diagnostics? 
# @ ess: Include estimated sample size diagnostics?
# @ trace_plot: Include visual trace plot diagnostics?
# @ rank_hist: include rank histogram diagnostics? 

convergence_check <- function(model_fit, 
                              params = c("kE", "kR", "a"), 
                              Rhat = TRUE, 
                              ess = TRUE,
                              trace_plot = TRUE, 
                              rank_hist = FALSE) {
  
  mu_params <- paste("mu", params, sep = "_")
  
  res <- list()
  
  
  if(Rhat | ess) {
    
    model_fit_sum <- cmdstanr::read_cmdstan_csv(
      model_fit$output_files(), variables = params
    )
    posterior_summarise_draws <- posterior::summarise_draws(model_fit_sum$post_warmup_draws)
    rm(model_fit_sum)
    
    if(Rhat) {
      res$Rhat <- summary(posterior_summarise_draws$rhat)
    }
    
    if(ess) {
      res$ess <- summary(posterior_summarise_draws$ess_bulk)
    }
    
    rm(posterior_summarise_draws)
  }
  
  
  if(trace_plot | rank_hist){
    
    model_fit_sum <- cmdstanr::read_cmdstan_csv(
      model_fit$output_files(), variables = mu_params
    )
    
    pal <- c("#cb793a", "#fecee9", "#e1bbc9", "#39a2ae", "#a93f55", "#386641")
    
    bayesplot::bayesplot_theme_set(
      cowplot::theme_half_open(
        font_family = "", font_size = 11
      )
    )
    
    bayesplot::color_scheme_set(pal)
    
    if(trace_plot) {
      res$trace_plot <- bayesplot::mcmc_trace(model_fit_sum$post_warmup_draws)
    }
    
    if(rank_hist) {
      res$rank_hist <- bayesplot::mcmc_rank_hist(model_fit_sum$post_warmup_draws)
    }
    
    rm(model_fit_sum)
    
  }
  
  return(res)
}

