################################################################################
###############----------- MODEL COMPARISON FUNCTION ------------###############
################################################################################

## Conducts model comparison in a memory saving manner -------------------------
# @ model_paths: character vector with directory paths to models
# @ loo_paths: character vector with directory paths to loo fits
# @ model_names: character vector with model names
# @ LOO_plot: visualize model comparison based on loo?
# @ ELPD_plot: visualize model comparison based on elpd?


model_comparison <- function(model_paths = NULL, 
                             loo_paths = NULL,
                             model_names = NULL,
                             LOO_plot = FALSE,
                             LOO_title = "Model comparision based on LOO",
                             LOO_col = "#5B9BD5",
                             ELPD_plot = FALSE, 
                             ELPD_title = "Model comparision based on ELPD",
                             ELPD_col = "#71AB48", 
                             LOO_ELPD_plot = TRUE,
                             LOO_ELPD_title = "LOO and ELPD model comparison",
                             LOO_ELPD_col = c("#5B9BD5", "#71AB48")){
  
  require(ggplot2)
  require(ggpubr)
  
  loo_output <- list()
  
  res <- list()
  
  if(!is.null(loo_paths)){
    for(i in seq_along(loo_paths)){
      loo_output[[model_names[i]]] <- readRDS(loo_paths[i])
    }
  } else {
    for(i in seq_along(model_paths)){
      model_fit <- readRDS(model_paths[i])
      if(is.null(model_names)){
        loo_output[[model_fit$metadata()$model_name]] <- model_fit$loo()
      } else {
        loo_output[[model_names[i]]] <- model_fit$loo()
      }
    }
  }
  res$loo_output <- map(loo_output, function(x) x$estimates)
  
  if(LOO_plot){
    model_comp_loo_dat <- data.frame("model" = names(loo_output),
                                     "loo" = as.numeric(lapply(loo_output, function(x) x$estimates[3,1])),
                                     "loo_sd" = as.numeric(lapply(loo_output, function(x) x$estimates[3,2])))
    
    model_comp_loo_plot <- ggplot(model_comp_loo_dat, aes(x=model, y=loo)) + 
      geom_bar(stat="identity", position=position_dodge(), fill = LOO_col) +
      geom_errorbar(aes(ymin=loo-loo_sd, ymax=loo+loo_sd), width=.2,
                    position=position_dodge(.9)) +
      #geom_text(aes(label = round(loo)), position = position_dodge(0.9), hjust=1.5, size = 3) +
      ylab("LOOIC (± SE)") + xlab("Model") +
      scale_y_continuous(breaks=(c(0, 4000)), limits = c(0, 4000)) +
      ggtitle(LOO_title) +
      coord_flip() +
      theme(legend.position = "none", 
            plot.title = element_text(size = 10),
            axis.title = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10)) 
    
    res$LOO_plot <- model_comp_loo_plot
  }
  
  if(ELPD_plot){
    model_comp_elpd_dat <- data.frame("model" = names(loo_output),
                                      "elpd" = as.numeric(lapply(loo_output, function(x) x$estimates[1,1])),
                                      "elpd_sd" = as.numeric(lapply(loo_output, function(x) x$estimates[1,2])))
    
    model_comp_elpd_plot <- ggplot(model_comp_elpd_dat, aes(x=model, y=elpd)) + 
      geom_bar(stat="identity", position=position_dodge(), fill = ELPD_col) +
      geom_errorbar(aes(ymin=elpd-elpd_sd, ymax=elpd+elpd_sd), width=.2,
                    position=position_dodge(.9)) +
      #geom_text(aes(label = round(elpd)), position = position_dodge(0.9), hjust=-0.5, size = 3) +
      ylab("Expected log posterior density (± SE)") + xlab("Model") +
      ggtitle(ELPD_title) +
      coord_flip() +
      theme(legend.position = "none", 
            plot.title = element_text(size = 10),
            axis.title = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10)) 
    
    res$ELPD_plot <- model_comp_elpd_plot
  }
  
  if(LOO_ELPD_plot){
    model_comp_loo_elpd_dat <- data.frame("model" = rep(names(loo_output), 2),
                                          "measure" = rep(c("loo", "elpd"), each = length(names(loo_output))),
                                          "loo_elpd" = c(as.numeric(lapply(loo_output, function(x) x$estimates[3,1])), 
                                                         as.numeric(lapply(loo_output, function(x) x$estimates[1,1]))),
                                          "loo_elpd_sd" = c(as.numeric(lapply(loo_output, function(x) x$estimates[3,2])), 
                                                            as.numeric(lapply(loo_output, function(x) x$estimates[1,2]))))
    
    model_comp_loo_elpd_plot <- ggplot(model_comp_loo_elpd_dat, aes(x=model, y=loo_elpd, fill = measure)) + 
      geom_bar(stat="identity", position="identity", alpha = 0.8) +
      scale_fill_manual(values = LOO_ELPD_col) +
      geom_errorbar(aes(ymin=loo_elpd-loo_elpd_sd, ymax=loo_elpd+loo_elpd_sd), width=.2,
                    position="identity") +
      ylab("Information criterion (± SE)") + xlab("") +
      ggtitle(LOO_ELPD_title) +
      coord_flip() +
      labs(fill = "Measure") +
      theme(legend.position = "bottom", 
            plot.title = element_text(size = 10),
            axis.title = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10), 
            legend.key.size = unit(0.25, 'cm'))
    
    indiv_plot <- annotate_figure(model_comp_loo_elpd_plot, 
                                  top = text_grob(LOO_ELPD_title, 
                                                  face = "bold", size = 10))

    res$LOO_ELPD_plot <- model_comp_loo_elpd_plot
  }
  
  return(res)
  
}













