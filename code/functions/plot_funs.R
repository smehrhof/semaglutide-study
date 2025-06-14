################################################################################
###############----------------- PLOT FUNCTIONS -----------------###############
################################################################################

### Plotting functions to visualize data

## Acceptance rate per effort and reward level ---------------------------------
# @ task_data: Task data to generate plot for
# @ main_title: Main title for the generated plot
# @ direction: Plot effort by reward or reward by effort

task_plot <- function(task_dat, 
                      main_title = "Acceptance proportion by effort and reward level",
                      direction = "effort_by_reward", 
                      arrange_cols = 2, arrange_rows = 2){
  
  require(ggplot2)
  require(ggpubr)
  
  choiceMean_subj <- aggregate(choice ~ subjID + effort_a + amount_a, 
                               FUN = function(x) {mean(x) * 100}, 
                               dat = task_dat)
  choiceMean <- aggregate(choice ~ effort_a + amount_a, 
                          FUN = mean, 
                          dat = choiceMean_subj)
  choiceSE_subj <- aggregate(choice ~ subjID + effort_a + amount_a, 
                        FUN = function(x) {(sd(x) / sqrt(length(x))) * 100}, 
                        dat = task_dat)
  choiceSE_subj$choice[is.na(choiceSE_subj$choice)] <- 0
  choiceSE <- aggregate(choice ~ effort_a + amount_a, 
                        FUN = mean,
                        dat = choiceSE_subj)
  
  choicesPlotDat <- cbind(choiceMean, choiceSE[,3])
  colnames(choicesPlotDat) <- c("effort_a", "amount_a", "mean", "se")      
  choicesPlotDat$effort_a <- as.factor(choicesPlotDat$effort_a)
  choicesPlotDat$amount_a <- as.factor(choicesPlotDat$amount_a)
  choicesPlotDat <- choicesPlotDat[choicesPlotDat$amount_a != 999,]
  
  choicePlot <- list()
  
  if(direction == "effort_by_reward"){
    rewards <- sort(unique(task_dat$amount_a))
    for(i in 1:4){
      
      choicePlot[[i]] <- ggplot(data=choicesPlotDat[choicesPlotDat$amount_a == rewards[i],], 
                                aes(x=effort_a, y=mean, fill=effort_a)) +
        geom_bar(stat="summary", fun = "mean", alpha = 0.8) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                      position=position_dodge(.9)) +
        scale_fill_manual(values=c('#E94D36','#5B9BD5','#71AB48','#FFBF00')) +
        theme(legend.position = "none", 
              title = element_text(size = 10),
              axis.title = element_text(size = 10),
              axis.text.x = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 10),
              axis.ticks.x = element_blank()) +
        ggtitle(paste("Reward Level ", i)) +
        scale_x_discrete(labels = c("1", "2", "3", "4")) +
        scale_y_continuous(breaks = c(0, 100), limits = c(0, 100)) #+
        #xlab("Effort Level") +
        #ylab("% Accepted (SE)") 
      
    }    
  } else {
    efforts <- sort(unique(task_dat$effort_a))
    for(i in 1:4){
      
      choicePlot[[i]] <- ggplot(data=choicesPlotDat[choicesPlotDat$effort_a == efforts[i],], 
                                aes(x=amount_a, y=mean, fill=amount_a)) +
        geom_bar(stat="summary", fun = "mean", alpha = 0.8) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                      position=position_dodge(.9)) +
        scale_fill_manual(values=c('#E94D36','#5B9BD5','#71AB48','#FFBF00')) +
        theme(legend.position = "none", 
              plot.title = element_text(size = 8),
              axis.title = element_text(size = 6),
              axis.text.x = element_text(size = 6),
              axis.text.y = element_text(size = 6)) +
        ggtitle(paste("Effort Level ", i)) +
        scale_x_discrete(labels = c("2", "3", "4", "5")) +
        scale_y_continuous(breaks = c(0, 50, 100), limits = c(0, 100)) +
        xlab("Reward Level") +
        ylab("% Accepted") 
    }
  }

  
  choicePlot <- ggarrange(choicePlot[[1]] + rremove("xlab"),
                          choicePlot[[2]] + rremove("ylab") + rremove("xlab"),
                          choicePlot[[3]], 
                          choicePlot[[4]] + rremove("ylab"),
                          ncol = arrange_cols, nrow = arrange_rows,
                          common.legend = FALSE)
  
  choicePlot <- annotate_figure(choicePlot, 
                                top = text_grob(main_title, 
                                                face = "bold", size = 18),
                                left = text_grob("% Accepted (SE)", size = 10, rot = 90),
                                bottom = text_grob("Effort level", size = 10)
                                )
  
  return(choicePlot)
} 



## Raincloud plot --------------------------------------------------------------
# @ dat: Data to plot
# @ title: Main title for the generated plot
# @ xlab: x axis label (when direction = horizontal this will be the y axis)
# @ ylab: y axis label (when direction = horizontal this will be the x axis)
# @ predictor_var: name of predictor variable
# @ outcome_var: name of outcome variable
# @ predictor_tick_lab: tick labels for predictor variable
# @ col: color(s) to use
# @ direction: plot horizontally or vertically?
# @ include_grouping: plot grouping variable?
# @ group_var: (if include_grouping is TRUE) name of grouping variable
# @ legendlab: (if include_grouping is TRUE) title of group legend

raincloud_plot <- function(dat, title, xlab, ylab,
                           predictor_var, outcome_var, 
                           predictor_tick_lab, col, direction = "vertical", 
                           include_grouping = FALSE, group_var, 
                           legendlab = "", scale_seq = c(0, 6, 2)){
  require(ggplot2)
  require(PupillometryR)
  
  if(include_grouping == FALSE){
    rain_plot <- ggplot(data = dat, 
                        aes_string(y = outcome_var, x = predictor_var, fill = predictor_var)) + 
      geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.6, trim = FALSE, color = NA) +
      geom_point(aes_string(y = outcome_var, color = predictor_var), position = position_jitter(width = 0.15, height = 0.025), 
                 size = 0.75, alpha = 0.6) +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x=.2)) + 
      guides(fill = "none", color = "none") +
      labs(title = title, x = xlab, y = ylab) 
  } else {
    rain_plot <- ggplot(data = dat, 
                        aes_string(y = outcome_var, x = predictor_var, fill = group_var)) + 
      geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.6, trim = FALSE, color = NA, bw = 10) +
      geom_point(aes_string(y = outcome_var, color = group_var), position = position_jitter(width = 0.15, height = 0.025), 
                 size = 0.75, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) + 
      guides(color = guide_legend(override.aes = list(size=10)), fill = "none") +
      labs(title = title, x = xlab, y = ylab, color = legendlab) 
  }
  rain_plot <- rain_plot +
    scale_fill_manual(values = col) +
    scale_color_manual(values = col) +
    scale_x_discrete(expand = c(0.1, 0.1), labels = predictor_tick_lab) +
    theme(plot.title = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "top")
  
  if(!is.null(scale_seq)){
    rain_plot <- rain_plot +
    scale_y_continuous(breaks = seq(scale_seq[1], scale_seq[2], scale_seq[3]), 
                       limits = c(scale_seq[1], scale_seq[2]))
  }

  if(direction == "horizontal"){
    rain_plot <- rain_plot + coord_flip()
  }
  rain_plot
}


## Cummulative lag plot --------------------------------------------------------
# @ dat: Data to plot
# @ var_labels: Variable names for legend
# @ title: Main title for the generated plot
# @ x_label: Label for x-axis
# @ y_label: Label for y-axis
# @ col: Color to use


cum_plot <- function(dat, var_labels, title,
                     x_label, y_label, y_lim, col, shape = 1){
  cum_plot <- ggplot(dat, aes(x = trial, y = mean_value, 
                              group = variable, color = variable)) +
    geom_line() +
    geom_point(size = 0.75, shape = shape) +
    geom_ribbon(aes(ymin = se_lower, 
                    ymax = se_upper, 
                    x = trial, group = variable, color = variable,
                    fill = variable), 
                outline.type = "both", alpha = 0.25) +
    scale_fill_manual(values = col, name = " ",
                      labels = var_labels) +
    scale_color_manual(values = col, name = " ",
                       labels = var_labels) +
    labs(title = title,
         x = x_label, y = y_label) +
    ylim(y_lim) +
    theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) 
  
  cum_plot
}

## Posterior predictive check plot ---------------------------------------------
# @ ppc_dat: PPC data to visualize (output from "posterior_predictions" function)
# @ group_plot: Group level plot (compares observed and predicted group-level 
# acceptance proportions per effort and reward level)?
# @ indiv_plot: Individual level plot (observed vs. predicted acceptance 
# proportions on a subject level for each effort and reward level)?

ppc_plots <- function(ppc_dat, 
                      group_plot = TRUE,
                      indiv_plot = TRUE, 
                      indiv_effort = TRUE, 
                      indiv_reward = TRUE,
                      group_plot_title = "Observed vs. model predicted acceptance proportions",
                      indiv_plot_title = "Subject wise observed vs. model predicted acceptance proportions "){
  
  require(ggplot2)
  require(ggpubr)
  
  amount <- sort(unique(ppc_dat$posterior_predictions_trial_type$amount_a))
  
  res_plots <- list()
  
  # Plot 1 - group wise
  if(group_plot == TRUE){
    obs_mean <- aggregate(ppc_dat$posterior_predictions_trial_type$observation ~ 
                            ppc_dat$posterior_predictions_trial_type$effort_a + 
                            ppc_dat$posterior_predictions_trial_type$amount_a, 
                          FUN = function(x) {mean(x) * 100})
    obs_se <- aggregate(ppc_dat$posterior_predictions_trial_type$observation ~
                          ppc_dat$posterior_predictions_trial_type$effort_a + 
                          ppc_dat$posterior_predictions_trial_type$amount_a, 
                        FUN = function(x) {(sd(x) / sqrt(length(x))) * 100})
    
    observed <- cbind(obs_mean, obs_se[,3])
    colnames(observed) <- c("effort_a", "amount_a", "mean", "se")      
    observed$effort_a <- as.factor(observed$effort_a)
    observed$amount_a <- as.factor(observed$amount_a)
    
    pred_mean <- aggregate(ppc_dat$posterior_predictions_trial_type$prediction_mean ~ 
                             ppc_dat$posterior_predictions_trial_type$effort_a + 
                             ppc_dat$posterior_predictions_trial_type$amount_a, 
                           FUN = function(x) {mean(x) * 100})
    pred_se <- aggregate(ppc_dat$posterior_predictions_trial_type$prediction_mean ~ 
                           ppc_dat$posterior_predictions_trial_type$effort_a + 
                           ppc_dat$posterior_predictions_trial_type$amount_a, 
                         FUN = function(x) {(sd(x) / sqrt(length(x))) * 100})
    predicted <- cbind(pred_mean, pred_se[,3])
    colnames(predicted) <- c("effort_a", "amount_a", "mean", "se")      
    predicted$effort_a <- as.factor(predicted$effort_a)
    predicted$amount_a <- as.factor(predicted$amount_a)
    
    group_plot_dat <- rbind(observed, predicted)
    group_plot_dat <- cbind(group_plot_dat, "group" = paste(rep(c("pred", "real"), each = 16), 
                                                            rep(1:4, 4), rep(1:4, each = 4), 
                                                            sep = ""))
    
    group_plot <- list()
    
    for(i in 1:4){
      group_plot[[i]] <- ggplot(group_plot_dat[group_plot_dat$amount_a == amount[i],], aes(x=effort_a, y=mean, fill=group)) +
        geom_bar(stat='identity', position=position_dodge(), alpha = 0.9) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3,
                      position=position_dodge(0.95)) +
        scale_fill_manual(values=c('#f29a8c', '#aecdea','#add293', '#ffdf80',
                                   '#E94D36', '#5B9BD5', '#71AB48', '#FFBF00'),
                          labels = c("Predicted", "Predicted", "Predicted", "Predicted",
                                     "Observed", "Observed", "Observed", "Observed")) + 
        ylab("% Accepted") + xlab("Effort level") +
        ggtitle(paste("Reward Level", i)) +
        guides(fill=guide_legend(title=" ")) + 
        scale_x_discrete(labels= 1:4) +
        theme(plot.title = element_text(size = 10),
              axis.title = element_text(size = 10),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 10)) 
    }
    
    group_plot <- ggarrange(group_plot[[1]],
                            group_plot[[2]],
                            group_plot[[3]],
                            group_plot[[4]],
                            ncol = 2, nrow = 2,
                            common.legend = TRUE,
                            legend="bottom")
    
    group_plot <- annotate_figure(group_plot, 
                                  top = text_grob(group_plot_title, 
                                                  face = "bold", size = 14))
    
    group_plot
    
    res_plots[["group_plot"]] <- group_plot
  }
  
  # Plot 2 - subject wise
  if(indiv_plot == TRUE){
    
    # plot by effort level
    indiv_plot_effort_dat <- ppc_dat$posterior_predictions_effort
    indiv_plot_effort_dat$effort_a <- as.factor(indiv_plot_effort_dat$effort_a)
    
    R_squared_effort <- round(cor(indiv_plot_effort_dat$observation, indiv_plot_effort_dat$prediction_mean)^2, 2)
    
    indiv_plot_effort <- ggplot(indiv_plot_effort_dat, aes(x=observation, y=prediction_mean, color=effort_a)) +
      geom_point(size=2, alpha=0.5) +
      geom_errorbar(aes(ymin=prediction_hdi_lower, ymax=prediction_hdi_higher), width=.025, alpha=0.25) +
      scale_color_manual(values=c('#E94D36', '#5B9BD5', '#71AB48', '#FFBF00'), 
                         labels=1:4) + 
      xlim(0,1) + ylim(0,1) +
      geom_abline(linetype = 3) +
      ylab("Predicted (± 95% HDI)") + xlab("Observed") +
      ggtitle(bquote("Across effort levels:"~R^{2}==.(R_squared_effort))) +
      guides(color=guide_legend(title="Effort/Reward level")) +
      theme(plot.title = element_text(size = 10),
            axis.title = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10)) +
      theme(panel.background = element_blank(), axis.line.x = element_blank()) +
      annotate(geom = "segment", y = 0, yend = 1, x = -Inf, xend = -Inf) +
      annotate(geom = "segment", y = -Inf, yend = -Inf, x = 0, xend = 1)
    
    # plot by reward level
    indiv_plot_reward_dat <- ppc_dat$posterior_predictions_reward
    indiv_plot_reward_dat$amount_a <- as.factor(indiv_plot_reward_dat$amount_a)
    
    R_squared_effort <- round(cor(indiv_plot_reward_dat$observation, indiv_plot_reward_dat$prediction_mean)^2, 2)
    
    indiv_plot_reward <- ggplot(indiv_plot_reward_dat, aes(x=observation, y=prediction_mean, color=amount_a)) +
      geom_point(size=2, alpha=0.5) +
      geom_errorbar(aes(ymin=prediction_hdi_lower, ymax=prediction_hdi_higher), width=.025, alpha=0.25) +
      scale_color_manual(values=c('#E94D36', '#5B9BD5', '#71AB48', '#FFBF00'), 
                         labels=1:4) + 
      xlim(0,1) + ylim(0,1) +
      geom_abline(linetype = 3) +
      ylab("Predicted (± 95% HDI)") + xlab("Observed") +
      ggtitle(bquote("Across reward levels:"~R^{2}==.(R_squared_effort))) +
      guides(color=guide_legend(title="Effort/Reward level")) +
      theme(plot.title = element_text(size = 10),
            axis.title = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10)) +
      theme(panel.background = element_blank(), axis.line.x = element_blank()) +
      annotate(geom = "segment", y = 0, yend = 1, x = -Inf, xend = -Inf) +
      annotate(geom = "segment", y = -Inf, yend = -Inf, x = 0, xend = 1) 
    
    
    if(indiv_effort == TRUE & indiv_reward == TRUE){
      indiv_plot <- ggarrange(indiv_plot_effort,
                              indiv_plot_reward,
                              ncol = 2, nrow = 1,
                              common.legend = TRUE,
                              legend="bottom") 
    } else if (indiv_effort == TRUE & indiv_reward == FALSE){
      indiv_plot <- indiv_plot_effort
    } else if (indiv_effort == FALSE & indiv_reward == TRUE){
      indiv_plot <- indiv_plot_reward
    }
    
    indiv_plot <- annotate_figure(indiv_plot, 
                                  top = text_grob(indiv_plot_title, 
                                                  face = "bold", size = 10))
    
    
    res_plots[["indiv_plot"]] <- indiv_plot
  }
  return(res_plots)
}


## Parameter recovery plot -----------------------------------------------------
# @ recovery_data: tibble with underlying and recovered parameter values 
# @ plot_title: title for plot
# @ col: color to use

params_recovery_plot <- function(recovery_data,
                                 plot_title,
                                 col){
  

  params_recovery_plot <- ggplot(recovery_data, 
                                 aes(x=real, y=recovered, 
                                     color=col)) +
    geom_point(size=2, alpha=0.5, color = col) + 
    geom_abline(linetype = 3) +
    ylab("Recovered parameter estimates") + xlab("Underlying parameters") +
    labs(title = plot_title) +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "none") 
  
  return(params_recovery_plot)
}

