#################################################################################################
######---------------- Comparing diabetes groups with non-diabetic controls ----------------#####
#################################################################################################

### In this script: 
# (1) Prepare data
# (2) Task data - model based
# (3) Questionnaire data

# Set working directory
here::i_am("github/semaglutide-study/code/analyses/6_non_diabetic_comparison.R")
setwd(here::here())

# source functions
source("github/semaglutide-study/code/functions/helper_funs.R")
source("github/semaglutide-study/code/functions/plot_funs.R")

# source dataset
main_data <- readRDS("github/semaglutide-study/data/processed_data/main_data.RDS")
non_diabetic_data <- readRDS("github/semaglutide-study/data/processed_data/non_diabetic_matched.RDS")
non_diabetic_normal_weight_data <- readRDS("github/semaglutide-study/data/processed_data/non_diabetic_normal_weight_matched.RDS")

# source parameter estimates
m3_para_treat_s1_params <- readRDS(here::here("github/semaglutide-study/data/model_fits/treatment_s1/m3_para_treat_s1_params.RDS"))

m3_para_control_params <- readRDS(here::here("github/semaglutide-study/data/model_fits/controls/m3_para_control_params.RDS"))

m3_para_non_diabetic_params <- readRDS(here::here("github/semaglutide-study/data/model_fits/non_diabetics/m3_para_non_diabetics_params.RDS"))

m3_para_normal_weight_params <- readRDS(here::here("github/semaglutide-study/data/model_fits/low_bmi/m3_para_normal_weight_params.RDS"))

# load required packages
librarian::shelf(ggplot2, ggpubr, tidyverse, dplyr, stringr, purrr, here, janitor, MatchIt, PupillometryR,
                 writexl, lubridate, magrittr, pushoverr, nlme, gridExtra, cmdstanr, rstanarm, bayestestR, rstatix, DescTools)

# Color pallet 
color_pal <- c("#E94D36", "#5B9BD5", "#71AB48", "#FDC219", "#8456B8", "#FF7236", "#1FD5B3", "#F781BE")

### (1) Prepare data  -----------------------------------------------

# Merge datasets for analyses
data <- bind_rows(
  # Diabetics
  main_data$demographic_data %>% 
  select(subj_id, group, age, gender, bmi, psych_neurdev, antidepressant) %>%
    mutate(antidepressant = ifelse(is.na(antidepressant), 0, antidepressant)) %>% 
  left_join(main_data$glp_data %>% 
              filter(session == 1) %>% 
              select(subj_id, start_date_glp, side_effects_glp, glp_dose_mg, hours_since_injection, local_testing_time), 
            by = "subj_id") %>% 
  left_join(main_data$questionnaire_data %>% 
              filter(session == 1) %>% 
              select(subj_id, aes_sumScore, findrisc_sumScore, mcq_discounting_rate, mctq_MSF_SC, meq_sumScore, 
                      shaps_sumScore, ipaq_sumScore),
            by = "subj_id") %>% 
  left_join(rbind(m3_para_treat_s1_params$individual_params %>% 
                    pivot_wider(id_cols = subj_id, names_from = parameter, 
                                values_from = c(estimate, hdi_lower, hdi_upper)), 
                  m3_para_control_params$individual_params %>% 
                    pivot_wider(id_cols = subj_id, names_from = parameter, 
                                values_from = c(estimate, hdi_lower, hdi_upper))), 
            by = "subj_id") %>% 
    # make MCTQ result numeric (minutes since 00:00)
    add_column(mctq_continuous = period_to_seconds(hm(.$mctq_MSF_SC))/60, 
               .before = "mctq_MSF_SC"), 
  
  # Non-diabetics, matched by BMI
  non_diabetic_data$demographic_data %>% 
    select(subj_id, age, gender, bmi, psych_neurdev, antidepressant) %>%
    mutate(antidepressant = ifelse(is.na(antidepressant), 0, antidepressant)) %>% 
    left_join(non_diabetic_data$questionnaire_data %>% 
                select(subj_id, aes_sumScore, findrisc_sumScore, mctq_MSF_SC, meq_sumScore, 
                       shaps_sumScore, ipaq_sumScore),
              by = "subj_id") %>% 
    left_join(rbind(m3_para_non_diabetic_params$individual_params %>% 
                      pivot_wider(id_cols = subj_id, names_from = parameter, 
                                  values_from = c(estimate, hdi_lower, hdi_upper))), 
              by = "subj_id") %>% 
    # make MCTQ result numeric (minutes since 00:00)
    add_column(mctq_continuous = period_to_seconds(hm(.$mctq_MSF_SC))/60, 
               .before = "mctq_MSF_SC") %>% 
    add_column(group = "non_diabetic", 
               .after = "subj_id"),
  
  # Non-diabetics, normal weight BMI
  non_diabetic_normal_weight_data$demographic_data %>% 
    select(subj_id, age, gender, bmi, psych_neurdev, antidepressant) %>%
    mutate(antidepressant = ifelse(is.na(antidepressant), 0, antidepressant)) %>% 
    left_join(non_diabetic_normal_weight_data$questionnaire_data %>% 
                select(subj_id, aes_sumScore, findrisc_sumScore, mctq_MSF_SC, meq_sumScore, 
                       shaps_sumScore, ipaq_sumScore),
              by = "subj_id") %>% 
    left_join(rbind(m3_para_normal_weight_params$individual_params %>% 
                      pivot_wider(id_cols = subj_id, names_from = parameter, 
                                  values_from = c(estimate, hdi_lower, hdi_upper))), 
              by = "subj_id") %>% 
    # make MCTQ result numeric (minutes since 00:00)
    add_column(mctq_continuous = period_to_seconds(hm(.$mctq_MSF_SC))/60, 
               .before = "mctq_MSF_SC") %>% 
    add_column(group = "normal_weight", 
               .after = "subj_id")
) %>% 
  mutate(chronotype = case_when(meq_sumScore > 58 & hm(mctq_MSF_SC) < hm("02:30") ~ "early",
                                meq_sumScore < 42 & hm(mctq_MSF_SC) > hm("05:30") ~ "late",
                                meq_sumScore > 58 & is.na(mctq_MSF_SC)  ~ "early",
                                meq_sumScore < 42 & is.na(mctq_MSF_SC)  ~ "late",
                                .default = "intermediate")) 

# Parameters
# Scale parameters to be between 0 and 1
data %<>% 
  ungroup %>% 
  mutate(across(c(estimate_kE, estimate_kR, estimate_a), rescale)) 

# Visualize distributions
ggplot(gather(data %>% select(c(estimate_kE:estimate_a))) %>% na.omit(), 
       aes(value)) + 
  geom_histogram(bins = 7) + 
  facet_wrap(~key, scales = 'free_x')

data %<>% 
  # positively skewed
  mutate_at(c("estimate_kE"), sqrt) %>% 
  # negatively skewed
  mutate_at(c("estimate_a"), norm_neg_skew)

# Questionnaires
# Scale parameters to be between 0 and 1
data %<>% 
  ungroup %>% 
  mutate_at(colnames(data)[c(5, 13:16, 19:20)], rescale)

# Visualize distributions
ggplot(gather(data %>% select(colnames(data)[c(5, 13:16, 19:20)])) %>% na.omit(), 
       aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

data %<>% 
  # positively skewed
  mutate_at(c("bmi", "ipaq_sumScore", "mctq_continuous", "shaps_sumScore"), sqrt) %>% 
  # negatively skewed
  mutate_at(c("aes_sumScore", "findrisc_sumScore"), norm_neg_skew)


### (2) Task data - model based -----------------------------------------------

data %<>% 
  mutate(group = factor(group, levels = c("normal_weight", "non_diabetic", "control", "treatment"))) %>% 
  mutate(antidepressant = as.factor(antidepressant)) %>% 
  mutate(psych_neurdev = as.factor(psych_neurdev))

### Effort sensitivity
# One-way ANOVA
kE_aov <- aov(estimate_kE ~ group, data = data)
summary(kE_aov)

### Reward sensitivity
# One-way ANOVA
kR_aov <- aov(estimate_kR ~ group, data = data)
summary(kR_aov)

### Choice bias
# One-way ANOVA
a_aov <- aov(estimate_a ~ group, data = data)
summary(a_aov)

# controlling for current or past psych disorder
a_aov_psych <- aov(estimate_a ~ group + psych_neurdev, data = data)
summary(a_aov_psych)

# controlling for antidepressant use
a_aov_antid <- aov(estimate_a ~ group + antidepressant, data = data)
summary(a_aov_antid)

# Mutiple comparisons
# HC, bmi matched & HC, low bmi
var.test(estimate_a ~ group, data = data %>% 
         filter(group == "non_diabetic" | group == "normal_weight"), 
       paired = FALSE)
t.test(estimate_a ~ group, data = data %>% 
         filter(group == "non_diabetic" | group == "normal_weight"), 
      var.equal = TRUE)

# T2D, off S & HC, low bmi
var.test(estimate_a ~ group, data = data %>% 
           filter(group == "control" | group == "normal_weight"), 
         paired = FALSE)
t.test(estimate_a ~ group, data = data %>% 
         filter(group == "control" | group == "normal_weight"), 
        var.equal = FALSE)

# T2D, on S & HC, low bmi
var.test(estimate_a ~ group, data = data %>% 
         filter(group == "treatment" | group == "normal_weight"), 
       paired = FALSE)
t.test(estimate_a ~ group, data = data %>% 
         filter(group == "treatment" | group == "normal_weight"), 
       var.equal = FALSE)

# T2D, off S & HC, matched bmi
var.test(estimate_a ~ group, data = data %>% 
         filter(group == "control" | group == "non_diabetic"), 
       paired = FALSE)
t.test(estimate_a ~ group, data = data %>% 
         filter(group == "control" | group == "non_diabetic"), 
       paired = FALSE, var.equal = FALSE)

# T2D, on S & HC, matched bmi
var.test(estimate_a ~ group, data = data %>% 
         filter(group == "treatment" | group == "non_diabetic"), 
       paired = FALSE)
t.test(estimate_a ~ group, data = data %>% 
         filter(group == "treatment" | group == "non_diabetic"), 
      var.equal = FALSE)

# T2D, on S & HC, T2D off S
t.test(estimate_a ~ group, data = data %>% 
         filter(group == "treatment" | group == "control"))


### Plot

a_plot <- raincloud_plot(dat = data, title = "", 
                         xlab = " ", ylab = "Acceptance bias", 
                         predictor_var = "group", outcome_var = "estimate_a", 
                         predictor_tick_lab = c("Non-diabetics\nBMI 18.5-25\n\n\n", "Non-diabetics\nBMI matched\n\n\n", "Type-2 diabetics\noff semaglutide\n\n\n", "Type-2 diabetics\non semaglutide\n\n\n"), 
                         col = c(color_pal[4], color_pal[3], color_pal[1], color_pal[2]), 
                         include_grouping = FALSE, direction = "horizontal", scale_seq = c(0,1, 0.5)) + 
  theme(legend.position = "none") +
  ggtitle("") +
 # theme(axis.text.y = element_text(angle = 25)) +
  theme(panel.background = element_blank(), axis.line.x = element_blank()) +
  theme(axis.ticks.y=element_blank()) +
  annotate(x=0.5, xend=0.5, y=-4, yend=8, colour="black", lwd=0.75, geom="segment")

a_plot

ggplot2::ggsave(filename = "fig_5.png", path = "/Users/saramehrhof/Desktop/", 
                dpi = 800, device = "png", width = 6, height = 4)


### FINDRISC predicting computational parameters
findrisc_glm <- glm(estimate_kE ~ findrisc_sumScore, data = data, family = "gaussian")
summary(findrisc_glm)

findrisc_glm <- glm(estimate_kR ~ findrisc_sumScore, data = data, family = "gaussian")
summary(findrisc_glm)

findrisc_glm <- glm(estimate_a ~ findrisc_sumScore, data = data, family = "gaussian")
summary(findrisc_glm)
EtaSq(findrisc_glm)  

# plot
a_findrisc_plot <- ggplot(data, aes(x = findrisc_sumScore, y = estimate_a)) + 
  geom_point(color = color_pal[5], alpha = 0.5, size = 1.25) +
  geom_smooth(method = glm, color = color_pal[5], fill = color_pal[5]) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0.13, 0.87), limits = c(0, 1), 
                     labels = c("low FINDRISC score", "high FINDRISC score")) +
  labs(y = "Acceptance bias") +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10)) +
        #axis.ticks.x = element_blank()) +
  theme(panel.background = element_blank(), axis.line.x = element_blank()) +
  ggtitle("") +
  annotate(geom = "segment", x = 0, xend = 1, y = -Inf, yend = -Inf) +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = 0, yend = 1)
  

a_findrisc_plot

ggplot2::ggsave(filename = "fig_6.png", path = "/Users/saramehrhof/Desktop/", 
                dpi = 800, device = "png", width = 4.2, height = 4)

pdf(file = here::here("output/figures/R_plots/results_all.pdf"),  
    width = 9, # The width of the plot in cm (transformed to inches)
    height = 6) # The height of the plot in cm (transformed to inches)


ggarrange(a_plot,
          a_findrisc_plot,
          ncol = 2, nrow = 1, 
          widths = c(1.6, 1.4))


dev.off()


# Controlling for current or past psych disorder
findrisc_antid_glm <- glm(estimate_a ~ findrisc_sumScore + psych_neurdev, data = data, family = "gaussian")
summary(findrisc_antid_glm)
EtaSq(findrisc_antid_glm) 

# Controlling for antidepressant use
findrisc_antid_glm <- glm(estimate_a ~ findrisc_sumScore + antidepressant, data = data, family = "gaussian")
summary(findrisc_antid_glm)
EtaSq(findrisc_antid_glm) 


### (3) Questionnaire data -----------------------------------------------

### Psychiatric questionnaires
data %>% 
  group_by(group) %>%
  summarise(mean_aes = mean(aes_sumScore), sd_aes = sd(aes_sumScore),
            mean_shaps = mean(shaps_sumScore), sd_shaps = sd(shaps_sumScore))

### AES
aes_aov <- aov(aes_sumScore ~ group, data = data)
summary(aes_aov)

### SHAPS
shaps_aov <- aov(shaps_sumScore ~ group, data = data)
summary(shaps_aov)


### Circadian questionnaires
### MCTQ
data %>% 
  group_by(group) %>%
  summarise(mean_mctq = mean(mctq_continuous, na.rm = TRUE), 
            sd_mctq = sd(mctq_continuous, na.rm = TRUE))

mctq_aov <- aov(mctq_continuous ~ group, data = data)
summary(mctq_aov)
TukeyHSD(mctq_aov)

mctq_findrisc_aov <- glm(mctq_continuous ~ findrisc_sumScore, data = data, family = "gaussian")
summary(mctq_findrisc_aov)

### MEQ
data %>% 
  group_by(group) %>%
  summarise(mean_meq = mean(meq_sumScore, na.rm = TRUE), 
            sd_meq = sd(meq_sumScore, na.rm = TRUE))
meq_aov <- aov(meq_sumScore ~ group, data = data)
summary(meq_aov)
TukeyHSD(meq_aov)

meq_findrisc_aov <- glm(meq_sumScore ~ findrisc_sumScore, data = data, family = "gaussian")
summary(meq_findrisc_aov)


# Both chronotype questionnaires indicate diabetic patients on treatment are have a higher evening preference
# now look at chronotypes
data %<>% 
  mutate(chronotype = case_when(meq_sumScore > 58 & hm(mctq_MSF_SC) < hm("02:30") ~ "early",
                                meq_sumScore < 42 & hm(mctq_MSF_SC) > hm("05:30") ~ "late",
                                is.na(mctq_MSF_SC) ~ "NA",
                                .default = "intermediate")) 

# Groups
data %>% 
  filter(!is.na(mctq_MSF_SC)) %>% 
  tabyl(group, chronotype) %>% 
  chisq.test()

# FINDRISC
data %<>% 
  filter(!is.na(mctq_MSF_SC)) %>%
  mutate(chronotype_num = case_when(chronotype == "early" ~ 0, 
                                    chronotype == "late" ~ 2, 
                                    chronotype == "intermediate" ~ 1))

chronotype_findrisc_aov <- glm(chronotype_num ~ findrisc_sumScore, data = data, family = "gaussian")
summary(chronotype_findrisc_aov)



