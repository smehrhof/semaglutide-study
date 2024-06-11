# A study on the effect of GLP-1 antagonis based injections of semaglutide on mental health and cognition

Sara Z. Mehrhof, Hugo Fleming & Camilla L. Nord

# Data

The cleaned data used in our study can be found in the `/data` directory.  
Data cleaning entailed subject anonymisation, processing of questionnaires according to questionnaire manuals, and applying pre-registered exclusion criteria.   

# Code

All code used to run analyses for our study can be found in the `/code` directory. Some scripts may call supporting functions and stan scripts found in the `/functions` and `/stan` directories.  

### Analyses

The `/analyses` directory contains all main study scripts, structured as follows:

#### 1 Screening

In this script, we parse screening data and identify participants eligible for the treatment and control group of the main study. We then randomize testing schedule to participants eligible for the treatment group and identify their testing days. Control participants are matched to treatment participants by age, gender, physical activity (as measured with the International Physical Activivty Questionnaire (IPAQ)) and BMI. 

#### 2 Task processing

In this script, we parse the data collected during main testing sessions. Data is merged with data collected during screening and exclusion criteria are applied.   
Control participants for the two non-diabetic groups are matched from an existing dataset. 

#### 3 Descriptives

In this script, we compute sample descriptives, including demographic data, psychiatric comorbidities and medications. We also confirm the replication of model agnostic effort discounting effects in all groups. 

#### 4 Model fitting

In this script, we fit all models included in our pre-registered model space to the data, by group. We perform convergence checks and model comparisons. The winning model is then validated by posterior predictive checks. 

#### 5 Primary analyses

In this script, we perform the main analyses described in our pre-registration (https://osf.io/rmz45/). This includes group comparisons between type-2 diabetic participants on and off semaglutide of model based task parameters and questionnaire variables.  

#### 6 Group comparison including non-diabetic controls

In this script, we perform an exploratory analysis described in our pre-registration (https://osf.io/rmz45/), including group comparisons of model based task parameters and questionnaire variables between type-2 diabetic participants on and off semaglutide, and non-diabetic controls with a matched BMI or a BMI restricted to 18.5-25. 

#### 7 Within subject comparison

In this script, we perform our pre-registered exploratory analysis of a within-subject analysis, comparing task and questionnaire measures of testing sessions one day vs. six days after treatment injection. 








