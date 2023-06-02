#Load packages-------------
library(tidymodels)  
tidymodels::tidymodels_prefer()
options(tidymodels.dark = TRUE)
# Helper packages
library(rpart)                 #needed to use using rpart.plot
library(rpart.plot)            # for visualizing a decision tree
library(vip)                   # for variable importance plots
library(bestNormalize)         #To use step_orderNorm
library(future)
library(foreach)
library(doFuture)            
library(rngtools)             #requires rng tools to call up doRNG
library(doRNG)
registerDoFuture()
library(furrr)
library(parallel)
library(iterators)
library(doParallel)
library(kernlab)
library(themis)
library(ranger)
library(kknn)
library(Matrix)
library(glmnet)
library(C50)
library(rules)
library(pROC)
library(kableExtra)
library(patchwork)
library(missForest)
library(rstatix)

#Step 1: Import dataset--------
predlum <- readRDS("Data/6_dta_symbol_remove.rds")

#Step 2.1Initial imputation-------------------------
#Ncite's Imputation script to impute OOR low and high values

#Use the "remove_symbols" function to classify missing data, add metadata variables to the data set and remove symbols.

remove_symbols <- function(dataName = "dataset_project"){
  #' @importFrom dplyr mutate
  #' @importFrom magrittr %>%
  dataName <- dataName %>%
    mutate(
      meta_truemissing=ifelse(
        obs_conc == ""
        , TRUE
        , FALSE
      )
      , meta_threestar=ifelse(
        grepl("***", obs_conc,fixed=TRUE)
        , TRUE
        , FALSE
      )
      , meta_onestar=ifelse(
        grepl("^\\*\\d+", obs_conc, perl=TRUE)
        , TRUE
        , FALSE
      )
      , meta_oorgt=ifelse(
        grepl("OOR >", obs_conc, perl=TRUE)
        , TRUE
        , FALSE
      )
      , meta_oorlt=ifelse(
        grepl("OOR <", obs_conc, perl=TRUE)
        , TRUE
        , FALSE
      )
      , obs_conc_numerics=as.numeric(
        gsub(
          "\\D*(\\d+(?:\\.\\d+)*(?:E[+-]\\d+)*)"
          ,"\\1"
          ,obs_conc,perl=T
        )
      )
    )
  
}
#Save in RDS object.
#saveRDS(dataName, "Data/rdsfold/predlum.rds")
predlum <- remove_symbols(dataName = predlum)


#Step 2.2: Imputation-----------
#Impute the OOR values, oor below (<) and oor above (>) values.
#Assuming the following variables and variable names are in your data set, e.g.,
#"analyte" containing the analyte names.
#"obs_conc_numerics" created in the previous step.

impute <- function(dataName = "dataset_project") {
  #' @importFrom dplyr arrange count group_by summarise
  #' @importFrom magrittr %>%
  #' @importFrom tidyr pivot_wider
  #' @importFrom utils write.csv
  df_rawminmaxv <- dataName %>%
    group_by(analyte) %>%
    summarise(
      Min = min(obs_conc_numerics, na.rm = TRUE),
      Max = max(obs_conc_numerics, na.rm = TRUE)
    )
  # Join experiment data and the table with the min and max analyte
  # concentrations.
  df_rawminmaxv <- left_join(dataName, df_rawminmaxv)
  # Impute OOR below (<) and OOR above (>) values.
  imp <- df_rawminmaxv %>%
    mutate(
      obs_conc_impute =
        ifelse(
          meta_oorlt == TRUE,
          Min - (Min * 0.0001),      #Min - (Min * 0.0001)
          ifelse(
            meta_oorgt == TRUE,
            Max + ((Max - Min) * 0.0001),        
            obs_conc_numerics
          )
        )
    ) %>%
    select(-Min,-Max)
  # Round Obs.Conc.Impute to two decimal places.
  imp$obs_conc_impute <- round(imp$obs_conc_impute, 2)
  return(imp)
}
predlum <- impute(dataName = predlum)


#Step 3: Select Observations only----------
#remove standards
predlum <- predlum %>% 
  filter(grepl(pattern = "X", type)) %>% 
  select(c("description", "analyte_simple", "obs_conc_impute")) |> 
  filter(description  != "ST01064373" & 
           description != "ST01061789" & 
           description  != "ST01061782"  & 
           description  != "ST01059102" & 
           description  != "ST01046593") 


#STEP 4: Create a vector with the PDs and TPs to be removed-------------
double_entry <- c("13031-W04", "13031-W08", "13031-W16", "13031-W24")
mislabeled <- c("130016-W24_7084169", "13009-W24", "13010-W08", "13010-W24", "13085-Dx", "13085-W04", 
                "13085-W08" , "13085-W16", "13085-W24", "13107-W16", "13108-W08")


#Step 5: Remove disputed observations--------------
predlum <- 
  predlum %>% 
  filter(!(description %in% double_entry )) |>      #remove concatenated double data points
  filter(!(description %in% mislabeled)) |>         #remove "mislabeled samples from the dataset
  mutate(description = ifelse(description == "15026 REC", "15026-W24", description)) |>  #Correct "15026 REC" to "15026-W24" as per luminex technician
  mutate(description = gsub('_.*', '', description)) |> 
  as.data.frame() %>% 
  replace(.=="NULL", NA)  |>   #Replace NULL characters in the dataset by coding them as NA; not enough sample volume for experiment
  rename("analyte" = "analyte_simple")


predlum <- predlum |>                                       
  unique() |> 
  pivot_wider(names_from = "analyte", values_from = "obs_conc_impute") 

predlum <- predlum |> 
  mutate(description1 = description) |> 
  select(description , description1, everything()) |> 
  separate(description1, into = c("PID", "Timepoint")) |> 
  arrange(PID) |> 
  mutate_at('PID', as.numeric)


#Step 6: Outcome from PETCT spreadsheet provided by Shawn (PredictTB)-----------
labs_petct <- readxl::read_xlsx("Data/Predict_luminex_Outcome_clinical_petct.xlsx") |> 
  rename(
    PID = SUBJID,
    HCT = LBORRES_HCT, 
    WBC = LBORRES_WBC, 
    HBA1C = LBORRES_HBA1C, 
    AST = LBORRES_AST, 
    HGB = LBORRES_HGB, 
    PLT = LBORRES_PLT, 
    CREAT = LBORRES_CREAT, 
    ALT = LBORRES_ALT, 
    bodymass = Weight, 
    Outcome = Cure_ConfRelapTF) |> 
  mutate(Outcome = ifelse(Outcome == "ConfRelapTF", "PoorOutcome", Outcome)) |>
  mutate_at("PID", as.character) |> 
  filter(PID != 13085) |> 
  select(1, 59, 3, 5, 6, 8, 11, 12, 14:23, 28, 29) |> 
  mutate(TBprev = ifelse(TBprev == 'otherwise', 'No', 'Yes')) |> 
  select(PID, Outcome, TBprev, CurrentSmoker, everything()) |> 
  mutate_at('PID', as.numeric)


#Step 6: Link observations to clinical outcome -----------
predlum <- predlum |>
  filter(Timepoint == "Dx") |> 
  left_join(labs_petct, by = "PID") 

predlum <- predlum |> 
  select(2, 54:56, everything()) |> 
  select(-c(description, Timepoint))


#Step 7: Baseline dataset MissForest Imputation for NA's--------
#The imputation is done on the dataset in the wide format
imputvars <- predlum |> 
  select(-(1:4)) |> 
  type.convert(as.is = FALSE) |>  #CONVERT CHARACTER VECTORS TO FACTORS, MISSForest rejects character vectors
  data.matrix() #mandatory transformation to a matrix for misForest


baseline_pid_nom <- predlum |>         #Store Outcome, PID, and nominal vars
  select(c(1:4)) 


set.seed(15440)
baseline_misF <- missForest(imputvars, verbose = TRUE)


data <- cbind(baseline_pid_nom, baseline_misF$ximp)
str(data)


#BASELINE MODEL-----
#STEP1: LOAD FINAL BASELINE DATASET -----------
data 


#Step2: Add BMI and TTD----------
#BMI and TTD from PETCT spreadsheet provided by Shawn (PredictTB)-----------
data <- data |> select(-PID)

vars_to_factor <- c('Outcome', 'TBprev', 'CurrentSmoker')
data[vars_to_factor] <- lapply(data[vars_to_factor], function(x) as.factor(x))
rm(vars_to_factor)
data$Outcome <- data$Outcome |> relevel(ref='PoorOutcome')

#Step 3: Corr, numeric transformation, recipe-----
##Step_corr-------
rec <- recipe(Outcome ~ ., data = data)
corr_filter <- rec |> 
  step_corr(all_numeric_predictors(), threshold = 0.9, method = "spearman") |> 
  prep(training = data) |> 
  bake(new_data = data)
data <- corr_filter


#Rank with Wilcoxin ---------
df1 <- data %>% group_by(Outcome) %>%
  summarise(across(apoa1:TLG_wk0, ~ mean(.x, na.rm = TRUE))) 
df2 <- data.frame(t(df1[-1]))
colnames(df2) <- c("Cured", "PoorOutcome")

data1 <- data |> 
  select(-c('TBprev', 'CurrentSmoker')) |> 
  select(Outcome, everything())
wilcxDx <- lapply(data1[c(-1)], function(x) wilcox.test(x ~ Outcome, data1, exact = FALSE))
wilcxDx <- map_df(wilcxDx, tidy) 
pvalDx <- df2 |> cbind(wilcxDx) |> 
  select(1, 2, 4) |> 
  arrange(p.value) 

#variables with p < 0.25--------------
#TTD's pval = 0.851, included initially because it features strongly in other models
#TTD not included in this model
#Nominal predictors: TBprev and CurrentSmoker included initially = models performed poorly
#Include TNFb, il5, and il9 since they have too many OOR low values

#variables with p < 0.25--------------
data <- data |> select(c(Outcome, BMI, tnfa, tnfb, il5,il9, mip1a, svegfr3, tnfri, CREAT, il1b, apoc3, 
                         ifng, XpertCT_wk0,il4ra, il6, CAVTOT_wk0, HCT, ip10, mmp2, il6ra, svegfr1, apoa1, 
                         il12, itac, HDV_wk0,c3, svegfr2, c3, HaveCav_wk0, c4, crp, TLG_wk0)
)

tranform_rec <- recipe(Outcome ~ ., data = data)  
transform <- step_orderNorm(tranform_rec, all_numeric_predictors())
prep_est <- prep(transform, training = data)
bake <- est <- bake(prep_est, data)
plot(density(data$tnfa), main = "before")
plot(density(bake$tnfa), main = "after")


##Initial recipe----
normalized_recipe <-  recipe(Outcome ~ ., data = data) |>    
  step_zv(all_numeric_predictors()) |> 
  step_normalize(all_numeric_predictors()) |>  
  themis::step_downsample(Outcome)


predictor_count <- sum(normalized_recipe$term_info$role == 'predictor')       
predictor_count

#Step4: Different Models--------------
C5.0_mod <-
  C5_rules(trees = tune(), min_n = tune()) |> 
  set_engine('C5.0') |> 
  set_mode("classification")
KNN_mod <-
  nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) |> 
  set_engine('kknn') |> 
  set_mode('classification')
Random_Forest_mod <-
  rand_forest(mtry = tune(), trees = tune(), min_n = tune()) |> 
  set_engine('ranger') |> 
  set_mode('classification')
Elastic_net_mod <-
  logistic_reg(penalty = tune(), mixture = tune()) |> 
  set_engine('glmnet') |> 
  set_mode("classification")


#Step5: CREATING THE WORKFLOW SET-----------
normalized_workflow <-
  workflowsets::workflow_set(
    preproc = list(normalised = normalized_recipe),
    models = list(C5.0 = C5.0_mod,
                  KNN = KNN_mod,
                  RF = Random_Forest_mod,
                  EN = Elastic_net_mod))

normalized_workflow <- normalized_workflow |> 
  mutate(wflow_id = gsub("(normalised_)|(numeric_)", "", wflow_id))           #to remove the normalised prefix


#Step6: Nested Crossvalidation---------------
set.seed(14193)
folds <- nested_cv(data,
                   outside = vfold_cv(v = 3, repeats = 3, strata = "Outcome"),   
                   inside = bootstraps(times = 20, strata = "Outcome"))


#Step5: Model parameters and workflows------------
C5.0_params <- parameters(trees(range = c(1,100)), min_n())
RF_params <- parameters(mtry(range=c(1, predictor_count)), trees(range = c(1,100)), min_n()) 
KNN_params <- parameters(neighbors(), weight_func(), dist_power()) |> 
  recipes::update(weight_func = weight_func(c("triangular", "biweight", "triweight", "cos", "gaussian", "rank", "optimal")))
EN_params <- parameters(penalty(), mixture())


#Step6: Workflows---------
workflows <- normalized_workflow |>
  option_add(param_info = C5.0_params, id = "C5.0") |> 
  option_add(param_info = KNN_params, id = "KNN") |> 
  option_add(param_info = RF_params, id = "RF") |> 
  option_add(param_info = EN_params, id = "EN")


#Step7: Tuning-------------
bayes_ctrl <- control_bayes(no_improve = 15L, 
                            save_pred = TRUE, 
                            parallel_over = "everything",    # replaced "everything" with "resamples" -  success when parallel =  "resamples" instead of parallel =  "everything" is used 
                            save_workflow = TRUE, 
                            allow_par = TRUE, 
                            verbose = TRUE)


#tune_results <- foreach(i=1:length(folds$splits)) %dorng% {
tune_results <- foreach(i=1:length(folds$splits)) %do% {
  library(rules)
  workflows %>% workflow_map(seed = 15440,
                             fn = "tune_bayes",
                             resamples = folds$inner_resamples[[i]],
                             metrics = metric_set(roc_auc),
                             objective = exp_improve(),
                             iter = 50,
                             control = bayes_ctrl)
}


saveRDS(tune_results, "Data/tune_results_baseline_petct2_transf.rds")


#Step8: Tune results object----------
tune_results <- readRDS("Data/tune_results_baseline_petct2_transf.rds")



#Step9: Updated Normalized recipe----------
normalised_training_recipe <- 
  recipe(Outcome ~ ., 
         data = data) |>
  step_zv(all_numeric_predictors()) |> 
  step_dummy(all_nominal_predictors())|> 
  step_normalize(all_numeric_predictors()) |>
  step_smote(Outcome) 


#Step10: Fit models function-------
fit_models <- function(model, grid, data, type){
  
  best_result <-  grid |> 
    show_best(n=1)
  
  model_fit <-
    grid |> extract_workflow(model) |>
    finalize_workflow(best_result) |>
    update_recipe(normalised_training_recipe) |> 
    fit(data=analysis(data))
  
  pred <-
    predict(model_fit, assessment(data)) |> 
    bind_cols(predict(model_fit, assessment(data), type = "prob")) |> 
    bind_cols(assessment(data) |> select(Outcome))
  
  roc_auc <- pred |> roc_auc(truth = Outcome, .pred_PoorOutcome )
  sens <- pred|> sensitivity(Outcome, .pred_class)
  spec <- pred |> specificity(Outcome, .pred_class)
  
  perf <- tibble(Model = model, data ="train") |> 
    mutate(auc = roc_auc$.estimate, sens =  sens$.estimate, spec = spec$.estimate)
  
  
  if (type==1){
    return(perf)
  }
  else if(type==2){
    pred <- pred |> mutate(Resamples = data$id)
    return(pred)
  }
  else if(type==3){
    return(best_result)
  }
}


#Step9:Fit------------
training_results <- foreach(x=1:length(folds$splits)) %do% {
  model_results <- foreach(y=1:length(workflows$wflow_id)) %do% {
    fit_models(tune_results[[x]]$wflow_id[[y]], tune_results[[x]]$result[[y]], folds$splits[[x]], 1)
  }
  bind_rows(model_results)
}

training_average <- foreach(x=1:length(workflows$wflow_id), .combine=rbind) %dorng% {
  bind_rows(training_results) |> filter(Model==workflows$wflow_id[[x]]) |>
    group_by(Model, data) |> summarise_if(is.numeric, list(mean = mean, sd = sd), na.rm=TRUE)
}

training_average[3:8] <- training_average[3:8] |> round(3)

new_training_average <- foreach(x=1:length(workflows$wflow_id)) %do% {
  tibble(Model = workflows$wflow_id[[x]], data = "train") |>
    mutate(auc = paste0(training_average$auc_mean[[x]], "\u00B1", training_average$auc_sd[[x]]),
           sens = paste0(training_average$sens_mean[[x]], "\u00B1", training_average$sens_sd[[x]]),
           spec = paste0(training_average$spec_mean[[x]], "\u00B1", training_average$spec_sd[[x]]))
} |> bind_rows()

arrange(new_training_average, desc(auc))|> 
  kable(align=rep('c')) |> 
  kable_classic(full_width = F)

#Step10: Predictions----------
prediction_results <- foreach(x=1:length(folds$splits)) %do% {
  model_perf <- foreach(y=1:length(workflows$wflow_id)) %do% {
    fit_models(tune_results[[x]]$wflow_id[[y]], tune_results[[x]]$result[[y]], folds$splits[[x]], 2)
  }
}

roc_curves_folds <- foreach(x=1:length(workflows$wflow_id)) %do% {
  model_pred <- foreach(y=1:length(folds$splits)) %do% {
    prediction_results[[y]][[x]] |> 
      unnest(cols = c(Resamples))|> 
      rename(Resamples=id)
  }
  bind_rows(model_pred) |> 
    group_by(Resamples) |> 
    roc_curve(Outcome, .pred_PoorOutcome) |> 
    autoplot() + 
    labs(title = workflows$wflow_id[[x]]) + 
    theme_update(plot.title = element_text(hjust = 0.5)) +
    theme_bw() +
    theme(legend.position="left")        #alternative legend = "none
}


roc_curves_folds[[1]] + roc_curves_folds[[2]] + roc_curves_folds[[3]] + roc_curves_folds[[4]] + 
  plot_layout(nrow = 2, byrow = FALSE)


roc_curves <- foreach(x=1:length(workflows$wflow_id)) %do% {
  model_pred <- foreach(y=1:length(folds$splits), .combine = rbind) %do% {
    prediction_results[[y]][[x]] |> 
      unnest(cols = c(Resamples))|> 
      rename(Resamples=id)
  }
  pROC::roc(Outcome ~ .pred_PoorOutcome, data=model_pred, auc=T,levels=c('Cured','PoorOutcome'), ci=TRUE, of = "se", ci.type = "bars", ci.method = "bootstrap", boot.n = 5000, parallel = TRUE, plot=FALSE)
}

par(mfrow=c(2,2))
for (i in 1:4) {
  plot(roc_curves[[i]], main=workflows$wflow_id[i], print.auc=T, grid=c(0.1, 0.2),
       print.thres=TRUE)
}

par(mfrow=c(1,1))       #if you want individual plots
for (i in 1:4) {
  plot(roc_curves[[i]], main=workflows$wflow_id[i], print.auc=T)
}


#?NB variables in the final model-----------BEST MODEL
RF <- tune_results[[3]]$result[[3]] |> 
  extract_workflow(tune_results[[3]]$wflow_id[[3]]) |> 
  finalize_workflow(show_best(tune_results[[3]]$result[[3]],n=1)) |> 
  fit(data=analysis(folds$splits[[1]]))  #?[[i]]


