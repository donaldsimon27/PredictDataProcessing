#Load packages----------
library(tidyverse)
library(magrittr)
library(utils)


#Set working directory---------
setwd("~/Data")


#Step 1: Import dataset--------
#Textfiles read using Jesse's pipeline
predlum <- readRDS("~/Data/6_dta_symbol_remove.rds")
dim(predlum)


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
#saveRDS(dataName, "~/Desktop/R.Projects/PredictLuminex/Data/rdsfold/predlum.rds")
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
  saveRDS(imp, "~/Data/predlum.rds")
  write.csv(imp, file = "~/Data/predlum.csv")
}
impute(dataName = predlum)


#Step 3: Select Observations only----------
#remove standards
predlum <- read.csv("~/Data/predlum.csv") %>% 
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
  arrange(PID)


#Step 6: Outcome from PETCT spreadsheet provided by Shawn (PredictTB)-----------
Outcome_PID <- readxl::read_xlsx("~/Data/Predict_luminex_Outcome_clinical_petct.xlsx") |> 
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
  select(PID, Outcome) |> 
  mutate(Outcome = ifelse(Outcome == "ConfRelapTF", "PoorOutcome", Outcome)) |> 
  filter(PID != 13085)
 


#Step 6: Link observations to clinical outcome -----------
predlum <- predlum |> 
  mutate_at("PID", as.numeric) |> 
  left_join(Outcome_PID, by = "PID") 

predlum <- predlum |> 
  select(1:3, 54, everything())


#Step 7: Baseline dataset MissForest Imputation for NA's--------
#The imputation is done on the dataset in the wide format
library(missForest)
set.seed(15440)
baseline_analytes <- predlum |> 
  filter(Timepoint == "Dx") |> 
  select(-c(1, 3))

baseline_outcome <- baseline_analytes |> 
  select(1, 2)

baseline_misF <- type.convert(baseline_analytes, as.is = FALSE) %>%   #CONVERT CHARACTER VECTORS TO FACTORS, MISSForest rejects character vectors
  data.matrix() %>%                                    #mandatory transformation to a matrix for misForest
  missForest(verbose = TRUE)

data <- baseline_misF$ximp |> as_tibble() |> 
  left_join(baseline_outcome, by = "PID")

data <- data |> 
  select(1, 2, 53, everything()) |> 
  select(-Outcome.x) |> 
  rename(Outcome  = "Outcome.y") |> 
  select(PID, Outcome, everything())

