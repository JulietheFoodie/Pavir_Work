
################################################################################################
#  Project Name :           SNF DID Analysis				 																		                                              	#
#  Principal Investigator : Burke																										                                                    #
#  Name of Program :        10_ControlVarCreation080724.R		      												                                              #
#  Programmer :             J Norman                                                                         							              #
#  Start Date :             08/07/24        	                                                               								            #
#  Related programs : 		/project/Burke_SNF_VBP/Analysis_2022/prog/DID_shared			  									                                #
#  Data Used:	   	  		 /project/Burke_SNF_VBP/Analysis_2022/prog/Heintz/Analytical Datasets/Final/BENE_TEMP.sas7bdat                             #
#  Data Used:	   	  		/project/Burke_SNF_VBP/Analysis_2022/prog/Norman/Crosswalk Data/ICDallxCCSHCC_cwAll.rds							            #
#  Data Used:	   	  		/project/Burke_SNF_VBP/Analysis_2022/prog/Norman/Crosswalk Data/ICD10xCCSedit.csv										            #
#  Final Output Dataset :				                                                                                                        #
#  Program Description :	Create control variables for model execution	    					                                                  #
#################################################################################################


#################################################################################################
### Step 0:  Upload Libraries and Cohort Data set ############################# 
#################################################################################################

library(tidyverse)
library(haven) #sas upload
library(lubridate) #year variable
library(lme4) #model
#library(reshape2) #hcc label cleaning
library(readxl) #crosswalk upload


coh_df <- read_sas("/project/Burke_SNF_VBP/Analysis_2022/prog/Heintz/Analytical Datasets/Final/BENE_FINAL")


#################################################################################################
### Step 1: Initial  Variable Creation #############################
#################################################################################################

coh_df1 <- coh_df %>% 
  mutate(
    acute90_cat = as.character(count90_PriorHospStays_cat),
    los_cat = as.character(Prox_LOS_cat),
    ICU_Ind = as.character(Prox_Intnsv_care),
    ICU_IndAll = as.character(anyICUdays_rev),
    covid_cat = as.character(covid_rf),
    # Disabled Eligibility
    Dis_Ind = case_when(
      OREC == 1 ~ 1,
      TRUE ~ 0
    ), # ESRD Eligibility
    ESRD_Ind = ESRD
    , #age/sex
    agesex_cat = case_when(
      m0_34 == 1 ~ "m0_34",
      m35_44 == 1 ~ "m35_44",
      m45_54 == 1 ~ "m45_54",
      m55_59 == 1 ~ "m55_59",
      m60_64 == 1 ~ "m60_64",
      m65_69 == 1 ~ "m65_69",
      m70_74 == 1 ~ "m70_74",
      m75_79 == 1 ~ "m75_79",
      m80_84 == 1 ~ "m80_84",
      m85_89 == 1 ~ "m85_89",
      m90_94 == 1 ~ "m90_94",
      m95_GT == 1 ~ "m95_GT",
      w0_34 == 1 ~ "w0_34",
      w35_44 == 1 ~ "w35_44",
      w45_54 == 1 ~ "w45_54",
      w55_59 == 1 ~ "w55_59",
      w60_64 == 1 ~ "w60_64",
      w65_69 == 1 ~ "w65_69",
      w70_74 == 1 ~ "w70_74",
      w75_79 == 1 ~ "w75_79",
      w80_84 == 1 ~ "w80_84",
      w85_89 == 1 ~ "w85_89",
      w90_94 == 1 ~ "w90_94",
      w95_GT == 1 ~ "w95_GT" )
  ) %>% 
  dplyr::select(-starts_with("m"), -starts_with("w"))

#################################################################################################
### Step 2: Create Procedure Variables #############################
#################################################################################################

#### Step 2A. Upload Crosswalk Data and Clean for Crosswalk Creation #####

proc_df <- coh_df1 %>% 
  select(contains("PROX_PRCDR_CD"), contains("PROX_SRGCL_PRCDR_VRSN_CD_"))

names(proc_df) <-toupper(names(proc_df))


# Create List of all ICD procedure codes in Cohort
proc_list <- list()

for (x in 1:25) {
  
  # create 0 leading number for 01-25   
  dig2 <- ifelse(nchar(as.character(x)) == 1, sprintf("%02d", x), x)
  # create variables
  icd_var <- paste0("PROX_PRCDR_CD", dig2)
  typ_var <- paste0("PROX_SRGCL_PRCDR_VRSN_CD_", dig2)
  
  proc_sel <- proc_df %>% 
    select(c(icd_var, typ_var)) %>% 
    filter(is.na(icd_var) == FALSE | is.na(typ_var) == FALSE) %>% 
    distinct() %>% 
    mutate(Diag_n = x)
  
  names(proc_sel) <- gsub(".{2}$", "", names(proc_sel))
  
  proc_list[[x]] <- proc_sel
  
}

procLng_df <- bind_rows(proc_list) %>% 
  select(-Diag_n) %>% 
  mutate(ICD = gsub("[^a-zA-Z0-9]", "", as.character(PROX_PRCDR_CD))) %>% 
  rename(ICD_typ = PROX_SRGCL_PRCDR_VRSN_CD_) %>% 
  select(ICD, ICD_typ) %>% 
  distinct() 

#### ICD to CCS Crosswalk ###

# ICD9
Proc_ICD9_CCS_raw <- read_csv("/project/Burke_SNF_VBP/Analysis_2022/prog/Norman/Crosswalk Data/Raw/CCSR_ICD9_Mapping.csv", skip = 1)

# remove apostrophes
Proc_ICD9_CCS_raw[] <-  lapply(Proc_ICD9_CCS_raw, gsub, pattern = "'", replacement = '')

names(Proc_ICD9_CCS_raw) = gsub(pattern = "'", replacement = "", x = names(Proc_ICD9_CCS_raw))

Proc_ICD9_CCS <- Proc_ICD9_CCS_raw %>% 
  mutate(ICD = gsub("[^a-zA-Z0-9]", "", as.character(`ICD-9-CM CODE`)),
         CCS9_num = as.numeric(`CCS CATEGORY`),
         ICD_typ = 9)  %>% 
  select(ICD, CCS9_num) %>% 
  distinct()

# Partial Match ICD9
Proc_ICD9ct_CCS <- Proc_ICD10_CCS_raw %>% 
  mutate(ICD = gsub("[^a-zA-Z0-9]", "", as.character(`ICD-10-PCS CODE`)),
         CCS9_num2 = as.numeric(`CCS CATEGORY`),
         ICD_typ = 9) %>% 
  mutate(ICD_ct = substr(ICD, 1, 4)) %>% 
  select(ICD_ct, CCS9_num2) %>% 
  distinct() %>% 
  mutate(is_dup = duplicated(ICD_ct) | duplicated(ICD_ct, fromLast = TRUE)) %>% 
  filter(is_dup == FALSE) %>% 
  select(-is_dup) %>% 
  distinct()

# ICD10

Proc_ICD10_CCS_raw <- read_excel("/project/Burke_SNF_VBP/FOIA/CCS proc codes Yale2020 PRA_2-15-2022_RTI fixs errors.xlsx")

Proc_ICD10_CCS <- Proc_ICD10_CCS_raw %>% 
  mutate(ICD = toupper(gsub("[^a-zA-Z0-9]", "", as.character(ICD10PCSCODE))),
         CCS10_num = as.numeric(CCSCATEGORY),
         ICD_typ = 0) %>% 
  select(ICD, CCS10_num) %>% 
  distinct()

# Partial Match ICD10
Proc_ICD10ct_CCS <- Proc_ICD10_CCS_raw %>% 
  mutate(ICD = gsub("[^a-zA-Z0-9]", "", as.character(ICD10PCSCODE)),
         CCS10_num2 = as.numeric(CCSCATEGORY),
         ICD_typ = 0) %>% 
  mutate(ICD_ct = substr(ICD, 1, 4)) %>% 
  select(ICD_ct, CCS10_num2) %>% 
  distinct() %>% 
  mutate(is_dup = duplicated(ICD_ct) | duplicated(ICD_ct, fromLast = TRUE)) %>% 
  filter(is_dup == FALSE) %>% 
  select(-is_dup) %>% 
  distinct

#### CCS to Procedure cateogry Crosswalk ###

ophtho_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                        sheet = "Ophtho") %>% 
  select(`Yale SAS Definition`, ccs)

ortho_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                       sheet = "Ortho") %>% 
  select(`Yale SAS Definition`, ccs)

vasc_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                      sheet = "Vascular") %>% 
  select(`Yale SAS Definition`, ccs)

gen_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                     sheet = "Gen") %>% 
  select(`Yale SAS Definition`, ccs)

neuro_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                       sheet = "Neuro") %>% 
  select(`Yale SAS Definition`, ccs)

plast_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                       sheet = "Plastic") %>% 
  select(`Yale SAS Definition`, ccs)

ent_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                     sheet = "Ent") %>% 
  select(`Yale SAS Definition`, ccs)

obgyn_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                       sheet = "Obgyn") %>% 
  select(`Yale SAS Definition`, ccs)

ct_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                    sheet = "CT_CCS") %>% 
  select(`Yale SAS Definition`, ccs)

uro_df <- read_excel("/project/Burke_SNF_VBP/FOIA/2021 Surgical_Cat_Specifications.xlsx", 
                     sheet = "Uro_CCS") %>% 
  select(`Yale SAS Definition`, ccs)

ccsProc_df <- rbind(ophtho_df, ortho_df) %>% 
  rbind(vasc_df) %>% 
  rbind(gen_df) %>% 
  rbind(neuro_df) %>% 
  rbind(plast_df) %>% 
  rbind(ent_df) %>% 
  rbind(obgyn_df) %>% 
  rbind(ct_df) %>% 
  rbind(uro_df) %>% 
  mutate(CCS_num = as.numeric(ccs )) %>% 
  rename(Proc_cat_ = `Yale SAS Definition`) %>% 
  select(!ccs)



#### Step 2B: Create ICD to CCS Crosswalk ####
procCW_df <- procLng_df  %>% 
  mutate(ICD_ct = substr(ICD, 1, 4)) %>% 
  left_join(Proc_ICD9_CCS) %>% #ICD9
  left_join(Proc_ICD10_CCS) %>% #ICD10 
  left_join(Proc_ICD9ct_CCS) %>% #ICD9 partial
  left_join(Proc_ICD10ct_CCS) %>% #ICD10 partial
  mutate(miss_ind = ifelse(is.na(CCS9_num) == TRUE & 
                             is.na(CCS10_num) == TRUE, 1, 0),
         doub_match = ifelse(is.na(CCS9_num) == FALSE & 
                               is.na(CCS10_num) == FALSE, 1, 0), 
         CCS_num = case_when( ICD_typ == 9 & is.na(CCS9_num) == FALSE ~ CCS9_num,
                              ICD_typ == 9 & is.na(CCS9_num) == TRUE ~ CCS9_num2,
                              ICD_typ == 0 & is.na(CCS10_num) == FALSE ~ CCS10_num,
                              ICD_typ == 0 & is.na(CCS10_num) == TRUE ~ CCS10_num2,
                              ICD_typ == 9 & is.na(CCS9_num) == TRUE &
                                is.na(CCS10_num) == FALSE ~ CCS10_num,
                              ICD_typ == 10 & is.na(CCS10_num) == TRUE &
                                is.na(CCS9_num) == FALSE ~ CCS9_num,
                              is.na(CCS9_num) == TRUE & 
                                is.na(CCS10_num) == TRUE &
                                is.na(CCS9_num2) == TRUE & 
                                is.na(CCS10_num2) == TRUE ~ NA, 
                              TRUE ~ NA)) %>% # 5% unmatched
  left_join(ccsProc_df) %>% 
  select(ICD, ICD_typ, Proc_cat_) %>% 
  filter(is.na(Proc_cat_) == FALSE) %>% 
  distinct() %>% 
  rename(PROX_PRCDR_CD = ICD,
         PROX_SRGCL_PRCDR_VRSN_CD_ = ICD_typ)

proc_comb <- coh_df1 %>% 
  select(BENE_ID, Prox_PRVDR_NUM, Prox_admsndt, 
         contains("PROX_PRCDR_CD"), contains("PROX_SRGCL_PRCDR_VRSN_CD_"))


for (x in 1:25) {
  
  # create 0 leading number   
  dig2 <- ifelse(nchar(as.character(x)) == 1, sprintf("%02d", x), x)
  
  # distinguish join by iteration
  ICD_join <- procCW_df 
  colnames(ICD_join) <- paste(colnames(ICD_join), dig2, sep = "")
  
  # join and create variable
  proc_comb <- proc_comb  %>% 
    left_join(ICD_join) 
  
}



#### Step 2C: Create Procedure Category Indicator ####


# list of surgical variables
proc_cols <- grep("^Proc_cat", colnames(proc_comb), value = TRUE)

# create surgery category vars
proc_comb2 <- proc_comb %>% 
  mutate(ophtho_ind = 
           as.integer(if_any(any_of(proc_cols), ~.x == 
                               "optho")),
         vasc_ind = 
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "vascular=")),
         ortho_ind = 
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "ortho")),
         gen_ind =
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "gen")),
         ent_ind = 
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "ent")),
         neuro_ind = 
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "neuro")),
         plast_ind = 
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "plastic")),
         obgyn_ind =
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "obgyn")), 
         ct_ind =
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "ct")),
         uro_ind = 
           as.integer(if_any(any_of(proc_cols), ~.x ==
                               "uro"))
  ) 

# sum procedure columns (to find empty columns)
ind_columns_sum <- proc_comb2 %>%
  dplyr::select(ends_with("_ind")) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

# Get the column names to remove (those with sum equal to 0)
columns_to_remove <- ind_columns_sum %>%
  dplyr::select(where(~ all(. == 0))) %>%
  names()


# Select columns that end with "_ind"
ind_cols <- grep("_ind$", colnames(proc_comb2), value = TRUE)

# create dataframe of procedure indicator columns
proc_comb_jn <- proc_comb2 %>%
  mutate_at(vars(all_of(ind_cols)), ~replace(., is.na(.), 0)) %>% # Replace NAs with 0
  dplyr::select(-one_of(columns_to_remove)) %>%  #remove 0 sum columns
  dplyr::select(BENE_ID, Prox_PRVDR_NUM, Prox_admsndt, contains("_ind"))

# join procedure indicator columns with cohort data

coh_df2 <- left_join(coh_df1, proc_comb_jn) %>% 
  dplyr::select(-contains("Proc_cat"))

#################################################################################################
### Step 3: Create CCS Primary Diagnosis Variable ####
#################################################################################################

#### Step 3A: Crosswalk Data Upload ####

##### ICD9 Raw Data Upload ####

CCS_ICD9_Map <- read_csv("/project/Burke_SNF_VBP/Analysis_2022/prog/Norman/Crosswalk Data/Raw/CCSR_ICD9_Mapping.csv", skip = 1) %>%
  select(c(1:4)) # Sourced from AHRQ

# remove apostrophes
CCS_ICD9_Map[] <-  lapply(CCS_ICD9_Map, gsub, pattern = "'", replacement = '')

names(CCS_ICD9_Map) = gsub(pattern = "'", replacement = "", x = names(CCS_ICD9_Map))

CCS_ICD9_join <- CCS_ICD9_Map %>% 
  mutate(ICD = toupper(gsub("[^a-zA-Z0-9]", "", `ICD-9-CM CODE`))) %>% 
  rename(CCS9 = `CCS CATEGORY`) %>% 
  select(ICD, CCS9)

##### ICD10 Raw Data Upload ####

CCS_ICD10_Map <- read_excel("/project/Burke_SNF_VBP/FOIA/CCS diag codes Yale2020 PRA_2-15-2022.xlsx")


CCS_ICD10_join <- CCS_ICD10_Map %>% 
  mutate(ICD = toupper(gsub("[^a-zA-Z0-9]", "", `ICD-10-ICD10CMCODE Code`))) %>%
  mutate(CCS10 = as.character(CCSCATEGORY)) %>% 
  select(ICD, CCS10)

CCS_ICD10_join2 <- CCS_ICD10_join %>% 
  mutate(ICD_ct = toupper(str_sub(ICD, 1, 4)), 
         CCS10B = CCS10) %>% 
  select(ICD_ct, CCS10B) %>% 
  distinct() %>% 
  mutate(is_dup = duplicated(ICD_ct) | duplicated(ICD_ct, fromLast = TRUE)) %>% 
  filter(is_dup == FALSE) %>% 
  select(-is_dup)



#### Step 3B: Crosswalk Data Cleaning ####


ICDxCCS_df <- coh_df2 %>% 
  mutate(ICD = toupper(gsub("[^a-zA-Z0-9]", "", as.character(Prox_DGNS_CD01))),
         ICD_typ = Prox_DGNS_VRSN_CD_01) %>% 
  select(ICD,  Prox_DGNS_VRSN_CD_01) %>% 
  distinct() %>% 
  left_join(CCS_ICD9_join) %>% 
  left_join(CCS_ICD10_join) %>% 
  mutate(ccs1 = case_when(Prox_DGNS_VRSN_CD_01 == 9 & is.na(CCS9) == FALSE ~ CCS9,
                          Prox_DGNS_VRSN_CD_01 == 9 & 
                            is.na(CCS9) == TRUE &
                            is.na(CCS10) == FALSE ~ CCS10,
                          Prox_DGNS_VRSN_CD_01 == 0 & is.na(CCS10) == FALSE ~ CCS10,
                          Prox_DGNS_VRSN_CD_01 == 0 & 
                            is.na(CCS9) == FALSE &
                            is.na(CCS10) == TRUE ~ CCS9,
                          is.na(CCS9) == TRUE &
                            is.na(CCS10) == TRUE ~ "0" ),
         doub_match = ifelse(is.na(CCS9) == FALSE &
                               is.na(CCS10) == FALSE, 1, 0), 
         ICD_ct = str_sub(ICD, 1, 4)) %>% 
  left_join(CCS_ICD10_join2) %>% 
  mutate(CCS_primaryDGNS = case_when(ccs1 ==  "0" & is.na(CCS10B) == FALSE ~ CCS10B, 
                                     ccs1 ==  "0" & is.na(CCS10B) == TRUE ~ "0",
                                     TRUE ~ ccs1)) %>% 
  select(ICD, Prox_DGNS_VRSN_CD_01, CCS_primaryDGNS)



#### Step 3C: Create CCS Variables ####

coh_df3 <- coh_df2 %>%
  left_join(ICDxCCS_df, by = c("Prox_DGNS_CD01" = "ICD", 
                               "Prox_DGNS_VRSN_CD_01" = "Prox_DGNS_VRSN_CD_01")) %>%
  mutate(HRRP_diag = case_when(CCS_primaryDGNS == "100" ~ "AMI", 
                               CCS_primaryDGNS == "127" ~ "COPD",
                               CCS_primaryDGNS == "10" ~ "HF",
                               CCS_primaryDGNS == "122" ~ "Pneu",
                               TRUE ~ "non_HRRP"),
         Covid_mainDGNS = ifelse(Prox_DGNS_CD01 == "U071" | 
                                   Prox_DGNS_CD01 == "U071", 1, 0))
#################################################################################################
### Step 4 Save and Export ####
#################################################################################################

write.csv(coh_df3, "/project/Burke_SNF_VBP/Analysis_2022/prog/Heintz/Analytical Datasets/Final/BeneFinal_controls.csv")

