#Libraries

library(data.table)
library(tidyverse)
library(lubridate)
library(broom)
library(survival)

#Variable definition

#RODCIS2 = unique personal identifier
#DATPRI = date of admission
#DATUKO = date of discharge
#VEK = age
#POHL = sex
#NBYDL = region of residence
#ZDG = primary diagnosis
#DAUMR = date of death

#Import all hospitalizations from 1994 to 2017

hospitalizations_1994_2015 <-  list.files(path = "/path/UZIS_data_raw", 
                                          full.names = TRUE) %>%
  .[str_detect(., "1[6-7].csv$", negate = TRUE)] %>%
  map_dfr(~ fread(.,
                  select = c("RODCIS2" = "character", 
                             "DATPRI" = "integer", 
                             "DATUKO" = "integer", 
                             "VEK"= "integer", 
                             "POHL" = "integer",
                             "NBYDL"= "character",
                             "ZDG" = "character"),
                  header = TRUE,
                  sep = ",",
                  dec = ".",
                  fill = TRUE,
                  encoding = "Latin-1",
                  nThread = 8))

hospitalizations_2016_2017 <- list.files(path = "/path/UZIS_data_raw", 
                                         full.names = TRUE) %>%
  .[str_detect(., "1[6-7].csv$", negate = FALSE)] %>%
  map_dfr(~ fread(.,
                  select = c("RODCIS" = "character", 
                             "DATPRI" = "integer", 
                             "DATUKO" = "integer", 
                             "VEK"= "integer", 
                             "POHL" = "integer",
                             "NBYDL"= "character",
                             "ZDG" = "character"),
                  header = TRUE,
                  sep = ";",
                  dec = ".",
                  fill = TRUE,
                  encoding = "Latin-1",
                  nThread = 8)) %>%
  rename(RODCIS2 = RODCIS)

#Combine all data 

hospitalizations_1994_2017 <- hospitalizations_1994_2015 %>%
  bind_rows(hospitalizations_2016_2017)

#Remove partial data

rm(hospitalizations_1994_2015)
rm(hospitalizations_2016_2017)

#Import data on mortality

deaths_1994_2013 <- fread(file = "/path/UZIS_mortality_raw/zem_1994_2013.csv")
deaths_2014 <- fread(file = "/path/UZIS_mortality_raw/zem_2014.csv")
deaths_2015 <- fread(file = "/path/UZIS_mortality_raw/zem_2015.csv")
deaths_2016 <- fread(file = "/path/UZIS_mortality_raw/zem_2016.csv")
deaths_2017 <- fread(file = "/path/UZIS_mortality_raw/zem_2017.csv")

#Unifying the format of mortality data

deaths_1994_2017 <- deaths_1994_2013 %>%
  transmute(RODCIS2 = RC,
            DAUMR = dmy(DAUMR),
            cause_of_death = trimws(DGP),
            external_cause_of_death = trimws(DGE)) %>%
  bind_rows(deaths_2014 %>%
              transmute(RODCIS2 = RC,
                        DAUMR = dmy(DAUMR),
                        cause_of_death = trimws(DGP),
                        external_cause_of_death = trimws(DGE)),
            deaths_2015 %>%
              transmute(RODCIS2 = RODCIS2,
                        DAUMR = ymd(DAUMR),
                        cause_of_death = trimws(DGP),
                        external_cause_of_death = trimws(DGE)),
            deaths_2016 %>%
              transmute(RODCIS2 = RCZEMAN2,
                        DAUMR = ymd(paste0(UMROK, str_pad(UMRMM, 2, pad = "0"), UMRDD)),
                        cause_of_death = trimws(DGUMR),
                        external_cause_of_death = trimws(DGUMR2)),
            deaths_2017 %>%
              transmute(RODCIS2 = RCZEMAN2,
                        DAUMR = ymd(paste0(UMROK, str_pad(UMRMM, 2, pad = "0"), UMRDD)),
                        cause_of_death = trimws(DGUMR),
                        external_cause_of_death = trimws(DGUMR2)))

#Remove partial data

rm(deaths_1994_2013)
rm(deaths_2014)
rm(deaths_2015)
rm(deaths_2016)
rm(deaths_2017)

#Excluding records with missing values on key variables
#Excluding records with invalid dates

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
  filter(rowSums(is.na(across(c(RODCIS2, DATPRI, DATUKO, VEK, POHL, NBYDL, ZDG)))) == 0) %>%
  filter(!is.na(ymd(DATPRI)) & !is.na(ymd(DATUKO)))

#Excluding individuals with more than one date of death
#Excluding individuals with hospitalizations after death

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
  anti_join(hospitalizations_1994_2017 %>%
              inner_join(deaths_1994_2017, 
                         by = c("RODCIS2" = "RODCIS2")) %>%
              mutate(DATUKO = ymd(DATUKO)) %>%
              group_by(RODCIS2) %>%
              filter(DAUMR < max(DATUKO) | n_distinct(DAUMR) > 1) %>%
              ungroup(),
            by = c("RODCIS2" = "RODCIS2"))

#Excluding records with invalid date overlaps (discharge date after the admission date of another record)

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
  anti_join(map_dfr(.x = hospitalizations_1994_2017 %>%
                      group_split(split_ID = frank(RODCIS2, ties.method = "dense") %/% 100000),
                    ~ .x %>% 
                      select(RODCIS2,
                             DATPRI_index = DATPRI,
                             DATUKO_index = DATUKO) %>%
                      inner_join(.x %>%
                                   select(RODCIS2,
                                          DATPRI_non_index = DATPRI,
                                          DATUKO_non_index = DATUKO), 
                                 by = c("RODCIS2" = "RODCIS2")) %>%
                      filter(DATPRI_non_index < DATPRI_index & DATUKO_non_index > DATPRI_index) %>%
                      pivot_longer(-RODCIS2,
                                   names_to = c(".value", "type"), 
                                   names_pattern = "([^_]+)_(.*)")),
            by = c("RODCIS2" = "RODCIS2",
                   "DATPRI" = "DATPRI",
                   "DATUKO" = "DATUKO"))

#Keeping hospitalizations that occurred between 1st January 1999 and 31st December 2017

hospitalizations_1999_2017 <- hospitalizations_1994_2017 %>%
  filter(year(ymd(DATPRI)) >= 1999 & year(ymd(DATUKO)) <= 2017) 

#Import data for the main analysis

load(file = "/path/Comorbidity_SMI/Data/data_main_analysis.RData")

#Establish mortality

data_models <- map(.x = data_main_analysis,
                   ~ .x %>%
                     mutate(mortality_natural_causes = as.numeric(int_overlaps(interval(ymd(DATUKO), ymd("2017-12-31")), interval(DAUMR, DAUMR)) * str_detect(external_cause_of_death, "^V|^X|^Y", negate = TRUE)),
                            mortality_natural_causes = replace(mortality_natural_causes, is.na(mortality_natural_causes), 0),
                            mortality_external_causes = as.numeric(int_overlaps(interval(ymd(DATUKO), ymd("2017-12-31")), interval(DAUMR, DAUMR)) * str_detect(external_cause_of_death, "^V|^X|^Y")),
                            mortality_external_causes = replace(mortality_external_causes, is.na(mortality_external_causes), 0),
                            years_until_death_natural_causes_or_censoring = case_when(mortality_natural_causes == 1 ~ years_until_death,
                                                                                       mortality_external_causes == 1 ~ years_until_death,
                                                                                       TRUE ~ as.numeric(years_diff_mortality_followup))))

#Stratified Cox proportional hazards models
#Fully adjusted

mfull <- imap_dfr(data_models,
                  ~ tidy(coxph(Surv(years_until_death_natural_causes_or_censoring, mortality_natural_causes) ~ exposure + VEK + POHL + discharge_year + strata(strata), 
                               data = .x),
                         conf.int = TRUE,
                         exponentiate = TRUE) %>%
                    mutate(cohort = .y) %>%
                    filter(term == "exposureexposed"))

#Export table with results

mfull %>%
  arrange(cohort) %>%
  transmute(cohort,
            `Natural causes of death` = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                                              paste0("(", 
                                                     formatC(round(conf.low, 2), format = "f", digits = 2),
                                                     "; ",
                                                     formatC(round(conf.high, 2), format = "f", digits = 2),
                                                     ")"))) %>%
  write.csv(file = "/path/Comorbidity_SMI/Results/HR_sens_analysis_only_natural_causes_of_death.csv",
            row.names = FALSE)
