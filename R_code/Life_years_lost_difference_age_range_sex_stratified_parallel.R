#Libraries

library(data.table)
library(tidyverse)
library(lillies)
library(furrr)

#Import data with cohorts

load(file = "/path/Comorbidity_SMI/Data/data_main_analysis.RData")

#Prepare the data for the procedure

data_lyl_males_no_smi <- map(.x = data_main_analysis,
                             ~ .x %>%
                               filter(POHL == 1 & exposure == "unexposed") %>%
                               select(mortality,
                                      age_death_or_censoring,
                                      age_first_hosp_comorbid_condition = VEK) %>%
                               as.data.frame())

data_lyl_males_smi <- map(.x = data_main_analysis,
                          ~ .x %>%
                            filter(POHL == 1 & exposure == "exposed") %>%
                            select(mortality,
                                   age_death_or_censoring,
                                   age_first_hosp_comorbid_condition = VEK) %>%
                            as.data.frame())

data_lyl_females_no_smi <- map(.x = data_main_analysis,
                               ~ .x %>%
                                 filter(POHL == 2 & exposure == "unexposed") %>%
                                 select(mortality,
                                        age_death_or_censoring,
                                        age_first_hosp_comorbid_condition = VEK) %>%
                                 as.data.frame())

data_lyl_females_smi <- map(.x = data_main_analysis,
                            ~ .x %>%
                              filter(POHL == 2 & exposure == "exposed") %>%
                              select(mortality,
                                     age_death_or_censoring,
                                     age_first_hosp_comorbid_condition = VEK) %>%
                              as.data.frame())

#Set up parallelization

future::plan(multisession, workers = 24)

#Estimate LYL
#Males
#Indices of cohorts with at least 10 at-risk individuals

ind_males_smi <- which(map_lgl(.x = data_lyl_males_smi,
                               ~ .x %>%
                                 summarise(n_at_risk = sum(age_first_hosp_comorbid_condition <= 80) >= 10) %>%
                                 pull(n_at_risk)))

lyl_smi_males <- future_map(.x = data_lyl_males_smi[ind_males_smi],
                            ~ lyl_ci(lyl_range(data = .x,
                                               t0 = age_first_hosp_comorbid_condition,
                                               t = age_death_or_censoring,
                                               status = mortality,
                                               age_begin = 0,
                                               age_end = pmin(80, max(.x$age_first_hosp_comorbid_condition)),
                                               tau = 81),
                                     niter = 10000),
                            .options = furrr_options(seed = 123))


lyl_no_smi_males <- future_map2(.x = data_lyl_males_no_smi[ind_males_smi],
                                .y = data_lyl_males_smi[ind_males_smi],
                                ~ lyl_ci(lyl_range(data = .x,
                                                   t0 = age_first_hosp_comorbid_condition,
                                                   t = age_death_or_censoring,
                                                   status = mortality,
                                                   age_begin = 0,
                                                   age_end = pmin(80, max(.y$age_first_hosp_comorbid_condition)),
                                                   tau = 81),
                                         niter = 10000),
                                .options = furrr_options(seed = 123))

#Females
#Indices of cohorts with at least 10 at-risk individuals

ind_females_smi <- which(map_lgl(.x = data_lyl_females_smi,
                                 ~ .x %>%
                                   summarise(n_at_risk = sum(age_first_hosp_comorbid_condition <= 80) >= 10) %>%
                                   pull(n_at_risk)))

lyl_smi_females <- future_map(.x = data_lyl_females_smi[ind_females_smi],
                              ~ lyl_ci(lyl_range(data = .x,
                                                 t0 = age_first_hosp_comorbid_condition,
                                                 t = age_death_or_censoring,
                                                 status = mortality,
                                                 age_begin = 0,
                                                 age_end = pmin(80, max(.x$age_first_hosp_comorbid_condition)),
                                                 tau = 81),
                                       niter = 10000),
                              .options = furrr_options(seed = 123))


lyl_no_smi_females <- future_map2(.x = data_lyl_females_no_smi[ind_females_smi],
                                  .y = data_lyl_females_smi[ind_females_smi],
                                  ~ lyl_ci(lyl_range(data = .x,
                                                     t0 = age_first_hosp_comorbid_condition,
                                                     t = age_death_or_censoring,
                                                     status = mortality,
                                                     age_begin = 0,
                                                     age_end = pmin(80, max(.y$age_first_hosp_comorbid_condition)),
                                                     tau = 81),
                                           niter = 10000),
                                  .options = furrr_options(seed = 123))


#Differences in life-years lost
#Males

lyl_diff_res_males <- imap_dfr(pmap(list(lyl_smi_males,
                                         lyl_no_smi_males,
                                         map(.x = data_lyl_males_smi[ind_males_smi],
                                             ~ .x %>%
                                               select(age_first_hosp_comorbid_condition) %>%
                                               pull())),
                                    function(x, y, z) {
                                      lyl_diff(x, 
                                               y,
                                               weights = z) %>%
                                        as.data.frame()}),
                               ~ .x %>%
                                 mutate(cohort = .y))

#Females

lyl_diff_res_females <- imap_dfr(pmap(list(lyl_smi_females,
                                           lyl_no_smi_females,
                                           map(.x = data_lyl_females_smi[ind_females_smi],
                                               ~ .x %>%
                                                 select(age_first_hosp_comorbid_condition) %>%
                                                 pull())),
                                      function(x, y, z) {
                                        lyl_diff(x, 
                                                 y,
                                                 weights = z) %>%
                                          as.data.frame()}),
                                 ~ .x %>%
                                   mutate(cohort = .y))

#Table

lyl_diff_res_males %>%
  mutate(sex = "males") %>%
  bind_rows(lyl_diff_res_females %>%
              mutate(sex = "females")) %>%
  transmute(cohort,
            sex,
            lyl_diff_estimate_with_ci = paste(formatC(round(lyl_estimate.TotalLYL, 2), format = "f", digits = 2),
                                              paste0("(", 
                                                     formatC(round(lyl_ci_left.TotalLYL, 2), format = "f", digits = 2),
                                                     "; ",
                                                     formatC(round(lyl_ci_right.TotalLYL, 2), format = "f", digits = 2),
                                                     ")"))) %>%
  pivot_wider(names_from = "sex",
              values_from = "lyl_diff_estimate_with_ci") %>%
  mutate(order = case_when(cohort == "Diseases of the circulatory system" ~ 1,
                           cohort == "Hypertension" ~ 2,
                           cohort == "Ischemic heart disease" ~ 3,
                           cohort == "Atrial fibrillation" ~ 4,
                           cohort == "Heart failure" ~ 5,
                           cohort == "Peripheral artery occlusive disease" ~ 6,
                           cohort == "Stroke" ~ 7,
                           cohort == "Diseases of the endocrine system" ~ 8,
                           cohort == "Diabetes mellitus" ~ 9,
                           cohort == "Thyroid disorder" ~ 10,
                           cohort == "Chronic pulmonary diseases" ~ 11,
                           cohort == "Diseases of the gastrointestinal system" ~ 12,
                           cohort == "Ulcer or chronic gastritis" ~ 13,
                           cohort == "Chronic liver disease" ~ 14,
                           cohort == "Inflammatory bowel disease" ~ 15,
                           cohort == "Diverticular disease of intestine" ~ 16,
                           cohort == "Diseases of the urogenital system" ~ 17,
                           cohort == "Chronic kidney disease" ~ 18,
                           cohort == "Prostate disorders" ~ 19,
                           cohort == "Connective tissue disorders" ~ 20,
                           cohort == "Cancers" ~ 21,
                           cohort == "Diseases of the neurological system" ~ 22,
                           cohort == "Epilepsy" ~ 23,
                           cohort == "Parkinson's disease" ~ 24,
                           cohort == "Multiple sclerosis" ~ 25,
                           cohort == "Infectious and parasitic diseases" ~ 26,
                           cohort == "Tuberculosis" ~ 27,
                           cohort == "Chronic viral hepatitis" ~ 28)) %>%
  arrange(order) %>%
  select(-order) %>%
  write.csv(file = "/path/Comorbidity_SMI/Results/Diff_life_years_lost_sex_stratified.csv",
            row.names = FALSE)
