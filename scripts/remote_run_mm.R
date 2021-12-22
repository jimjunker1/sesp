#!/usr/bin/env Rscript --vanilla
## Code to clean and run SESP isotope and fatty acid data
setwd("/home/jrjunker/projects/sesp")
# load library calls
library(magrittr)
library(dplyr)
library(fastinR)
# load and clean SESP tracer & resource tracer data
mixing_files = list.files("./data/", "*.rds", full.names = TRUE)
mixing_objects = lapply(mixing_files, readRDS) 
# setNames and load into global
mixing_objects = setNames(mixing_objects, nm = gsub(".rds","", sapply(strsplit(mixing_files,"/"),"[",11))) %>%
  list2env(.,.GlobalEnv)

### 2013 ----------
## Let's run fastinR on just 2013 SESP liver
sesp_liver_2013 = sesp_liver %>%
  dplyr::filter(year == "2013") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_liver_2013 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_liver_2013 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.1,
                               Rnot = 0.2,
                               Rnot_SI = 0.2,
                               plott = FALSE,
                               spawn = FALSE)

saveRDS(fastin_run, file = "./derived-data/sesp_2013_liver_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2013_liver_fastinR.Rdata")

## Let's run fastinR on just 2013 SESP muscle
sesp_muscle_2013 = sesp_muscle %>%
  dplyr::filter(year == "2013") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_muscle_2013 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_muscle_2013 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.1,
                               Rnot = 0.3,
                               Rnot_SI = 0.3,
                               plott = FALSE,
                               spawn = FALSE)

saveRDS(fastin_run, file = "./derived-data/sesp_2013_muscle_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2013_muscle_fastinR.Rdata")

### 2014 ----------
## Let's run fastinR on just 2014 SESP liver
sesp_liver_2014 = sesp_liver %>%
  dplyr::filter(year == "2014") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_liver_2014 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_liver_2014 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.1,
                               Rnot = 0.2,
                               Rnot_SI = 0.2,
                               plott = FALSE,
                               spawn = FALSE)
saveRDS(fastin_run, file = "./derived-data/sesp_2014_liver_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2014_liver_fastinR.Rdata")

## Let's run fastinR on just 2014 SESP muscle
sesp_muscle_2014 = sesp_muscle %>%
  dplyr::filter(year == "2014") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_muscle_2014 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_muscle_2014 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.1,
                               Rnot = 0.3,
                               Rnot_SI = 0.3,
                               plott = FALSE,
                               spawn = FALSE)
saveRDS(fastin_run, file = "./derived-data/sesp_2014_muscle_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2014_muscle_fastinR.Rdata")

### 2015 ----------
## Let's run fastinR on just 2015 SESP liver
sesp_liver_2015 = sesp_liver %>%
  dplyr::filter(year == "2015") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_liver_2015 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_liver_2015 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
dataplot(dats)
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.1,
                               Rnot = 0.2,
                               Rnot_SI = 0.2,
                               plott = FALSE,
                               spawn = FALSE#,
                               # file = here::here("sub-project/SESP/data/models/sesp_2015_liver.Rdata"),
                               # n.adapt = 5000
                               )
saveRDS(fastin_run, file = "./derived-data/sesp_2015_liver_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2015_liver_fastinR.Rdata")


## Let's run fastinR on just 2015 SESP muscle
sesp_muscle_2015 = sesp_muscle %>%
  dplyr::filter(year == "2015") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_muscle_2015 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_muscle_2015 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame
  
FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.5,
                               Rnot = 0.8,
                               Rnot_SI = 0.9,
                               plott = FALSE,
                               spawn = FALSE)

saveRDS(fastin_run, file = "./derived-data/sesp_2015_muscle_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2015_muscle_fastinR.Rdata")

### 2016 ----------
## Let's run fastinR on just 2016 SESP liver
sesp_liver_2016 = sesp_liver %>%
  dplyr::filter(year == "2016") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_liver_2016 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_liver_2016 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.3,
                               Rnot = 0.5,
                               Rnot_SI = 0.2,
                               plott = FALSE,
                               spawn = FALSE)
saveRDS(fastin_run, file = "./derived-data/sesp_2016_liver_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2016_liver_fastinR.Rdata")

## Let's run fastinR on just 2016 SESP muscle
sesp_muscle_2016 = sesp_muscle %>%
  dplyr::filter(year == "2016") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_muscle_2016 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_muscle_2016 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.3,
                               Rnot = 0.3,
                               Rnot_SI = 0.3,
                               plott = FALSE,
                               spawn = FALSE)
saveRDS(fastin_run, file = "./derived-data/sesp_2016_muscle_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2016_muscle_fastinR.Rdata")


#### 2017 -----
## Let's run fastinR on just 2017 SESP liver
sesp_liver_2017 = sesp_liver %>%
  dplyr::filter(year == "2017") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_liver_2017 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_liver_2017 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.5,
                               Rnot = 0.8,
                               Rnot_SI = 0.9,
                               plott = FALSE,
                               spawn = FALSE)
saveRDS(fastin_run, file = "./derived-data/sesp_2017_liver_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2017_liver_fastinR.Rdata")

## Let's run fastinR on just 2017 SESP muscle
sesp_muscle_2017 = sesp_muscle %>%
  dplyr::filter(year == "2017") %>%
  select(-tissue, `%c`,`%n`,cn)

mixing_sesp_si = sesp_muscle_2017 %>%
  dplyr::select(samp_id, d13c, d15n) %>%
  column_to_rownames('samp_id')

# set up the stable isotope data
# debugonce(add_SI)
dats = fastinR::add_SI(SI.predators = mixing_sesp_si,
                       SI.preys = mixing_basal_si,
                       Frac.Coeffs.mean = basal_frac_mean,
                       Frac.Coeffs.var = basal_frac_var)

# add in the FA dat

mixing_sesp_fa = sesp_muscle_2017 %>%
  dplyr::select(samp_id, `14_0`:`22_6n3`) %>%
  column_to_rownames('samp_id') %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::mutate(across(matches("\\d{2}_.*"), ~ case_when(
    is.na(.x) ~ 0.001,
    .x <= 0 ~ 0.001,
    grepl("measured", .x, ignore.case = TRUE) ~ NA_real_,
    TRUE ~ .x
  ))) %>%
  data.frame

FA_list = mixing_sesp_fa %>% 
  names %>% unlist
mixing_basal_fa = mixing_basal_fa %>%
  dplyr::select(any_of(c('node_code', FA_list)))

# debugonce(add_FA)
dats = fastinR::add_FA(
  FA.predators = mixing_sesp_fa,
  FA.preys = mixing_basal_fa,
  CC.mean = 1, 
  CC.var = 0.1,
  datas = dats
)

## check the FA list
dats_subset = fastinR::select_vars(dats)
1 2 3 4 5 6 7 8 9
## run the model

set.seed(42)
# debugonce(fastinR::run_MCMC)
fastin_run = fastinR::run_MCMC(datas = dats_subset,
                               Covs = NULL,
                               nIter = 100000,
                               nBurnin = 90000,
                               nChains = 4,
                               nThin = 100,
                               Data.Type = 'Combined.Analysis',
                               Analysis.Type = 'Individual.proportions',
                               even = 0.7,
                               Rnot = 0.8,
                               Rnot_SI = 0.9,
                               plott = FALSE,
                               spawn = FALSE)
saveRDS(fastin_run, file = "./derived-data/sesp_2017_muscle_fastinR.rds")
file.rename("jagsdata.Rdata", "./data/models/sesp_2017_muscle_fastinR.Rdata")
