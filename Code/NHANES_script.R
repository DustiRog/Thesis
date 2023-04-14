library(here)
library(fastGFPCA)
library(mgcv)
library(ggplot2)
library(splines)
library(tidyverse)
library(refund)
library(lme4)
library(parallel)
library(R.utils)

source(here::here("NHANES_funks.R"))

nhanes_df=readRDS(here("NHANES_df.rds")) %>%
  select(SEQN, Z, sind, age_years_interview, gender) %>%
  mutate("id" = SEQN)

cores = 8

#N = 8763
#J = 1440

params = data.frame("run" = "NHANES")

############ Process le data
print(Sys.time())
fgfosr_start_t = Sys.time() #Time that sucka

# if run on windows remove cores (default is 1) or set to 1 to remove parallel running
mod_list = nhanes_get_fgfosr_mod(nhanes_df , 20, cores)

fgfosr_end_t = Sys.time() #The cake is done!
fgfosr_time_diff =  as.double(difftime(fgfosr_end_t, fgfosr_start_t, units="mins"))

saveRDS(mod_list, "NHANES_mod_list.rds")
writeLines(as.character(fgfosr_time_diff), "runtime.txt")     


