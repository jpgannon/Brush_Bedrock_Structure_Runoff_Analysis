library(tidyverse)
library(lubridate)

source("Ls_to_mmd.R")

#Beast
#Q = 1.510428e-05 * Stage^4.2727  R^2= 0.94
#MDH
#Q = 1.737801e-21 * Stage^14.4394   R^2 = 0.98

mcd <- read_csv("McDonald_Hollow_MonStat_Data_Sept21.csv") %>%
  mutate(Q = 1.737801e-21 * (Stage_m_pt * 100)^14.4394) %>%
  mutate(Qmmd = Ls_to_mmd(Q, 62.84)) %>%
  select(datetime, Q, Qmmd)

beast <- read_csv("Qdat_beast20_share.csv") %>%
  mutate(Qmmd = Ls_to_mmd(Q, 82.85)) %>%
  select(datetime = dt.clean, Q, Qmmd) %>%
  mutate(datetime = mdy_hm(datetime))

mcd_beast <- inner_join(mcd, beast, by = "datetime", suffix = c("_mcd", "_beast"))

write_csv(mcd_beast, "mcd_beast_Q.csv")
