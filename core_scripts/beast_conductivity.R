library(tidyverse)
library(lubridate)

cond <- read_csv("20210404_beast_cond_150m_bowie.csv", skip = 1)%>%
  rename(DateTime = `Date Time, GMT-05:00`,
         Cond = `Full Range, μS/cm (LGR S/N: 20059128, SEN S/N: 20059128)`)%>%
  mutate(DateTime = mdy_hms(DateTime))

cond2 <- read_csv("20211022_beast_cond_150m_David.csv", skip = 1) %>%
  rename(DateTime = `Date Time, GMT-04:00`,
         Cond = `Full Range, μS/cm (LGR S/N: 20195878, SEN S/N: 20195878)`)%>%
  mutate(DateTime = mdy_hms(DateTime))
  
  
cond2 %>%
  ggplot(aes(x = DateTime, y = Cond))+
  geom_line()+
  geom_line(data = cond,mapping = aes(x = DateTime, y = Cond))
