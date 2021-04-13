library(ggplot2)
library(dplyr)

## Read in data of number of infections aggregated by week
weekly_cases <- read.csv("weekly_cases.csv")

## Specify the health zones we are plotting, change date to a week
weekly_cases2  <- weekly_cases %>%
  group_by(terr) %>%
  mutate(cum_confirm = cumsum(case)) %>%
  filter(terr %in% c("beni", "butembo", 'kalunguta', 'katwa', 'mabalako', 'mandima')) %>%
  mutate(Date = as.Date(week)) %>%
  mutate(SCOnset = terr)

## Change from dates to number of weeks from 1 to n=total number of weeks
week_number <- weekly_cases2 %>%
  arrange(week)
weeks <- as.data.frame(unique(week_number$week))
weeks$weeknum <- seq(1:61)
colnames(weeks) <- c("Date", "Week")
weeks$Date <- as.Date(weeks$Date)

weekly_cases2 <- left_join(weekly_cases2,weeks,by="Date")

## Plot a stacked bar chart of aggregated number of weeks
ggplot(weekly_cases2, aes(fill=SCOnset, x=Week, y=case)) + labs(fill = "Health Zone") +
  geom_bar(position = "stack", stat = "identity") + 
  scale_x_continuous(breaks = seq(1,61,by=3)) +
  scale_y_continuous(breaks = seq(0,180, by=10)) +
  theme(legend.position = c(0.15,0.85), legend.title = element_text(size=26), 
        legend.text = element_text(size=20)) + 
  ylab("Cases by Week") + xlab("Week") +
  theme(axis.title.y = element_text(size = rel(1.4), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.4)))

