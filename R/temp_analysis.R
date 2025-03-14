#R script for Orrock et al. 2025 Canadian Journal of Forest Research
#Freeze-thaw events differently affect survival of seeds of two native
#and two invasive woody species

##Code to plot temperature during the freeze-thaw treatment and within the growth chamber

#load packages
library(tidyverse)
library(here)
library(lubridate)

#source files/functions
'%!in%' <- function(x,y)!('%in%'(x,y))

## Custom ggplot theme
clrs <- MetBrewer::met.brewer("Tiepolo")

theme_nice <- function() {
  theme_minimal(base_family = "Calibri") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(family = "Calibri", face = "bold"),
          axis.title = element_text(family = "Calibri"),
          strip.text = element_text(family = "Calibri", face = "bold",
                                    size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          axis.line = element_line(size=.7, color="black"))
}

#load data
temps <- read.csv(here::here("data", "ft_temp_data.csv"))
logger.id <- read.csv(here::here("data", "logger_id_freeze_thaw.csv"))

#drop extra loggers
temps <- temps[, !names(temps) %in% c("D20000003F5F6521.1", "D20000003F5F6521.2")]

#convert column to date/time
temps$Time <- ymd_hms(temps$Time,tz=Sys.timezone())

ft.temp <- temps %>% 
  filter(between(Time, ymd_hms("2023-01-02 18:00:00"), ymd_hms("2023-01-09 08:00:00"))) %>% 
  pivot_longer(!Time, names_to = "DeviceAddress", values_to = "temp") %>% 
  left_join(., logger.id) %>% 
  group_by(Time, FreezeThaw) %>% 
  summarize(mean = mean(temp), se=(sd(temp)/sqrt(length(temp))))

png("ft_treatment_temp.png", height=1500, width=1800, res=300)
ggplot(ft.temp, aes(Time, mean, color=FreezeThaw))+
  geom_line(size=1)+
  geom_ribbon(aes(ymax=mean+se, ymin=mean-se, fill=FreezeThaw), alpha=0.25)+
  scale_color_manual(name="Freeze-thaw", values = c(clrs[2], clrs[7]))+
  scale_fill_manual(name="Freeze-thaw", values = c(clrs[2], clrs[7]))+
  labs(y="Mean temperature (C)")+
  theme_nice()+
  theme(text = element_text(size=16), legend.position = c(0.8,0.8))
dev.off()

##  
gc.temp<-temps %>% 
  filter(Time>ymd_hms("2023-01-09 18:00:00")) %>% 
  pivot_longer(!Time, names_to = "DeviceAddress", values_to = "temp") %>% 
  left_join(., logger.id) %>% 
  group_by(Time) %>% 
  summarize(mean=mean(temp), se=(sd(temp)/sqrt(length(temp))))

png("growth_chamber_temp.png", height=1500, width=1800, res=300)
ggplot(gc.temp, aes(Time, mean))+
  geom_line(size=1)+
  geom_ribbon(aes(ymax=mean+se, ymin=mean-se), fill="grey", alpha=0.5)+
  #scale_color_manual(name="Freeze-thaw", values = c(clrs[2], clrs[7]))+
  #scale_fill_manual(name="Freeze-thaw", values = c(clrs[2], clrs[7]))+
  labs(y="Mean temperature (C)")+
  theme_nice()+
  theme(text = element_text(size=16), legend.position = "none")
dev.off()
