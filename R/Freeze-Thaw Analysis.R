#R script for Orrock et al. 2025 Canadian Journal of Forest Research
#Freeze-thaw events differently affect survival of seeds of two native
#and two invasive woody species

#R version 4.2.1

#loading packages
library(tidyverse)
library(here)
library(lubridate)
library(blme)
library(lme4)
library(DHARMa)
library(car)
library(emmeans)

#adding different colors
clrs <- MetBrewer::met.brewer("Tiepolo")

#custom ggplot theme
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

##Species labels for figures
sp.label<-c("ACESAC"="A. saccharum", "LONMAC"="L. maackii",
            "PINSTR"="P. strobus", "RHACAT"="R. cathartica")

#custom functions
'%!in%' <- function(x,y)!('%in%'(x,y))
std.er <- function(x) sd(x)/sqrt(length(x))

### planned contrast tables ---------------------------------------------
#Freeze thaw vs. Fungicide
ft.fung.con<-list(FT.Treatment.vs.Control.NoFungicide=c(1,-1,0,0),
                  FT.Treatment.vs.Control.YesFungicide=c(0,0,1,-1),
                  Fungicide.Treatment.vs.Control.Freeze=c(1,0,-1,0),
                  Fungicide.Treatment.vs.Control.FreezeThaw=c(0,1,0,-1))

#Freeze thaw by Species
ft.sp.con<-list(FreezeThaw.Treatment.vs.Control.Acer=c(1,-1,0,0,0,0,0,0),
                FreezeThaw.Treatment.vs.Control.Lonicera=c(0,0,1,-1,0,0,0,0),
                FreezeThaw.Treatment.vs.Control.Pinus=c(0,0,0,0,1,-1,0,0),
                FreezeThaw.Treatment.vs.Control.Rhamnus=c(0,0,0,0,0,0,1,-1),
                
                Acer.vs.Lonicera.Freeze=c(1,0,-1,0,0,0,0,0),
                Acer.vs.Lonicera.FreezeThaw=c(0,1,0,-1,0,0,0,0),
                Acer.vs.Pinus.Freeze=c(1,0,0,0,-1,0,0,0),
                Acer.vs.Pinus.FreezeThaw=c(0,1,0,0,0,-1,0,0),
                Acer.vs.Rhamnus.Freeze=c(1,0,0,0,0,0,-1,0),
                Acer.vs.Rhamnus.FreezeThaw=c(0,1,0,0,0,0,0,-1),
                
                Lonicera.vs.Pinus.Freeze=c(0,0,1,0,-1,0,0,0),
                Lonicera.vs.Pinus.FreezeThaw=c(0,0,0,1,0,-1,0,0),
                Lonicera.vs.Rhamnus.Freeze=c(0,0,1,0,-1,0,-1,0),
                Lonicera.vs.Rhamnus.FreezeThaw=c(0,0,0,1,0,0,0,-1),
                
                Pinus.vs.Rhamnus.Freeze=c(0,0,0,0,1,0,-1,0),
                Pinus.vs.Rhamnus.FreezeThaw=c(0,0,0,0,0,1,0,-1))

#Freeze Thaw by Fungicide by Species
ft.fung.sp.con<-list(FT.Treatment.vs.Control.NoFungicide.Acer=c(1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     FT.Treatment.vs.Control.YesFungicide.Acer=c(0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0),
                     Fungicide.Treatment.vs.Control.Freeze.Acer=c(1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                     Fungicide.Treatment.vs.Control.FreezeThaw.Acer=c(0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0),
                     
                     FT.Treatment.vs.Control.NoFungicide.Lonicera=c(0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0),
                     FT.Treatment.vs.Control.YesFungicide.Lonicera=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0),
                     Fungicide.Treatment.vs.Control.Freeze.Lonicera=c(0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0),
                     Fungicide.Treatment.vs.Control.FreezeThaw.Lonicera=c(0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0),
                     
                     FT.Treatment.vs.Control.NoFungicide.Pinus=c(0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0),
                     FT.Treatment.vs.Control.YesFungicide.Pinus=c(0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0),
                     Fungicide.Treatment.vs.Control.Freeze.Pinus=c(0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0),
                     Fungicide.Treatment.vs.Control.FreezeThaw.Pinus=c(0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0),
                     
                     FT.Treatment.vs.Control.NoFungicide.Rhamnus=c(0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0),
                     FT.Treatment.vs.Control.YesFungicide.Rhamnus=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1),
                     Fungicide.Treatment.vs.Control.Freeze.Rhamnus=c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0),
                     Fungicide.Treatment.vs.Control.FreezeThaw.Rhamnus=c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1),
                     
                     Acer.vs.Lonicera.Freeze.NoFungicide=c(1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0),
                     Acer.vs.Lonicera.FreezeThaw.NoFungicide=c(0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0),
                     Acer.vs.Lonicera.Freeze.YesFungicide=c(0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0),
                     Acer.vs.Lonicera.FreezeThaw.YesFungicide=c(0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0),
                     
                     Acer.vs.Pinus.Freeze.NoFungicide=c(1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0),
                     Acer.vs.Pinus.FreezeThaw.NoFungicide=c(0,1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0),
                     Acer.vs.Pinus.Freeze.YesFungicide=c(0,0,1,0,0,0,0,0,0,0,-1,0,0,0,0,0),
                     Acer.vs.Pinus.FreezeThaw.YesFungicide=c(0,0,0,1,0,0,0,0,0,0,0,-1,0,0,0,0),
                     
                     Acer.vs.Rhamnus.Freeze.NoFungicide=c(1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0),
                     Acer.vs.Rhamnus.FreezeThaw.NoFungicide=c(0,1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0),
                     Acer.vs.Rhamnus.Freeze.YesFungicide=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,-1,0),
                     Acer.vs.Rhamnus.FreezeThaw.YesFungicide=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,-1),
                     
                     Lonicera.vs.Pinus.Freeze.NoFungicide=c(0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0),
                     Lonicera.vs.Pinus.FreezeThaw.NoFungicide=c(0,0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0),
                     Lonicera.vs.Pinus.Freeze.YesFungicide=c(0,0,0,0,0,0,1,0,0,0,-1,0,0,0,0,0),
                     Lonicera.vs.Pinus.FreezeThaw.YesFungicide=c(0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,0),
                     
                     Lonicera.vs.Rhamnus.Freeze.NoFungicide=c(0,0,0,0,1,0,0,0,0,0,0,0,-1,0,0,0),
                     Lonicera.vs.Rhamnus.FreezeThaw.NoFungicide=c(0,0,0,0,0,1,0,0,0,0,0,0,0,-1,0,0),
                     Lonicera.vs.Rhamnus.Freeze.YesFungicide=c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1,0),
                     Lonicera.vs.Rhamnus.FreezeThaw.YesFungicide=c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1),
                     
                     Pinus.vs.Rhamnus.Freeze.NoFungicide=c(0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0),
                     Pinus.vs.Rhamnus.FreezeThaw.NoFungicide=c(0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0),
                     Pinus.vs.Rhamnus.Freeze.YesFungicide=c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0),
                     Pinus.vs.Rhamnus.FreezeThaw.YesFungicide=c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1))

### Load and clean data ----------------------------------------------------
#load data
ft.data <- read.csv(here::here("data", "FreezeThawGerminationData021323.csv"))

#adding invasive species status
ft.data$InvasionStatus <- ifelse(ft.data$Species %in% "RHACAT" | ft.data$Species %in% "LONMAC", 
                               "Invasive", "Native")

#calculate seeds not germinated
ft.data$not.germinated <- ft.data$Num.seeds - ft.data$germinated

#calculate number of dead seeds
ft.data$died <- ft.data$Num.seeds - ft.data$total

#summary of germinated seeds
ft.data %>% 
  group_by(Species) %>% 
  summarize(germ.sum = sum(germinated), seed.sum=sum(Num.seeds))

ft.data %>% 
  select(FreezeThaw, Species, Fungicide, germinated, Num.seeds) %>% 
  group_by(Species) %>% 
  summarise(mean=mean(germinated/Num.seeds),sd=sd(germinated/Num.seeds),
            se=(sd(germinated/Num.seeds)/sqrt(length(germinated/Num.seeds))),
            ci=(sd(germinated/Num.seeds)/sqrt(length(germinated/Num.seeds))) * qt((0.95/2 +0.5), (length(germinated/Num.seeds)-1)))

ft.data %>% 
  select(FreezeThaw, Species, Fungicide, germinated, viable) %>% 
  group_by(Species) %>% 
  mutate(prop.germ = germinated/(germinated + viable)) %>% 
  filter(is.finite(prop.germ), !is.na(prop.germ)) %>% 
  summarise(mean=mean(prop.germ),sd=sd(prop.germ),
            se=(sd(prop.germ)/sqrt(length(prop.germ))),
            ci=(sd(prop.germ)/sqrt(length(prop.germ))) * qt((0.95/2 +0.5), (length(prop.germ)-1)))

###Survival rate----------------------------------------------------------------
#summary of surviving seeds
ft.data %>% 
  group_by(Species) %>% 
  summarize(survive.sum = sum(total), seed.sum = sum(Num.seeds))

#Change the optimizer and number of iterations
control <- glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))

#creating a version of the Species variable with "Sum to zero" contrasts
ft.data$Species.s2z <- as.factor(ft.data$Species)
contrasts(ft.data$Species.s2z) = contr.sum(4) 

#running glmm with a weakly informative prior to account for complete separation
bayes.surv.mod.sp <- bglmer(cbind(total, died) ~ FreezeThaw*Fungicide*Species.s2z
                       +(1|Tray),
                     family=binomial, ft.data,
                     control = control,
                     fixef.prior = normal(cov=diag(6,16)))

summary(bayes.surv.mod.sp)
Anova(bayes.surv.mod.sp, type=3)

#diagnostics
plot(DHARMa::simulateResiduals(bayes.surv.mod.sp, 1000))
DHARMa::testDispersion(DHARMa::simulateResiduals(bayes.surv.mod.sp, 1000))

##contrasts
#pairwise species contrasts
pairs(emmeans(bayes.surv.mod.sp, "Species.s2z", type="response"))

#FT by Fun
emmeans(bayes.surv.mod.sp, list(~FreezeThaw + Fungicide), contr = ft.fung.con,
        adjust = "mvt", type="response")

#FT x Sp
emmeans(bayes.surv.mod.sp, list(~FreezeThaw + Species.s2z), contr = ft.sp.con,
        adjust = "mvt", type="response")

##FT x Fun x Sp
#options(scipen=999)
#options(scipen=0)
emmeans(bayes.surv.mod.sp, list(~FreezeThaw + Fungicide + Species.s2z), contr = ft.fung.sp.con,
        adjust = "mvt", type="response")

##Plot
surv.sum.sp <- ft.data %>% 
  select(FreezeThaw, Fungicide, Species, total, Num.seeds) %>% 
  group_by(FreezeThaw, Fungicide, Species) %>% 
  summarise(mean=mean(total/Num.seeds),sd=sd(total/Num.seeds),
            se=(sd(total/Num.seeds)/sqrt(length(total/Num.seeds))),
            ci=(sd(total/Num.seeds)/sqrt(length(total/Num.seeds))) * qt((0.95/2 +0.5), (length(total/Num.seeds)-1)))

png("survival-sp.png", height=2000, width=2000, res=400)
ggplot(surv.sum.sp, aes(FreezeThaw, mean, color=Fungicide))+
  geom_point(size=2.5, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size=1, width=0.2, position=position_dodge(width=0.5))+
  facet_wrap(~Species, nrow=2, labeller= as_labeller(sp.label))+
  labs(x="Thermal treatment", y="Proportion of seeds surviving")+
  scale_color_manual(values = c(clrs[2], clrs[7]))+
  theme_nice()+
  theme(text = element_text(size=16),
        strip.text= element_text(face= "bold.italic"),
        legend.position=c(0.1, 0.9))
dev.off()

### Fungal attack rate ---------------------------------------------------------
#fungus model
ft.data$NoFungus<-ft.data$Num.seeds-ft.data$fungus

##fungal mortality model
fung.mod <- ft.data %>% 
            mutate(prop.fung=fungus/Num.seeds) %>% 
            glmmTMB(cbind(died, total) ~ prop.fung + FreezeThaw + Fungicide +
              (1|Species) + (1|Tray),
            family=binomial, data=.)

#diagnostics
plot(DHARMa::simulateResiduals(mod, 1000))
DHARMa::testDispersion(DHARMa::simulateResiduals(mod, 1000))

#plotting
fun.sum<-ft.data %>% 
  select(Fungicide, Species, fungus, Num.seeds) %>% 
  group_by(Fungicide, Species ) %>% 
  summarise(mean=mean(fungus/Num.seeds),sd=sd(fungus/Num.seeds),
            se=(sd(fungus/Num.seeds)/sqrt(length(fungus/Num.seeds))),
            ci=(sd(fungus/Num.seeds)/sqrt(length(fungus/Num.seeds))) * qt((0.95/2 +0.5), (length(fungus/Num.seeds)-1)))

png("fung-trt.png", height=2000, width=3000, res=400)
ggplot(fun.sum, aes(Fungicide, mean))+
  geom_point(size=2.5, position=position_dodge(width=0.5))+
  facet_grid(cols=vars(Species), labeller= as_labeller(sp.label))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size=1, width=0.2, position=position_dodge(width=0.5))+
  labs(x="Fungicide treatment", y="Proportion of seeds with fungal growth")+
  theme_classic()+
  theme(text = element_text(size=16),
        strip.text= element_text(face= "bold.italic"),
        legend.position="none")
dev.off()

##Correlation between fungal attack and seed mortality
#All species
ft.data %>% 
  group_by(Species, FreezeThaw, Fungicide) %>% 
  summarize(prop.fung=mean(fungus/Num.seeds), prop.dead = mean(died/Num.seeds)) %>% 
  with(cor.test(prop.fung, prop.dead, use="complete.obs"))

#No Pinus
ft.data %>% 
  group_by(Species, FreezeThaw, Fungicide) %>% 
  filter(Species %!in% "PINSTR") %>% 
  summarize(prop.fung = mean(fungus/Num.seeds), prop.dead = mean(died/Num.seeds)) %>% 
  with(cor.test(prop.fung, prop.dead, use="complete.obs"))

#Plot fungus x mortality correlation
png("fung-mort.png", height=2000, width=3000, res=400)
ft.data %>% 
  group_by(Species, FreezeThaw, Fungicide) %>% 
  summarize(prop.fung=mean(fungus/Num.seeds), prop.dead=mean(died/Num.seeds)) %>% 
  mutate(Treatment= ifelse(FreezeThaw=="Control" & Fungicide=="Control", "Freeze/Fungicide -",
                           ifelse(FreezeThaw=="Control" & Fungicide=="Treatment", "Freeze/Fungicide +",
                                  ifelse(FreezeThaw=="Treatment" & Fungicide=="Control", "Freeze-Thaw/Fungicide -",
                                         "Freeze-Thaw/Fungicide +")))) %>% 
ggplot(aes(prop.fung, prop.dead))+
  geom_point(data= . %>% filter(Treatment=="Freeze/Fungicide -"), size=2.5, shape=21, stroke=1.5, aes(color=Species, fill=Species))+
  geom_point(data= . %>% filter(Treatment=="Freeze/Fungicide +"), size=2.5, shape=1,stroke=1.5, aes(color=Species, fill=Species))+
  geom_point(data= . %>% filter(Treatment=="Freeze-Thaw/Fungicide -"), size=2.5, shape=22, stroke=1.5, aes(color=Species, fill=Species))+
  geom_point(data= . %>% filter(Treatment=="Freeze-Thaw/Fungicide +"), size=2.5, shape=0, stroke=1.5, aes(color=Species, fill=Species))+
  stat_smooth(method="lm", se=F, color="black")+
  stat_smooth(data=. %>% filter(Species%!in%"PINSTR"),
              method="lm", se=F, color="black", linetype="dashed")+
  labs(x="Proportion of seeds with fungal growth", y="Proportion of seeds killed")+
  scale_fill_manual(labels=sp.label, values = c(clrs[1], clrs[4], clrs[5], clrs[8]))+
  scale_color_manual(guide=F, values = c(clrs[1], clrs[4], clrs[5], clrs[8]))+
  theme_nice()+
  theme(text = element_text(size=16),
        legend.text = element_text(face = "italic"))
dev.off()
