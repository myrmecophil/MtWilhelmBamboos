##### R Script [presentation BES] Phil Nov 2022
## Climate data visualization

package_list <- 
  c("rlang",
    "dplyr",
    "tidyr",
    "stringr",
    "plotrix",
    "ggpubr",
    "ggplot2",    
    "ggrepel",
    "ggthemes",
    "lmerTest",
    "lme4") 

# install all packages
#sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# load data

# three canopy data logger files in NP1, but only one ground
# but one of them, NP1B top, is same as NP2B top 
# therefore, NP1B bottom (data_94231383_2022_09_23_0.csv) missing - now added

KP1.top.raw<-read.csv(file="./Datalogger/KP1/data_91221507_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
KP1.bottom.raw<-read.csv(file="./Datalogger/KP1/data_94231387_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

KP2.top.raw<-read.csv(file="./Datalogger/KP2/data_91221508_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
KP2.bottom.raw<-read.csv(file="./Datalogger/KP2/data_94231388_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

KP4.top.raw<-read.csv(file="./Datalogger/KP4/data_91221509_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
KP4.bottom.raw<-read.csv(file="./Datalogger/KP4/data_94231390_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

KP6.top.raw<-read.csv(file="./Datalogger/KP6/data_91221510_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
KP6.bottom.raw<-read.csv(file="./Datalogger/KP6/data_94231389_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

NP1A.top.raw<-read.csv(file="./Datalogger/NP1/data_91221502_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
#NP1B.top.raw<-read.csv(file="./Datalogger/NP1/data_91221505_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP1B.top.raw<-read.csv(file="./Datalogger/NP1/data_91221512_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

NP1A.bottom.raw<-read.csv(file="./Datalogger/NP1/data_94213981_2022_04_20_0.csv", header=F, sep = ';', dec = ",")
NP1B.bottom.raw<-read.csv(file="./Datalogger/NP1/data_94231383_2022_11_30_0.csv", header=F, sep = ';', dec = ",")

NP2A.top.raw<-read.csv(file="./Datalogger/NP2/data_91221501_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP2B.top.raw<-read.csv(file="./Datalogger/NP2/data_91221505_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP2A.bottom.raw<-read.csv(file="./Datalogger/NP2/data_94213982_2022_04_20_0.csv", header=F, sep = ';', dec = ",")
NP2B.bottom.raw<-read.csv(file="./Datalogger/NP2/data_94231385_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

NP4A.top.raw<-read.csv(file="./Datalogger/NP4/data_91221504_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP4B.top.raw<-read.csv(file="./Datalogger/NP4/data_91221511_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP4A.bottom.raw<-read.csv(file="./Datalogger/NP4/data_94231386_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP4B.bottom.raw<-read.csv(file="./Datalogger/NP4/data_94231392_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

NP5A.top.raw<-read.csv(file="./Datalogger/NP5/data_91221503_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP5B.top.raw<-read.csv(file="./Datalogger/NP5/data_91221506_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP5A.bottom.raw<-read.csv(file="./Datalogger/NP5/data_94231384_2022_09_23_0.csv", header=F, sep = ';', dec = ",")
NP5B.bottom.raw<-read.csv(file="./Datalogger/NP5/data_94231391_2022_09_23_0.csv", header=F, sep = ';', dec = ",")

## merge climate data

# KP 1
KP1.top<- KP1.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
KP1.bottom<- KP1.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
KP1.logger<-merge(KP1.top, KP1.bottom)
KP1.logger$Plot<-"KP1"
KP1.logger$Location<-"Kausi"

# date subset- approximate time of data loggers: 12. june - 2022.09.20 23:30 (need to be confirmed!)
KP1.logger<-KP1.logger[c(11201:20867),]

# KP 2
KP2.top<- KP2.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
KP2.bottom<- KP2.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
KP2.logger<-merge(KP2.top, KP2.bottom)
KP2.logger$Plot<-"KP2"
KP2.logger$Location<-"Kausi"

# date subset / 2022.06.12 06:30 - 2022.09.21 20:30
KP2.logger<-KP2.logger[c(11195:20947),]

# KP 4
KP4.top<- KP4.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
KP4.bottom<- KP4.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
KP4.logger<-merge(KP4.top, KP4.bottom)
KP4.logger$Plot<-"KP4"
KP4.logger$Location<-"Kausi"
# date subset / 2022.06.15 04:15 - 2022.09.20 23:15
KP4.logger<-KP4.logger[c(11473:20861),]

# KP 6
KP6.top<- KP6.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
KP6.bottom<- KP6.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
KP6.logger<-merge(KP6.top, KP6.bottom)
KP6.logger$Plot<-"KP6"
KP6.logger$Location<-"Kausi"
# date 2022.06.15 02:30 - 2022.09.21 01:45
KP6.logger<-KP6.logger[c(11464:20869),]

### Numba
# NP 1 - using B since A has wrong dates
NP1.top<- NP1B.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
NP1.bottom<- NP1B.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
NP1.logger<-merge(NP1.top, NP1.bottom)
NP1.logger$Plot<-"NP1"
NP1.logger$Location<-"Numba"

# date subset- approximate time of data loggers: 2022.06.23 03:30 -2022.09.17 01:30
NP1.logger<-NP1.logger[c(12237:20485),]

# NP2 - using B since A has date wrong
NP2.top<- NP2B.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
NP2.bottom<- NP2B.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
NP2.logger<-merge(NP2.top, NP2.bottom)
NP2.logger$Plot<-"NP2"
NP2.logger$Location<-"Numba"

# date subset- approximate time of data loggers: 2022.06.22 01:45 - 2022.09.17 01:45
NP2.logger<-NP2.logger[c(12141:20493),]

# NP 4
NP4.top<- NP4A.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
NP4.bottom<- NP4A.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
NP4.logger<-merge(NP4.top, NP4.bottom)
NP4.logger$Plot<-"NP4"
NP4.logger$Location<-"Numba"

# date subset- approximate time of data loggers: 2022.06.07 00:45 - 2022.09.17 01:30
NP4.logger<-NP4.logger[c(10697:20492),]

# NP5
NP5.top<- NP5A.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
NP5.bottom<- NP5A.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
NP5.logger<-merge(NP5.top, NP5.bottom)
NP5.logger$Plot<-"NP5"
NP5.logger$Location<-"Numba"

# date subset- approximate time of data loggers: 2022.06.23 02:00 - 2022.09.17 00:45
NP5.logger<-NP5.logger[c(12234:20485),]

##
## bottom1 = deepest, bottom 3= highest position
datalogger<-bind_rows(KP1.logger, KP2.logger, KP4.logger, KP6.logger, NP1.logger, NP2.logger, NP4.logger, NP5.logger)
boxplot(datalogger$temp.top~datalogger$Location)
boxplot(datalogger$temp.bottom1~datalogger$Location)

boxplot(datalogger$humidity.bottom~datalogger$Location)

boxplot(datalogger$temp.top~datalogger$Plot)
boxplot(datalogger$temp.bottom1~datalogger$Plot)
boxplot(datalogger$temp.bottom2~datalogger$Plot)
boxplot(datalogger$temp.bottom3~datalogger$Plot)
boxplot(datalogger$humidity.bottom~datalogger$Plot)

kausi.datalogger <-subset(datalogger, Location == "Kausi")
numba.datalogger <-subset(datalogger, Location == "Numba")

toptemp<-ggplot(datalogger, aes(x=Plot, y=temp.top, fill=Location)) +
  ggtitle("Canopy Temperature") +
  geom_boxplot()+
  ylim(10,35)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Temperature [?C]")+
  xlab("")+
  theme_hc()+
  theme(axis.text.x = element_text(size=12))
toptemp

soiltemp<-ggplot(datalogger, aes(x=Plot, y=temp.bottom1, fill=Location)) +
  ggtitle("Soil Temperature") +
  geom_boxplot()+
  ylim(10,26)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Temperature [?C]")+
  xlab("")+
  theme_hc()+
  theme(axis.text.x = element_text(size=12))
soiltemp

humidity<-ggplot(datalogger, aes(x=Plot, y=humidity.bottom, fill=Location)) +
  ggtitle("Soil Humidity") +
  geom_boxplot()+
  #ylim(10,26)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("raw TDT")+
  xlab("")+
  theme_hc()+
  theme(axis.text.x = element_text(size=12))
humidity
