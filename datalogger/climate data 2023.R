##### R Script Phil March 2024
## Climate data visualization of data collected in 2023


# Two pairs (canopy + udnerstorey) were located close to the ant collection points in 2023.

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

# Kausi data
KP23.A.top.raw<-read.csv(file="./Datalogger/Kausi2023/top/data_91221501_2024_03_22_0.csv", header=F, sep = ';', dec = ",")
KP23.A.bottom.raw<-read.csv(file="./Datalogger/Kausi2023/bottom/data_94231387_2024_03_22_0.csv", header=F, sep = ';', dec = ",")

KP23.B.top.raw<-read.csv(file=".//Datalogger/Kausi2023/top/data_91221507_2024_03_22_0.csv", header=F, sep = ';', dec = ",")
KP23.B.bottom.raw<-read.csv(file=".//Datalogger/Kausi2023/bottom/data_94231388_2024_03_22_0.csv", header=F, sep = ';', dec = ",")

# Numba data
NP23.A.top.raw<-read.csv(file=".//Datalogger/Numba2023/top/data_91221503_2024_03_22_0.csv", header=F, sep = ';', dec = ",")
NP23.A.bottom.raw<-read.csv(file=".//Datalogger/Numba2023/bottom/data_94231389_2024_03_22_0.csv", header=F, sep = ';', dec = ",")

NP23.B.top.raw<-read.csv(file=".//Datalogger/Numba2023//top/data_91221512_2024_03_22_0.csv", header=F, sep = ';', dec = ",")
NP23.B.bottom.raw<-read.csv(file=".//Datalogger/Numba2023//bottom/data_94231386_2024_03_22_0.csv", header=F, sep = ';', dec = ",")

## merge climate data

# Kausi
KP23A.top<- KP23.A.top.raw %>% dplyr::select(V2,V4)%>% rename(time = V2,
         temp.top = V4)
KP23A.bottom<- KP23.A.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
KP23A.logger<-merge(KP23A.top, KP23A.bottom)
KP23A.logger$Plot<-"Kausi.A"
KP23A.logger$Location<-"Kausi"

# date subset- time of data loggers: 04. sep 2023 - 15:00 - 09 sep 2023 15:00
KP23A.logger<-KP23A.logger[c(54331:54811),]

# KP23 B
KP23B.top<- KP23.B.top.raw %>% dplyr::select(V2,V4)%>% rename(time = V2,
                                                              temp.top = V4)
KP23B.bottom<- KP23.B.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
KP23B.logger<-merge(KP23A.top, KP23A.bottom)
KP23B.logger$Plot<-"Kausi.B"
KP23B.logger$Location<-"Kausi"

# date subset- time of data loggers: 05. sep 2023 - 00:00 - 08 sep 2023 00:00
KP23B.logger<-KP23B.logger[c(54331:54811),]

# Numba
NP23A.top<- NP23.A.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
NP23A.bottom<- NP23.A.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
NP23A.logger<-merge(NP23A.bottom, NP23A.top)
NP23A.logger$Plot<-"Numba.A"
NP23A.logger$Location<-"Numba"

# date subset- time of data loggers:  19 aug 2023 00:00 - 30 Aug 2023 00:00
NP23A.logger<-NP23A.logger[c(53512:54568),]
# 
NP23B.top<- NP23.B.top.raw %>% dplyr::select(V2,V4)%>%
  rename(time = V2,
         temp.top = V4)
NP23B.bottom<- NP23.B.bottom.raw %>% dplyr::select(V2,V4:V7)%>%
  rename(time = V2,
         temp.bottom1 = V4,
         temp.bottom2 = V5,
         temp.bottom3 = V6,
         humidity.bottom = V7)
NP23B.logger<-merge(NP23B.bottom, NP23B.top)
NP23B.logger$Plot<-"Numba.B"
NP23B.logger$Location<-"Numba"

# date subset- time of data loggers:  19 aug 2023 00:00 - 30 Aug 2023 00:00
NP23B.logger<-NP23B.logger[c(53231:54287),]
##
## bottom1 = deepest, bottom 3= highest position
datalogger<-bind_rows(KP23A.logger, KP23B.logger,NP23A.logger, NP23B.logger)


datalogger$temp.top<-as.numeric(datalogger$temp.top)
datalogger$temp.bottom1<-as.numeric(datalogger$temp.bottom1)
datalogger$temp.bottom2<-as.numeric(datalogger$temp.bottom2)
datalogger$temp.bottom3<-as.numeric(datalogger$temp.bottom3)


boxplot(datalogger$temp.top~datalogger$Location)
boxplot(datalogger$temp.bottom1~datalogger$Location)#inside soil
boxplot(datalogger$temp.bottom2~datalogger$Location) #on the soil surface (more realistic what ants experience)

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
  ylab("Temperature [C]")+
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

groundtemp<-ggplot(datalogger, aes(x=Plot, y=temp.bottom2, fill=Location)) +
  ggtitle("Ground level Temperature") +
  geom_boxplot()+
  ylim(15,35)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Temperature [?C]")+
  xlab("")+
  theme_hc()+
  theme(axis.text.x = element_text(size=12))
groundtemp

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

# long format
datalogger.temp<-datalogger[,c(1,2,4,7,8)]
datalogger.templong <- tidyr::gather(datalogger.temp, key = "Stratum", value = "temperature", -c(time, Location, Plot))

labs <- expression("lowland", "midelevation")


# stratum 
temp_stratum<-ggplot(datalogger.templong, aes(x=Location, y=temperature, fill=Stratum)) +
  ggtitle("Temperature [°C]") +
  geom_violin(lwd=1)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(labels=c("Understory", "Canopy"), values=c("#E69F00", "#0072B2"))+
  ylab("")+
  #ylim(0,35)+
  xlab("")+
  theme_minimal(15)
temp_stratum

# average temp
means_by_location <- datalogger.templong %>%
  group_by(Location, Stratum) %>%
  summarise(T_mean = mean(temperature), sd = sd(temperature))
means_by_location
