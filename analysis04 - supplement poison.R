## Supplement Poison Analysis: Numba-Kausi Baits & nests - R script by Phil, April 2024


# This part plots the effects of the poisoning on baiting and bamboos

### Associated csv files:

# BaitsData.csv: data on the baiting experiment
# NestsData.csv: data on bamboon nesting
# Plots_2022_AntProject.csv: Plot attributes
# TreeAttributes.csv: Tree attributes

#----------------------------------------------------------#
### List of R-packages
#----------------------------------------------------------#

package_list <- 
  c("dplyr",
    "tidyr",
    "stringr",
    "reshape2",
    "ggplot2",
    "ggthemes",
    'ggpubr',
    "glmmTMB",
    "corrplot",
    "plotrix",
    "DHARMa")

# install all packages
#sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# Citations
#sapply(package_list, citation)

set.seed(1234)
setwd("~/GitHub/MtWilhelmBamboos")


---------------------------------------------------------#
  # Baits analysis for poison plots -----
#----------------------------------------------------------#

# Get data
baiting.raw<-read.csv(file="raw data/BaitsData.csv", header=T)

# Define a function to parse the abundances
parse_abundance <- function(x) {
  # If the abundance contains a ">" symbol, remove it
  x <- gsub(">", "", x)
  # If the abundance contains "to" or "-", split it into range
  if (grepl("to|-", x)) {
    range <- strsplit(x, "to|-")[[1]]
    # Calculate the mean of the range, represented as whole number
    round(mean(as.numeric(range)))
  } else {
    # If the abundance is a single value, convert it to numeric
    as.numeric(x)
  }
}

### Raw Baiting data
# Import raw occurrence data

# calculate abundance
baiting.raw$Abundance <- apply(baiting.raw["Abundance"], 1, parse_abundance)

# define empty baits as 0 abundance
baiting.raw$Abundance[is.na(baiting.raw$Abundance)] <- 0

# remove unused data
baiting.raw<-subset(baiting.raw, Season!="Baiting")
baiting.raw<-subset(baiting.raw, Season!="Kausi September Baiting")

# select relevant columns
baiting.data <- baiting.raw %>%
  dplyr::select(Block, Plot, Subplot, Stratum, Code, Forest, Season, AntSpCODE, Abundance, Bait.missing, Season)

# remove missing baits
baiting.data<-subset(baiting.data, Bait.missing==FALSE)

# remove doubles in the data (if two species were at one bait)
baiting.incidence<-baiting.data%>%
  group_by(Block, Plot, Subplot, Stratum, Code, Forest, Season) %>%
  summarize(Abundance = sum(Abundance))

# make new variable which counts bait occupancy as binary variable
baiting.incidence<- baiting.incidence %>%
  mutate(occupancy = case_when(Abundance == '0' ~ 0,
                               Abundance != '0' ~ 1)) # if ant was is there it counts it as 1, if not 0
baiting.incidence<-as.data.frame(baiting.incidence)
# sort
baiting.incidence$Season<-as.factor(baiting.incidence$Season)
baiting.incidence$Season <- relevel(baiting.incidence$Season, "Kausi Pre-Poison Baiting")

labs=c('pre-poison', 'post-poison')

# Bait abundance plot
bait.abundance.poison<-ggplot(baiting.incidence, aes(x=Season, y=log(Abundance+1), fill = Stratum)) +
  ggtitle("Bait abundance") +
  ylab("ln(number of ants+1)")+
  xlab("")+ 
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
bait.abundance.poison

# average poison
mean(subset(baiting.incidence, Season=="Kausi Pre-Poison Baiting")$Abundance)
std.error(subset(baiting.incidence, Season=="Kausi Pre-Poison Baiting")$Abundance)

mean(subset(baiting.incidence, Season=="Kausi Post-Poison Baiting")$Abundance)
std.error(subset(baiting.incidence, Season=="Kausi Post-Poison Baiting")$Abundance)

---------------------------------------------------------#
  # Bamboo analysis for poison plots -----
#----------------------------------------------------------#

# Get data
phase.nests <- read.csv(file="phase.nests.csv", header=T)

# define poison/no poisoned plots
phase.nest.poison <-phase.nests %>%
  mutate(poison = case_when(Plot == "KP1-S3" ~ "poison",
                            Plot == "KP2-S3" ~ "poison",
                            Plot == "KP4-S3" ~ "poison",
                            Plot == "KP6-S3" ~ "poison",
                            Plot == "KP1-S2" ~ "no poison",
                            Plot == "KP2-S2" ~ "no poison",
                            Plot == "KP4-S2" ~ "no poison",
                            Plot == "KP6-S2" ~ "no poison",
                            Plot == "KP1-S1" ~ "no poison",
                            Plot == "KP2-S1" ~ "no poison",
                            Plot == "KP4-S1" ~ "no poison",
                            Plot == "KP6-S1" ~ "no poison"))

### Poison plot effects plot, as average proportion per plot
phase1.n<-subset(phase.nest.poison, phase=="phase 1")
phase2.n<-subset(phase.nest.poison, phase=="phase 2")

#count total bamboos per plot per phase
phase1.c<-count(phase1.n, Plot)
phase2.c<-count(phase2.n, Plot)

# summarize occupied per plot per phase
phase1 <- phase1.n %>%
  group_by(poison, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase1$phase <- "phase 1"

# merge counts
phase1<-merge(phase1, phase1.c)

# same for phase 2
phase2<- phase2.n %>%
  group_by(poison, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase2$phase<-"phase 2"
phase2<-merge(phase2, phase2.c)

# merge
phase.final<-rbind(phase1, phase2)
phase.final <- na.omit(phase.final)

# proportion
phase.final$proportion<-phase.final$occupancy/phase.final$n*100
phase.final$poison<-as.factor(phase.final$poison)

phase2<-subset(phase.final, phase=="phase 2")

# average poison
mean(subset(phase.final, poison=="poison")$proportion)
std.error(subset(phase.final, poison=="poison")$proportion)

mean(subset(phase.final, poison=="no poison")$proportion)
std.error(subset(phase.final, poison=="no poison")$proportion)

mean(subset(phase2, poison=="poison")$proportion)
std.error(subset(phase.final, poison=="poison")$proportion)

mean(subset(phase2, poison=="no poison")$proportion)
std.error(subset(phase.final, poison=="poison")$proportion)

# Bamboo nest occupancy by poison
nest.occupancy.poison.plot<-ggplot(phase.final, aes(x=poison, y=proportion, fill=phase)) +
  ggtitle("Nest occupancy in lowland") +
  labs(fill = "Phase")+  
  ylab("%")+
  xlab("")+ 
  ylim(0,50)+
  geom_point(aes(fill=phase), size=3, shape=21, colour="grey20", position=position_jitterdodge(0.1))+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete()+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.occupancy.poison.plot



