
## Numba-Kausi Baits & nests - R script by Phil 18 April 2023

# to do: 
# optimize all models
# include soil humidity in environmental data
# delete survivals and see if sth changes
# nesting abundance estimate only occupied nests

rm(list=ls()) 

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
    "vegan",
    "ggplot2",
    "ggpubr",
    "glmmTMB",
    "lme4",
    "lmerTest",
    "DHARMa",
    "corrplot",
    "car",
    "emmeans",
    "effects",
    "broom.mixed",
    "openxlsx")

# install all packages
#sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# Citations
#sapply(package_list, citation)

#----------------------------------------------------------#
# raw data transformation -----
#----------------------------------------------------------#

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
baiting.raw =read.csv(file="raw data/BaitsData.csv", header=T)

# remove post poison data
baiting.raw<-subset(baiting.raw, Season!="Kausi September Baiting")
baiting.raw<-subset(baiting.raw, Season!="Kausi Post-Poison Baiting")

# calculate abundance
baiting.raw$Abundance <- apply(baiting.raw["Abundance"], 1, parse_abundance)

# define empty baits as 0 abundance
baiting.raw$Abundance[is.na(baiting.raw$Abundance)] <- 0

# select relevant columns
baiting.data <- baiting.raw %>%
  dplyr::select(Block, Plot, Subplot, Stratum, Code, Forest, Season, AntSpCODE, Abundance, Bait.missing)

# remove missing baits
baiting.data<-subset(baiting.data, Bait.missing==FALSE)

# remove doubles in the data (if two species were at one bait)
baiting.incidence<-distinct(baiting.data, Code,  .keep_all = TRUE)

# make new variable which counts bait occupancy as binary variable
baiting.incidence<- baiting.incidence %>%
  mutate(occupancy = case_when(AntSpCODE == '' ~ 0,
                               AntSpCODE != '' ~ 1)) # if ant was is there it counts it as 1, if not 0


### Raw nesting data

# Import bamboo nesting data
nest.raw =read.csv(file="raw data/NestsData.csv", header=T)

# split column 'moved to' into two separate columns
nest.raw$Moved.to.block <- str_split(nest.raw$Moved.to.code, "-", simplify = TRUE)[,1]
nest.raw$Moved.to.plot <- str_split(nest.raw$Moved.to.code, "-", simplify = TRUE)[,1:2] %>% apply(1, paste, collapse = "-")

nest.raw$Moved.to.plot<-as.factor(nest.raw$Moved.to.plot)
nest.raw$Moved.to.block<-as.factor(nest.raw$Moved.to.block)

# remove bamboo data from poisoned plots (KX-S3)
nest.raw2<-subset(nest.raw, Moved.to.plot!= "KP1-S3" & Moved.to.plot!= "KP2-S3" & Moved.to.plot!= "KP4-S3"& Moved.to.plot!= "KP6-S3")

# parse number of nesting ants
nest.raw2$nesting.estimate <- apply(nest.raw2["N.ants..camera."],1, parse_abundance)

# Count a bamboo nest as occupied if there were >2 worker ants inside or 1 queen
# add 'moved to' forest and block
nesters <-nest.raw2 %>%
  mutate(phase1.occupancy = case_when(nesting.estimate > 2 ~ 1,
                                      TRUE ~ 0),
         phase2.occupancy = case_when(N.ants > 2 ~ 1,
                                      N.queens > 0 ~ 1,
                                      TRUE ~ 0), # if single queen, it counts as occupied
         moved.to.forest  = case_when(Moved.to.plot == 'NP1-S1' | Moved.to.plot == 'NP1-S2'| Moved.to.plot == 'NP1-S3'| Moved.to.plot == 'NP1-S4'| Moved.to.plot == 'NP1-S5' ~ 'Numba Primary',
                                      Moved.to.plot == 'NP2-S1' | Moved.to.plot == 'NP2-S2'| Moved.to.plot == 'NP2-S3'| Moved.to.plot == 'NP2-S4'| Moved.to.plot == 'NP2-S5' ~ 'Numba Primary',
                                      Moved.to.plot == 'NP4-S1' | Moved.to.plot == 'NP4-S2'| Moved.to.plot == 'NP4-S3'| Moved.to.plot == 'NP4-S4'| Moved.to.plot == 'NP4-S5' ~ 'Numba Primary',
                                      Moved.to.plot == 'NP5-S1' | Moved.to.plot == 'NP5-S2'| Moved.to.plot == 'NP5-S3'| Moved.to.plot == 'NP5-S4'| Moved.to.plot == 'NP5-S5' ~ 'Numba Primary',
                                      Moved.to.plot == 'KP1-S1' | Moved.to.plot == 'KP1-S2'| Moved.to.plot == 'KP1-S3'| Moved.to.plot == 'KP1-S4'| Moved.to.plot == 'KP1-S5' ~ 'Kausi Primary',
                                      Moved.to.plot == 'KP2-S1' | Moved.to.plot == 'KP2-S2'| Moved.to.plot == 'KP2-S3'| Moved.to.plot == 'KP2-S4'| Moved.to.plot == 'KP2-S5' ~ 'Kausi Primary',
                                      Moved.to.plot == 'KP4-S1' | Moved.to.plot == 'KP4-S2'| Moved.to.plot == 'KP4-S3'| Moved.to.plot == 'KP4-S4'| Moved.to.plot == 'KP4-S5' ~ 'Kausi Primary',
                                      Moved.to.plot == 'KP6-S1' | Moved.to.plot == 'KP6-S2'| Moved.to.plot == 'KP6-S3'| Moved.to.plot == 'KP6-S4'| Moved.to.plot == 'KP6-S5' ~ 'Kausi Primary'))


## Separate data sets, then put them back together so that counts and occupancies of ants from phase 1 and 2 are in the same column

phase1.nests <- nesters %>%
  dplyr::select(Block, Plot, Stratum, Code, Forest, Treatment, AntSpCode, nesting.estimate, phase1.occupancy,LOST.REPLACED)%>%
  rename(occupancy = phase1.occupancy, lost = LOST.REPLACED)

phase2.nests <- nesters %>%
  dplyr::select(Moved.to.block, Moved.to.plot, Stratum, Code, moved.to.forest, Treatment, AntSpCODE.final, N.ants, phase2.occupancy, LOST.FINAL) %>%
  rename(Block = Moved.to.block, Plot = Moved.to.plot, Forest = moved.to.forest, occupancy = phase2.occupancy,lost = LOST.FINAL, AntSpCode = AntSpCODE.final, nesting.estimate=N.ants)

phase1.nests$phase <- "phase 1"
phase2.nests$phase <- "phase 2"

phase.nests <-rbind(phase1.nests, phase2.nests)
phase.nests$phase<-as.factor(phase.nests$phase)

# remove lost bamboos
phase.nests<-subset(phase.nests,lost != 1)

# define NA abundance as true 0
phase.nests$nesting.estimate[is.na(phase.nests$nesting.estimate)] <- 0

#----------------------------------------------------------#
#  Plot attributes -----
#----------------------------------------------------------#

# upload raw plot data
plot.metada.raw<-read.csv(file="raw data/Plots_2022_AntProject.csv", header=T)
plot.metada.raw$Block<-as.factor(plot.metada.raw$Block)

# calculate slope variation as mean of the 4 absolute slope values
plot.metada.raw$slope.var <- apply(abs(plot.metada.raw[,c(6:9)]), 1, mean)

# upload raw tree dat
tree.metada.raw<-read.csv(file="raw data/TreeAttributes.csv", header=T)

# turn number of lianas into numeric
tree.metada.raw$Lianas.n <- as.numeric(gsub('>', '', tree.metada.raw$Lianas))

# define plot.stratum 
tree.metada.raw$plot.stratum<-paste(tree.metada.raw$Plot,tree.metada.raw$Stratum)

# Prepare tree-level metadata for later models
tree.meta.sub<-tree.metada.raw %>%
  dplyr::select(Plot, plot.stratum, mastercode, Trunk.perimenter..cm., tree.ID, Lianas.n, deadwood.., deadwood., Bait.Bamboo.Height.m.)%>%
  rename(trunk =Trunk.perimenter..cm., Code = mastercode, height = Bait.Bamboo.Height.m., dw.percent = deadwood., dw.number=deadwood..)

plot.metada.sub<-  plot.metada.raw %>%
  dplyr::select(Plot_code, Block, forest_type, Canopy.cover..Caco., elevation, slope.var)%>%
  rename(Plot=Plot_code, Forest = forest_type, Caco =Canopy.cover..Caco.)

## calculate tree metadata as means on plot level
# first get rid of tree duplicates measures
tree.meta.sub.ID<-distinct(tree.meta.sub, tree.ID, .keep_all = TRUE)

tree.meta.plot<-tree.meta.sub.ID %>%
  group_by(Plot) %>%  summarise(across(
    .cols = where(is.numeric), 
    .fns = list(mean = ~mean(., na.rm = TRUE)),
    .names = "{col}_{fn}"))

# tree attributes incl. plot attributes
tree.meta<-merge(tree.meta.sub,plot.metada.sub)

# plot attributes incl. tree level means
plot.meta<-merge(plot.metada.sub,tree.meta.plot)

# Prepare stratum-level metadata for later models

## calculate tree metadata as means on stratum level
tree.meta.plot2<-tree.meta.sub.ID %>%
  group_by(Plot, plot.stratum) %>%  summarise(across(
    .cols = where(is.numeric), 
    .fns = list(mean = ~mean(., na.rm = TRUE)),
    .names = "{col}_{fn}"))

# stratum level plot attributes (averages)
plot.meta2<-merge(tree.meta.plot2,plot.metada.sub)

#----------------------------------------------------------#
# 1.1 Bait occupancy -----
#----------------------------------------------------------#

## count the double occupation of baits (number of species per bait)
bait.double <- baiting.data %>%
  group_by(Forest, Block, Plot, Code) %>%
  summarize(species.per.bait.raw = n())

# to remove baits with no ants, merge with incidence data
bait.double<-merge(bait.double, baiting.incidence)

# if incidence is 0, put also 0 species, if 2 species on bait, put 2 otherwise 1 
bait.double <-bait.double %>%
  mutate(species.per.bait = case_when(occupancy == 0 ~ 0,
                                      species.per.bait.raw > 1 ~ 2,
                                       .default = 1)) #


### Plot bait occupancy as average proportion per plot-stratum

# define plot-stratum
baiting.incidence$plot.stratum <-paste(baiting.incidence$Plot,baiting.incidence$Stratum)

#count total baits per stratum per plot
baits.c<-count(baiting.incidence, plot.stratum)

# summarize occupied per stratum per plot
bait1 <- baiting.incidence %>%
  group_by(Forest, Plot, Stratum, plot.stratum) %>%
  summarize(occupancy = sum(occupancy))

# merge with total baits
baits.final<-merge(bait1, baits.c)

# proportion
baits.final$proportion<-baits.final$occupancy/baits.final$n*100

# plot it
bait.occupancy.stratum<-ggplot(baits.final, aes(x=Forest, y=proportion, fill = Stratum)) +
  ggtitle("Bait occupancy") +
  ylab("Occupancy [%]")+
  xlab("")+ 
  #ylim(0,50)+
  geom_boxplot()+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  theme_bw()
bait.occupancy.stratum

### Plot bait occupancy as average proportion per plot

#count total baits per plot
baits.c<-count(baiting.incidence, Plot)

# summarize occupied per plot
bait1 <- baiting.incidence %>%
  group_by(Forest, Plot) %>%
  summarize(occupancy = sum(occupancy))

# merge with total baits
baits.final<-merge(bait1, baits.c)

# proportion
baits.final$proportion<-baits.final$occupancy/baits.final$n*100

# plot it
bait.occupancy.plot<-ggplot(baits.final, aes(x=Forest, y=proportion, fill = Forest)) +
  ggtitle("Bait occupancy") +
  ylab("Occupancy [%]")+
  xlab("")+ 
  #ylim(0,50)+
  geom_boxplot()+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  theme_bw()
bait.occupancy.plot

# plot as glm prediction against elevation

ggplot(baiting.incidence.e, aes(x=elevation, y=occupancy)) + geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)

#----------------------------------------------------------#
# 1.2 Bait diversity -----
#----------------------------------------------------------#

## On plot level
# plot x species incidence
ant_m<-dcast(baiting.data, formula = Plot ~ AntSpCODE, length)
rownames(ant_m)<-ant_m[,1] # Plot as rownames
ant_m<-ant_m[,-c(1,2)]

# extract diversity values
data2 <- data.frame(row.names(ant_m))     # create a new file with plot numbers
names(data2) <- "Plot" # rename variable "plot"
data2$Shannon <- diversity(ant_m, index = "shannon", base = exp(1)) # Shannon diversity per plot
data2$Richness <- specnumber(ant_m)                          # Richness per plot
data2$expH <- exp(data2$Shannon)                                         # exponential Shannon diversity per plot
data2$expH[data2$Richness == 0] <- 0                                     # define exp(Shannon) = 0
data2$evenness <- data2$Shannon/log(specnumber(ant_m)) # evenness per plot
str(data2)

# add Plot location
plot_m<-distinct(baiting.data, Plot, .keep_all = TRUE)
diversity<-merge(data2, plot_m)

## Plot it
#Labs
labs <- expression("Kausi", "Numba")

baitdiversity.plot<-ggplot(diversity, aes(x=Forest, y=expH, fill=Forest)) +
  ggtitle("Bait species diversity") +
  geom_boxplot()+
  ylim(0,20)+
  scale_fill_manual(labels=labs, values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=labs)+
  ylab("Species diversity [expH]")+
  xlab("")+ 
  guides(fill="none")+
  theme_bw()
baitdiversity.plot

## On stratum level
baiting.data$plot.stratum <-paste(baiting.data$Plot,baiting.data$Stratum)

# plot.stratum x species incidence
ant_m<-dcast(baiting.data, formula = plot.stratum ~ AntSpCODE, length)
rownames(ant_m)<-ant_m[,1] # Plot as rownames
ant_m<-ant_m[,-c(1,2)]

# extract diversity values
data3 <- data.frame(row.names(ant_m))     # create a new file with plot numbers
names(data3) <- "plot.stratum" # rename variable "plot.stratum"
data3$Shannon <- diversity(ant_m, index = "shannon", base = exp(1)) # Shannon diversity per plot
data3$Richness <- specnumber(ant_m)                          # Richness per plot
data3$expH <- exp(data3$Shannon)                                         # exponential Shannon diversity per plot
data3$expH[data3$Richness == 0] <- 0                                     # define exp(Shannon) = 0
data3$evenness <- data3$Shannon/log(specnumber(ant_m)) # evenness per plot
str(data3)

# add Plot location
plot_m2<-distinct(baiting.data, plot.stratum, .keep_all = TRUE)
diversity.stratum<-merge(data3, plot_m2)
# 

## Plot it
#Labs
labs <- expression("Kausi", "Numba")

baitdiversity.stratum<-ggplot(diversity.stratum, aes(x=Forest, y=expH, fill = Stratum)) +
  ggtitle("Bait diversity") +
  ylab("Shannon diversity")+
  xlab("")+ 
  #ylim(0,50)+
  geom_boxplot()+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  theme_bw()
baitdiversity.stratum

#----------------------------------------------------------#
# 1.3 Bait composition -----
#----------------------------------------------------------#

# summarize bait incidence counts
bait.abundance <- baiting.data %>%
  group_by(Forest, AntSpCODE) %>%
  summarize(Abundance = sum(Abundance))

# Remove empty baits
bait.abundance <- bait.abundance[bait.abundance$AntSpCODE != "", ]

# species with NA abundance have at least 1 worker, so set to 1
bait.abundance$Abundance[is.na(bait.abundance$Abundance)] <-1

n.kausi<-sum(bait.abundance$Abundance[bait.abundance$Forest=="Kausi Primary"])
n.numba<-sum(bait.abundance$Abundance[bait.abundance$Forest=="Numba Primary"])

bait.abundance$percentage<-0

bait.abundance$percentage[bait.abundance$Forest=="Kausi Primary"] <-bait.abundance$Abundance[bait.abundance$Forest=="Kausi Primary"]/n.kausi*100
bait.abundance$percentage[bait.abundance$Forest=="Numba Primary"] <-bait.abundance$Abundance[bait.abundance$Forest=="Numba Primary"]/n.numba*100

# check sums
sum(bait.abundance$percentage[bait.abundance$Forest=="Kausi Primary"])
sum(bait.abundance$percentage[bait.abundance$Forest=="Numba Primary"])

# Define "Rest" group at <3% of total abundance 
bait.abundance$AntSpCODE[bait.abundance$percentage < 3] <- "other species"

bait.abundance$AntSpCODE<-as.factor(bait.abundance$AntSpCODE)
bait.abundance$Forest<-as.factor(bait.abundance$Forest)

# summarize bait percentage counts
bait.abundance2<- bait.abundance%>% 
  group_by(Forest, AntSpCODE) %>% 
  summarise(percentage = sum(percentage))

# Plot it

# Change order of levels
bait.abundance2$AntSpCODE <- relevel(bait.abundance2$AntSpCODE, "CREM 014")
bait.abundance2$AntSpCODE <- relevel(bait.abundance2$AntSpCODE, "CREM 003")
bait.abundance2$AntSpCODE <- relevel(bait.abundance2$AntSpCODE, "other species")

abundance.plot<-ggplot(bait.abundance2, aes(x = Forest, y = percentage, fill = AntSpCODE)) +
  geom_bar(stat = "identity",position= position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative abundance [%]") +
  ggtitle("bait species composition") +
  theme_minimal()
abundance.plot

 # summarize bait incidence counts
bait.incidence <- baiting.data %>%
  group_by(Forest, AntSpCODE,) %>%
  summarize(count = n())

# Remove empty baits
#bait.abundance <- bait.abundance[bait.abundance$AntSpCODE != "", ]

# include empty baits
bait.incidence$AntSpCODE[bait.incidence$AntSpCODE == ""]<-'baits empty'

n.kausi<-sum(bait.incidence$count[bait.incidence$Forest=="Kausi Primary"])
n.numba<-sum(bait.incidence$count[bait.incidence$Forest=="Numba Primary"])

bait.incidence$percentage<-0

bait.incidence$percentage[bait.incidence$Forest=="Kausi Primary"] <-bait.incidence$count[bait.incidence$Forest=="Kausi Primary"]/n.kausi*100
bait.incidence$percentage[bait.incidence$Forest=="Numba Primary"] <-bait.incidence$count[bait.incidence$Forest=="Numba Primary"]/n.numba*100

sum(bait.incidence$percentage[bait.incidence$Forest=="Kausi Primary"])
sum(bait.incidence$percentage[bait.incidence$Forest=="Numba Primary"])

# Define "Rest" group at <3% of total incidence 
bait.incidence$AntSpCODE[bait.incidence$percentage < 3] <- "other species"

# summarize bait incidence counts
bait.incidence2<- bait.incidence%>% 
  group_by(Forest, AntSpCODE) %>% 
  summarise(percentage = sum(percentage))

# Plot it

# Change order of levels
bait.incidence2$AntSpCODE<-as.factor(bait.incidence2$AntSpCODE)
bait.incidence2$AntSpCODE <- relevel(bait.incidence2$AntSpCODE, "CREM 014")
bait.incidence2$AntSpCODE <- relevel(bait.incidence2$AntSpCODE, "CREM 003")
bait.incidence2$AntSpCODE <- relevel(bait.incidence2$AntSpCODE, "other species")
bait.incidence2$AntSpCODE <- relevel(bait.incidence2$AntSpCODE, "baits empty")

bait.incidence2$AntSpCODE<-as.factor(bait.incidence2$AntSpCODE)
bait.incidence2$Forest<-as.factor(bait.incidence2$Forest)

#
bait.incidence.plot<-ggplot(bait.incidence2, aes(x = Forest, y = percentage, fill = AntSpCODE)) +
  geom_bar(stat = "identity",position= position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative incidence [%]") +
  ggtitle("bait species composition") +
  #scale_fill_brewer(palette="Dark2")+
  theme_minimal()
bait.incidence.plot


## Bait rank abundance (incidence) curves
library(ggrepel)
library("BiodiversityR")
data(dune.env)
data(dune)

# forest x species incidence
ant_m<-dcast(baiting.data, formula = Forest ~ AntSpCODE, length)
rownames(ant_m)<-ant_m[,1] # Plot as rownames
ant_m<-ant_m[,c(-2)]

ant_m.env<-data.frame(forest  = c("Kausi Primary", "Numba Primary"))
ant_m.env$forest<-as.factor(ant_m.env$forest)

RA.data <- rankabuncomp(ant_m[,-1], y=ant_m.env, factor='forest', 
                        return.data=TRUE, specnames=c(1:2), legend=FALSE)
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

uncorrected.rank.incidence <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(aes(colour=Grouping, shape=Grouping), size=1, alpha=0.7) +
  geom_text_repel(data=subset(RA.data, labelit == TRUE), 
                  aes(colour=Grouping, label=species), 
                  angle=0, nudge_x=1, nudge_y=1, show.legend=FALSE) +
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  labs(x = "rank", y = "incidence", colour = "Forest", shape = "Forest")
uncorrected.rank.incidence

#----------------------------------------------------------#
# 2.1 Nest occupancy -----
#----------------------------------------------------------#

# Survivor Plots

# create extra variable
nest.raw$Forest.Treatment<-paste(nest.raw$Forest, nest.raw$Treatment)

nest.raw3 <-nest.raw %>%
  mutate(survived.d = case_when(Survived == "TRUE" ~ "survived",
                                Reoccupied == "TRUE" ~ "not survived",
                                Abandoned == "TRUE" ~ "not survived",
                                TRUE ~ "empty"))

ggplot(nest.raw3, aes(fill=survived.d, y=Forest.Treatment, x=Forest.Treatment)) + 
  geom_bar(position="fill", stat="identity")

nest.raw3 %>%
  count(Forest.Treatment, survived.d) %>%
  ggplot(aes(Forest.Treatment, n, fill = survived.d)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("Treatment")+
  scale_x_discrete(labels=labs1)+
  ylab("% Bamboos translocated") +
  ggtitle("Survived as proportion of total bamboos") +
  scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(0, .25))
  #geom_text(aes(label = n), vjust=1, colour = "black", size = 5)

    
labs1 <- expression("low elevation control", "high elevation control", "other elevation transfer", "same elevation transfer")

# Survived as proportion of occupied nests
nest.raw %>% filter(Occupied == TRUE)%>%
  count(Forest.Treatment, Survived) %>%
  ggplot(aes(Forest.Treatment, n, fill = Survived)) + 
  xlab("Treatment") +
  scale_x_discrete(labels=labs1)+
  ylab("occupied nests") +
  ggtitle("Survived as proportion of occupied nests") +
  scale_fill_brewer(palette="Dark2")+
  geom_col() + 
  geom_text(aes(label = n), 
            position = position_stack(vjust=.5),colour = "black", size = 5)

# Survived as proportion of total nests
nest.raw %>%
  count(Forest.Treatment, Survived) %>%
  ggplot(aes(Forest.Treatment, n, fill = Survived)) + 
  xlab("Treatment") +
  scale_x_discrete(labels=labs1)+
  ylab("Bamboos translocated") +
  ggtitle("Survived as proportion of total bamboos") +
  scale_fill_brewer(palette="Dark2")+
  geom_col() + 
  geom_text(aes(label = n), 
            position = position_stack(vjust=.5),colour = "black", size = 5)


# define treatments.

phase.nests <-phase.nests %>%
  mutate(new.Treatment = case_when(phase == "phase 1" ~ "first.placement",
                                   Treatment == "other elevation" ~ "other elevation",
                                   Treatment == "same elevation" ~ "same elevation", 
                                   Treatment == "control" ~ "control"))

phase.nests <-phase.nests %>%
  mutate(new.Treatment = case_when(phase == "phase 1" ~ "first.placement",
                                   Treatment == "other elevation" ~ "elevation",
                                   Treatment == "same elevation" ~ "elevation", 
                                   Treatment == "control" ~ "control"))

# Get environmental metadata
bamboo.incidence.e<-merge(phase.nests, tree.meta)


### Plot nest occupancy as average proportion per plot
phase1.n<-subset(phase.nests, phase=="phase 1")
phase2.n<-subset(phase.nests, phase=="phase 2")

#count total bamboos per plot per phase
phase1.c<-count(phase1.n, Plot)
phase2.c<-count(phase2.n, Plot)
counts<-rbind(phase1.c,phase2.c)

# summarize occupied per plot per phase
phase1 <- phase1.n %>%
  group_by(Forest, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase1$phase <- "phase 1"

phase2<- phase2.n %>%
  group_by(Forest, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase2$phase<-"phase 2"

# merge
phase.t<-rbind(phase1, phase2)

# merge with total nests
phase.final<-merge(phase.t, counts)

# proportion
phase.final$proportion<-phase.final$occupancy/phase.final$n*100

# plot it
nest.occupancy.phase.plot<-ggplot(phase.final, aes(x=Forest, y=proportion, fill=phase)) +
  ggtitle("Bamboo nest occupancy") +
  labs(fill = "Phase")+  
  ylab("Occupancy [%]")+
  xlab("")+ 
  ylim(0,50)+
  geom_boxplot()+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  theme_bw()
nest.occupancy.phase.plot

### plot nest occupancy as proportion per forest
nest.proportion <- phase.nests %>%
  group_by(Forest, AntSpCode, phase) %>%
  summarize(count = n())

nest.proportion$AntSpCode[nest.proportion$AntSpCode==""] <- "empty"
nest.proportion$percentage<-0

n.kausi.phase1<-sum(nest.proportion$count[nest.proportion$Forest=="Kausi Primary" & nest.proportion$phase == "phase 1"])
n.kausi.phase2<-sum(nest.proportion$count[nest.proportion$Forest=="Kausi Primary" & nest.proportion$phase == "phase 2"])

n.numba.phase1<-sum(nest.proportion$count[nest.proportion$Forest=="Numba Primary"& nest.proportion$phase == "phase 1"])
n.numba.phase2<-sum(nest.proportion$count[nest.proportion$Forest=="Numba Primary"& nest.proportion$phase == "phase 2"])


nest.proportion$percentage[nest.proportion$Forest=="Kausi Primary"& nest.proportion$phase == "phase 1"] <- nest.proportion$count[nest.proportion$Forest=="Kausi Primary"& nest.proportion$phase == "phase 1"]/n.kausi.phase1*100
nest.proportion$percentage[nest.proportion$Forest=="Kausi Primary"& nest.proportion$phase == "phase 2"] <- nest.proportion$count[nest.proportion$Forest=="Kausi Primary"& nest.proportion$phase == "phase 2"]/n.kausi.phase2*100
nest.proportion$percentage[nest.proportion$Forest=="Numba Primary"& nest.proportion$phase == "phase 1"] <- nest.proportion$count[nest.proportion$Forest=="Numba Primary"& nest.proportion$phase == "phase 1"]/n.numba.phase1*100
nest.proportion$percentage[nest.proportion$Forest=="Numba Primary"& nest.proportion$phase == "phase 2"] <- nest.proportion$count[nest.proportion$Forest=="Numba Primary"& nest.proportion$phase == "phase 2"]/n.numba.phase2*100

# check if sums add up
sum(nest.proportion$percentage[nest.proportion$Forest=="Kausi Primary" & nest.proportion$phase == "phase 2"])
sum(nest.proportion$percentage[nest.proportion$Forest=="Numba Primary" & nest.proportion$phase == "phase 1"])

# Define "occupied" group 
nest.proportion$AntSpCode[nest.proportion$percentage < 1] <- "Rest"

# summarize 
nest.proportion2<- nest.proportion%>% 
  group_by(Forest, AntSpCode, phase) %>% 
  summarise(percentage = sum(percentage))

# Plot it
nest.phase1.plot<-ggplot(subset(nest.proportion2, phase == "phase 1"), aes(x = Forest, y = percentage, fill=AntSpCode)) +
  geom_bar(stat = "identity", color='black') +
  xlab("") +
  labs(fill = "")+
  scale_x_discrete(labels=labs)+
  ylab("relative nest occupation [%]") +
  ggtitle("Phase 1 occupied nests") +
  scale_fill_brewer(palette="Dark2")+
  theme_minimal()
nest.phase1.plot

nest.phase2.plot<-ggplot(subset(nest.proportion2, phase == "phase 2"), aes(x = Forest, y = percentage, fill=AntSpCode)) +
  geom_bar(stat = "identity", color='black') +
  xlab("") +
  labs(fill = "")+
  scale_x_discrete(labels=labs)+
  ylab("relative nest occupation [%]") +
  ggtitle("Phase 2 occupied nests") +
  scale_fill_brewer(palette="Dark2")+
  theme_minimal()
nest.phase2.plot

##### ignoring phase and empty bamboos: Proportion of occupancy

### Plot nest occupancy as average proportion per plot
phase1.n<-subset(phase.nests, phase=="phase 1")
phase2.n<-subset(phase.nests, phase=="phase 2")

# define plot-stratum
phase.nests$plot.stratum<-paste(phase.nests$Plot, phase.nests$Stratum)

#count total nests per plot and stratum
nests.c<-count(phase.nests, plot.stratum)

# summarize occupied per plot
nest1 <- phase.nests %>%
  group_by(Forest, Plot, Stratum, plot.stratum) %>%
  summarize(occupancy = sum(occupancy))

# merge with total baits
nests.final<-merge(nest1, nests.c)

# proportion
nests.final$proportion<-nests.final$occupancy/nests.final$n*100

# plot it
nests.final.occupancy.plot<-ggplot(nests.final, aes(x=Forest, y=proportion, fill = Stratum)) +
  ggtitle("Nest occupancy") +
  ylab("Occupancy [%]")+
  xlab("")+ 
  #ylim(0,50)+
  geom_boxplot()+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  theme_bw()
nests.final.occupancy.plot

# plot nest occupancy as proportion
nest.proportion3 <- phase.nests %>%
  group_by(Forest, AntSpCode) %>%
  summarize(count = n())

# Remove empty values
nest.proportion3 <- nest.proportion3[nest.proportion3$AntSpCode != "", ]

# sum of occupancies
n.nests.kausi<-sum(nest.proportion3$count[nest.proportion3$Forest=="Kausi Primary"])
n.nests.numba<-sum(nest.proportion3$count[nest.proportion3$Forest=="Numba Primary"])

nest.proportion3$percentage<-0

nest.proportion3$percentage[nest.proportion3$Forest=="Kausi Primary"] <-nest.proportion3$count[nest.proportion3$Forest=="Kausi Primary"]/n.nests.kausi*100
nest.proportion3$percentage[nest.proportion3$Forest=="Numba Primary"] <-nest.proportion3$count[nest.proportion3$Forest=="Numba Primary"]/n.nests.numba*100

# sum check
sum(nest.proportion3$percentage[nest.proportion3$Forest=="Kausi Primary"])
sum(nest.proportion3$percentage[nest.proportion3$Forest=="Numba Primary"])

# Define "Rest" group at <5% of total incidence 
nest.proportion3$AntSpCode[nest.proportion3$percentage < 5] <- "other species"

# summarize bait incidence counts
nest.proportion3<- nest.proportion3%>% 
  group_by(Forest, AntSpCode) %>% 
  summarise(percentage = sum(percentage))

# Plot it
# relevel
nest.proportion3$AntSpCode <- as.factor(nest.proportion3$AntSpCode)
nest.proportion3$AntSpCode <- relevel(nest.proportion3$AntSpCode, "other species")

nesting.proportion<-ggplot(nest.proportion3, aes(x = Forest, y = percentage, fill = AntSpCode)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative incidence [%]") +
  ggtitle("bamboo species composition") +
  #scale_fill_brewer(palette="Dark2")+
  theme_minimal()
nesting.proportion

##### ignoring phase and empty bamboos: Proportion of occupancy per stratum per plot
phase.nests

#count total bamboos per plot per stratum
ca.c<-count(phase.nests, Plot)

counts<-rbind(phase1.c,phase2.c)

# summarize occupied per plot per phase
phase1 <- phase1.n %>%
  group_by(Forest, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase1$phase <- "phase 1"

phase2<- phase2.n %>%
  group_by(Forest, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase2$phase<-"phase 2"

# merge
phase.t<-rbind(phase1, phase2)

# merge with total nests
phase.final<-merge(phase.t, counts)

# proportion
phase.final$proportion<-phase.final$occupancy/phase.final$n*100

# plot it
nest.occupancy.phase.plot<-ggplot(phase.final, aes(x=Forest, y=proportion, fill=phase)) +
  ggtitle("Bamboo nest occupancy") +
  labs(fill = "Phase")+  
  ylab("Occupancy [%]")+
  xlab("")+ 
  ylim(0,50)+
  geom_boxplot()+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  theme_bw()
nest.occupancy.phase.plot


##### ignoring phase and empty bamboos: Proportion of abundance

# plot nest abundance as proportion
nest.abundance <- phase.nests %>%
  group_by(Forest, AntSpCode) %>%
  summarize(Abundance = sum(nesting.estimate))

# Remove empty values
nest.abundance <- nest.abundance[nest.abundance$AntSpCode != "", ]

# sum of occupancies
n.nests.kausi<-sum(nest.abundance$Abundance[nest.abundance$Forest=="Kausi Primary"])
n.nests.numba<-sum(nest.abundance$Abundance[nest.abundance$Forest=="Numba Primary"])

nest.abundance$percentage<-0

nest.abundance$percentage[nest.abundance$Forest=="Kausi Primary"] <-nest.abundance$Abundance[nest.abundance$Forest=="Kausi Primary"]/n.nests.kausi*100
nest.abundance$percentage[nest.abundance$Forest=="Numba Primary"] <-nest.abundance$Abundance[nest.abundance$Forest=="Numba Primary"]/n.nests.numba*100

# sum check
sum(nest.abundance$percentage[nest.abundance$Forest=="Kausi Primary"])
sum(nest.abundance$percentage[nest.abundance$Forest=="Numba Primary"])

# Define "Rest" group at <5% of total incidence 
nest.abundance$AntSpCode[nest.abundance$percentage < 5] <- "other species"

# summarize bait incidence counts
nest.abundance<- nest.abundance%>% 
  group_by(Forest, AntSpCode) %>% 
  summarise(percentage = sum(percentage))

# Plot it
# relevel
nest.abundance$AntSpCode <- as.factor(nest.abundance$AntSpCode)
nest.abundance$AntSpCode <- relevel(nest.abundance$AntSpCode, "other species")

abundance.proportion<-ggplot(nest.abundance, aes(x = Forest, y = percentage, fill = AntSpCode)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative abundance [%]") +
  ggtitle("bamboo species composition") +
  #scale_fill_brewer(palette="Dark2")+
  theme_minimal()
abundance.proportion

# plot as glm prediction against elevation

ggplot(bamboo.incidence.e, aes(x=elevation, y=occupancy)) + geom_point() + 
  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)

#----------------------------------------------------------#
# 2.2 Nest diversity -----
#----------------------------------------------------------#

# remove species from unoccupied bamboos (i.e., foragers)
nesters<-phase.nests
nesters$AntSpCode[nesters$occupancy==0]  <- ""

# define plot.stratum
nesters$plot.stratum<-paste(nesters$Plot,nesters$Stratum)

# subset phase
phase1<-subset(nesters, phase=="phase 1")
phase2<-subset(nesters, phase=="phase 2")

# phase 1 plot.stratum x species incidence
ant_n.phase1<-dcast(phase1, formula = plot.stratum ~ AntSpCode, length)
rownames(ant_n.phase1)<-ant_n.phase1[,1] # plot.stratum as rownames
ant_n.phase1<-ant_n.phase1[,-c(1,2)]

# phase 2 plot.stratum x species incidence
ant_n.phase2<-dcast(phase2, formula = plot.stratum ~ AntSpCode, length)
rownames(ant_n.phase2)<-ant_n.phase2[,1] # plot.stratum as rownames
ant_n.phase2<-ant_n.phase2[,-c(1,2)]

# extract diversity values for phase 1
data2 <- data.frame(row.names(ant_n.phase1))     # create a new file with plot numbers
names(data2) <- "plot.stratum" # rename variable "plot.stratum"
data2$Shannon <- diversity(ant_n.phase1, index = "shannon", base = exp(1)) # Shannon diversity per plot
data2$Richness <- specnumber(ant_n.phase1)                          # Richness per plot
data2$expH <- exp(data2$Shannon)                                         # exponential Shannon diversity per plot
data2$expH[data2$Richness == 0] <- 0                                  # define exp(Shannon) = 0
data2$evenness <- data2$Shannon/log(specnumber(ant_n.phase1)) # evenness per plot
str(data2)

# add Plot location
plot_m3<-subset(nesters, phase=="phase 1")
plot_m3<-distinct(plot_m3, plot.stratum, .keep_all = TRUE)
phase1.diversity<-merge(data2, plot_m3)

# extract diversity values for phase 2
data2 <- data.frame(row.names(ant_n.phase2))     # create a new file with plot numbers
names(data2) <- "plot.stratum" # rename variable "plot.stratum"
data2$Shannon <- diversity(ant_n.phase2, index = "shannon", base = exp(1)) # Shannon diversity per plot
data2$Richness <- specnumber(ant_n.phase2)                          # Richness per plot
data2$expH <- exp(data2$Shannon)                                         # exponential Shannon diversity per plot
data2$expH[data2$Richness == 0] <- 0                                     # define exp(Shannon) = 0
data2$evenness <- data2$Shannon/log(specnumber(ant_n.phase2)) # evenness per plot
str(data2)

# add Plot location
plot_m4<-subset(nesters, phase=="phase 2")
plot_m4<-distinct(plot_m4, plot.stratum, .keep_all = TRUE)
phase2.diversity<-merge(data2, plot_m4)

## gg plot it
nests.phase1<-ggplot(phase1.diversity, aes(x=Forest, y=expH, fill=Stratum)) +
  ggtitle("Phase 1 Bamboo nester diversity") +
  geom_boxplot()+
  ylim(0,5)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Species diversity")+
  xlab("")+ 
  theme(axis.text.x = element_text(size=12))
nests.phase1

nests.phase2 <-ggplot(phase2.diversity, aes(x=Forest, y=expH, fill=Stratum)) +
  ggtitle("Phase 2 Bamboo nester diversity") +
  geom_boxplot()+
  #ylim(0,5)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Species diversity")+
  xlab("")+ 
  theme(axis.text.x = element_text(size=12))
nests.phase2
#
diversity.nester <-rbind(phase1.diversity, phase2.diversity)

# Plot nester diversity, phase dependence
nests.diversity.phase.plot<-ggplot(diversity.nester, aes(x=Forest, y=expH, fill=phase)) +
  ggtitle("Bamboo nester diversity") +
  labs(fill = "Phase")+  
  ylab("Species diversity [expH]")+
  xlab("")+ 
  geom_boxplot()+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  theme_bw()
nests.diversity.phase.plot

#----------------------------------------------------------#
# 3.1 Statistic: Environment  -----
#----------------------------------------------------------#

# test for overdispersion
overdisp_fun<-function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Plot attribute transformations: Are environmental factors skewed?
hist(tree.meta$trunk) #right skewed, log+1 transformation fine
hist(log(tree.meta$trunk+1))

tree.meta$trunk.log<-log(tree.meta$trunk+1)
#
hist(tree.meta$Lianas.n)#right skewed, log+1 transformation fine
hist(log(tree.meta$Lianas.n+1))
tree.meta$Lianas.n.log<-log(tree.meta$Lianas.n+1)
#
hist(tree.meta$height)#right skewed,  sqrt transformation
hist(sqrt(tree.meta$height))
tree.meta$height.sqrt<-sqrt(tree.meta$height)
#  
hist(tree.meta$Caco)#looks ok
#
hist(tree.meta$dw.percent)#right skewed, log+1 transformation
hist(log(tree.meta$dw.percent+1))
tree.meta$dw.percent.log<-log(tree.meta$dw.percent+1)
#
hist(tree.meta$dw.number)#right skewed, log+1 transformation
hist(log(tree.meta$dw.number+1))
tree.meta$dw.number.log<-log(tree.meta$dw.number+1) 
#
hist(tree.meta$slope.var)#slightly right skewed, sqrt transformation
hist(sqrt(tree.meta$slope.var))
tree.meta$slope.var.log<-log(tree.meta$slope.var+1) 

# Same for plot means
hist(plot.meta$Caco)
plot.meta$Caco.log<-log(plot.meta$Caco+1)
hist(plot.meta$Caco.log)

# this one is ok
hist(plot.meta$slope.var)

# skewed
hist(plot.meta$trunk_mean)
plot.meta$trunk_mean.log<-log(plot.meta$trunk_mean+1)
hist(plot.meta$trunk_mean.log)


hist(plot.meta$Lianas.n_mean)
plot.meta$Lianas.n_mean.log<-log(plot.meta$Lianas.n_mean+1)
hist(plot.meta$Lianas.n_mean.log)


hist(plot.meta$dw.percent_mean)
plot.meta$dw.percent_mean.log<-log(plot.meta$dw.percent_mean+1)
hist(plot.meta$dw.percent_mean.log)

hist(plot.meta$dw.number_mean)
plot.meta$dw.number_mean.log<-log(plot.meta$dw.number_mean+1)
hist(plot.meta$dw.number_mean.log)

hist(plot.meta$height_mean)
plot.meta$height_mean.log<-log(plot.meta$height_mean+1)
hist(plot.meta$height_mean.log)

# for plot.meta2 same
plot.meta2$Caco.log<-log(plot.meta2$Caco+1)
plot.meta2$trunk_mean.log<-log(plot.meta2$trunk_mean+1)
plot.meta2$Lianas.n_mean.log<-log(plot.meta2$Lianas.n_mean+1)
plot.meta2$dw.percent_mean.log<-log(plot.meta2$dw.percent_mean+1)
plot.meta2$dw.number_mean.log<-log(plot.meta2$dw.number_mean+1)
plot.meta2$height_mean.log<-log(plot.meta2$height_mean+1)

# Are any plot attributes different between elevations?

# Are canopy covers different between forests?
caco_model <- glmmTMB(Caco.log~  Forest +(1|Block), data = plot.meta, family = gaussian)
summary(caco_model) #ns

# Is slope different?
slope_model <- glmmTMB(slope.var~  Forest +(1|Block), data = plot.meta, family = gaussian)
summary(caco_model) # ns

# Is dbh different between forests?
model_dbh<- glmmTMB(trunk.log ~  Forest + (1|Block/Plot), data = tree.meta, family = gaussian)
summary(model_dbh)# ns

# Is number of deadwood pieces different?
model_dw.n<- glmmTMB(dw.number.log ~  Forest + (1|Block/Plot), data= tree.meta, family = gaussian)
summary(model_dw.n)# ns

# is % deadwood cover different?
model_dw.p<- glmmTMB(dw.percent.log ~  Forest + (1|Block/Plot), data= tree.meta, family = gaussian)
summary(model_dw.p)# ns

# Is number of lianas different?
model_lianas<- glmmTMB(Lianas.n.log ~  Forest + (1|Block/Plot), data= tree.meta, family = gaussian)
summary(model_lianas)
# =>  more lianas in higher elevation

# bamboo/bait height different?
model_height<- glmmTMB(height.sqrt ~  Forest + (1|Block/Plot), data= tree.meta, family = gaussian)
summary(model_height) # ns

# correlation plot of environmental variables
names(tree.meta)
mydata.cor <- cor(tree.meta[, c(12,13, 15:20)], method = c("spearman"), use = "pairwise.complete.obs")
corrplot(mydata.cor,
         method = "color", type = "lower", order = "AOE", diag = F,
         tl.col = "black", outline = T, addCoef.col = "black", number.cex = 0.8,
         tl.cex = 1.1, cl.cex = 0.9
)

#----------------------------------------------------------#
# 3.2 Bait Statistics -----
#----------------------------------------------------------#

# Are there more species per bait in higher elevation?

# get environmental data
bait.double.e<-merge(bait.double, tree.meta)
bait.double.e<-subset(bait.double.e, species.per.bait!=0) # remove empty baits
bait.double.e$species.per.bait[bait.double.e$species.per.bait == 1] <- 0 # if there was one species on a bait, make it zero
bait.double.e$species.per.bait[bait.double.e$species.per.bait == 2] <- 1 # if there were 2, make it 1

#
bait.double.model1 <- glmmTMB(species.per.bait~Forest*Stratum+(1|Block/Plot),
                              data=bait.double.e, 
                              family = binomial)

bait.double.model2 <- glmmTMB(species.per.bait~Forest+Stratum+(1|Block/Plot), 
                              data=bait.double.e,
                              family = binomial)
# 
anova(bait.double.model1, bait.double.model2) # ns, better without interactions
summary(bait.double.model2)
overdisp_fun(bait.double.model2)
#
testDispersion(bait.double.model2) # ok
simulateResiduals(bait.double.model2, plot = T) # good
testZeroInflation(simulateResiduals(fittedModel = bait.double.model2)) # ok

# with environmental factors
bait.double.model.e1 <- glmmTMB(species.per.bait~Forest+Lianas.n.log+Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var+(1|Block/Plot),
                                 data=bait.double.e, family= binomial)
bait.double.model.e2 <- glmmTMB(species.per.bait~Forest*Lianas.n.log*Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var+(1|Block/Plot),
                                data=bait.double.e, family= binomial)
anova(bait.double.model.e1, bait.double.model.e2) # ns, without interaction

summary(bait.double.model.e1)
overdisp_fun(bait.double.model.e1)
#
testDispersion(bait.double.model.e1) # ok
simulateResiduals(bait.double.model.e1, plot = T) # good
testZeroInflation(simulateResiduals(fittedModel = bait.double.model.e1)) # ok

## Bait occupancy
# get environmental data
baiting.incidence.e<-merge(baiting.incidence, tree.meta)

# Binominal regression
baitoccupancy.model1 <- glmmTMB(occupancy~Forest*Stratum+(1|Block/Plot),
                               data=baiting.incidence.e,
                               family=binomial)
baitoccupancy.model2 <- glmmTMB(occupancy~Forest+Stratum+(1|Block/Plot),
                                data=baiting.incidence.e,
                                family=binomial)
anova(baitoccupancy.model1, baitoccupancy.model2) # interaction model better

summary(baitoccupancy.model1)
#
overdisp_fun(baitoccupancy.model1)
testDispersion(baitoccupancy.model1) # ok
simulateResiduals(fittedModel = baitoccupancy.model1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model1)) # ok
plot(allEffects(baitoccupancy.model1)) # model visualization

# with environmental factors
baitoccupancy.model.e1 <- glmmTMB(occupancy~Forest*Lianas.n.log*Stratum+dw.percent.log+dw.number.log*trunk.log+Caco+slope.var.log+(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=binomial)
baitoccupancy.model.e2 <- glmmTMB(occupancy~Forest+Lianas.n.log+Stratum+dw.percent.log+dw.number.log*trunk.log+Caco+slope.var.log+(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=binomial)

anova(baitoccupancy.model.e1, baitoccupancy.model.e2) # interaction model better

summary(baitoccupancy.model.e1)
overdisp_fun(baitoccupancy.model.e1)
#
testDispersion(baitoccupancy.model.e1) # ok
simulateResiduals(baitoccupancy.model.e1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model.e1)) # ok
plot(allEffects(baitoccupancy.model.e1)) # model visualization

# with abundance 
baitabundance.model1 <- glmmTMB(Abundance~Forest+Stratum+(1|Block/Plot),
                                zi=~Forest*Stratum,
                                data=baiting.incidence.e,
                                family=nbinom1)
baitabundance.model2 <- glmmTMB(Abundance~Forest*Stratum+(1|Block/Plot),
                                zi=~Forest*Stratum,
                                data=baiting.incidence.e,
                                family=nbinom1)
anova(baitabundance.model1, baitabundance.model2) # interaction better

summary(baitabundance.model2)
overdisp_fun(baitabundance.model2)
#
testDispersion(baitabundance.model2) # ok
simulateResiduals(baitabundance.model2, plot = T) # ok-ish
testZeroInflation(simulateResiduals(fittedModel = baitabundance.model2)) # ok 
plot(allEffects(baitabundance.model2)) # model visualization

# 

# with abundance and environmental
baitabundance.model.e1 <- glmmTMB(Abundance~Forest+Lianas.n.log+Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+(1|Block/Plot),
                                  zi=~Forest*Stratum,
                                  data=baiting.incidence.e,
                                  family=nbinom1)
baitabundance.model.e2 <- glmmTMB(Abundance~Forest*Lianas.n.log*Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+(1|Block/Plot),
                                  zi=~Forest*Stratum,
                                  data=baiting.incidence.e,
                                  family=nbinom1)
anova(baitabundance.model.e1, baitabundance.model.e2) # better without interaction

summary(baitabundance.model.e1)
overdisp_fun(baitabundance.model.e1)
#
testDispersion(baitabundance.model.e1) # ok
simulateResiduals(baitabundance.model.e1, plot = T) #  maybe ok
testZeroInflation(simulateResiduals(fittedModel = baitabundance.model.e1)) # ok

#

#####  Baiting diversity 
# NOTE: Here we test plot-level diversity, which aggregates strata

# get environmental data: Plot.meta has averages for variables on Plot-level
baiting.diversity.e<-merge(diversity, plot.meta)

# Fit a LMM for shannon diversity
baitdiversity.model <- glmmTMB(expH ~ Forest + (1|Block), data = baiting.diversity.e, family = gaussian)

summary(baitdiversity.model)
#
testDispersion(baitdiversity.model) # ok
simulateResiduals(baitdiversity.model, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model)) # ok

# add environmental factors
baitdiversity.model.e1 <- glmmTMB(expH ~ Forest+Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log + (1|Block),
                              data = baiting.diversity.e, family = gaussian)
baitdiversity.model.e2 <- glmmTMB(expH ~ Forest*Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log + (1|Block),
                                  data = baiting.diversity.e, family = gaussian)
anova(baitdiversity.model.e1, baitdiversity.model.e2) # better without interaction

summary(baitdiversity.model.e1)
overdisp_fun(baitdiversity.model.e1) # 
#
testDispersion(baitdiversity.model.e1) # ok
simulateResiduals(baitdiversity.model.e1, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.e1)) # ok

# evenness (will maybe be removed?)
baitdiversity.model.eve <- glmmTMB(evenness ~ Forest + (1|Block), data = baiting.diversity.e, family = gaussian)
summary(baitdiversity.model.eve)
overdisp_fun(baitdiversity.model.eve)
#
testDispersion(baitdiversity.model.eve) # ok
simulateResiduals(baitdiversity.model.eve, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.eve)) # ok

# add environmental factors
baitdiversity.model.eve.e1 <- glmmTMB(evenness ~Forest+slope.var+Caco.log+trunk_mean.log+Lianas.n_mean.log+dw.number_mean.log+dw.percent_mean.log + (1|Block),
                                  data = baiting.diversity.e, family = gaussian)
baitdiversity.model.eve.e2 <- glmmTMB(evenness~Forest*Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log+ (1|Block),
                                  data = baiting.diversity.e, family = gaussian)
anova(baitdiversity.model.eve.e1, baitdiversity.model.eve.e2) # ns, no interaction

summary(baitdiversity.model.eve.e1)
overdisp_fun(baitdiversity.model.eve.e1)
#
testDispersion(baitdiversity.model.eve.e1) # ok
simulateResiduals(baitdiversity.model.eve.e1, plot = T) # ok-ish
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.eve.e1)) # ok

### Bait diversity incl. stratum
# NOTE: Here, we look at stratum-level diversity, i.e. twice for each plot (understory+canopy)

# get environmental data: plot.meta2 has averages for variables on Stratum-level (2 strata per plot)
diversity.stratum.e<-merge(diversity.stratum, plot.meta2, by.x = 'plot.stratum', by.y = 'plot.stratum', all.x = T)

# Fit a LMM for Shannon diversity. Needs also plot as random factor since strata are not independent
baitdiversity.stratum.model1 <- glmmTMB(expH ~ Forest.x + Stratum+ (1|Block.x/Plot.x), data = diversity.stratum.e, family = gaussian)
baitdiversity.stratum.model2 <- glmmTMB(expH ~ Forest.x * Stratum+ (1|Block.x/Plot.x), data = diversity.stratum.e, family = gaussian)

anova(baitdiversity.stratum.model1, baitdiversity.stratum.model2) # ns, both are similar, probably best w/o interaction

summary(baitdiversity.stratum.model1)
overdisp_fun(baitdiversity.stratum.model1)  
#
testDispersion(baitdiversity.stratum.model1) # ok
simulateResiduals(baitdiversity.stratum.model1, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.stratum.model1)) # ok

# add environmental factors
baitdiversity.stratum.model.e1 <- glmmTMB(expH ~Forest.x+Stratum+slope.var+Caco.log+trunk_mean.log+Lianas.n_mean.log+dw.number_mean.log+dw.percent_mean.log + (1|Block.x/Plot.x),
                              data = diversity.stratum.e, family = gaussian)
baitdiversity.stratum.model.e2 <- glmmTMB(expH ~Forest.x*Stratum*Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log + (1|Block.x/Plot.x),
                                         data = diversity.stratum.e, family = gaussian)
anova(baitdiversity.stratum.model.e1, baitdiversity.stratum.model.e2) # interaction is better

summary(baitdiversity.stratum.model.e2)
overdisp_fun(baitdiversity.stratum.model.e2) #
#
testDispersion(baitdiversity.stratum.model.e2) # ok
simulateResiduals(baitdiversity.stratum.model.e2, plot = T) # ok-ish
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.stratum.model.e2)) # ok


#----------------------------------------------------------#
# 3.3 Bamboo Nest Statistics -----
#----------------------------------------------------------#

# bamboo occupancy models

# Get environmental metadata
tree<-tree.meta %>%
  dplyr::select(Code:height, Caco:slope.var.log)
#
bamboo.incidence.e<-merge(phase.nests, tree, by = "Code")
bamboo.incidence.e$new.Treatment<-as.factor(bamboo.incidence.e$new.Treatment)


# binominal model of bamboo nesting using phase
bamboo.occupancy.model1 <- glmmTMB(occupancy~Forest*Stratum+new.Treatment+(1|Block/Plot),
                                  data=phase.nests,
                                  family=binomial)
bamboo.occupancy.model2 <- glmmTMB(occupancy~Forest+Stratum+new.Treatment+(1|Block/Plot),
                                   data=phase.nests,
                                   family=binomial)
anova(bamboo.occupancy.model2, bamboo.occupancy.model1)
# ns, no interaction

summary(bamboo.occupancy.model2)
overdisp_fun(bamboo.occupancy.model2)
#
testDispersion(bamboo.occupancy.model2) # ok
simulateResiduals(bamboo.occupancy.model2, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboo.occupancy.model2)) # ok 

# binominal model of bamboo nesting with environment
bamboooccupancy.model.e1 <- glmmTMB(occupancy~Forest*Lianas.n.log*Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+new.Treatment+(1|Block/Plot),
                                   data=bamboo.incidence.e,
                                   family=binomial)
bamboooccupancy.model.e2 <- glmmTMB(occupancy~Forest+Lianas.n.log+Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+new.Treatment+(1|Block/Plot),
                                    data=bamboo.incidence.e,
                                    family=binomial)
anova(bamboooccupancy.model.e1, bamboooccupancy.model.e2) # ns -  better using no interactions

summary(bamboooccupancy.model.e2)
overdisp_fun(bamboooccupancy.model.e2)
#
testDispersion(bamboooccupancy.model.e2) # ok
simulateResiduals(bamboooccupancy.model.e2, plot = T) # all good
testZeroInflation(simulateResiduals(fittedModel = bamboooccupancy.model.e2)) # ok 

# with abundance 
bamboo.abundance.model3 <- glmmTMB(nesting.estimate~Forest*Stratum+new.Treatment+(1|Block/Plot),
                                   data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                   family=nbinom2)
bamboo.abundance.model4 <- glmmTMB(nesting.estimate~Forest+Stratum+new.Treatment+(1|Block/Plot),
                                   data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                   family=nbinom2)
anova(bamboo.abundance.model3,bamboo.abundance.model4) # ns, no interaction
summary(bamboo.abundance.model4)
overdisp_fun(bamboo.abundance.model4)
#
testDispersion(bamboo.abundance.model4) # ok
simulateResiduals(bamboo.abundance.model4, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboo.abundance.model4)) # ok 
#  residuals better in nbinom2

# with abundance and environmental
bamboo.abundance.model.e4 <- glmmTMB(nesting.estimate~Forest+Lianas.n.log+Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+new.Treatment+(1|Block/Plot),
                                     data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                     family=nbinom2) #convergence troubles
bamboo.abundance.model.e5 <- glmmTMB(nesting.estimate~Forest*Lianas.n.log*Stratum+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+new.Treatment+(1|Block/Plot),
                                     data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                     family=nbinom2)  #convergence troubles
anova(bamboo.abundance.model.e4, bamboo.abundance.model.e5) 

summary(bamboo.abundance.model.e4)
overdisp_fun(bamboo.abundance.model.e4)
#
testDispersion(bamboo.abundance.model.e4) # 
simulateResiduals(bamboo.abundance.model.e4, plot = T) # 
testZeroInflation(simulateResiduals(fittedModel = bamboo.abundance.model.e4)) # 

#####   Bamboo Species diversity
# get environment
plot<- plot.meta2%>% 
  dplyr::select(plot.stratum:height_mean, Caco:height_mean.log)

diversity.nester.e<-merge(diversity.nester,plot, by = "plot.stratum", all.x = TRUE)
#
bamboo.diversity.model1 <- glmmTMB(expH ~ Forest+new.Treatment+Stratum+(1|Block/Plot), zi=~Forest*new.Treatment, data = diversity.nester.e, family=gaussian)
bamboo.diversity.model2 <- glmmTMB(expH ~ Forest*Stratum*new.Treatment+(1|Block/Plot), zi=~Forest*Stratum, data = diversity.nester.e, family=gaussian)

anova(bamboo.diversity.model2,bamboo.diversity.model1) # ns, no interaction
#
summary(bamboo.diversity.model1)

overdisp_fun(bamboo.diversity.model1)#ok
#
testDispersion(bamboo.diversity.model1) # ok
simulateResiduals(bamboo.diversity.model1, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model1)) # slightly zero inflated (p=0.048)
plot(allEffects(bamboo.diversity.model1)) # model visualization

# bamboo diversity changes with environment metadata
bamboo.diversity.model.e2 <- glmmTMB(expH ~Forest*Stratum*Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log+new.Treatment+(1|Block/Plot),
                                     zi=~Forest,
                                     data = diversity.nester.e, family=gaussian)

bamboo.diversity.model.e1 <- glmmTMB(expH ~ Forest+Stratum+Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log+new.Treatment+(1|Block/Plot), 
                                     zi=~Forest,
                                     data = diversity.nester.e, family=gaussian) # zi gives convergence troubles, not included

bamboo.diversity.model.e1 <- glmmTMB(expH ~Forest*Stratum*Lianas.n_mean.log+ trunk_mean.log+dw.number_mean.log+dw.percent_mean.log+new.Treatment+(1|Block/Plot),
                                     zi=~Forest,
                                     data = diversity.nester.e, family=gaussian) # alternative without slope and caco?

anova(bamboo.diversity.model.e1, bamboo.diversity.model.e2) # 

bamboo.diversity.model3 <- glmmTMB((1+expH) ~ Forest+new.Treatment+Stratum+Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log+(1|Block/Plot), data = diversity.nester.e, family=gaussian(link=log))
bamboo.diversity.model4 <- glmmTMB((1+expH) ~ Forest*new.Treatment*Stratum*Lianas.n_mean.log+slope.var+Caco.log+trunk_mean.log+dw.number_mean.log+dw.percent_mean.log+(1|Block/Plot), data = diversity.nester.e, family=gaussian(link=log))
anova(bamboo.diversity.model3, bamboo.diversity.model4) # ns, no interaction better

summary(bamboo.diversity.model3)
overdisp_fun(bamboo.diversity.model.e2) 
#
testDispersion(bamboo.diversity.model3) # ok
simulateResiduals(bamboo.diversity.model3, plot = T) # weird 
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model3)) # slightly zero inflated (p=0.024)

plot(allEffects(bamboo.diversity.model3)) # model visualization


#----------------------------------------------------------#
# 3.4 Summary Statistics -----
#----------------------------------------------------------#

# This part is just a selection of the most important models and the questions they address, and their interpretation

# 1) Does the number of species per bait change?
summary(bait.double.model2)

# 2) Does bait occupancy change?
summary(baitoccupancy.model1)
summary(baitoccupancy.model.e2)

# - higher baits -> higher occupancy
# - interaction: lower elevation has higher occupancy in canopy
# bigger trees, higher chance of occupancy
# interaction Lianas.n.log:StratumUN: higher chance of occupancy on understory trees w. many lianas, negative effects in the canopy in kausi but no effects in canopy in numba
# more lianas, lower chance of occupancy?

# 3) Does bait abundance change?
summary(baitabundance.model2)
summary(baitabundance.model.e1)

# - lower baits, higher abundance 
# - lower elevation, higher abundance

# 4) Does bait diversity change?
summary(baitdiversity.stratum.model1)
summary(baitdiversity.stratum.model.e1)

# - higher elevation, higher diversity
# - higher diversity in understory
 
# 5) Does bamboo occupancy change?
summary(bamboo.occupancy.model2)
summary(bamboooccupancy.model.e1)

# - higher elevation lower occupancy
# - higher occupancy in canopy
# - first phase treatment slightly lower occupancy

# 6) Does bamboo diversity change?
summary(bamboo.diversity.model1)
summary(bamboo.diversity.model.e2)
# - higher elevation, lower nesting diversity
# - first phase lower diversity
# - bigger trees, higher diversity


# abundance of ants inside bamboo is not as interesting as bait abundance, so I would not present it in the main manuscript
# evenness is probably not so interesting since its partly accounted for in Shannon diversity (expH)

#----------------------------------------------------------#
# Model summary export -----
#----------------------------------------------------------#

# Create list of models
models_list <- list()
models_list[[1]] <- bait.double.model2
models_list[[2]] <- baitoccupancy.model1
models_list[[3]] <- baitoccupancy.model.e2
models_list[[4]] <- baitdiversity.stratum.model1
models_list[[5]] <- baitdiversity.model.eve.e1
models_list[[6]] <- baitabundance.model2
models_list[[7]] <- baitabundance.model.e1
models_list[[8]] <- bamboo.occupancy.model2
models_list[[9]] <- bamboooccupancy.model.e1
models_list[[10]] <- bamboo.diversity.model1
models_list[[11]] <- bamboo.diversity.model.e2

# Create Excel workbook
wb <- createWorkbook()

# Loop through models and generate model summaries
for (i in 1:length(models_list)) {
  # Generate model summary
  model_summary <- broom.mixed::tidy(models_list[[i]], effects = "fixed")
  
  # Extract model formula and add to header
  formula_text <- as.character(formula(models_list[[i]]))
  header <- c(paste0("Model ", i, ": ", formula_text), rep("", ncol(model_summary)-1))
  model_summary <- rbind(header, model_summary)
  
  # Write model summary to Excel sheet
  addWorksheet(wb, sheetName = paste0("Model ", i))
  writeData(wb, sheet = i, x = model_summary)
}

# Save Excel file
saveWorkbook(wb, "model_summaries.xlsx", overwrite = TRUE)

#----------------------------------------------------------#
# Figure summaries -----
#----------------------------------------------------------#

#Bamboo diversity
nests.diversity.phase.plot

#Bait diversity
baitdiversity.stratum

#bamboo occupied
nest.occupancy.phase.plot

#Bamboo 
nesting.proportion

# Bait composition
# as proportion of incidence
bait.incidence.plot
# as proportion of abundance
abundance.plot


figure_baitsbamboos <- ggarrange(nests.diversity.phase.plot, baitdiversity.stratum, nesting.proportion, abundance.plot,
                               labels = c("A", "B", "C", "D"),
                               ncol = 2, nrow = 2, common.legend = F
)
figure_baitsbamboos
