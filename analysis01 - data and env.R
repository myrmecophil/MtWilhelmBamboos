## Numba-Kausi Baits & nests - R script by Phil, April 2024


# This part of the script does the raw data conversion and saves it as tables that are being used in the other scripts.
# Further, it includes the statistical analysis of the environmental data.

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
baiting.raw<-read.csv(file="raw data/BaitsData.csv", header=T)


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
baiting.incidence<-baiting.data%>%
  group_by(Block, Plot, Subplot, Stratum, Code, Forest, Season) %>%
  summarize(Abundance = sum(Abundance))

# make new variable which counts bait occupancy as binary variable
baiting.incidence<- baiting.incidence %>%
  mutate(occupancy = case_when(Abundance == '0' ~ 0,
                               Abundance != '0' ~ 1)) # if ant was is there it counts it as 1, if not 0
baiting.incidence<-as.data.frame(baiting.incidence)

# export
write.csv(baiting.incidence, "baiting.incidence.csv", row.names=FALSE)
write.csv(baiting.data, "baiting.data.csv", row.names=FALSE)

### Raw nesting data

# Import bamboo nesting data
nest.raw =read.csv(file="raw data/NestsData.csv", header=T)

# split column 'moved to' into two separate columns
nest.raw$Moved.to.block <- str_split(nest.raw$Moved.to.code, "-", simplify = TRUE)[,1]
nest.raw$Moved.to.plot <- str_split(nest.raw$Moved.to.code, "-", simplify = TRUE)[,1:2] %>% apply(1, paste, collapse = "-")
nest.raw$Moved.to.plot<-as.factor(nest.raw$Moved.to.plot)
nest.raw$Moved.to.block<-as.factor(nest.raw$Moved.to.block)

nest.raw2<-nest.raw

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

# export
write.csv(phase.nests, "phase.nests.csv", row.names=FALSE)
write.csv(nest.raw, "nest.raw.csv", row.names=FALSE)

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

# Plot attribute transformations: Are environmental factors skewed? 
# Most of them are, so for models we will use log-transformation and normalization

# tree trunk circumference 
hist(tree.meta$trunk) 
hist(scale(log(tree.meta$trunk+1)))

# number of lianas
hist(tree.meta$Lianas.n)
hist(scale(tree.meta$Lianas.n))
hist(scale(log(tree.meta$Lianas.n+1)))

# height of the hanged bamboo. Note this is not used in the models since its equivalent to stratum
hist(tree.meta$height)#
hist(scale(log(tree.meta$height+1)))

# deadwood percentage on the ground in 1 m circumference around the tree base 
hist(tree.meta$dw.percent)
hist(scale(log(tree.meta$dw.percent+1)))

# deadwood pieces on the ground in 1 m circumference around the tree base 
hist(tree.meta$dw.number)
hist(scale(log(tree.meta$dw.number+1)))

## Same for plot variables 

# canopy cover index
hist(plot.meta$Caco)
hist(scale(log(plot.meta$Caco+1)))

# slope inclination angle
hist(plot.meta$slope.var)
hist(scale(log(plot.meta$slope.var+1)))

# plot means of tree variables
hist(plot.meta$trunk_mean)
hist(scale(log(plot.meta$trunk_mean+1)))
#
hist(plot.meta$Lianas.n_mean)
hist(scale(log(plot.meta$Lianas.n_mean+1)))
#
hist(plot.meta$dw.percent_mean)
hist(scale(log(plot.meta$dw.percent_mean+1)))
#
hist(plot.meta$dw.number_mean)
hist(scale(log(plot.meta$dw.number_mean+1)))
#
hist(plot.meta$height_mean)
hist(scale(log(plot.meta$height_mean+1)))

# export
write.csv(plot.meta2, "plot.meta2.csv", row.names=FALSE)
write.csv(plot.meta, "plot.meta.csv", row.names=FALSE)
write.csv(tree.meta, "tree.meta.csv", row.names=FALSE)

#----------------------------------------------------------#
# 3.1 Statistic: Environment  -----
#----------------------------------------------------------#

### 1. Are any plot attributes different between elevations?

# Are canopy cover different between forests? - no
caco_model <- glmmTMB(scale(log(Caco+1))~  Forest +(1|Block), data = plot.meta, family = gaussian)
summary(caco_model) #ns
#
testDispersion(caco_model) # ok
simulateResiduals(caco_model, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = caco_model)) # ok

# Is slope different?
slope_model <- glmmTMB(scale(log(slope.var+1))~  Forest +(1|Block), data = plot.meta, family = gaussian)
summary(slope_model) # ns
#
testDispersion(slope_model) # ok
simulateResiduals(slope_model, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = slope_model)) # ok

### 2. Are any tree attributes different between elevations? - none of them are

# make tree level data with removed tree duplicates
tree.data<-merge(tree.meta.sub.ID, baiting.incidence)

# Is dbh different between forests?
model_dbh<- glmmTMB(scale(log(trunk+1)) ~  Forest + (1|Block/Plot), data = tree.data, family = gaussian)
summary(model_dbh)# ns
#
testDispersion(model_dbh) # ok
simulateResiduals(model_dbh, plot = T) # slight deviation but ok
testZeroInflation(simulateResiduals(fittedModel = model_dbh)) # ok

# Is number of deadwood pieces different?
model_dw.n<- glmmTMB(scale(log(dw.number+1)) ~ Forest + (1|Block/Plot), data= tree.data, family = gaussian)
summary(model_dw.n)# ns
#
testDispersion(model_dw.n) # ok
simulateResiduals(model_dw.n, plot = T) # some troubles
testZeroInflation(simulateResiduals(fittedModel = model_dw.n)) # ok

# is % deadwood cover different?
model_dw.p<- glmmTMB(scale(log(dw.percent+1)) ~  Forest + (1|Block/Plot), data= tree.data, family = gaussian)
summary(model_dw.p)# ns
#
testDispersion(model_dw.p) # ok
simulateResiduals(model_dw.p, plot = T) # lots of troubles
testZeroInflation(simulateResiduals(fittedModel = model_dw.p)) # ok

# Is number of lianas different?
model_lianas1<- glmmTMB(Lianas.n ~  Forest+(1|Block/Plot), data= tree.data, family = nbinom2())
summary(model_lianas1) # ns
#
testDispersion(model_lianas1) # ok
simulateResiduals(model_lianas1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = model_lianas1)) # ok

# correlation plot of environmental variables. Note: does not include tree duplicates
names(tree.meta)
# remove duplicate tree IDs
tree.meta.sub2<-distinct(tree.meta, tree.ID, .keep_all = TRUE)

mydata.cor <- cor(tree.meta.sub2[, c(4, 6:9,12:14)], method = c("spearman"), use = "pairwise.complete.obs")
corrplot(mydata.cor,
         method = "color", type = "lower", order = "AOE", diag = F,
         tl.col = "black", outline = T, addCoef.col = "black", number.cex = 0.8,
         tl.cex = 1.1, cl.cex = 0.9
)

# plot lianas (without duplicates)
labs <- expression("lowland", "midelevation")

lianas.plot<-ggplot(tree.data, aes(x=Forest, y=log(Lianas.n+1), fill = Stratum)) +
  ggtitle("Number of Lianas") +
  ylab("ln(n+1)")+
  xlab("")+ 
  #ylim(0,4)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
lianas.plot

# plot tree circumference
trunk.plot<-ggplot(tree.data, aes(x=Forest, y=log(trunk+1), fill = Stratum)) +
  ggtitle("Tree circumference [ln(n+1)]") +
  ylab("")+
  xlab("")+ 
  #ylim(0,4)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
trunk.plot

# plot %deadwood
dw.percent.plot<-ggplot(tree.data, aes(x=Forest, y=dw.percent, fill = Stratum)) +
  ggtitle("Deadwood on ground (%)") +
  ylab("")+
  xlab("")+ 
  ylim(0,100)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
lianas.plot

# plot deadwood pieces
dw.number.plot<-ggplot(tree.data, aes(x=Forest, y=log(dw.number+1), fill = Stratum)) +
  ggtitle("Deadwood on ground [ln(n+1)]") +
  ylab("")+
  xlab("")+ 
  #ylim(0,100)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
dw.number.plot

# plot caco
caco.plot<-ggplot(plot.meta, aes(x=Forest, y=Caco, fill = Forest)) +
  ggtitle("Canopy cover index") +
  ylab("")+
  xlab("")+ 
  #ylim(0,100)+
  geom_point(aes(fill=Forest), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
caco.plot

# plot slope
slope.plot<-ggplot(plot.meta, aes(x=Forest, y=slope.var, fill = Forest)) +
  ggtitle("Slope inclination [°]") +
  ylab("")+
  xlab("")+ 
  #ylim(0,100)+
  geom_point(aes(fill=Forest), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
slope.plot


