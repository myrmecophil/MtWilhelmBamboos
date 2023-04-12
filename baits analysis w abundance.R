## Numba-Kausi Baits & nests - R script by Phil 10 march 2023

# to do: 
# check all statiscal models

# Associated csv files:

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
    "DHARMa")

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

# test for overdispersion
overdisp_fun<-function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
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
nest.raw<-subset(nest.raw, Moved.to.plot!= "KP1-S3" & Moved.to.plot!= "KP2-S3" & Moved.to.plot!= "KP4-S3"& Moved.to.plot!= "KP6-S3")

# parse number of nesting ants
nest.raw$nesting.estimate <- apply(nest.raw["N.ants..camera."],1, parse_abundance)

# Count a bamboo nest as occupied if there were >2 worker ants inside or 1 queen
# add 'moved to' forest and block
nesters <-nest.raw %>%
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

# Prepare tree-level metadata for later models
tree.meta.sub<-tree.metada.raw %>%
  dplyr::select(Plot, mastercode, Trunk.perimenter..cm., tree.ID, Lianas.n, deadwood.., deadwood., Bait.Bamboo.Height.m.)%>%
  rename(trunk =Trunk.perimenter..cm., Code = mastercode, height = Bait.Bamboo.Height.m., dw.percent = deadwood., dw.number=deadwood..)

plot.metada.sub<-  plot.metada.raw %>%
  dplyr::select(Plot_code, Block, Canopy.cover..Caco., elevation, slope.var)%>%
  rename(Plot=Plot_code, Caco =Canopy.cover..Caco.)

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

#----------------------------------------------------------#
# 1.2 Bait diversity -----
#----------------------------------------------------------#

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

# 
baiting.diversity.e<-merge(diversity, plot.meta)

## Plot it
#Labs
labs <- expression("Kausi", "Numba")

baitdiversity<-ggplot(diversity, aes(x=Forest, y=expH, fill=Forest)) +
  ggtitle("Bait species diversity") +
  geom_boxplot()+
  ylim(0,20)+
  scale_fill_manual(labels=labs, values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=labs)+
  ylab("Species diversity [expH]")+
  xlab("")+ 
  guides(fill="none")+
  theme_bw()
baitdiversity

#----------------------------------------------------------#
#  1.3 Bait composition -----
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

#----------------------------------------------------------#
# 2.1 Nest occupancy -----
#----------------------------------------------------------#
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

#----------------------------------------------------------#
# 2.2 Nest diversity -----
#----------------------------------------------------------#

# remove species from unoccupied bamboos (i.e., foragers)
nesters<-phase.nests
nesters$AntSpCode[nesters$occupancy==0]  <- ""

# subset phase
phase1<-subset(nesters, phase=="phase 1")
phase2<-subset(nesters, phase=="phase 2")

# phase 1 plot x species incidence
ant_n.phase1<-dcast(phase1, formula = Plot ~ AntSpCode, length)
rownames(ant_n.phase1)<-ant_n.phase1[,1] # Plot as rownames
ant_n.phase1<-ant_n.phase1[,-c(1,2)]

# phase 2 plot x species incidence
ant_n.phase2<-dcast(phase2, formula = Plot ~ AntSpCode, length)
rownames(ant_n.phase2)<-ant_n.phase2[,1] # Plot as rownames
ant_n.phase2<-ant_n.phase2[,-c(1,2)]

# extract diversity values for phase 1
data2 <- data.frame(row.names(ant_n.phase1))     # create a new file with plot numbers
names(data2) <- "Plot" # rename variable "plot"
data2$Shannon <- diversity(ant_n.phase1, index = "shannon", base = exp(1)) # Shannon diversity per plot
data2$Richness <- specnumber(ant_n.phase1)                          # Richness per plot
data2$expH <- exp(data2$Shannon)                                         # exponential Shannon diversity per plot
data2$expH[data2$Richness == 0] <- 0                                  # define exp(Shannon) = 0
data2$evenness <- data2$Shannon/log(specnumber(ant_n.phase1)) # evenness per plot
str(data2)

# add Plot location
phase1.diversity<-merge(data2, plot_m)
phase1.diversity$Plot<-as.factor(phase1.diversity$Plot)

# extract diversity values for phase 2
data2 <- data.frame(row.names(ant_n.phase2))     # create a new file with plot numbers
names(data2) <- "Plot" # rename variable "plot"
data2$Shannon <- diversity(ant_n.phase2, index = "shannon", base = exp(1)) # Shannon diversity per plot
data2$Richness <- specnumber(ant_n.phase2)                          # Richness per plot
data2$expH <- exp(data2$Shannon)                                         # exponential Shannon diversity per plot
data2$expH[data2$Richness == 0] <- 0                                     # define exp(Shannon) = 0
data2$evenness <- data2$Shannon/log(specnumber(ant_n.phase2)) # evenness per plot
str(data2)

# add Plot location
phase2.diversity<-merge(data2, plot_m)
phase2.diversity$Plot<-as.factor(phase2.diversity$Plot)

## gg plot it
nests.phase1<-ggplot(phase1.diversity, aes(x=Forest, y=expH, fill=Forest)) +
  ggtitle("Phase 1 Bamboo nester diversity") +
  geom_boxplot()+
  ylim(0,5)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Species diversity")+
  xlab("")+ 
  theme(axis.text.x = element_text(size=12))
nests.phase1

nests.phase2 <-ggplot(phase2.diversity, aes(x=Forest, y=expH, fill=Forest)) +
  ggtitle("Phase 2 Bamboo nester diversity") +
  geom_boxplot()+
  #ylim(0,5)+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Species diversity")+
  xlab("")+ 
  theme(axis.text.x = element_text(size=12))
nests.phase2
#
phase1.diversity$phase<- "phase 1"
phase2.diversity$phase<- "phase 2"

diversity.nester <-rbind(phase1.diversity, phase2.diversity)

#
diversity.nester.e<-merge(diversity.nester, plot.meta)

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
# Statistics -----
#----------------------------------------------------------#

# Q1: Should plot means of environmental variables be transformed as well? Such as dw.percent_mean
# Q2: Can/should 'tree ID' be included in the models? For a couple baits/bamboos, the trees are the same (and thus, also dbh and liana counts). 
# One possibility would be to include them as another crossed random factor, but not sure if it makes sense since its not always the case: (1|Block/Plot/tree.ID)

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

# Same for plot means? maybe not
hist(plot.meta$Caco)
hist(plot.meta$slope.var)
hist(plot.meta$trunk_mean)
hist(plot.meta$Lianas.n_mean)
hist(plot.meta$dw.percent_mean)
hist(plot.meta$dw.number_mean)
hist(plot.meta$height_mean)

# Are any plot attributes different between elevations?

# Are canopy covers different between forests?
caco_model <- lmer(Caco~  elevation +(1|Block), data = plot.meta)
summary(caco_model)

# Is slope different?
slope_model <- lmer(slope.var~  elevation +(1|Block), data = plot.meta)
summary(caco_model) # ns

# Is dbh different between forests?
model_dbh<- lmer(trunk.log ~  elevation + (1|Block/Plot), data = tree.meta)
summary(model_dbh)# ns

# Is number of deadwood pieces different?
model_dw.n<- lmer(dw.number.log ~  elevation + (1|Block/Plot), data= tree.meta)
summary(model_dw.n)# ns

# is % deadwood cover different?
model_dw.p<- lmer(dw.percent.log ~  elevation + (1|Block/Plot), data= tree.meta)
summary(model_dw.p)# ns

# Is number of lianas different?
model_lianas<- lmer(Lianas.n.log ~  elevation + (1|Block/Plot), data= tree.meta)
summary(model_lianas)
# => slighly more lianas in higher elevation

# bamboo/bait height different?
model_height<- lmer(height.sqrt ~  elevation + (1|Block/Plot), data= tree.meta)
summary(model_height) # ns

#----------------------------------------------------------#
## Bait statistics
#----------------------------------------------------------#

# Are there more species per bait in higher elevation?

# get environmental data
bait.double.e<-merge(bait.double, tree.meta)

# Binominal regression
bait.double.model <- glmmTMB(species.per.bait~elevation+(1|Block/Plot),
                               data=bait.double.e,
                               family=nbinom2)
summary(bait.double.model)
overdisp_fun(bait.double.model)

# model does not run properly, what to do? same for nbinom1

# with environmental factors
bait.double.model.e <- glmmTMB(species.per.bait~elevation+Lianas.n+height+dw.percent+dw.number+trunk+Caco+slope.var+(1|Block/Plot),
                                 data=bait.double.e,
                                 family=nbinom2)
summary(bait.double.model.e)
overdisp_fun(bait.double.model.e)
# model does not run properly, what to do? same for nbinom1

## Bait occupancy
# get environmental data
baiting.incidence.e<-merge(baiting.incidence, tree.meta)

# Binominal regression
baitoccupancy.model <- glmmTMB(occupancy~elevation+(1|Block/Plot),
                               data=baiting.incidence.e,
                               family=binomial)
summary(baitoccupancy.model)
overdisp_fun(baitoccupancy.model)
testDispersion(baitoccupancy.model)
simulateResiduals(fittedModel = baitoccupancy.model, plot = T) # within group deviation from normality - probably no big problem?
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model))

# with environmental factors
baitoccupancy.model.e <- glmmTMB(occupancy~elevation+Lianas.n.log+height.sqrt+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+(1|Block/Plot),
                                 data=baiting.incidence.e,
                                 family=binomial)
summary(baitoccupancy.model.e)
overdisp_fun(baitoccupancy.model.e)
#
testDispersion(baitoccupancy.model.e) # ok
simulateResiduals(baitoccupancy.model.e, plot = T) # lower quantiles a bit weird  - probably not a big problem?
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model.e)) # ok

# with abundance 
baitoccupancy.model1 <- glmmTMB(Abundance~elevation+(1|Block/Plot),
                                data=baiting.incidence.e,
                                family=nbinom1)
summary(baitoccupancy.model1)
overdisp_fun(baitoccupancy.model1)

#
testDispersion(baitoccupancy.model1) # for nbinom1 overdispersion ok, for nbinom2 not good
simulateResiduals(baitoccupancy.model1, plot = T) # doesnt look good for either family
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model1)) # ok for binom1, not ok for binom2

# with abundance and environmental
baitoccupancy.model.e1 <- glmmTMB(Abundance~elevation+Lianas.n.log+height.sqrt+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=nbinom1)
summary(baitoccupancy.model.e1)
overdisp_fun(baitoccupancy.model.e1)
#
testDispersion(baitoccupancy.model.e1) # for nbinom1 overdispersion ok, for nbinom2 not good
simulateResiduals(baitoccupancy.model.e1, plot = T) # doesnt look good for either family
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model.e1)) # ok for binom1, not ok for binom2

#----------------------------------------------------------#
## Bamboo statistics
#----------------------------------------------------------#

# bamboo occupancy models

# Get environmental metadata
bamboo.incidence.e<-merge(phase.nests, tree.meta)

# binominal model of bamboo nesting using phase
bamboo.occupancy.model <- glmmTMB(occupancy~elevation+phase+(1|Block/Plot),
                                  data=bamboo.incidence.e,
                                  family=binomial)
summary(bamboo.occupancy.model)
overdisp_fun(bamboo.occupancy.model)
#
testDispersion(bamboo.occupancy.model) # ok
simulateResiduals(bamboo.occupancy.model, plot = T) # slight deviance in lower quantiles, probably fine
testZeroInflation(simulateResiduals(fittedModel = bamboo.occupancy.model)) # ok 


# binominal model of bamboo nesting using only controls
bamboo.occupancy.model <- glmmTMB(occupancy~elevation+phase+(1|Block/Plot),
                                  data=subset(bamboo.incidence.e, Treatment=="control"),
                                  family=binomial)
summary(bamboo.occupancy.model)
overdisp_fun(bamboo.occupancy.model)
#
testDispersion(bamboo.occupancy.model) # ok
simulateResiduals(bamboo.occupancy.model, plot = T) # slight deviance in lower quantiles, probably fine
testZeroInflation(simulateResiduals(fittedModel = bamboo.occupancy.model)) # ok 

# or alternatively, splitting phase 2 into different treatments. Works only as fixed factor since its levels are not well replicated across forest types
bamboo.occupancy.model2 <- glmmTMB(occupancy~elevation+new.Treatment+(1|Block/Plot),
                                   data=bamboo.incidence.e,
                                   family=binomial)
summary(bamboo.occupancy.model2)
overdisp_fun(bamboo.occupancy.model2)
#
testDispersion(bamboo.occupancy.model2) # ok
simulateResiduals(bamboo.occupancy.model2, plot = T) # slight deviance in lower quantiles, probably fine
testZeroInflation(simulateResiduals(fittedModel = bamboo.occupancy.model2)) # ok 


# binominal model of bamboo nesting using phase with environment
bamboooccupancy.model.e <- glmmTMB(occupancy~elevation+Lianas.n.log+height.sqrt+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+phase+(1|Block/Plot),
                                   data=bamboo.incidence.e,
                                   family=binomial)
summary(bamboooccupancy.model.e)
overdisp_fun(bamboooccupancy.model.e)
#
testDispersion(bamboooccupancy.model.e) # ok
simulateResiduals(bamboooccupancy.model.e, plot = T) # all good
testZeroInflation(simulateResiduals(fittedModel = bamboooccupancy.model.e)) # ok 

# with abundance 
bamboo.occupancy.model3 <- glmmTMB(nesting.estimate~elevation+phase+(1|Block/Plot),
                                   data=bamboo.incidence.e,
                                   family=nbinom1)
summary(bamboo.occupancy.model3)
overdisp_fun(bamboo.occupancy.model3)
#
testDispersion(bamboo.occupancy.model3) # ok
simulateResiduals(bamboo.occupancy.model3, plot = T) # ok for nbinom1 / not good for nbinom2
testZeroInflation(simulateResiduals(fittedModel = bamboo.occupancy.model3)) # ok 

# with abundance and environmental
bamboo.occupancy.model.e4 <- glmmTMB(nesting.estimate~elevation+Lianas.n.log+height.sqrt+dw.percent.log+dw.number.log+trunk.log+Caco+slope.var.log+phase+(1|Block/Plot),
                                     data=bamboo.incidence.e,
                                     family=nbinom1)
summary(bamboo.occupancy.model.e4)
overdisp_fun(bamboo.occupancy.model.e4)
#
testDispersion(bamboo.occupancy.model.e4) # ok
simulateResiduals(bamboo.occupancy.model.e4, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = bamboo.occupancy.model.e4)) # ok 

# Species diversity in bamboo model
bamboo.diversity.model <- lmer(expH ~ elevation * phase+ (1|Block), data = diversity.nester.e)
summary(bamboo.diversity.model)
overdisp_fun(bamboo.diversity.model)
#
testDispersion(bamboo.diversity.model) # ok
simulateResiduals(bamboo.diversity.model, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model)) # zero inflated - troubles! 

# diversity changes with environment metadata
bamboo.diversity.model.e <- lmer(expH ~ elevation+slope.var+Caco+trunk_mean+Lianas.n_mean+dw.number_mean+dw.percent_mean+phase+ (1|Block), data = diversity.nester.e)
summary(bamboo.diversity.model.e)
overdisp_fun(bamboo.diversity.model.e) # overdispersion?
#
testDispersion(bamboo.diversity.model.e) # ok
simulateResiduals(bamboo.diversity.model.e, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model.e)) # zero inflated - troubles! 

# Baiting diversity 
# Fit a LMM for shannon diversity
baitdiversity.model <- lmer(expH ~ elevation + (1|Block), data = baiting.diversity.e)
summary(baitdiversity.model)
overdisp_fun(baitdiversity.model)
#
testDispersion(baitdiversity.model) # ok
simulateResiduals(baitdiversity.model, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model)) # ?

# add environmental factors
baitdiversity.model.e <- lmer(expH ~elevation+slope.var+Caco+trunk_mean+Lianas.n_mean+dw.number_mean+dw.percent_mean + (1|Block),
                              data = baiting.diversity.e)
summary(baitdiversity.model.e)
overdisp_fun(baitdiversity.model.e) # overdipsersion?
#
testDispersion(baitdiversity.model.e) # ok
simulateResiduals(baitdiversity.model.e, plot = T) # weird 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.e)) # weird

# Fit a LMM for evenness
baitdiversity.model.eve <- lmer(evenness ~ elevation + (1|Block), data = baiting.diversity.e)
summary(baitdiversity.model.eve)
overdisp_fun(baitdiversity.model.eve)
#
testDispersion(baitdiversity.model.eve) # ok
simulateResiduals(baitdiversity.model.eve, plot = T) # weird 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.eve)) # weird

# add environmental factors
baitdiversity.model.eve.e <- lmer(evenness ~elevation+slope.var+Caco+trunk_mean+Lianas.n_mean+dw.number_mean+dw.percent_mean + (1|Block),
                              data = baiting.diversity.e)
summary(baitdiversity.model.eve.e)
overdisp_fun(baitdiversity.model.eve.e)
#
testDispersion(baitdiversity.model.eve.e) # ok
simulateResiduals(baitdiversity.model.eve.e, plot = T) # weird 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.eve.e)) # weird

#----------------------------------------------------------#
# Figure summaries -----
#----------------------------------------------------------#
## with incidences

#Bamboo diversity
nests.diversity.phase.plot

#Bait diversity
baitdiversity

#bamboo occupied
nest.occupancy.phase.plot

#Bamboo 
nesting.proportion

### Bait composition
# as proportion of incidence
bait.incidence.plot
# as proportion of abundance
abundance.plot

figure_baitsbamboos <- ggarrange(nests.diversity.phase.plot, baitdiversity, nesting.proportion, abundance.plot,
                               labels = c("A", "B", "C", "D"),
                               ncol = 2, nrow = 2, common.legend = F
)
figure_baitsbamboos

