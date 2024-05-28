# Bamboo subscript

#----------------------------------------------------------#
### List of R-packages
#----------------------------------------------------------#

package_list <- 
  c("dplyr",
    "tidyr",
    "reshape2",
    "vegan",
    "betapart",
    "ggplot2",
    "ggthemes",
    "ggpubr",
    "glmmTMB",
    "DHARMa",
    "car",
    "emmeans",
    "effects",
    "plotrix")

# install all packages
#sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# Citations
#sapply(package_list, citation)

setwd("~/GitHub/MtWilhelmBamboos")

set.seed(1234)

# load data
nest.raw <- read.csv(file="nest.raw.csv", header=T)
baiting.data <-read.csv(file="baiting.data.csv", header=T)
phase.nests <- read.csv(file="phase.nests.csv", header=T)
tree.meta <- read.csv(file="tree.meta.csv", header=T)
plot.meta<- read.csv(file="plot.meta.csv", header=T)
plot.meta2 <- read.csv(file="plot.meta2.csv", header=T)

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

labs1 <- expression("low elevation control", "midelevation control", "other elevation transfer", "same elevation transfer")

nest.r<-nest.raw3%>%
  count(Forest.Treatment, survived.d)

nest.r$percent<-nest.r$n/128*100
nest.r$percent[nest.r$Forest.Treatment== "Numba Primary other elevation"] <- nest.r$n[nest.r$Forest.Treatment== "Numba Primary other elevation"]/253*100 # more total bamboos in other elevation treatment (3 of 256 lost)

# Survived as proportion of total nests
nest.r%>% filter(survived.d != 'empty')%>%
  ggplot(aes(Forest.Treatment, percent, fill = survived.d)) + 
  xlab("Treatment") +
  scale_x_discrete(labels=labs1)+
  ylab("% of Bamboos translocated") +
  labs(fill = "Occupied bamboo nest")+
  ggtitle("Survived as proportion of total bamboos") +
  scale_fill_brewer(palette="Dark2")+
  geom_col() + 
  geom_text(aes(label = n), 
            position = position_stack(vjust=.5),colour = "black", size = 5)+
  theme_hc()

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

# this step removes poison plots KPX-S3
phase.nests<-subset(phase.nests, Plot!= "KP1-S3" & Plot!= "KP2-S3" & Plot!= "KP4-S3"& Plot!= "KP6-S3")


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

#PK: phase.final is odd: it contains 52 plots and occupancy for each, which fits 
#28+(28-4)= all our  established plots minus eradicated, so even those plots not having bamboos at a phase, are included (duplicated?)  
#e.g. KP-6-S1 is twice in data (example of duplicate)
#I am not sure where data duplication happen
#The phase.final  data are used perhaps only for Nest occupancy figure, not models, so correction is needed for the Fig. but shall not touch modelling
#it may still affect the mean and offsets points I think in the figure though

# Bamboo nest occupancy by phase
nest.occupancy.phase.plot<-ggplot(phase.final, aes(x=Forest, y=proportion, fill=phase)) +
  ggtitle("Nest occupancy [%]") +
  labs(fill = "Phase")+  
  ylab("")+
  xlab("")+ 
  ylim(0,50)+
  geom_point(aes(fill=phase), size=3, shape=21, colour="grey20", position=position_jitterdodge(0.1))+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.occupancy.phase.plot


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

#PK: Nest.final (plotting only): 56 replicates by plots and strata, but some are pooled across two phases (16 vs 32 bamboos)
#We shall be using for plotting rather plot x  time  x stratum replicate (each about equal effort 16 bamboos in all as a "sample") I think
# as this would be our sample if we would compare plot-wise proportions as mimics to our binomial tests we do on bamboos
# (=if we would use actual plot-proportions for quasi-binomial test with phase included, you would not pool phases per plot)
# the update will not change main pattern I think but will be a bit more data points

# proportion
nests.final$proportion<-nests.final$occupancy/nests.final$n*100

# plot it
nests.final.occupancy.plot<-ggplot(nests.final, aes(x=Forest, y=proportion, fill = Stratum)) +
  ggtitle("Nest occupancy [%]") +
  ylab("")+
  xlab("")+ 
  ylim(0,50)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nests.final.occupancy.plot


##### ignoring phase and empty bamboos: Proportion of abundance
# nest abundance (=size) plot

#PK: "nest.abundance1": We claim in paper that we analyse and plot mean abundance of ants per occupied nest, so
#we shall have each point as a nest, as we say in legends and method, and how also GLMM is done,
# (not as a pool of  indiv. per plot across nests and phases as now, which is incorrect approach), 
#to correct the  nest.abundance1 data calculation and the figure

# summarize bamboo nest abundance
nest.abundance1 <- phase.nests %>%
  group_by(Forest, Plot,Stratum) %>%
  summarize(Abundance = sum(nesting.estimate))

# Remove empty values
nest.abundance1 <- subset(nest.abundance1, Abundance!=0)

nest.abundance.stratum<-ggplot(nest.abundance1, aes(x=Forest, y=log(Abundance), fill = Stratum)) +
  ggtitle("Nest size [log]") +
  ylab("")+
  xlab("")+ 
  #ylim(0,50)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1))+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.abundance.stratum

##### Phase, ignoring stratum: abundance (nest size)

# summarize bamboo nest abundance
nest.abundance2 <- phase.nests %>%
  group_by(Forest, Plot, phase) %>%
  summarize(Abundance = sum(nesting.estimate))

# Remove empty values and plot
nest.abundance2 <- subset(nest.abundance2, Abundance!=0)

nest.abundance.phase<-ggplot(nest.abundance2, aes(x=Forest, y=log(Abundance), fill = phase)) +
  ggtitle("Nest size [log]") +
  ylab("")+
  xlab("")+ 
  #ylim(0,50)+
  geom_point(aes(fill=phase), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1))+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.abundance.phase

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
# PK: Similar to baits you are summing not nest incidences but all workers counts across the nests 
#This  differs how  it is it now in the paper and is not usual for ants (I think  composition by n of nests would be more logical) 
# PK: in case of  baits it was 3% (should not we use same threshold in same figure for both, eihter 3 or 5?)
nest.abundance$AntSpCode[nest.abundance$percentage < 5] <- "other species"

# summarize bait incidence counts
nest.abundance<- nest.abundance%>% 
  group_by(Forest, AntSpCode) %>% 
  summarise(percentage = sum(percentage))

# Plot it
# relevel
nest.abundance$AntSpCode <- as.factor(nest.abundance$AntSpCode)
nest.abundance$AntSpCode <- relevel(nest.abundance$AntSpCode, "other species")

abundance.proportion.nest<-ggplot(nest.abundance, aes(x = Forest, y = percentage, fill = AntSpCode)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative abundance [%]") +
  ggtitle("bamboo species composition") +
  #scale_fill_brewer(palette="Dark2")+
  theme_minimal()
abundance.proportion.nest

#----------------------------------------------------------#
# 2.2 Nest diversity -----
#----------------------------------------------------------#

# remove species from unoccupied bamboos (i.e., foragers)
nesters<-phase.nests
#PK: This line does not remove any data? (occupied and not occupied still in matrix, but seems not matter for the code so perhaps redundant)
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
data2$evenness <- data2$expH/specnumber(ant_n.phase1) # evenness per plot
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
data2$evenness <- data2$expH/specnumber(ant_n.phase2) # evenness per plot
str(data2)

# add Plot location
plot_m4<-subset(nesters, phase=="phase 2")
plot_m4<-distinct(plot_m4, plot.stratum, .keep_all = TRUE)
phase2.diversity<-merge(data2, plot_m4)

diversity.nester <-rbind(phase1.diversity, phase2.diversity)

# Plot nester diversity, phase dependence
nest.diversity.phase.plot<-ggplot(diversity.nester, aes(x=Forest, y=expH, fill=phase)) +
  ggtitle("Nester diversity [expH]") +
  labs(fill = "Phase")+  
  ylab("")+
  xlab("")+ 
  geom_point(aes(fill=phase), size=3, shape=21, colour="grey20", position=position_jitterdodge(0.1))+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.diversity.phase.plot

# Plot nester diversity, stratum dependent
nest.diversity.stratum.plot<-ggplot(diversity.nester, aes(x=Forest, y=expH, fill=Stratum)) +
  ggtitle("Nest species diversity [expH]") +
  labs(fill = "Stratum")+  
  ylab("")+
  xlab("")+ 
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20", position=position_jitterdodge(0.1), alpha=0.8)+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.diversity.stratum.plot

#----------------------------------------------------------#
# 2.2 Nest composition -----
#----------------------------------------------------------#

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

# sum check (should be 100)
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

#PK: I think this is correct chart by n of nests and shall be used (and the one on baits done yet similarly)
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


#### beta partitioning, abundance based

## Species overlap, canopy only
# make matrices
bamboo.matrix.ca <- dcast(subset(nesters, Stratum=="CA"), formula = Plot ~ AntSpCode, length)
bamboo.matrix.ca

# set rownames
rownames(bamboo.matrix.ca) <- bamboo.matrix.ca[, 1]
bamboo.matrix.ca <- bamboo.matrix.ca[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.bamboo <- vegdist(bamboo.matrix.ca, method = "bray")

# global matrix of plot dissimiliarities
dist.bamboo <- as.matrix(dist.bamboo)

## select subsets
# midelevation vs lowland
b.m_vs_l.ca <- dist.bamboo[c(1:8), c(9:28)]

distance.bamboo.elev.ca <- as.data.frame(rowMeans(b.m_vs_l.ca))
mean(colMeans(b.m_vs_l.ca))
std.error(b.m_vs_l.ca)

## Species overlap, understorey only
# make matrices
bamboo.matrix.un <- dcast(subset(nesters, Stratum=="UN"), formula = Plot ~ AntSpCode, length)
bamboo.matrix.un

# set rownames
rownames(bamboo.matrix.un) <- bamboo.matrix.un[, 1]
bamboo.matrix.un <- bamboo.matrix.un[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.bamboo <- vegdist(bamboo.matrix.un, method = "bray")

# global matrix of plot dissimiliarities
dist.bamboo <- as.matrix(dist.bamboo)

## select subsets
# midelevation vs lowland
b.m_vs_l.un <- dist.bamboo[c(1:8), c(9:28)]

distance.bamboo.elev.un <- as.data.frame(rowMeans(b.m_vs_l.un))
mean(colMeans(b.m_vs_l.un))

## Species overlap, total
# make matrices
bamboo.matrix <- dcast(nesters, formula = Plot ~ AntSpCode, length)
bamboo.matrix

# set rownames
rownames(bamboo.matrix) <- bamboo.matrix[, 1]
bamboo.matrix <- bamboo.matrix[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.bamboo <- vegdist(bamboo.matrix, method = "bray")

# global matrix of plot dissimiliarities
dist.bamboo <- as.matrix(dist.bamboo)

## select subsets
# midelevation vs lowland
b.m_vs_l <- dist.bamboo[c(1:8), c(9:28)]

distance.bamboo.elev.un <- as.data.frame(rowMeans(b.m_vs_l))
mean(colMeans(b.m_vs_l))

#----------------------------------------------------------#
# 3.3 Bamboo Nest Statistics -----
#----------------------------------------------------------#

# bamboo occupancy models

# Get environmental metadata
bamboo.incidence.e<-merge(phase.nests, tree.meta, by = "Code")
bamboo.incidence.e$new.Treatment<-as.factor(bamboo.incidence.e$new.Treatment)

# re-level for models
bamboo.incidence.e$new.Treatment <- relevel(bamboo.incidence.e$new.Treatment, ref = "first.placement")

# binominal model of bamboo nesting using phase
bamboo.occupancy.model1 <- glmmTMB(occupancy~Forest.x
                                   +Stratum
                                   +new.Treatment
                                   +(1|Block.x/Plot.x),
                                   data=bamboo.incidence.e,
                                   family=binomial)
bamboo.occupancy.model2 <- glmmTMB(occupancy~Forest.x
                                   *Stratum
                                   +new.Treatment+(1|Block.x/Plot.x),
                                   data=bamboo.incidence.e,
                                   family=binomial)

anova(bamboo.occupancy.model1, bamboo.occupancy.model2)
# ns, no interaction

summary(bamboo.occupancy.model1)
#
testDispersion(bamboo.occupancy.model1) # ok
simulateResiduals(bamboo.occupancy.model1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboo.occupancy.model1)) # ok 

# binominal model of bamboo nesting with environment
bamboooccupancy.model.e1 <- glmmTMB(occupancy~Forest.x
                                    +Stratum
                                    +scale(log(Lianas.n+1))
                                    +scale(log(dw.percent+1))
                                    +scale(log(dw.number+1))
                                    +scale(log(trunk+1))
                                    +scale(log(Caco+1))
                                    +scale(log(slope.var+1))
                                    +new.Treatment
                                    +(1|Block.x/Plot.x),
                                    data=bamboo.incidence.e,
                                    family=binomial)

bamboooccupancy.model.e2 <- glmmTMB(occupancy~Forest.x
                                    *Stratum
                                    *scale(log(Lianas.n+1))
                                    +scale(log(dw.percent+1))
                                    +scale(log(dw.number+1))
                                    +scale(log(trunk+1))
                                    +scale(log(Caco+1))
                                    +scale(log(slope.var+1))
                                    +new.Treatment
                                    +(1|Block.x/Plot.x),
                                    data=bamboo.incidence.e,
                                    family=binomial)
anova(bamboooccupancy.model.e1, bamboooccupancy.model.e2) # ns -  better using no interactions

summary(bamboooccupancy.model.e1)
#
testDispersion(bamboooccupancy.model.e1) # ok
simulateResiduals(bamboooccupancy.model.e1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboooccupancy.model.e1)) # ok 

# with nest size (abundance) 
bamboo.abundance.model1 <- glmmTMB(nesting.estimate~Forest.x
                                   +Stratum
                                   +new.Treatment
                                   +(1|Block.x/Plot.x),
                                   data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                   family=nbinom2)
bamboo.abundance.model2 <- glmmTMB(nesting.estimate~Forest.x
                                   *Stratum
                                   +new.Treatment
                                   +(1|Block.x/Plot.x),
                                   data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                   family=nbinom2)
anova(bamboo.abundance.model1,bamboo.abundance.model2) # ns, no interaction
summary(bamboo.abundance.model1)
#
testDispersion(bamboo.abundance.model1) # ok
simulateResiduals(bamboo.abundance.model1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboo.abundance.model1)) # ok 

# with nest size (abundance) and environment
bamboo.abundance.model.e4 <- glmmTMB(nesting.estimate~Forest.x
                                     +Stratum
                                     +scale(log(Lianas.n+1))
                                     +scale(log(dw.percent+1))
                                     +scale(log(dw.number+1))
                                     +scale(log(trunk+1))
                                     +scale(log(Caco+1))
                                     +scale(log(slope.var+1))
                                     +new.Treatment
                                     +(1|Block.x/Plot.x),
                                     data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                     family=nbinom2)

bamboo.abundance.model.e5 <- glmmTMB(nesting.estimate~Forest.x
                                     *Stratum
                                     *scale(log(Lianas.n+1))
                                     +scale(log(dw.percent+1))
                                     +scale(log(dw.number+1))
                                     +scale(log(trunk+1))
                                     +scale(log(Caco+1))
                                     +scale(log(slope.var+1))
                                     +new.Treatment
                                     +(1|Block.x/Plot.x),
                                     data=subset(bamboo.incidence.e, nesting.estimate!=0),
                                     family=nbinom2)

anova(bamboo.abundance.model.e4, bamboo.abundance.model.e5) # ns, no interaction 

summary(bamboo.abundance.model.e4)
#
testDispersion(bamboo.abundance.model.e4) # ok
simulateResiduals(bamboo.abundance.model.e4, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboo.abundance.model.e4)) # ok
# 
#####   Bamboo Species diversity
# get environment
diversity.nester.e<-merge(diversity.nester,plot.meta2, by = "plot.stratum", all.x = TRUE)

# re-level for models
diversity.nester.e$new.Treatment<-as.factor(diversity.nester.e$new.Treatment)
diversity.nester.e$new.Treatment <- relevel(diversity.nester.e$new.Treatment, ref = "first.placement")
#
bamboo.diversity.model1 <- glmmTMB((expH+1) ~ Forest.x
                                   +Stratum
                                   +new.Treatment
                                   +(1|Block.x/Plot.x), 
                                   #zi=~Forest+new.Treatment, 
                                   data = diversity.nester.e, 
                                   family=gaussian(link="log"))

bamboo.diversity.model2 <- glmmTMB((expH+1) ~ Forest.x
                                   *Stratum
                                   +new.Treatment
                                   +(1|Block.x/Plot.x), 
                                   #zi=~Forest+new.Treatment, 
                                   data = diversity.nester.e, 
                                   family=gaussian(link="log"))


anova(bamboo.diversity.model1,bamboo.diversity.model2) # ns, no interaction
#
summary(bamboo.diversity.model1)

#
testDispersion(bamboo.diversity.model1) # ok
simulateResiduals(bamboo.diversity.model1, plot = T) # ok-ish
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model1)) # ok

# bamboo diversity changes with environment metadata
bamboo.diversity.model.e1 <- glmmTMB((1+expH) ~ Forest.x
                                     +Stratum
                                     +scale(log(Lianas.n_mean))
                                     +scale(log(dw.percent_mean))
                                     +scale(log(dw.number_mean))
                                     +scale(log(trunk_mean))
                                     +scale(log(Caco+1))
                                     +scale(log(slope.var+1))
                                     +new.Treatment
                                     +(1|Block.x/Plot.x),
                                     zi=~Forest.x+new.Treatment,
                                     data = diversity.nester.e, family=gaussian(link=log))

bamboo.diversity.model.e2 <- glmmTMB((1+expH) ~ Forest.x
                                     *Stratum
                                     *scale(log(Lianas.n_mean))
                                     +scale(log(dw.percent_mean))
                                     +scale(log(dw.number_mean))
                                     +scale(log(trunk_mean))
                                     +scale(log(Caco+1))
                                     +scale(log(slope.var+1))
                                     +new.Treatment
                                     +(1|Block.x/Plot.x),
                                     zi=~Forest.x+new.Treatment,
                                     data = diversity.nester.e, family=gaussian(link=log))

anova(bamboo.diversity.model.e1, bamboo.diversity.model.e2) # no interaction better

summary(bamboo.diversity.model.e1)
#
testDispersion(bamboo.diversity.model.e1) # ok
simulateResiduals(bamboo.diversity.model.e1, plot = T) # ok-ish
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model.e1)) # ok

