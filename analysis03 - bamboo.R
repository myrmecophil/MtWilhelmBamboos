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

# set labels
labs <- expression("lowland", "mid-elevation")

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

labs1 <- expression("low elevation control", "mid-elevation control", "other elevation transfer", "same elevation transfer")

nest.r<-nest.raw3%>%
  count(Forest.Treatment, survived.d)

nest.r$percent<-nest.r$n/128*100
nest.r$percent[nest.r$Forest.Treatment== "Numba Primary other elevation"] <- nest.r$n[nest.r$Forest.Treatment== "Numba Primary other elevation"]/255*100 # more total bamboos in other elevation treatment

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

# define unique code plot.stratum
phase.nests$plot.stratum<-paste(phase.nests$Plot,phase.nests$Stratum)

# Get environmental metadata
bamboo.incidence.e<-merge(phase.nests, tree.meta)

### Plot nest occupancy as average proportion per plot
phase1.n<-subset(phase.nests, phase=="phase 1")
phase2.n<-subset(phase.nests, phase=="phase 2")

#count total bamboos per plot per phase
phase1.c<-count(phase1.n, Plot)
phase2.c<-count(phase2.n, Plot)

# summarize occupied per plot per phase
phase1 <- phase1.n %>%
  group_by(Forest, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase1$phase <- "phase 1"

# merge counts
phase1<-merge(phase1, phase1.c)

# same for phase 2
phase2<- phase2.n %>%
  group_by(Forest, Plot) %>%
  summarize(occupancy = sum(occupancy))

phase2$phase<-"phase 2"
phase2<-merge(phase2, phase2.c)

# merge
phase.final<-rbind(phase1, phase2)

# proportion
phase.final$proportion<-phase.final$occupancy/phase.final$n*100

# Bamboo nest occupancy by phase
nest.occupancy.phase.plot<-ggplot(phase.final, aes(x=Forest, y=proportion, fill=phase)) +
  ggtitle("Nest occupancy") +
  labs(fill = "Phase")+  
  ylab("%")+
  xlab("")+ 
  ylim(0,50)+
  geom_point(aes(fill=phase), size=3, shape=21, colour="grey20", position=position_jitterdodge(0.1))+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.occupancy.phase.plot

##### ignoring empty bamboos: Proportion of occupancy

### Plot nest occupancy as average proportion per plot
phase1.n<-subset(phase.nests, phase=="phase 1")
phase2.n<-subset(phase.nests, phase=="phase 2")

#count total nests per plot and stratum
nests.p1.c<-count(phase1.n, plot.stratum)
nests.p2.c<-count(phase2.n, plot.stratum)

# summarize occupied per plot
nests.p1 <- phase1.n %>%
  group_by(Forest, Plot, Stratum, plot.stratum) %>%
  summarize(occupancy = sum(occupancy))

# merge with nest counts
nests.p1$phase<-"phase 1"
nests.p1<-merge(nests.p1, nests.p1.c)

# same with phase 2
nests.p2 <- phase2.n %>%
  group_by(Forest, Plot, Stratum, plot.stratum) %>%
  summarize(occupancy = sum(occupancy))

# merge with nest counts
nests.p2$phase<-"phase 2"
nests.p2<-merge(nests.p2, nests.p2.c)

# merge phases
nests.final<-rbind(nests.p1, nests.p2)

# proportion
nests.final$proportion<-nests.final$occupancy/nests.final$n*100

# plot it
nests.final.occupancy.plot<-ggplot(nests.final, aes(x=Forest, y=proportion, fill = Stratum)) +
  ggtitle("Nest occupancy") +
  ylab("%")+
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

#  bamboo nest abundance
# Remove empty values
nest.abundance1 <- subset(phase.nests, occupancy!=0)

nest.abundance.stratum<-ggplot(nest.abundance1, aes(x=Forest, y=log(nesting.estimate), fill = Stratum)) +
  ggtitle("Nest size") +
  ylab("ln(number of ants)")+
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
nest.abundance.phase<-ggplot(nest.abundance1, aes(x=Forest, y=log(nesting.estimate), fill = phase)) +
  ggtitle("Nest size") +
  ylab("ln(number of ants)")+
  xlab("")+ 
  #ylim(0,50)+
  geom_point(aes(fill=phase), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1))+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.abundance.phase

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
  ggtitle("Nester diversity") +
  labs(fill = "Phase")+  
  ylab("expH")+
  xlab("")+ 
  geom_point(aes(fill=phase), size=3, shape=21, colour="grey20", position=position_jitterdodge(0.1))+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
nest.diversity.phase.plot

# Plot nester diversity, stratum dependent
nest.diversity.stratum.plot<-ggplot(diversity.nester, aes(x=Forest, y=expH, fill=Stratum)) +
  ggtitle("Nest species diversity") +
  labs(fill = "Stratum")+  
  ylab("expH")+
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

incidence.proportion.nest<-ggplot(nest.proportion3, aes(x = Forest, y = percentage, fill = AntSpCode)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative occupancy [%]") +
  ggtitle("bamboo species composition") +
  #scale_fill_brewer(palette="Dark2")+
  theme_minimal()
incidence.proportion.nest


#### beta partitioning, abundance based

## Species overlap: midelevation vs lowland
# make matrices
bamboo.matrix <- dcast(nesters, formula = Plot ~ AntSpCode, length)

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
std.error(colMeans(b.m_vs_l))

#----------------------------------------------------------#
# 3.3 Bamboo Nest Statistics -----
#----------------------------------------------------------#

# bamboo occupancy models

### Does bamboo cavity size affect occupancy?
cav.mod1<- glmmTMB(occupancy ~  cav.size + (1|Block/Plot),
                   data = phase.nests, 
                   family=binomial)
summary(cav.mod1) # ns
#
testDispersion(cav.mod1) # ok
simulateResiduals(cav.mod1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = cav.mod1)) # ok

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
anova(bamboooccupancy.model.e1, bamboooccupancy.model.e2) # ns  no interactions

summary(bamboooccupancy.model.e1)
#
testDispersion(bamboooccupancy.model.e1) # ok
simulateResiduals(bamboooccupancy.model.e1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboooccupancy.model.e1)) # ok 

#  nest size (abundance)
## Does bamboo cavity size affect nest size?
# 
cav.mod2<- glmmTMB(nesting.estimate ~  cav.size + (1|Block/Plot),
                   data=subset(phase.nests, nesting.estimate!=0),
                   family=nbinom2)
summary(cav.mod2) # ns
#
testDispersion(cav.mod2) # ok
simulateResiduals(cav.mod2, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = cav.mod2)) # ok


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

anova(bamboo.abundance.model.e4, bamboo.abundance.model.e5) # ns

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
simulateResiduals(bamboo.diversity.model1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model1)) # ok

# bamboo diversity changes with environment metadata
bamboo.diversity.model.e1 <- glmmTMB((expH) ~ Forest.x
                                     +Stratum
                                     +scale(log(Lianas.n_mean))
                                     +scale(log(dw.percent_mean))
                                     +scale(log(dw.number_mean))
                                     +scale(log(trunk_mean))
                                     +scale(log(Caco+1))
                                     +scale(log(slope.var+1))
                                     +new.Treatment
                                     +(1|Block.x/Plot.x),
                                     zi=~Forest.x,
                                     data = diversity.nester.e, family=gaussian(link=identity))

bamboo.diversity.model.e2 <- glmmTMB((expH) ~ Forest.x
                                     *Stratum
                                     +scale(log(Lianas.n_mean))
                                     +scale(log(dw.percent_mean))
                                     +scale(log(dw.number_mean))
                                     +scale(log(trunk_mean))
                                     +scale(log(Caco+1))
                                     +scale(log(slope.var+1))
                                     +new.Treatment
                                     +(1|Block.x/Plot.x),
                                     zi=~Forest.x,
                                     data = diversity.nester.e, family=gaussian(link=identity))

anova(bamboo.diversity.model.e1, bamboo.diversity.model.e2) # no interaction better

summary(bamboo.diversity.model.e1)
#
testDispersion(bamboo.diversity.model.e1) # ok
simulateResiduals(bamboo.diversity.model.e1, plot = T) # ok-ish
testZeroInflation(simulateResiduals(fittedModel = bamboo.diversity.model.e1)) # ok

