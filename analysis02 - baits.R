# Bait subscript

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
baiting.data <- read.csv(file="baiting.data.csv", header=T)
baiting.incidence <- read.csv(file="baiting.incidence.csv", header=T)
tree.meta <- read.csv(file="tree.meta.csv", header=T)
plot.meta<- read.csv(file="plot.meta.csv", header=T)
plot.meta2 <- read.csv(file="plot.meta2.csv", header=T)

#----------------------------------------------------------#
# 1.1 Bait occupancy -----
#----------------------------------------------------------#

## count the double occupation of baits (number of species per bait)
bait.double <- baiting.data %>%
  group_by(Forest, Block, Plot, Code) %>%
  summarize(species.per.bait.raw = n())

# to remove baits with no ants, merge with incidence data
bait.double<-merge(bait.double, baiting.incidence)

bait.double <-bait.double %>%
  mutate(species.per.bait = case_when(occupancy == 0 ~ 0,
                                      species.per.bait.raw > 1 ~ 2,
                                      .default = 1)) #

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

# labels
labs <- expression("lowland", "mid-elevation")

# plot it
bait.occupancy.stratum<-ggplot(baits.final, aes(x=Forest, y=proportion, fill = Stratum)) +
  ggtitle("Bait occupancy") +
  ylab("%")+
  xlab("")+ 
  ylim(0,100)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
bait.occupancy.stratum

#----------------------------------------------------------#
# 1.2 Bait diversity -----
#----------------------------------------------------------#

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
data3$evenness <- data3$expH/specnumber(ant_m) # evenness per plot
str(data3)

# add Plot location
plot_m2<-distinct(baiting.data, plot.stratum, .keep_all = TRUE)
diversity.stratum<-merge(data3, plot_m2)
# 

## Plot it
bait.diversity.stratum<-ggplot(diversity.stratum, aes(x=Forest, y=expH, fill = Stratum)) +
  ggtitle("Bait species diversity") +
  ylab("expH")+
  xlab("")+ 
  ylim(0,15)+
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
bait.diversity.stratum

#----------------------------------------------------------#
# 1.3 Bait composition -----
#----------------------------------------------------------#

# Plot scale species overlap
# make occurance matrices
baits.matrix <- dcast(baiting.data, formula = Plot ~ AntSpCODE, length)

# set rownames
rownames(baits.matrix) <- baits.matrix[, 1]
baits.matrix <- baits.matrix[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.baits <- vegdist(baits.matrix, method = "bray")

# global matrix of plot dissimiliarities
dist.baits <- as.matrix(dist.baits)

## select plot distances: midelevation vs lowland
m_vs_l <- dist.baits[c(1:12), c(13:32)]

distance.baits.elev <- as.data.frame(rowMeans(m_vs_l))
mean(colMeans(m_vs_l))
std.error(colMeans(m_vs_l))

#### Species identity figure

# summarize bait incidence counts
bait.incidence<- baiting.data %>%
  group_by(Forest, AntSpCODE) %>%
  summarize(incidence = n())

# Remove empty baits
bait.incidence <- bait.incidence[bait.incidence$AntSpCODE != "", ]

# species with NA abundance have at least 1 worker, so set to 1
bait.incidence$incidence[is.na(bait.incidence$incidence)] <-1

n.kausi<-sum(bait.incidence$incidence[bait.incidence$Forest=="Kausi Primary"])
n.numba<-sum(bait.incidence$incidence[bait.incidence$Forest=="Numba Primary"])

bait.incidence$percentage<-0

bait.incidence$percentage[bait.incidence$Forest=="Kausi Primary"] <-bait.incidence$incidence[bait.incidence$Forest=="Kausi Primary"]/n.kausi*100
bait.incidence$percentage[bait.incidence$Forest=="Numba Primary"] <-bait.incidence$incidence[bait.incidence$Forest=="Numba Primary"]/n.numba*100

# check sums
sum(bait.incidence$percentage[bait.incidence$Forest=="Kausi Primary"])
sum(bait.incidence$percentage[bait.incidence$Forest=="Numba Primary"])

# Define "Rest" group at <5% of total abundance 
bait.incidence$AntSpCODE[bait.incidence$percentage < 3] <- "other species"

bait.incidence$AntSpCODE<-as.factor(bait.incidence$AntSpCODE)
bait.incidence$Forest<-as.factor(bait.incidence$Forest)

# summarize bait percentage counts
bait.incidence<- bait.incidence%>% 
  group_by(Forest, AntSpCODE) %>% 
  summarise(percentage = sum(percentage))

# Plot it

# Change order of levels
bait.incidence$AntSpCODE <- relevel(bait.incidence$AntSpCODE, "CREM 014")
bait.incidence$AntSpCODE <- relevel(bait.incidence$AntSpCODE, "CREM 003")
bait.incidence$AntSpCODE <- relevel(bait.incidence$AntSpCODE, "other species")

incidence.proportion.bait<-ggplot(bait.incidence, aes(x = Forest, y = percentage, fill = AntSpCODE)) +
  geom_bar(stat = "identity",position= position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative incidence [%]") +
  ggtitle("bait species composition") +
  theme_minimal()
incidence.proportion.bait

# summarize bait incidence counts
bait.incidence <- baiting.data %>%
  group_by(Forest, AntSpCODE,) %>%
  summarize(count = n())

# Bait abundance plot
bait.abundance.stratum<-ggplot(baiting.incidence, aes(x=Forest, y=log(Abundance+1), fill = Stratum)) +
  ggtitle("Bait abundance") +
  ylab("ln(number of ants+1)")+
  xlab("")+ 
  geom_point(aes(fill=Stratum), size=3, shape=21, colour="grey20",
             position=position_jitterdodge(0.1), alpha=0.8)+
  geom_boxplot(outlier.color=NA, lwd=1, alpha=0.6)+
  scale_x_discrete(labels=labs)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_minimal(15)
bait.abundance.stratum

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
#
testDispersion(bait.double.model2) # ok
simulateResiduals(bait.double.model2, plot = T) # good
testZeroInflation(simulateResiduals(fittedModel = bait.double.model2)) # ok

# with environmental factors
bait.double.model.e1 <- glmmTMB(species.per.bait~Forest
                                +Stratum
                                +scale(log(Lianas.n+1))
                                +scale(log(dw.percent+1))
                                +scale(log(dw.number+1))
                                +scale(log(trunk+1))
                                +scale(log(Caco+1))
                                +scale(log(slope.var+1))
                                +(1|Block/Plot),
                                data=bait.double.e, family= binomial)

bait.double.model.e2 <- glmmTMB(species.per.bait~Forest
                                *Stratum
                                +scale(log(Lianas.n+1))
                                +scale(log(dw.percent+1))
                                +scale(log(dw.number+1))
                                +scale(log(trunk+1))
                                +scale(log(Caco+1))
                                +scale(log(slope.var+1))
                                +(1|Block/Plot),
                                data=bait.double.e, family= binomial)

anova(bait.double.model.e1, bait.double.model.e2) # ns, without interaction

summary(bait.double.model.e1)
#
testDispersion(bait.double.model.e1) # ok
simulateResiduals(bait.double.model.e1, plot = T) # good
testZeroInflation(simulateResiduals(fittedModel = bait.double.model.e1)) # ok

## Bait occupancy
# get environmental data
baiting.incidence.e<-merge(baiting.incidence, tree.meta)

# Binominal regression
baitoccupancy.model1 <- glmmTMB(occupancy~Forest+Stratum+(1|Block/Plot),
                                data=baiting.incidence.e,
                                family=binomial)
baitoccupancy.model2 <- glmmTMB(occupancy~Forest*Stratum+(1|Block/Plot),
                                data=baiting.incidence.e,
                                family=binomial)
anova(baitoccupancy.model1, baitoccupancy.model2) # interaction model better

summary(baitoccupancy.model2)
#
testDispersion(baitoccupancy.model2) # ok
simulateResiduals(fittedModel = baitoccupancy.model2, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model2)) # ok
plot(allEffects(baitoccupancy.model2)) # model visualization

# with environmental factors
baitoccupancy.model.e1 <- glmmTMB(occupancy~Forest
                                  +Stratum
                                  +scale(log(Lianas.n+1))
                                  +scale(log(dw.percent+1))
                                  +scale(log(dw.number+1))
                                  +scale(log(trunk+1))
                                  +scale(log(Caco+1))
                                  +scale(log(slope.var+1))
                                  +(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=binomial)

baitoccupancy.model.e2 <- glmmTMB(occupancy~Forest
                                  *Stratum
                                  +scale(log(Lianas.n+1))
                                  +scale(log(dw.percent+1))
                                  +scale(log(dw.number+1))
                                  +scale(log(trunk+1))
                                  +scale(log(Caco+1))
                                  +scale(log(slope.var+1))
                                  +(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=binomial)

anova(baitoccupancy.model.e1, baitoccupancy.model.e2) # interaction model better

summary(baitoccupancy.model.e2)
#
testDispersion(baitoccupancy.model.e2) # ok
simulateResiduals(baitoccupancy.model.e2, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = baitoccupancy.model.e2)) # ok
plot(allEffects(baitoccupancy.model.e2)) # model visualization

# with abundance 
baitabundance.model1 <- glmmTMB(Abundance~Forest+Stratum+(1|Block/Plot),
                                zi=~Forest*Stratum,
                                data=baiting.incidence.e,
                                family=nbinom1)
baitabundance.model2 <- glmmTMB(Abundance~Forest*Stratum+(1|Block/Plot),
                                zi=~Forest*Stratum,
                                data=baiting.incidence.e,
                                family=nbinom1)
anova(baitabundance.model1, baitabundance.model2) # interaction not better

summary(baitabundance.model1)
#
testDispersion(baitabundance.model1) # ok
simulateResiduals(baitabundance.model1, plot = T) # ok-ish
testZeroInflation(simulateResiduals(fittedModel = baitabundance.model1)) # ok 
plot(allEffects(baitabundance.model1)) # model visualization
# 
# with abundance and environmental
baitabundance.model.e1 <- glmmTMB(Abundance~Forest
                                  +Stratum
                                  +scale(log(Lianas.n+1))
                                  +scale(log(dw.percent+1))
                                  +scale(log(dw.number+1))
                                  +scale(log(trunk+1))
                                  +scale(log(Caco+1))
                                  +scale(log(slope.var+1))
                                  +(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=nbinom1)

baitabundance.model.e2 <- glmmTMB(Abundance~Forest
                                  *Stratum
                                  +scale(log(Lianas.n+1))
                                  +scale(log(dw.percent+1))
                                  +scale(log(dw.number+1))
                                  +scale(log(trunk+1))
                                  +scale(log(Caco+1))
                                  +scale(log(slope.var+1))
                                  +(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=nbinom1)

anova(baitabundance.model.e1, baitabundance.model.e2) # no interaction

summary(baitabundance.model.e1)
#
testDispersion(baitabundance.model.e1) # ok
simulateResiduals(baitabundance.model.e1, plot = T) #   ok-ish
testZeroInflation(simulateResiduals(fittedModel = baitabundance.model.e1)) # ok

#####  Baiting diversity 

### Bait diversity incl. stratum
# NOTE: Here, we look at stratum-level diversity, i.e. twice for each plot (understorey+canopy)

# get environmental data: plot.meta2 has averages for variables on Stratum-level (2 strata per plot)
diversity.stratum.e<-merge(diversity.stratum, plot.meta2, by.x = 'plot.stratum', by.y = 'plot.stratum', all.x = T)

# Fit a LMM for Shannon diversity. Needs also plot as random factor since strata are not independent
baitdiversity.stratum.model1 <- glmmTMB(expH ~ Forest.x + Stratum+ (1|Block.x/Plot.x), data = diversity.stratum.e, family = gaussian)
baitdiversity.stratum.model2 <- glmmTMB(expH ~ Forest.x * Stratum+ (1|Block.x/Plot.x), data = diversity.stratum.e, family = gaussian)

anova(baitdiversity.stratum.model1, baitdiversity.stratum.model2) # ns, no interaction

summary(baitdiversity.stratum.model1)
#
testDispersion(baitdiversity.stratum.model1) # ok
simulateResiduals(baitdiversity.stratum.model1, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.stratum.model1)) # ok

# add environmental factors
baitdiversity.stratum.model.e1 <- glmmTMB(expH ~Forest.x
                                          +Stratum
                                          +scale(Lianas.n_mean)
                                          +scale(log(slope.var+1))
                                          +scale(log(Caco+1))
                                          +scale(trunk_mean)
                                          +scale(dw.number_mean)
                                          +scale(dw.percent_mean) 
                                          +(1|Block.x/Plot.x),
                                          data = diversity.stratum.e, family = gaussian(link="log"))

baitdiversity.stratum.model.e2 <- glmmTMB(expH ~Forest.x
                                          *Stratum
                                          +scale(Lianas.n_mean)
                                          +scale(log(slope.var+1))
                                          +scale(log(Caco+1))
                                          +scale(trunk_mean)
                                          +scale(dw.number_mean)
                                          +scale(dw.percent_mean) 
                                          +(1|Block.x/Plot.x),
                                          data = diversity.stratum.e, family = gaussian(link="log"))

anova(baitdiversity.stratum.model.e1, baitdiversity.stratum.model.e2) # no interaction

summary(baitdiversity.stratum.model.e1)
#
testDispersion(baitdiversity.stratum.model.e1) # ok
simulateResiduals(baitdiversity.stratum.model.e1, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.stratum.model.e1)) # ok

# evenness
baitdiversity.model.eve1 <- glmmTMB(evenness ~ Forest.x + Stratum+ (1|Block.x/Plot.x), data = diversity.stratum.e, family = gaussian)
baitdiversity.model.eve2 <- glmmTMB(evenness ~ Forest.x * Stratum+ (1|Block.x/Plot.x), data = diversity.stratum.e, family = gaussian)
anova(baitdiversity.model.eve1, baitdiversity.model.eve2) # no interaction

summary(baitdiversity.model.eve1)
#
testDispersion(baitdiversity.model.eve1) # ok
simulateResiduals(baitdiversity.model.eve1, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.eve1)) # ok

# evenness with environmental factors
baitdiversity.model.eve.e1 <- glmmTMB(evenness ~Forest.x
                                      +Stratum
                                      +scale(Lianas.n_mean)
                                      +scale(log(slope.var+1))
                                      +scale(log(Caco+1))
                                      +scale(trunk_mean)
                                      +scale(dw.number_mean)
                                      +scale(dw.percent_mean) 
                                      +(1|Block.x/Plot.x),
                                      data = diversity.stratum.e, family = gaussian)

baitdiversity.model.eve.e2 <- glmmTMB(evenness ~Forest.x
                                      *Stratum
                                      +scale(Lianas.n_mean)
                                      +scale(log(slope.var+1))
                                      +scale(log(Caco+1))
                                      +scale(trunk_mean)
                                      +scale(dw.number_mean)
                                      +scale(dw.percent_mean) 
                                      +(1|Block.x/Plot.x),
                                      data = diversity.stratum.e, family = gaussian)

anova(baitdiversity.model.eve.e1, baitdiversity.model.eve.e2) #  interaction

summary(baitdiversity.model.eve.e1)
#
testDispersion(baitdiversity.model.eve.e1) # ok
simulateResiduals(baitdiversity.model.eve.e1, plot = T) # ok 
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.model.eve.e1)) # ok

### Bait species beta diversity, abundance-based

# NOTE: Here, we look at stratum-level, i.e. twice for each plot (understory+canopy)

# get environmental data: plot.meta2 has averages for variables on Stratum-level (2 strata per plot)
beta.stratum.a<-merge(beta.all.abund, plot.meta2, by.x = 'plot.stratum', by.y = 'plot.stratum', all.x = T)

# Fit a model for species turnover (Beta sim). 
bait.sim.model1 <- glmmTMB(sim ~ Stratum+ (1|Block/Plot), data = beta.stratum.a, family = gaussian)

summary(bait.sim.model1) #
#
testDispersion(bait.sim.model1) # ok
simulateResiduals(bait.sim.model1, plot = T) # some heterogeneity in variance, but ok
testZeroInflation(simulateResiduals(fittedModel = bait.sim.model1)) # ok

# Fit a model for nestedness
bait.sne.model1 <- glmmTMB((sne+1) ~ Stratum+(1|Block/Plot), data = beta.stratum.a, family = gaussian(link="log"))

summary(bait.sne.model1) 
#
testDispersion(bait.sne.model1) # ok
simulateResiduals(bait.sne.model1, plot = T) # some heterogeneity in variance, but ok
testZeroInflation(simulateResiduals(fittedModel = bait.sne.model1)) # ok

# Fit a model for total turnover (Soerensen)
bait.total.model1 <- glmmTMB(total ~ Stratum+ (1|Block/Plot), data = beta.stratum.a, family=gaussian)

summary(bait.total.model1) # 

testDispersion(bait.total.model1) # ok
simulateResiduals(bait.total.model1, plot = T) #some heterogeneity in variance, but ok
testZeroInflation(simulateResiduals(fittedModel = bait.total.model1)) # ok
