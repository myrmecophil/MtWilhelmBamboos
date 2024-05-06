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

# if incidence is 0, put also 0 species, if 2 species on bait, put 2 otherwise 1 
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

labs <- expression("lowland", "midelevation")

# plot it
bait.occupancy.stratum<-ggplot(baits.final, aes(x=Forest, y=proportion, fill = Stratum)) +
  ggtitle("Bait occupancy [%]") +
  ylab("")+
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
  ggtitle("Bait species diversity [expH]") +
  ylab("")+
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

# Plot scale species overlap, all
# make matrices
baits.matrix <- dcast(baiting.data, formula = Plot ~ AntSpCODE, length)
baits.matrix

# set rownames
rownames(baits.matrix) <- baits.matrix[, 1]
baits.matrix <- baits.matrix[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.baits <- vegdist(baits.matrix, method = "bray")

# global matrix of plot dissimiliarities
dist.baits <- as.matrix(dist.baits)

## select subsets
# midelevation vs lowland
m_vs_l <- dist.baits[c(1:12), c(13:32)]

distance.baits.elev <- as.data.frame(rowMeans(m_vs_l))
mean(colMeans(m_vs_l))
std.error(colMeans(m_vs_l))

## Species overlap, canopy only
# make matrices
baits.matrix.ca <- dcast(subset(baiting.data, Stratum=="CA"), formula = Plot ~ AntSpCODE, length)
baits.matrix.ca

# set rownames
rownames(baits.matrix.ca) <- baits.matrix.ca[, 1]
baits.matrix.ca <- baits.matrix.ca[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.baits <- vegdist(baits.matrix.ca, method = "bray")

# global matrix of plot dissimiliarities
dist.baits <- as.matrix(dist.baits)

## select subsets
# midelevation vs lowland
m_vs_l.ca <- dist.baits[c(1:12), c(13:32)]

distance.baits.elev.ca <- as.data.frame(rowMeans(m_vs_l.ca))
mean(colMeans(m_vs_l.ca))
std.error(colMeans(m_vs_l.ca))


## Species overlap, understorey only
# make matrices
baits.matrix.un <- dcast(subset(baiting.data, Stratum=="UN"), formula = Plot ~ AntSpCODE, length)
baits.matrix.un

# set rownames
rownames(baits.matrix.un) <- baits.matrix.un[, 1]
baits.matrix.un <- baits.matrix.un[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.baits <- vegdist(baits.matrix.un, method = "bray")

# global matrix of plot dissimiliarities
dist.baits <- as.matrix(dist.baits)

## select subsets
# midelevation vs lowland
m_vs_l.un <- dist.baits[c(1:12), c(13:32)]

distance.baits.elev.un <- as.data.frame(rowMeans(m_vs_l.un))
mean(colMeans(m_vs_l.un))
std.error(colMeans(m_vs_l.un))

## Species overlap, understorey vs canopy
# in midelevation



#### beta partitioning, abundance based

#  abundance matrix
baits.matrix.un.oc <-baits.matrix.un
baits.matrix.ca.oc <-baits.matrix.ca

beta.sim_un <- as.matrix(beta.pair.abund(baits.matrix.un.oc, index.family = "bray")$beta.bray.bal)
beta.sne_un <- as.matrix(beta.pair.abund(baits.matrix.un.oc, index.family = "bray")$beta.bray.gra)
beta.total_un <- as.matrix(beta.pair.abund(baits.matrix.un.oc, index.family = "bray")$beta.bray)


beta.sim_ca <- as.matrix(beta.pair.abund(baits.matrix.ca.oc, index.family = "bray")$beta.bray.bal)
beta.sne_ca <- as.matrix(beta.pair.abund(baits.matrix.ca.oc, index.family = "bray")$beta.bray.gra)
beta.total_ca <- as.matrix(beta.pair.abund(baits.matrix.ca.oc, index.family = "bray")$beta.bray)

# species replacement midelevation vs lowland, UN
beta.sim_un <- beta.sim_un[c(1:12), c(13:32)]

beta_un.sim.means <- as.data.frame(rowMeans(beta.sim_un))
colnames(beta_un.sim.means)<-"sim"
beta_un.sim.means$Stratum<-'UN'
beta_un.sim.means$plot<-rownames(beta_un.sim.means)
beta_un.sim.means$plot.stratum<-paste(beta_un.sim.means$plot,beta_un.sim.means$Stratum)
mean(beta_un.sim.means$sim)

# species replacement midelevation vs lowland, CA
beta.sim_ca <- beta.sim_ca[c(1:12), c(13:32)]

beta_ca.sim.means <- as.data.frame(rowMeans(beta.sim_ca))
colnames(beta_ca.sim.means)<-"sim"
beta_ca.sim.means$Stratum<-"CA"
beta_ca.sim.means$plot<-rownames(beta_ca.sim.means)
beta_ca.sim.means$plot.stratum<-paste(beta_ca.sim.means$plot,beta_ca.sim.means$Stratum)
mean(beta_ca.sim.means$sim)

# nestedness midelevation vs lowland, UN
beta.sne_un <- beta.sne_un[c(1:12), c(13:32)]

beta_un.sne.means <- as.data.frame(rowMeans(beta.sne_un))
colnames(beta_un.sne.means)<-"sne"
beta_un.sne.means$Stratum<-"UN"
beta_un.sne.means$plot<-rownames(beta_un.sne.means)
beta_un.sne.means$plot.stratum<-paste(beta_un.sne.means$plot,beta_un.sne.means$Stratum)
mean(beta_un.sne.means$sne)

# nestedness midelevation vs lowland, CA
beta.sne_ca <- beta.sne_ca[c(1:12), c(13:32)]

beta_ca.sne.means <- as.data.frame(rowMeans(beta.sne_ca))
colnames(beta_ca.sne.means)<-"sne"
beta_ca.sne.means$Stratum<-"CA"
beta_ca.sne.means$plot<-rownames(beta_ca.sne.means)
beta_ca.sne.means$plot.stratum<-paste(beta_ca.sne.means$plot,beta_ca.sne.means$Stratum)
mean(beta_ca.sne.means$sne)

# total turnover, midelevation vs lowland, UN
beta.total_un <- beta.total_un[c(1:12), c(13:32)]

beta_un.total.means <- as.data.frame(rowMeans(beta.total_un))
colnames(beta_un.total.means)<-"total"
beta_un.total.means$Stratum<-"UN"
beta_un.total.means$plot<-rownames(beta_un.total.means)
beta_un.total.means$plot.stratum<-paste(beta_un.total.means$plot,beta_un.total.means$Stratum)
mean(beta_un.total.means$total)

# total turnover, midelevation vs lowland, CA
beta.total_ca <- beta.total_ca[c(1:12), c(13:32)]

beta_ca.total.means <- as.data.frame(rowMeans(beta.total_ca))
colnames(beta_ca.total.means)<-"total"
beta_ca.total.means$Stratum<-"CA"
beta_ca.total.means$plot<-rownames(beta_ca.total.means)
beta_ca.total.means$plot.stratum<-paste(beta_ca.total.means$plot,beta_ca.total.means$Stratum)
mean(beta_ca.total.means$total)

# combine all values
df<-rbind(beta_ca.sne.means, beta_un.sne.means)
df1<-rbind(beta_ca.sim.means, beta_un.sim.means)
df2<-rbind(beta_ca.total.means, beta_un.total.means)

beta.all<-merge(df, df1, by="plot.stratum")
beta.all<-merge(beta.all, df2, by="plot.stratum")

beta.all.abund<-beta.all[,-c(3,4,6,7)]

# plot abundance based
bray_plot<-ggplot(beta.all.abund, aes(x=Stratum, y=total, colour = Stratum)) +
  ggtitle("Elevational Species Turnover") +
  ylab("Average Bray Curits Distance")+
  xlab("Stratum")+ 
  ylim(0.5,1.01)+
  geom_boxplot(outlier.colour=NA, lwd=1)+
  geom_jitter(data = beta.all.abund, aes(x = Stratum, y = total, colour = Stratum), alpha = 0.5, size = 3)+
  scale_color_manual(values=c("#0072B2", "#E69F00"))+
  guides(colour = "none")+
  theme_minimal(20)
bray_plot

#### Species identity figure

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

abundance.proportion.bait<-ggplot(bait.abundance2, aes(x = Forest, y = percentage, fill = AntSpCODE)) +
  geom_bar(stat = "identity",position= position_fill(reverse = TRUE), color='black') +
  xlab("") +
  labs(fill = "Ant species")+
  scale_x_discrete(labels=labs)+
  ylab("relative abundance [%]") +
  ggtitle("bait species composition") +
  theme_minimal()
abundance.proportion.bait

# summarize bait incidence counts
bait.incidence <- baiting.data %>%
  group_by(Forest, AntSpCODE,) %>%
  summarize(count = n())

# Bait abundance plot
bait.abundance.stratum<-ggplot(baiting.incidence, aes(x=Forest, y=log(Abundance+1), fill = Stratum)) +
  ggtitle("Bait abundance [log+1]") +
  ylab("")+
  xlab("")+ 
  #ylim(0,50)+
  geom_violin(lwd=1)+
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
                                *scale(log(Lianas.n+1))
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
                                  *scale(log(Lianas.n+1))
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
                                  *scale(log(Lianas.n+1))
                                  +scale(log(dw.percent+1))
                                  +scale(log(dw.number+1))
                                  +scale(log(trunk+1))
                                  +scale(log(Caco+1))
                                  +scale(log(slope.var+1))
                                  +(1|Block/Plot),
                                  data=baiting.incidence.e,
                                  family=nbinom1)

anova(baitabundance.model.e1, baitabundance.model.e2) # better with interaction

summary(baitabundance.model.e2)
#
testDispersion(baitabundance.model.e2) # ok
simulateResiduals(baitabundance.model.e2, plot = T) #   ok-ish
testZeroInflation(simulateResiduals(fittedModel = baitabundance.model.e2)) # ok
plot(allEffects(baitabundance.model.e2)) # model visualization

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
                                          *scale(Lianas.n_mean)
                                          +scale(log(slope.var+1))
                                          +scale(log(Caco+1))
                                          +scale(trunk_mean)
                                          +scale(dw.number_mean)
                                          +scale(dw.percent_mean) 
                                          +(1|Block.x/Plot.x),
                                          data = diversity.stratum.e, family = gaussian(link="log"))

anova(baitdiversity.stratum.model.e1, baitdiversity.stratum.model.e2) # no interaction

summary(baitdiversity.stratum.model.e2)
#
testDispersion(baitdiversity.stratum.model.e2) # ok
simulateResiduals(baitdiversity.stratum.model.e2, plot = T) # ok
testZeroInflation(simulateResiduals(fittedModel = baitdiversity.stratum.model.e2)) # ok

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
                                      *scale(Lianas.n_mean)
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
