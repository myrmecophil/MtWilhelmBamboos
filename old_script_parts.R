### OLD unused Script parts



# 1.1 Bait occupancy -----
### Plot bait occupancy as average proportion per plot, not on stratum

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
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  theme_bw()
bait.occupancy.plot

# plot as glm prediction against elevation

#ggplot(baiting.incidence.e, aes(x=elevation, y=occupancy)) + geom_point() + 
#  stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE)


# 1.2 Bait diversity -----

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
data2$evenness <- data2$expH/specnumber(ant_m) # evenness per plot - changed to Hill evenness, earlier defined H/log(specnumber)
str(data2)

# add Plot location
plot_m<-distinct(baiting.data, Plot, .keep_all = TRUE)
diversity<-merge(data2, plot_m)

## Plot it
#Labs
labs <- expression("lowland", "midelevation")

baitdiversity.plot<-ggplot(diversity, aes(x=Forest, y=expH, fill=Forest)) +
  ggtitle("Bait species diversity") +
  geom_boxplot()+
  ylim(0,20)+
  scale_fill_manual(labels=labs, values=c("#0072B2", "#E69F00"))+
  scale_x_discrete(labels=labs)+
  ylab("Species diversity [expH]")+
  xlab("")+ 
  guides(fill="none")+
  theme_bw()
baitdiversity.plot



# 1.3 Bait composition -----

#### Beta partitioning, occurrence based

# make occurrence matrix
baits.matrix.un.oc <-ifelse(baits.matrix.un >0, 1, 0)
baits.matrix.ca.oc <-ifelse(baits.matrix.ca >0, 1, 0)

beta.sim_un <- as.matrix(beta.pair.abund(baits.matrix.un.oc, index.family = "bray")$beta.bray.bal)

beta.sim_un <- as.matrix(beta.pair(baits.matrix.un.oc, index.family = "sorensen")$beta.sim)
beta.sne_un <- as.matrix(beta.pair(baits.matrix.un.oc, index.family = "sorensen")$beta.sne)
beta.total_un <- as.matrix(beta.pair(baits.matrix.un.oc, index.family = "sorensen")$beta.sor)


beta.sim_ca <- as.matrix(beta.pair(baits.matrix.ca.oc, index.family = "sorensen")$beta.sim)
beta.sne_ca <- as.matrix(beta.pair(baits.matrix.ca.oc, index.family = "sorensen")$beta.sne)
beta.total_ca <- as.matrix(beta.pair(baits.matrix.ca.oc, index.family = "sorensen")$beta.sor)

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

beta.all.occurrence<-beta.all[,-c(3,4,6,7)]


# plot occurrence based

soerensen_plot<-ggplot(beta.all.occurrence, aes(x=Stratum, y=total, colour = Stratum)) +
  ggtitle("Elevational Species Turnover") +
  ylab("Soerensen Dissimilarity")+
  xlab("Stratum")+ 
  ylim(0.5,1.01)+
  geom_boxplot(outlier.colour=NA, lwd=1)+
  geom_jitter(data = beta.all.occurrence, aes(x = Stratum, y = total, colour = Stratum), alpha = 0.5, size = 3)+
  scale_color_manual(values=c("#0072B2", "#E69F00"))+
  guides(colour = "none")+
  theme_minimal(20)
soerensen_plot


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

uncorrected.rank.incidence.baits <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
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
uncorrected.rank.incidence.baits

# 2.1 Nest occupancy -----

# Survived as proportion of total amount of nests
nest.raw3 %>%
  count(Forest.Treatment, survived.d) %>%
  ggplot(aes(Forest.Treatment, n, fill = survived.d)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("Treatment")+
  scale_x_discrete(labels=labs1)+
  ylab("% Bamboos translocated") +
  ggtitle("Survived as proportion of total bamboos") +
  scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(0, .25))+
  geom_text(aes(label = n), vjust=1, colour = "black", size = 5)

# Survived as proportion of occupied nests
nest.raw3 %>% filter(Occupied == TRUE)%>%
  count(Forest.Treatment, survived.d) %>%
  ggplot(aes(Forest.Treatment, n, fill = survived.d)) + 
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

### plot species occupancy as proportion per forest and phase
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


## Bait rank abundance (incidence) curves
# forest x species incidence
ant_m<-dcast(phase.nests, formula = Forest ~ AntSpCode, length)
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

uncorrected.rank.incidence.bamboo <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
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
uncorrected.rank.incidence.bamboo


# 2.2 Nest diversity -----

## gg plot it
nest.phase1<-ggplot(phase1.diversity, aes(x=Forest, y=expH, fill=Stratum)) +
  ggtitle("Phase 1 Bamboo nester diversity") +
  geom_boxplot()+
  ylim(0,5)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  ylab("Species diversity")+
  xlab("")+ 
  theme(axis.text.x = element_text(size=12))
nest.phase1

nest.phase2 <-ggplot(phase2.diversity, aes(x=Forest, y=expH, fill=Stratum)) +
  ggtitle("Phase 2 Bamboo nester diversity") +
  geom_boxplot()+
  #ylim(0,5)+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  ylab("Species diversity")+
  xlab("")+ 
  theme(axis.text.x = element_text(size=12))
nest.phase2
#