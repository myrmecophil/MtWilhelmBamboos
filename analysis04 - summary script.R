# Summary Script
# (only run if all other scripts have been run)


## Species overlap between baits and bamboos  

# Bait species overlap by Forest
baits.matrix <- dcast(baiting.data, formula = Forest ~ AntSpCODE, length)
df1<-baits.matrix

# bamboo species overlap by forest
bamboo.matrix <- dcast(nesters, formula = Forest ~ AntSpCode, length)
df2<-bamboo.matrix

# fill in non-overlapping columns with NAs
df1[setdiff(names(df2), names(df1))] <- NA
df2[setdiff(names(df1), names(df2))] <- NA

#merge
bamboo.baits.matrix<-rbind(df1,df2)

# set NAs as 0
bamboo.baits.matrix[is.na(bamboo.baits.matrix)] <- 0

# set rownames
rownames(bamboo.baits.matrix) <- c("Kausi Baits", "Numba Baits", "Kausi Bamboo", "Numba Bamboo")
bamboo.baits.matrix <- bamboo.baits.matrix[, -c(1,2)]

# Bray-Curtis Dissimilarity
dist.kausi.numba.bray <- vegdist(bamboo.baits.matrix, method = "bray")

# global matrix of forest dissimilarities
dist.kausi.numba.bray <- as.matrix(dist.kausi.numba.bray)
dist.kausi.numba.bray

#----------------------------------------------------------#
# Figure summaries -----
#----------------------------------------------------------#

# Main figures
# diversity
nest.diversity.stratum.plot
bait.diversity.stratum

# occupancy
nests.final.occupancy.plot
bait.occupancy.stratum

# log abundance
bait.abundance.stratum
nest.abundance.stratum

# proportional composition
abundance.proportion.bait
abundance.proportion.nest


main_figure <- ggarrange(bait.occupancy.stratum, nests.final.occupancy.plot, bait.abundance.stratum, nest.abundance.stratum, bait.diversity.stratum, nest.diversity.stratum.plot,
                         labels = c("A", "B", "C", "D", "E", "F"),
                         ncol = 2, nrow = 3, common.legend = T
)
main_figure

species_composition <- ggarrange(abundance.proportion.bait, abundance.proportion.nest,
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1, common.legend =F
)
species_composition 



# Supplement figure: Phase divided

supplement_phase_fig <- ggarrange(nest.occupancy.phase.plot, nest.diversity.phase.plot, nest.abundance.phase,
                                  labels = c("A", "B", "C"),
                                  ncol = 3, nrow = 1, common.legend = T
)
supplement_phase_fig




# maybe include in supplement..
obs.sim.model = ggplot()+
  geom_boxplot(data = beta.all.abund, aes(x = Stratum, y = total, colour = Stratum))+
  geom_jitter(data = beta.all.abund, aes(x = Stratum, y = total, colour = Stratum), stroke = 1.2, alpha = 0.6, size = 2)+
  scale_color_manual(values = c("#56B4E9", "#E69F00"))+
  ylab(expression(Species~turnover~-~beta[total]))+ xlab("Stratum")+labs(shape = "")+
  ylim(0.5,1.01)+
  theme_classic(20)+
  theme(legend.position = "none")
obs.sim.model

