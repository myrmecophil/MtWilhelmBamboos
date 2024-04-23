## Numba-Kausi Nest accumulation curves - R script by Phil 01 April 2024


rm(list=ls()) 

#----------------------------------------------------------#
### List of R-packages
#----------------------------------------------------------#


package_list <- 
  c("dplyr",
    "tidyr",
    "reshape2",
    "vegan",
    "ggplot2",
    "ggthemes",
    "ggpubr",
    "iNEXT")

# install all packages
#sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# Citations
#sapply(package_list, citation)

set.seed(1234)

#----------------------------------------------------------#
### Load data from other script
#----------------------------------------------------------#

# work with phase.nests, which contains both phases separately
phase.nests<-read.csv(file="phase.nests.csv", header=T)

# this step removes poison plots KPX-S3
phase.nests<-subset(phase.nests, Plot!= "KP1-S3" & Plot!= "KP2-S3" & Plot!= "KP4-S3"& Plot!= "KP6-S3")

#----------------------------------------------------------#
### bamboo level, ignoring phase accumulation curve
#----------------------------------------------------------#

acc.forest<-dcast(phase.nests,formula = Code ~ AntSpCode, length)
rownames(acc.forest)<-acc.forest[, 1]
acc.forest<-acc.forest[, -c(1,2)]

# transform to incidence
acc.forest <- ifelse(acc.forest > 0, 1, 0)

list.acc <-list() # empty list

list.acc$Kausi<-t(acc.forest[1:128,])
list.acc$Numba<-t(acc.forest[129:639,])

nrow(acc.forest[1:128,])
nrow(acc.forest[129:639,])

out.stage <- iNEXT(list.acc, datatype="incidence_raw", endpoint=511, se=TRUE, conf=0.95, nboot=1000) # endpoint  is highest number on the X axis
rarecurve1<-ggiNEXT(out.stage,se=TRUE)+theme_bw(base_size = 18)+theme(legend.position="right", legend.title=element_blank())+
  xlab("number of samping units")
rarecurve1

#----------------------------------------------------------#
### bamboo level, comparing phases accumulation curve
#----------------------------------------------------------#

acc.phase1<-dcast(subset(phase.nests, phase == "phase 1"),formula = Code ~ AntSpCode, length)
rownames(acc.phase1)<-acc.phase1[, 1]
acc.phase1<-acc.phase1[, -c(1,2)]

acc.phase2<-dcast(subset(phase.nests, phase == "phase 2"),formula = Code ~ AntSpCode, length)
rownames(acc.phase2)<-acc.phase2[, 1]
acc.phase2<-acc.phase2[, -c(1,2)]

list.acc <-list() # empty list

list.acc$Kausi.phase1<-t(acc.phase1[1:127,])
list.acc$Numba.phase1<-t(acc.phase1[128:634,])

list.acc$Kausi.phase2<-t(acc.phase2[1:126,])
list.acc$Numba.phase2<-t(acc.phase2[127:507,])

nrow(acc.phase1[1:127,])
nrow(acc.phase1[128:634,])
nrow(acc.phase2[127:507,])
nrow(acc.phase2[1:126,])

out.stage <- iNEXT(list.acc, datatype="incidence_raw", endpoint=507, se=TRUE, conf=0.95, nboot=1000) # endpoint  is highest number on the X axis
rarecurve2<-ggiNEXT(out.stage,se=TRUE)+theme_bw(base_size = 18)+theme(legend.position="right", legend.title=element_blank())+
  xlab("number of samping units")
rarecurve2


#----------------------------------------------------------#
### bamboo level, accumulation curves with removed unoccupied bamboos
#----------------------------------------------------------#

phase.nests1<-subset(phase.nests, occupancy !="0")

# ignoring phase accumulation curve

acc.forest<-dcast(phase.nests1,formula = Code ~ AntSpCode, length)
rownames(acc.forest)<-acc.forest[, 1]
acc.forest<-acc.forest[, -c(1,2)]

# transform to incidence
acc.forest <- ifelse(acc.forest > 0, 1, 0)

list.acc <-list() # empty list

list.acc$Kausi<-t(acc.forest[1:34,])
list.acc$Numba<-t(acc.forest[35:96,])

nrow(acc.forest[1:34,])
nrow(acc.forest[35:96,])

out.stage <- iNEXT(list.acc, datatype="incidence_raw", endpoint=62, se=TRUE, conf=0.95, nboot=1000) # endpoint  is highest number on the X axis
rarecurve1<-ggiNEXT(out.stage,se=TRUE)+theme_bw(base_size = 18)+theme(legend.position="right", legend.title=element_blank())+
  xlab("number of samping units")
rarecurve1

#----------------------------------------------------------#
### plot level, accumulation curves
#----------------------------------------------------------#

# ignoring phases 
acc.forest<-dcast(subset(phase.nests),formula = Plot ~ AntSpCode, length)
rownames(acc.forest)<-acc.forest[, 1]
acc.forest<-acc.forest[, -c(1,2)]
# transform to incidence
acc.forest <- ifelse(acc.forest > 0, 1, 0)

list.acc <-list() # empty list

list.acc$Kausi<-t(acc.forest[1:8,])
list.acc$Numba<-t(acc.forest[9:28,])

nrow(acc.forest[1:8,])
nrow(acc.forest[9:28,])


out.stage <- iNEXT(list.acc, datatype="incidence_raw", endpoint=20, se=TRUE, conf=0.95, nboot=1000) # endpoint  is highest number on the X axis
rare.plot<-ggiNEXT(out.stage,se=TRUE)+theme_bw(base_size = 18)+theme(legend.position="right", legend.title=element_blank())+
  xlab("number of samping units")
rare.plot

# comparing phase accumulation curve
acc.phase1<-dcast(subset(phase.nests, phase == "phase 1"),formula = Plot ~ AntSpCode, length)
rownames(acc.phase1)<-acc.phase1[, 1]
acc.phase1<-acc.phase1[, -c(1,2)]

acc.phase2<-dcast(subset(phase.nests, phase == "phase 2"),formula = Plot ~ AntSpCode, length)
rownames(acc.phase2)<-acc.phase2[, 1]
acc.phase2<-acc.phase2[, -c(1,2)]

list.acc <-list() # empty list

list.acc$Kausi.phase1<-t(acc.phase1[1:4,])
list.acc$Numba.phase1<-t(acc.phase1[5:20,])

list.acc$Kausi.phase2<-t(acc.phase2[1:8,])
list.acc$Numba.phase2<-t(acc.phase2[9:16,])

nrow(acc.phase1[1:4,])
nrow(acc.phase1[5:20,])
nrow(acc.phase2[1:8,])
nrow(acc.phase2[9:16,])

out.stage <- iNEXT(list.acc, datatype="incidence_raw", endpoint=16, se=TRUE, conf=0.95, nboot=1000) # endpoint  is highest number on the X axis
rarecurve2<-ggiNEXT(out.stage,se=TRUE)+theme_bw(base_size = 18)+theme(legend.position="right", legend.title=element_blank())+
  xlab("number of samping units")
rarecurve2


