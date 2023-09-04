##################################################################
#            modeling INTERpopulation spatial structures
#               Phan Biogeo 2019 Project
#  #(THIS SCRIPT FOLLOWS the script called Rarefy_Shannon_dissimilarities_final.R)
##################################################################
###### CHECKING MAY 30, 2022 #######
# Description: This script performs distance-based Moran's eigenvector map models. Specifically, it
# separately examines potential interpopulation structure separately in P. vindex, P. difformis, 
# and all beetles (i.e. both species) together. Importantly, we only ended up using the interpopulation
# spatial structure for analyses using ALL beetles together in the manuscript.

# The lat/long of each population is the centroid
# of all lat/longs for each trap where beetles were caught (calculated in get_centroid_data.R). However,
# for purposes of protecting landowners' privacy, we only reported latitude and longitude to the nearest
# tenth of a degree (available in Table S1 in manuscripts's Supplementary Materials). This can be used 
# for re-creating the results in the publication. 

# Most of the work here was done following examples in Borcard, Gillet, and Legendre's Numerical Ecology 
# with R, chapter 7. 

setwd("/Users/clairewinfrey/Desktop/UTK_Research/copyOf_All_pops")

set.seed(93)

# LOAD REQUIRED PACKAGES
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(phyloseq)
library(SoDA)
library(lme4)

# LOAD OBJECTS MADE IN R DOCUMENT "setUp_RarefyShannon.R"
load("savedObjectsMarch27_2023/ps_rare_diss_shann_.RData")
#load("PV_vars.RData")
#load("gut_coords.RData")

# FILES CREATED IN THIS SCRIPT #
# load("among.pop.allguts.dbMEM.RData") # this RData file was created in this script
# load(file = "among.pop.PV.dbMEM.RData") # this RData file was created in this script
#load(file = "among.pop.PD.dbMEM.RData") # this RData file was created in this script

################## GET COORDINATES BY CENTROID OF POPULATION ###################
##### P. VINDEX
# 1. LOAD CENTROID DATA (ALREADY MADE GEODETIC CARTESIAN WITH GEOXY() )
allgutscoordxy.file <- read.csv(file = "allgutscoordxy.csv") #MADE IS get_centroid_data.R
allgutscoordxy.file

# FILTER BASED ON P. VINDEX VERSUS P. DIFFORMIS
#P. vindex
PV_coord <- dplyr::filter(allgutscoordxy.file, Species =="P.vindex")
dim(PV_coord)
class(PV_coord)

PV_centr_coord <- PV_coord[,7:8]
rownames(PV_centr_coord) <-PV_coord[,1]
colnames(PV_centr_coord)

##### P. DIFFORMIS
PD_coord <- dplyr::filter(allgutscoordxy.file, Species =="P.difformis")
dim(PD_coord)
colnames(PD_coord)

PD_centr_coord <- PD_coord[,9:10]
rownames(PD_centr_coord) <-PD_coord[,1]

##### ALL GUTS 

guts_centr_coord <- allgutscoordxy.file[,5:6] # this is the centroid for each population
rownames(guts_centr_coord) <- allgutscoordxy.file[,1]

########################## P. VINDEX ##########################

######## CONSTRUCT DBMEMS ON CENTROID COORDS OF P. VINDEX #########

###### POSITIVE DBMEMS ######
# CONSTRUCT POSITIVE ONLY DBMEM MODEL ON COORDS 
PV.AmPop.pos.dbmem <- dbmem(PV_centr_coord, silent = FALSE, MEM.autocor = "positive") # non-null displays both positive and negative

# DISPLAY EIGENVALUES
attributes(PV.AmPop.pos.dbmem)$values # only 3 positive: ***** 0.18887325 0.10566015 0.01454909

##### NEGATIVE DBMEMS ######
# CONSTRUCT NEGATIVE ONLY DBMEM MODEL ON COORDS 
PV.AmPop.neg.dbmem <- dbmem(PV_centr_coord, silent = FALSE, MEM.autocor = "negative") 
PV.AmPop.dbmem.thresh <- give.thresh(dist(PV_centr_coord))
PV.AmPop.dbmem.thresh #337.2224

# DISPLAY EIGENVALUES
attributes(PV.AmPop.neg.dbmem)$values 
length(PV.AmPop.neg.dbmem) #85


########## JACCARD #########

# 1. CHECK FOR LINEAR TREND AND DETREND DATA
set.seed(93)
PV.jacc.AmPop.detrend.mod <-dbrda(PV.jacc.dist~PV_centr_coord[,1]+PV_centr_coord[,2])
PV.jacc.AmPop.detrend.mod.results <- anova.cca(PV.jacc.AmPop.detrend.mod, permutations=9999) # significant
# March 2, 2023:  1e-04 ***
# ***** April 3, 2023: 1.9153  1e-04 *** (same as before, I think)

PV.jacc.AmPop.dt.dat<-resid(PV.jacc.AmPop.detrend.mod) #save residuals
PV.jacc.AmPop.dt.dat

# RUN POS DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
PV.AmPop.jacc.pos.dbmem.dbrda <- dbrda(PV.jacc.AmPop.dt.dat ~ as.matrix(PV.AmPop.pos.dbmem))
PV.AmPop.jacc.pos.dbmem.dbrda.results <- anova.cca(PV.AmPop.jacc.pos.dbmem.dbrda, permutations=9999)
PV.AmPop.jacc.pos.dbmem.dbrda.results # NOT SIGNIFICANT, so we don't continue with this
# March 2, 2023: 0.9912
# *****APRIL 3, 2023: 0.7651 0.9912 NOT SIGNIFICANT

# RUN NEG DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
PV.AmPop.jacc.neg.dbmem.dbrda <- dbrda(PV.jacc.AmPop.dt.dat ~ as.matrix(PV.AmPop.neg.dbmem))
PV.AmPop.jacc.neg.dbmem.dbrda.results <- anova.cca(PV.AmPop.jacc.neg.dbmem.dbrda, permutations=9999)
PV.AmPop.jacc.neg.dbmem.dbrda.results # ** QUITE SIGNIFICANT!
# March 2, 2023: 0.0089 **
# *****APRIL 3, 2023: .307 0.0089 **

# FORWARD MODEL SELECTION USING NEGATIVE EIGENVECTORS

# FIRST MAKE A MODEL WITH INTERCEPT ONLY
set.seed(93)
PV.AmPop.jacc.neg.dbmem.dbrda.0 <- dbrda(PV.jacc.AmPop.dt.dat ~ 1, data = PV.AmPop.neg.dbmem)

# MAKE A MODEL WITH ALL POTENTIAL MEMS.
# If I didn't make this super explicit using ~. (i.e. if I fed PV.AmPop.neg.dbmem.dbrda into "scope" of ordiR2step, it wouldn't
# do model sel across each MEM individually)
set.seed(93)
PV.AmPop.jacc.neg.dbmem.dbrda.all <- dbrda(PV.jacc.AmPop.dt.dat ~ ., data = PV.AmPop.neg.dbmem)


# MODEL SELECTION
set.seed(93)
PV.AmPop.jacc.neg.forsel <- ordiR2step(PV.AmPop.jacc.neg.dbmem.dbrda.0, scope = PV.AmPop.jacc.neg.dbmem.dbrda.all, permutations =9999)
PV.AmPop.jacc.neg.forsel.results <- PV.AmPop.jacc.neg.forsel$anova # significant MEMs and thus will go into big P. vindex Jaccard model
# # *****April 3, 2023:                 R2.adj Df    AIC      F  Pr(>F)   
# R2.adj Df    AIC      F Pr(>F)    
# + MEM85         0.010093  1 317.83 1.8972 0.0029 ** 
#   + MEM10         0.020151  1 317.89 1.8930 0.0008 ***
#   + MEM52         0.029885  1 317.96 1.8629 0.0015 ** 
#   + MEM29         0.037963  1 318.17 1.7137 0.0051 ** 
#   + MEM59         0.044920  1 318.45 1.6119 0.0094 ** 
#   + MEM7          0.051455  1 318.76 1.5718 0.0093 ** 
#   + MEM30         0.057360  1 319.12 1.5137 0.0171 *  
#   + MEM84         0.063255  1 319.45 1.5097 0.0157 *  
#   + MEM53         0.068466  1 319.84 1.4475 0.0302 *  
#   + MEM51         0.073725  1 320.20 1.4485 0.0245 *  
#   + MEM70         0.079068  1 320.54 1.4525 0.0267 *  
#   <All variables> 0.228707                            
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                    

# APPLY SIDAK CORRECTION
library(mutoss)
PV.jacc.AmPop.final <- sidak(pValues= PV.AmPop.jacc.neg.forsel$anova$`Pr(>F)`)
#  0.005539739 0.041035538
# #***** April 3, 2023 :  [1] 0.034250271 0.009557872 0.017852240 0.059512191 0.107147161 0.106064970 0.186959797 0.172953703 0.307872417 0.257447295
# [11] 0.277295630          NA


########## (PV) wUniFrac #########

# 1. CHECK FOR LINEAR TREND AND DETREND DATA
set.seed(93)
PV.wUF.AmPop.detrend.mod <-dbrda(PV.wUF.dist~PV_centr_coord[,1]+PV_centr_coord[,2])
PV.wUF.AmPop.detrend.mod.results <- anova.cca(PV.wUF.AmPop.detrend.mod, permutations=9999) #3e-04 ***
# April 3,2023 *****:3.4452  3e-04 ***

PV.wUF.AmPop.dt.dat<-resid(PV.wUF.AmPop.detrend.mod)
PV.wUF.AmPop.dt.dat


# RUN POS DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
PV.AmPop.wUF.pos.dbmem.dbrda <- dbrda(PV.wUF.AmPop.dt.dat ~ as.matrix(PV.AmPop.pos.dbmem))
PV.AmPop.wUF.pos.dbmem.dbrda.results <- anova.cca(PV.AmPop.wUF.pos.dbmem.dbrda, permutations=9999)
PV.AmPop.wUF.pos.dbmem.dbrda.results # NOT SIGNIFICANT, so we don't continue with this
# APRIL 3, 2023: 
# ***** 0.8702 0.6025

# RUN NEG DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
PV.AmPop.wUF.neg.dbmem.dbrda <- dbrda(PV.wUF.AmPop.dt.dat ~ as.matrix(PV.AmPop.neg.dbmem))
PV.AmPop.wUF.neg.dbmem.dbrda.results <- anova.cca(PV.AmPop.wUF.neg.dbmem.dbrda, permutations=9999)
PV.AmPop.wUF.neg.dbmem.dbrda.results # Pr(>F)= 3e-04 *** QUITE SIGNIFICANT!
# ***** April 3, 2023 : NOT SIGNIFICANT: 0.3978

# ***** # ALL OF THIS IS GRAYED OUT SINCE NEGATIVE NO LONGER SIGNIFICANT!! *****
# # FORWARD MODEL SELECTION USING NEGATIVE EIGENVECTORS

# 
# # FIRST MAKE A MODEL WITH INTERCEPT ONLY
# set.seed(93)
# PV.AmPop.wUF.neg.dbmem.dbrda.0 <- dbrda(PV.wUF.AmPop.dt.dat ~ 1, data = PV.AmPop.neg.dbmem)
# 
# # MAKE A MODEL WITH ALL POTENTIAL MEMS.
# # If I didn't make this super explicit using ~. (i.e. if I fed PV.AmPop.neg.dbmem.dbrda into "scope" of ordiR2step, it wouldn't
# # do model sel across each MEM individually)
# set.seed(93)
# PV.AmPop.wUF.neg.dbmem.dbrda.all <- dbrda(PV.wUF.AmPop.dt.dat ~ ., data = PV.AmPop.neg.dbmem)
# 
# 
# # MODEL SELECTION
# set.seed(93)
# PV.AmPop.wUF.neg.forsel <- ordiR2step(PV.AmPop.wUF.neg.dbmem.dbrda.0, scope = PV.AmPop.wUF.neg.dbmem.dbrda.all, permutations =9999)
# PV.AmPop.wUF.neg.forsel.results <- PV.AmPop.wUF.neg.forsel$anova # MEM7, MEM4, and MEM8 are significant KEEP THESE?
# # March 2, 2023 :
# # > PV.AmPop.wUF.neg.forsel.results
# # R2.adj Df    AIC      F Pr(>F)  
# # + MEM84         0.019758  1 135.71 2.7737 0.0145 *
# #   + MEM51         0.034220  1 135.36 2.3028 0.0354 *
# #   + MEM52         0.048994  1 134.95 2.3360 0.0358 *
# #   + MEM59         0.063866  1 134.49 2.3504 0.0279 *
# #   + MEM50         0.077603  1 134.11 2.2510 0.0378 *
# #   <All variables> 0.125928                 

# # APPLY SIDAK CORRECTION
# library(mutoss)
# sidak(pValues= PV.AmPop.neg.forsel$anova$`Pr(>F)`)
# PV.wUF.AmPop.final <-sidak(pValues= PV.AmPop.wUF.neg.forsel$anova$`Pr(>F)`) # After corrrection: 0.1176970 0.1168227 0.1757102 NONE ARE SIGNIFICANT
# # March 3, 2023:
# # > PV.wUF.AmPop.final
# # $adjPValues
# # [1] 0.08390656 0.19446661 0.19646877 0.15614921 0.20641744         NA
# 
# # ***** # ALL OF THIS IS GRAYED OUT SINCE NEGATIVE NO LONGER SIGNIFICANT!! *****
# # ***** CAN DELETE ABOVE THIS LINE

########################## P. DIFFORMIS ##########################

######## CONSTRUCT DBMEMS ON CENTROID COORDS OF P. DIFFORMIS #########

###### POSITIVE DBMEMS ######
# CONSTRUCT POSITIVE ONLY DBMEM MODEL ON COORDS 
PD.AmPop.pos.dbmem <- dbmem(PD_centr_coord, silent = FALSE, MEM.autocor = "positive") # non-null displays both positive and negative
PD.AmPop.pos.dbmem.thresh <- give.thresh(dist(PD_centr_coord))
PD.AmPop.pos.dbmem.thresh #211.6707

# DISPLAY EIGENVALUES
attributes(PD.AmPop.pos.dbmem)$values

##### NEGATIVE DBMEMS ######
# CONSTRUCT NEGATIVE ONLY DBMEM MODEL ON COORDS 
PD.AmPop.neg.dbmem <- dbmem(PD_centr_coord, silent = FALSE, MEM.autocor = "negative") 

# DISPLAY EIGENVALUES
attributes(PD.AmPop.neg.dbmem)$values 

########## JACCARD #########

# 1. CHECK FOR LINEAR TREND AND DETREND DATA
set.seed(93)
PD.jacc.AmPop.detrend.mod <- dbrda(PD.jacc.dist~PD_centr_coord[,1]+PD_centr_coord[,2])
PD.jacc.AmPop.detrend.mod.results <- anova.cca(PD.jacc.AmPop.detrend.mod, permutations=9999) # Pr(>F) = 0.00042 ***
# ****** April 3, 2023:  3e-04 ***

PD.jacc.AmPop.dt.dat<-resid(PD.jacc.AmPop.detrend.mod)
PD.jacc.AmPop.dt.dat

# RUN POS DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
PD.AmPop.jacc.pos.dbmem.dbrda <- dbrda(PD.jacc.AmPop.dt.dat ~ as.matrix(PD.AmPop.pos.dbmem))
PD.AmPop.jacc.pos.dbmem.dbrda.results <- anova.cca(PD.AmPop.jacc.pos.dbmem.dbrda, permutations=9999)
PD.AmPop.jacc.pos.dbmem.dbrda.results #  Pr(>F)= 0.2014 NOT SIGNIFICANT
# ***** April 3, 2023: SIGNIFICANT 0.0062 **

# RUN NEG DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
PD.AmPop.jacc.neg.dbmem.dbrda <- dbrda(PD.jacc.AmPop.dt.dat ~ as.matrix(PD.AmPop.neg.dbmem))
PD.AmPop.jacc.neg.dbmem.dbrda.results <- anova.cca(PD.AmPop.jacc.neg.dbmem.dbrda, permutations=9999)
PD.AmPop.jacc.neg.dbmem.dbrda.results # OLD RESULTS: Pr(>F)= 1e-05 *** VERY SIGNIFICANT!
# ****** And april 4, 2023 MARCH 3, 2023: NOT SIGNIFICANT: 0.9939

# ***** NEW MOD SELECTION WITH POSITIVE
# FIRST MAKE A MODEL WITH INTERCEPT ONLY
set.seed(93)
PD.AmPop.jacc.pos.dbmem.dbrda.0 <- dbrda(PD.jacc.AmPop.dt.dat ~ 1, data = PD.AmPop.pos.dbmem)

# MAKE A MODEL WITH ALL POTENTIAL MEMS.
# If I didn't make this super explicit using ~. (i.e. if I fed PD.AmPop.neg.dbmem.dbrda into "scope" of ordiR2step, it wouldn't
# do model sel across each MEM individually)
set.seed(93)
PD.AmPop.jacc.pos.dbmem.dbrda.all <- dbrda(PD.jacc.AmPop.dt.dat ~ ., data = PD.AmPop.pos.dbmem)

# MODEL SELECTION
set.seed(93)
PD.AmPop.jacc.pos.forsel <- ordiR2step(PD.AmPop.jacc.pos.dbmem.dbrda.0, scope = PD.AmPop.jacc.pos.dbmem.dbrda.all, permutations =9999)
PD.AmPop.jacc.pos.forsel.results <- PD.AmPop.jacc.pos.forsel$anova 
PD.AmPop.jacc.pos.forsel.results #*****  R2.adj Df    AIC      F Pr(>F)   
# + MEM5          0.0088786  1 413.78 1.9764 0.0012 **
#   <All variables> 0.0126989                           

# APPLY SIDAK CORRECTION
library(mutoss)
PD.jacc.AmPop.final <- sidak(pValues= PD.AmPop.jacc.neg.forsel$anova$`Pr(>F)`)
PD.jacc.AmPop.final # three out of the four are significant

# # FORWARD MODEL SELECTION USING NEGATIVE EIGENVECTORS
# #***** GRAYED OUT SINCE NO LONGER SIGNIFICANT
# # FIRST MAKE A MODEL WITH INTERCEPT ONLY
# set.seed(93)
# PD.AmPop.jacc.neg.dbmem.dbrda.0 <- dbrda(PD.jacc.AmPop.dt.dat ~ 1, data = PD.AmPop.neg.dbmem)
# 
# # MAKE A MODEL WITH ALL POTENTIAL MEMS.
# # If I didn't make this super explicit using ~. (i.e. if I fed PD.AmPop.neg.dbmem.dbrda into "scope" of ordiR2step, it wouldn't
# # do model sel across each MEM individually)
# set.seed(93)
# PD.AmPop.jacc.neg.dbmem.dbrda.all <- dbrda(PD.jacc.AmPop.dt.dat ~ ., data = PD.AmPop.neg.dbmem)
# 
# # MODEL SELECTION
# set.seed(93)
# PD.AmPop.jacc.neg.forsel <- ordiR2step(PD.AmPop.jacc.neg.dbmem.dbrda.0, scope = PD.AmPop.jacc.neg.dbmem.dbrda.all, permutations =9999)
# PD.AmPop.jacc.neg.forsel.results <- PD.AmPop.jacc.neg.forsel$anova # MEM2, 4, 6, and 3. r2 for all is 0.0265738
# 
# # APPLY SIDAK CORRECTION
# library(mutoss)
# PD.jacc.AmPop.final <- sidak(pValues= PD.AmPop.jacc.neg.forsel$anova$`Pr(>F)`)
# PD.jacc.AmPop.final # three out of the four are significant
# ***** GRAYED OUT SINCE NO LONGER SIGNIFICANT

########## (PD) wUniFrac #########

# 1. CHECK FOR LINEAR TREND AND DETREND DATA
set.seed(93)
PD.wUF.AmPop.detrend.mod <-dbrda(PD.wUF.dist~PD_centr_coord[,1]+PD_centr_coord[,2])
PD.wUF.AmPop.detrend.mod.results <- anova.cca(PD.wUF.AmPop.detrend.mod, permutations=9999) # Pr(>F) = 0.1516 NOT SIGNIFICANT
# April 3, 2023: STILL NOT SIGNFICANT: 1.4219 0.1498

##### RAN WITH DETRENDED BELOW, HASEHD IT OUT AND RAN ON REGULAR DATA BECAUSE LINEAR TREND NOT SIGNIFICANT
#PD.wUF.AmPop.dt.dat<-resid(PD.wUF.AmPop.detrend.mod)
#PD.wUF.AmPop.dt.dat

## RUN POS DBMEM ANALYSIS ON DETRENDED DATA
#set.seed(93)
#PD.AmPop.wUF.pos.dbmem.dbrda.dt <- dbrda(PD.wUF.AmPop.dt.dat ~ as.matrix(PD.AmPop.pos.dbmem))
#PD.AmPop.wUF.pos.dbmem.dbrda.dt.results <- anova.cca(PD.AmPop.wUF.pos.dbmem.dbrda.dt, permutations=9999)
#PD.AmPop.wUF.pos.dbmem.dbrda.dt.results #  Pr(>F)= 0.064 NOT SIGNIFICANT, BUT CLOSE

# RUN NEG DBMEM ANALYSIS ON DETRENDED DATA
#set.seed(93)
#PD.AmPop.wUF.neg.dbmem.dbrda.dt <- dbrda(PD.wUF.AmPop.dt.dat ~ as.matrix(PD.AmPop.neg.dbmem))
#PD.AmPop.wUF.neg.dbmem.dbrda.dt.results <- anova.cca(PD.AmPop.wUF.neg.dbmem.dbrda.dt, permutations=9999)
#PD.AmPop.wUF.neg.dbmem.dbrda.dt.results # Pr(>F)= 0.3284, NOT SIGNIFICANT


# FORWARD MODEL SELECTION USING NEGATIVE EIGENVECTORS

# FIRST MAKE A MODEL WITH INTERCEPT ONLY
#set.seed(93)
#PD.AmPop.wUF.neg.dbmem.dbrda.dt.0 <- dbrda(PD.wUF.AmPop.dt.dat ~ 1, data = PD.AmPop.neg.dbmem)

# MAKE A MODEL WITH ALL POTENTIAL MEMS.
# If I didn't make this super explicit using ~. (i.e. if I fed PD.AmPop.neg.dbmem.dbrda into "scope" of ordiR2step, it wouldn't
# do model sel across each MEM individually)
#set.seed(93)
#PD.AmPop.wUF.neg.dbmem.dbrda.dt.all <- dbrda(PD.wUF.AmPop.dt.dat ~ ., data = PD.AmPop.neg.dbmem)

# MODEL SELECTION
#set.seed(93)
#PD.AmPop.wUF.dt.neg.forsel <- ordiR2step(PD.AmPop.wUF.neg.dbmem.dbrda.0, scope = PD.AmPop.wUF.neg.dbmem.dbrda.all, permutations =9999)
#PD.AmPop.wUF.dt.neg.forsel.results <- PD.AmPop.wUF.neg.forsel$anova # NULL SO DID NOT APPLY SIDAK

# RUN POS DBMEM ANALYSIS ON full DATA
set.seed(93)
PD.AmPop.wUF.pos.dbmem.Ndet.dbrda <- dbrda(PD.wUF.dist ~ as.matrix(PD.AmPop.pos.dbmem))
PD.AmPop.wUF.pos.dbmem.Ndet.dbrda.results <- anova.cca(PD.AmPop.wUF.pos.dbmem.Ndet.dbrda, permutations=9999)
PD.AmPop.wUF.pos.dbmem.Ndet.dbrda.results #  Pr(>F)= 0.00375 ** SIGNIFICANT
# MARCH 3: same result, 0.0036 **

# RUN NEG DBMEM ANALYSIS ON full DATA
set.seed(93)
PD.AmPop.wUF.neg.dbmem.Ndet.dbrda <- dbrda(PD.wUF.dist ~ as.matrix(PD.AmPop.neg.dbmem))
PD.AmPop.wUF.neg.dbmem.Ndet.dbrda.results <- anova.cca(PD.AmPop.wUF.neg.dbmem.Ndet.dbrda, permutations=9999)
PD.AmPop.wUF.neg.dbmem.Ndet.dbrda.results # Pr(>F)=  0.25 INSIGNIFICANT
# MARCH 3, 2023: STILL INSIGNFICANT: 0.9965

# FORWARD MODEL SELECTION USING POSITIVE EIGENVECTORS

# FIRST MAKE A MODEL WITH INTERCEPT ONLY
set.seed(93)
PD.AmPop.wUF.pos.Ndet.dbmem.dbrda.0 <- dbrda(PD.wUF.dist ~ 1, data = PD.AmPop.pos.dbmem)

# MAKE A MODEL WITH ALL POTENTIAL pos MEMS.
set.seed(93)
PD.AmPop.wUF.pos.Ndet.dbmem.dbrda.all <- dbrda(PD.wUF.dist ~ ., data = PD.AmPop.pos.dbmem)

# MODEL SELECTION
set.seed(93)
PD.AmPop.wUF.pos.Ndet.forsel <- ordiR2step(PD.AmPop.wUF.pos.Ndet.dbmem.dbrda.0, scope = PD.AmPop.wUF.pos.Ndet.dbmem.dbrda.all, permutations =9999)
PD.AmPop.wUF.pos.Ndet.forsel.results <- PD.AmPop.wUF.pos.Ndet.forsel$anova # 
PD.AmPop.wUF.pos.Ndet.forsel.results #ONLY MEM 2 Pr(>F): 0.00035 ***
# march 3, 2023: same result: 9e-04 *** # ***** same result April 3, 2023

# APPLY SIDAK CORRECTION
PD.wUF.AmPop.Ndet.final <- sidak(pValues= PD.AmPop.wUF.pos.Ndet.forsel$anova$`Pr(>F)`)
PD.wUF.AmPop.Ndet.final #0.0006998775 

########################## ALL GUT SAMPLES ##########################

######## CONSTRUCT DBMEMS ON CENTROID COORDS OF ALL GUT SAMPLES #########

###### POSITIVE DBMEMS ######
# CONSTRUCT POSITIVE ONLY DBMEM MODEL ON COORDS 
guts.AmPop.pos.dbmem <- dbmem(guts_centr_coord, silent = FALSE, MEM.autocor = "positive") #Orthobasis with 199 rows and 7 columns
guts.AmPop.pos.dbmem.thresh <- give.thresh(dist(guts_centr_coord))
guts.AmPop.pos.dbmem.thresh #337.2224

##### THIS LITTLE BIT IS JUST CHECKING CALCULATIONS #####
set.seed(1)
guts.AmPop.pos.dbmem1 <- dbmem(guts_centr_coord, silent = FALSE, MEM.autocor = "positive")
set.seed(93)
guts.AmPop.pos.dbmem93 <- dbmem(guts_centr_coord, silent = FALSE, MEM.autocor = "positive")
set.seed(2)
guts.AmPop.pos.dbmem2 <- dbmem(guts_centr_coord, silent = FALSE, MEM.autocor = "positive")

unique(guts.AmPop.pos.dbmem1 == guts.AmPop.pos.dbmem2) #all true 
unique(guts.AmPop.pos.dbmem93 == guts.AmPop.pos.dbmem2) #all true 

unique(guts_centr_coord)

########################

# DISPLAY EIGENVALUES
attributes(guts.AmPop.pos.dbmem)$values # 

##### NEGATIVE DBMEMS ######
# CONSTRUCT NEGATIVE ONLY DBMEM MODEL ON COORDS 
guts.AmPop.neg.dbmem <- dbmem(guts_centr_coord, silent = FALSE, MEM.autocor = "negative") 
#truncation level = 337.2224 

# DISPLAY EIGENVALUES
attributes(guts.AmPop.neg.dbmem)$values 

######### JACCARD ##########

# 1. CHECK FOR LINEAR TREND AND DETREND DATA
set.seed(93)
guts.jacc.AmPop.detrend.mod <-dbrda(gut.jacc.dist~guts_centr_coord[,1]+guts_centr_coord[,2])
guts.jacc.AmPop.detrend.mod.results <- anova.cca(guts.jacc.AmPop.detrend.mod, permutations=9999) # very significant linear trend
# March 3, 2023 AND APRIL 4, 2023 *****: SAME OVERALL 1e-04 *** 

# detrend data
guts.jacc.AmPop.dt.dat<-resid(guts.jacc.AmPop.detrend.mod) #save residuals 
guts.jacc.AmPop.dt.dat

# RUN POS DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
guts.AmPop.jacc.pos.dbmem.dbrda <- dbrda(guts.jacc.AmPop.dt.dat ~ as.matrix(guts.AmPop.pos.dbmem))
guts.AmPop.jacc.pos.dbmem.dbrda.results <- anova.cca(guts.AmPop.jacc.pos.dbmem.dbrda, permutations=9999)
guts.AmPop.jacc.pos.dbmem.dbrda.results #  Pr(>F)= 0.6997 NOT SIGNIFICANT **May 30, 2022- 0.00241 **
# April and ***** March 3, 2023: SIGNIFICANT: 0.0019 **

# RUN NEG DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
guts.AmPop.jacc.neg.dbmem.dbrda <- dbrda(guts.jacc.AmPop.dt.dat ~ as.matrix(guts.AmPop.neg.dbmem))
guts.AmPop.jacc.neg.dbmem.dbrda.results <- anova.cca(guts.AmPop.jacc.neg.dbmem.dbrda, permutations=9999)
guts.AmPop.jacc.neg.dbmem.dbrda.results # Pr(>F)= 1e-05 *** SIGNIFICANT ** May 30, 2022 - 0.9976
# and april 4, 2023 ***** March 3, 2023: NOT SIGNIFICANT 0.9982
############################################################################################################
#### ADD MAY 30, 2022 TRYING MODEL SELECTION WITH POSITIVE EIGENVECTORS

# FIRST MAKE A MODEL WITH INTERCEPT ONLY
set.seed(93)
guts.AmPop.jacc.pos.dbmem.dbrda.0 <- dbrda(guts.jacc.AmPop.dt.dat ~ 1, data = guts.AmPop.pos.dbmem)

# MAKE A MODEL WITH ALL POTENTIAL MEMS.
set.seed(93)
guts.AmPop.jacc.pos.dbmem.dbrda.all <- dbrda(guts.jacc.AmPop.dt.dat ~ ., data = guts.AmPop.pos.dbmem)

# MODEL SELECTION
set.seed(93) #CHANGE TO MORE PERMUTATIONS IF NECESSARY
guts.AmPop.jacc.pos.forsel <- ordiR2step(guts.AmPop.jacc.pos.dbmem.dbrda.0, scope = guts.AmPop.jacc.pos.dbmem.dbrda.all, permutations =9999)
guts.AmPop.jacc.pos.forsel.results <- guts.AmPop.jacc.pos.forsel$anova 
guts.AmPop.jacc.pos.forsel$adjust #1

# APPLY SIDAK CORRECTION
library(mutoss)
guts.jacc.AmPop.pos.final <- sidak(pValues= guts.AmPop.jacc.pos.forsel$anova$`Pr(>F)`)
guts.jacc.AmPop.pos.final #[1] 0.01787886 0.03044608 0.14604346         NA (these are mems7, 5, and 4)
###########END MAY 30 additions #################################################################################################

# ***** APRIL 4, 2023 GRAYED OUT SINCE NEGATIVE NO LONGER SIGNIFICANT!
# FORWARD MODEL SELECTION USING NEGATIVE EIGENVECTORS

# # FIRST MAKE A MODEL WITH INTERCEPT ONLY
# set.seed(93)
# guts.AmPop.jacc.neg.dbmem.dbrda.0 <- dbrda(guts.jacc.AmPop.dt.dat ~ 1, data = guts.AmPop.neg.dbmem)
# 
# # MAKE A MODEL WITH ALL POTENTIAL MEMS.
# set.seed(93)
# guts.AmPop.jacc.neg.dbmem.dbrda.all <- dbrda(guts.jacc.AmPop.dt.dat ~ ., data = guts.AmPop.neg.dbmem)
# 
# # MODEL SELECTION
# set.seed(93)
# guts.AmPop.jacc.neg.forsel <- ordiR2step(guts.AmPop.jacc.neg.dbmem.dbrda.0, scope = guts.AmPop.jacc.neg.dbmem.dbrda.all, permutations =9999)
# guts.AmPop.jacc.neg.forsel.results <- guts.AmPop.jacc.neg.forsel$anova # 10 significant ones. But, adjusted P-values are only 0.040210
# guts.AmPop.jacc.neg.forsel$adjust
# 
# # APPLY SIDAK CORRECTION
# library(mutoss)
# guts.jacc.AmPop.final <- sidak(pValues= guts.AmPop.jacc.neg.forsel$anova$`Pr(>F)`)
# guts.jacc.AmPop.final
# guts.jacc.AmPop.final$adjPValues #0.000659802 0.004829366 0.015075690 0.012033669 0.038151857 0.037408350 0.061981623 0.091643835 0.191323778 0.295573459

########### wUF ###############
# 1. CHECK FOR LINEAR TREND AND DETREND DATA
set.seed(93)
guts.wUF.AmPop.detrend.mod <-dbrda(guts.wUF.dist~guts_centr_coord[,1]+guts_centr_coord[,2])
guts.wUF.AmPop.detrend.mod.results <- anova.cca(guts.wUF.AmPop.detrend.mod, permutations=9999) # Pr(>F) = 0.00195 **
# MARCH 2, 2023: SAME 0.0021 **

guts.wUF.AmPop.dt.dat<-resid(guts.wUF.AmPop.detrend.mod)
guts.wUF.AmPop.dt.dat

# RUN POS DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
guts.AmPop.wUF.pos.dbmem.dbrda <- dbrda(guts.wUF.AmPop.dt.dat ~ as.matrix(guts.AmPop.pos.dbmem))
guts.AmPop.wUF.pos.dbmem.dbrda.results <- anova.cca(guts.AmPop.wUF.pos.dbmem.dbrda, permutations=9999)
guts.AmPop.wUF.pos.dbmem.dbrda.results #  Pr(>F)= 0.1257 NOT SIGNIFICANT **May 30, 2022: 0.04894 *
# when ran May 30, 2022, I got that it was barely significant 0.04894...
# MARCH 2, 2023 DIFFERENT: SIGNIFICANT 0.0473 *

# New ***** April 4, 2023
# FORWARD MODEL SELECTION USING POSITIVE EIGENVECTORS

# FIRST MAKE A MODEL WITH INTERCEPT ONLY
set.seed(93)
guts.AmPop.wUF.pos.dbmem.dbrda.0 <- dbrda(guts.wUF.AmPop.dt.dat ~ 1, data = guts.AmPop.pos.dbmem)

# MAKE A MODEL WITH ALL POTENTIAL MEMS.
set.seed(93)
guts.AmPop.wUF.pos.dbmem.dbrda.all <- dbrda(guts.wUF.AmPop.dt.dat ~ ., data = guts.AmPop.pos.dbmem)

# MODEL SELECTION
set.seed(93)
guts.AmPop.wUF.pos.forsel <- ordiR2step(guts.AmPop.wUF.pos.dbmem.dbrda.0, scope = guts.AmPop.wUF.pos.dbmem.dbrda.all, permutations =9999)
guts.AmPop.wUF.pos.forsel.results <- guts.AmPop.wUF.pos.forsel$anova # NO SIGNIFICANT ONES 

# APPLY SIDAK CORRECTION
guts.wUF.AmPop.pos.final <- sidak(pValues= guts.AmPop.wUF.pos.forsel$anova$`Pr(>F)`)
guts.wUF.AmPop.pos.final #  0.03786317 0.06797273 0.11256030 0.19921546   SO ONLY TWO SIGNIFICANT REMAIN

# RUN NEG DBMEM ANALYSIS ON DETRENDED DATA
set.seed(93)
guts.AmPop.wUF.neg.dbmem.dbrda <- dbrda(guts.wUF.AmPop.dt.dat ~ as.matrix(guts.AmPop.neg.dbmem))
guts.AmPop.wUF.neg.dbmem.dbrda.results <- anova.cca(guts.AmPop.wUF.neg.dbmem.dbrda, permutations=9999)
guts.AmPop.wUF.neg.dbmem.dbrda.results # Pr(>F)=  9e-05 *** VERY SIGNIFICANT **May 30, 2022: 0.9511
# MARCH 2, 2023 DIFFERENT: NOT SIGNIFICANT 0.9528
# ***** APRIL 4, 2023 : 0.7065 0.9528

# # FORWARD MODEL SELECTION USING NEGATIVE EIGENVECTORS
# ***** GRAYED OUT SINCE NO LONGER SIGNIFICANT
# # FIRST MAKE A MODEL WITH INTERCEPT ONLY
# set.seed(93)
# guts.AmPop.wUF.neg.dbmem.dbrda.0 <- dbrda(guts.wUF.AmPop.dt.dat ~ 1, data = guts.AmPop.neg.dbmem)
# 
# # MAKE A MODEL WITH ALL POTENTIAL MEMS.
# set.seed(93)
# guts.AmPop.wUF.neg.dbmem.dbrda.all <- dbrda(guts.wUF.AmPop.dt.dat ~ ., data = guts.AmPop.neg.dbmem)
# 
# # MODEL SELECTION
# set.seed(93)
# guts.AmPop.wUF.neg.forsel <- ordiR2step(guts.AmPop.wUF.neg.dbmem.dbrda.0, scope = guts.AmPop.wUF.neg.dbmem.dbrda.all, permutations =9999)
# guts.AmPop.wUF.neg.forsel.results <- guts.AmPop.wUF.neg.forsel$anova # 4 significant ones. But, adjusted P-values are only 0.049836 
# 
# # APPLY SIDAK CORRECTION
# guts.wUF.AmPop.neg.final <- sidak(pValues= guts.AmPop.wUF.neg.forsel$anova$`Pr(>F)`)
# guts.wUF.AmPop.neg.final #  0.03786317 0.06797273 0.11256030 0.19921546   SO ONLY TWO SIGNIFICANT REMAIN

guts.AmPop.pos.dbmem
guts.AmPop.pos.dbmem.thresh
guts.AmPop.neg.dbmem
guts.jacc.AmPop.detrend.mod
guts.jacc.AmPop.detrend.mod.results
guts.jacc.AmPop.dt.dat
guts.AmPop.jacc.pos.dbmem.dbrda
guts.AmPop.jacc.pos.dbmem.dbrda.results
guts.AmPop.jacc.neg.dbmem.dbrda
guts.AmPop.jacc.neg.dbmem.dbrda.results
guts.AmPop.jacc.neg.forsel
guts.AmPop.jacc.neg.forsel.results
guts.jacc.AmPop.final
guts.wUF.AmPop.detrend.mod
guts.wUF.AmPop.detrend.mod.results
guts.wUF.AmPop.dt.dat
guts.AmPop.wUF.pos.dbmem.dbrda
guts.AmPop.wUF.pos.dbmem.dbrda.results
guts.AmPop.wUF.neg.dbmem.dbrda
guts.AmPop.wUF.neg.dbmem.dbrda.results
guts.AmPop.wUF.neg.forsel
guts.AmPop.wUF.neg.forsel.results
guts.wUF.AmPop.neg.final

#####################################################
# This below was all saved March 2, 2023 in "copy" folder
# save(PV.AmPop.pos.dbmem, PV.AmPop.neg.dbmem, PV.AmPop.dbmem.thresh, PV.AmPop.jacc.pos.dbmem.dbrda, PV.AmPop.jacc.pos.dbmem.dbrda.results, PV.AmPop.jacc.neg.dbmem.dbrda, PV.AmPop.jacc.neg.dbmem.dbrda.results, PV.AmPop.jacc.neg.forsel, PV.AmPop.jacc.neg.forsel.results, PV.jacc.AmPop.final, PV.wUF.AmPop.detrend.mod, PV.wUF.AmPop.detrend.mod.results, PV.wUF.AmPop.dt.dat, PV.AmPop.wUF.pos.dbmem.dbrda, PV.AmPop.wUF.pos.dbmem.dbrda.results, PV.AmPop.wUF.neg.dbmem.dbrda, PV.AmPop.wUF.neg.dbmem.dbrda.results, PV.AmPop.wUF.neg.forsel, PV.AmPop.wUF.neg.forsel.results, PV.wUF.AmPop.final, file = "among.pop.PV.dbMEM.RData" )
# save(PD.AmPop.pos.dbmem, PD.AmPop.pos.dbmem.thresh, PD.AmPop.neg.dbmem, PD.jacc.AmPop.detrend.mod, PD.jacc.AmPop.detrend.mod.results, PD.jacc.AmPop.dt.dat, PD.AmPop.jacc.pos.dbmem.dbrda, PD.AmPop.jacc.pos.dbmem.dbrda.results, PD.AmPop.jacc.neg.dbmem.dbrda, PD.AmPop.jacc.neg.dbmem.dbrda.results, PD.AmPop.jacc.neg.forsel, PD.AmPop.jacc.neg.forsel.results, PD.jacc.AmPop.final, PD.wUF.AmPop.detrend.mod, PD.wUF.AmPop.detrend.mod.results, PD.AmPop.wUF.pos.dbmem.Ndet.dbrda, PD.AmPop.wUF.pos.dbmem.Ndet.dbrda.results, PD.AmPop.wUF.neg.dbmem.Ndet.dbrda, PD.AmPop.wUF.neg.dbmem.Ndet.dbrda.results, PD.AmPop.wUF.pos.Ndet.forsel, PD.AmPop.wUF.pos.Ndet.forsel.results, PD.wUF.AmPop.Ndet.final, file = "among.pop.PD.dbMEM.RData")
# save(guts.AmPop.pos.dbmem, guts.AmPop.pos.dbmem.thresh, guts.AmPop.neg.dbmem, guts.jacc.AmPop.detrend.mod, guts.jacc.AmPop.detrend.mod.results, guts.jacc.AmPop.dt.dat, guts.AmPop.jacc.pos.dbmem.dbrda, guts.AmPop.jacc.pos.dbmem.dbrda.results, guts.AmPop.jacc.neg.dbmem.dbrda, guts.AmPop.jacc.neg.dbmem.dbrda.results, guts.AmPop.jacc.neg.forsel, guts.AmPop.jacc.neg.forsel.results, guts.jacc.AmPop.final, guts.wUF.AmPop.detrend.mod, guts.wUF.AmPop.detrend.mod.results, guts.wUF.AmPop.dt.dat, guts.AmPop.wUF.pos.dbmem.dbrda, guts.AmPop.wUF.pos.dbmem.dbrda.results, guts.AmPop.wUF.neg.dbmem.dbrda, guts.AmPop.wUF.neg.dbmem.dbrda.results, guts.AmPop.wUF.neg.forsel, guts.AmPop.wUF.neg.forsel.results, guts.wUF.AmPop.neg.final, file = "among.pop.allguts.dbMEM.RData")
# 
# #  saved March 2, 2023 in "copy" folder
# save(guts.AmPop.pos.dbmem, guts.AmPop.pos.dbmem.thresh, guts.AmPop.neg.dbmem, guts.jacc.AmPop.detrend.mod, guts.jacc.AmPop.detrend.mod.results, guts.jacc.AmPop.dt.dat, guts.AmPop.jacc.pos.dbmem.dbrda, guts.AmPop.jacc.pos.dbmem.dbrda.results, guts.AmPop.jacc.neg.dbmem.dbrda, guts.AmPop.jacc.neg.dbmem.dbrda.results, guts.AmPop.jacc.neg.forsel, guts.AmPop.jacc.neg.forsel.results, guts.jacc.AmPop.final, guts.wUF.AmPop.detrend.mod, guts.wUF.AmPop.detrend.mod.results, guts.wUF.AmPop.dt.dat, guts.AmPop.wUF.pos.dbmem.dbrda, guts.AmPop.wUF.pos.dbmem.dbrda.results, guts.AmPop.wUF.neg.dbmem.dbrda, guts.AmPop.wUF.neg.dbmem.dbrda.results, guts.AmPop.wUF.neg.forsel, guts.AmPop.wUF.neg.forsel.results, guts.wUF.AmPop.neg.final, file = "among.pop.allguts.dbMEM_MAY302022.RData")
# #  saved APRIL 4, 2023 in "copy" folder
save(guts.AmPop.pos.dbmem, guts.AmPop.pos.dbmem.thresh, guts.AmPop.neg.dbmem, guts.jacc.AmPop.detrend.mod, guts.jacc.AmPop.detrend.mod.results, guts.jacc.AmPop.dt.dat, guts.AmPop.jacc.pos.dbmem.dbrda, guts.AmPop.jacc.pos.dbmem.dbrda.results, guts.AmPop.jacc.neg.dbmem.dbrda, guts.AmPop.jacc.neg.dbmem.dbrda.results, guts.wUF.AmPop.detrend.mod, guts.wUF.AmPop.detrend.mod.results, guts.wUF.AmPop.dt.dat, guts.AmPop.wUF.pos.dbmem.dbrda, guts.AmPop.wUF.pos.dbmem.dbrda.results, guts.AmPop.wUF.neg.dbmem.dbrda, guts.AmPop.wUF.neg.dbmem.dbrda.results, file = "among.pop.allguts.dbMEM_april4_2023.RData")
