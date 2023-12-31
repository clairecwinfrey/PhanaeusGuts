#############################################################################################
#                 ENV MODEL SEL WITH NO SOIL AND ONLY AUTO-DUMMY CODING
#                            JUNE 26, 2020
# (THIS SCRIPT FOLLOWS the script called inter_pop_dbMEMUpdated.R. Earlier version was called
# inter_pop_dbMEM.R)
#############################################################################################

setwd("~/Desktop/UTK_Research/copyOf_All_pops")

library(vegan)
library(phyloseq)

load("ps_rare_diss_shann_.RData") # this has all the distances, phyloseq objects, etc. Made in Rarefy_Shannon_dissimilarities_final.R
# load("among.pop.allguts.dbMEM_april4_2023.RData") # this is where dbMEMs that were created and selected were saved. From "inter_pop_dbMEMUpadted.R"
load("inter_pop_dbMEMs_Nov28_2023.RData")
load("gut_coords.RData")

gut_r_veg_final_sampdat #Ord_Cattle_PresenceRank is correct
set.seed(93)

# RData files saved in this document:
#### *** NOT UPDATED YET NOV. 2023 ####
# 1. "env_geo_for_varpart_April4_2023.RData -- all model sel before varpart (i.e. compiling significant geo vars and mod sel on environ.)
# 2. "varpart_April4_2023.RData" -- variation partitioning 
# 3. "dbRDAsVarPart_allGuts.RData" -- model selection, VIFs, dbRDAs, and ANOVAs for all of the gut samples considered together
# 4."PV_PDseparate_dbRDAs.RData") -- db-RDAs using same vars as before on PV and PD samples separated
# 5. "alloSymp_separateDbRDAs.RData" - db-RDAs using same vars as before on allopatric and sympatric samples separated

################################################################ ALL GUTS ################################################################ 
###################### JACCARD ########################

####### 1.  MAKE A CATEGORY CALLED "GEO" THAT HAS X,Y COORDS AND SIGNIFICANT DBMEMS
# which MEMs were significant? Only positive
#guts.AmPop.pos.dbmem
guts.AmPop.jacc.pos.forsel.results #3 were significant . 
#Now afterr Sidak correction: #MEMs 7 and 5 only 
guts.jacc.AmPop.pos.final #Now afterr Sidak correction: #MEMs 7 and 5 only 

# MEMS 7 and 5
guts.jacc.geo <- cbind(guts_centr_coord[,1], guts_centr_coord[,2], guts.AmPop.pos.dbmem$MEM7, guts.AmPop.pos.dbmem$MEM5)  
colnames(guts.jacc.geo) <- c("centr_X", "centr_Y", "posMEM_7", "posMEM_5")
rownames(guts.jacc.geo) <- rownames(guts_centr_coord)

################ MODEL SELECTION WITH ENVIRONMENTAL VARIABLES ##################
set.seed(93)
colnames(gut_r_veg_final_sampdat)
gut_env_vars <- gut_r_veg_final_sampdat[, c(13:15,20:21)] #annual temp, temp. variability, shannon soil, ord cattle presence rank
gut_env_vars$Ord_Cattle_PresenceRank #still ordinal
colnames(gut_env_vars)

#1. make model with only intercept 
guts.jacc.env.dbrda.0 <- dbrda(gut.jacc.dist ~1, data = gut_env_vars) 

#2. make mod with all possible env vars
guts.jacc.env.dbrda.all <- dbrda(gut.jacc.dist ~ ., data = gut_env_vars) 

#3. double-stopping forward model sel
set.seed(93)
guts.jacc.env.for.sel <- ordiR2step(guts.jacc.env.dbrda.0, scope = guts.jacc.env.dbrda.all, permutations = 9999)
guts.jacc.env.for.sel.results <- guts.jacc.env.for.sel$anova
guts.jacc.env.for.sel.results 
#cattle, soil shann, temp, precip, temp. var

# R2.adj Df    AIC      F Pr(>F)    
# + Ord_Cattle_PresenceRank 0.010721  2 877.64 2.0728 0.0001 ***
#   + Soil.Shann.Avg          0.015704  1 877.61 1.9923 0.0007 ***
#   + Annual.Temp             0.020615  1 877.59 1.9780 0.0004 ***
#   + Annual.Precip           0.025039  1 877.67 1.8803 0.0011 ** 
#   + Temp.Variability        0.027239  1 878.18 1.4363 0.0326 *  
#   <All variables>           0.027239     

# MAKE DATAFRAME OF ENV VAR TO USE
colnames(gut_env_vars)
# guts.jacc.env.final <- gut_env_vars[,c(1,3:5)] #now all are significant
guts.jacc.env.final <- gut_env_vars
rownames(guts.jacc.env.final) 
colnames(guts.jacc.env.final) 
class(guts.jacc.env.final) # dataframe so....LOOKS GOOD AND READY FOR VARPART!
dim(guts.jacc.env.final)


#objects saved below were created April 4, 2023
#save(guts.jacc.env.for.sel.results, guts.jacc.env.final, guts.jacc.geo, guts.AmPop.jacc.pos.forsel.results, guts.jacc.AmPop.pos.final, guts.AmPop.wUF.pos.forsel.results, guts.wUF.env.for.sel.results, guts.wUF.env.final, file = "env_geo_for_varpart_April4_2023.RData")

################### VARIATION PARTITIONING ###################

####### JACCARD #######
set.seed(93)
gut.jacc.varpart <- varpart(gut.jacc.dist, as.data.frame(guts.jacc.geo), guts.jacc.env.final, as.data.frame(gut_r_veg_final_sampdat$Species), as.data.frame(gut_r_veg_final_sampdat$Patry))  
colnames(guts.jacc.env.final)

par(mfrow = c(1,2))
quartz()
gut.jacc.varpart.plot <-plot(gut.jacc.varpart, digits = 3, bg = c("yellow", "navyblue", "red", "lightblue"), Xnames = c("Geo.", "Envir.", "Species", "Patry"), id.size = 1.2) 



#objects saved below were created April 4, 2023
# save(gut.jacc.varpart, gut.jacc.varpart.plot, gut.wUF.varpart, gut.wUF.varpart.plot, file = "varpart_April4_2023.RData")

################### DBRDAS ###################
####### JACCARD #########
######## 1. MODEL SELECTION USING ALL PRE-SELECTED ENV AND GEO VARIABLES, AS WELL AS BIOTIC VARS OF INTEREST ##########

colnames(gut_r_veg_final_sampdat) #col ten is beetle mass, col 2 is beetle species, and col 4 is patry


guts.jacc.geo.df <- as.data.frame(guts.jacc.geo) #huh, I can't coerce it unless I save it as an object.
class(guts.jacc.geo.df)
colnames(guts.jacc.geo.df)

set.seed(93)
# included beetle mass too
gut.jacc.dbrda.full <- dbrda(gut.jacc.dist ~ guts.jacc.geo.df$centr_X + guts.jacc.geo.df$centr_Y + guts.jacc.geo.df$posMEM_7 + guts.jacc.geo.df$posMEM_5 + guts.jacc.env.final$Annual.Temp + guts.jacc.env.final$Annual.Precip + guts.jacc.env.final$Soil.Shann.Avg + guts.jacc.env.final$Temp.Variability + guts.jacc.env.final$Ord_Cattle_PresenceRank +gut_r_veg_final_sampdat[,10] + gut_r_veg_final_sampdat$Species + gut_r_veg_final_sampdat$Patry)
gut.jacc.dbrda.0 <- dbrda(gut.jacc.dist ~1)
# split up guts.jacc.geo in gut.jacc.dbrda.full. This allowed me to take out a variable!

set.seed(93)
gut.jacc.mod.forsel <- ordiR2step(gut.jacc.dbrda.0, scope = gut.jacc.dbrda.full, permutations = 9999)
gut.jacc.mod.forsel_results <- gut.jacc.mod.forsel$anova # 
gut.jacc.mod.forsel_results

# R2.adj Df    AIC      F Pr(>F)    
# + guts.jacc.env.final$Ord_Cattle_PresenceRank 0.010721  2 877.64 2.0728 0.0001 ***
#   + gut_r_veg_final_sampdat$Species             0.019945  1 876.75 2.8447 0.0001 ***
#   + gut_r_veg_final_sampdat$Patry               0.026399  1 876.42 2.2928 0.0001 ***
#   + guts.jacc.geo.df$centr_Y                    0.031405  1 876.36 2.0025 0.0003 ***
#   + guts.jacc.env.final$Annual.Precip           0.036306  1 876.32 1.9816 0.0011 ** 
#   + guts.jacc.geo.df$posMEM_5                   0.039490  1 876.62 1.6364 0.0074 ** 
#   + guts.jacc.geo.df$centr_X                    0.041466  1 877.17 1.3938 0.0384 *  
#   + guts.jacc.env.final$Soil.Shann.Avg          0.043605  1 877.67 1.4249 0.0327 *  
#   + guts.jacc.env.final$Annual.Temp             0.045587  1 878.20 1.3925 0.0396 *  
#   + guts.jacc.env.final$Temp.Variability        0.047758  1 878.69 1.4286 0.0322 *  
#   <All variables>                               0.049899                            

##### MAKE NEW GEO WITH ONLY THOSE THAT WE JUST FOUND TO BE SIGNIFICANT 
# (i.e. x, y, and MEM 5)
guts.jacc.geo.postsel <- guts.jacc.geo.df[,c(1,2,4)]
colnames(guts.jacc.geo.postsel)
class(guts.jacc.geo.postsel)
guts.jacc.geo.postsel.mat <- as.matrix(guts.jacc.geo.postsel)
class(guts.jacc.geo.postsel.mat)

###### MAKE DBRDA MODEL ######
set.seed(93) 
##### removed precip from this below because it had highest vif index of 22.91 (assumed overlap with X). After removal, all vif under 10 (except interaction effect and consituents)
# gut_r_veg_final_sampdat[ ,2] = Species
# gut_r_veg_final_sampdat[ ,4] = Patry

gut.jacc.fin.dbrda <- dbrda(gut.jacc.dist ~ guts.jacc.geo.postsel.mat + guts.jacc.env.final$Soil.Shann.Avg + guts.jacc.env.final$Ord_Cattle_PresenceRank + gut_r_veg_final_sampdat[ ,2] * gut_r_veg_final_sampdat[ ,4] + 
                              gut_r_veg_final_sampdat$Annual.Precip + gut_r_veg_final_sampdat$Annual.Temp + guts.jacc.env.final$Temp.Variability)
gut.jacc.fin.dbrda.results <- anova.cca(gut.jacc.fin.dbrda, permutations = 9999, by = "margin")
gut.jacc.fin.dbrda.results 
summary(gut.jacc.fin.dbrda)
# Model: dbrda(formula = gut.jacc.dist ~ guts.jacc.geo.postsel.mat + guts.jacc.env.final$Soil.Shann.Avg + guts.jacc.env.final$Ord_Cattle_PresenceRank + gut_r_veg_final_sampdat[, 2] * gut_r_veg_final_sampdat[, 4] + gut_r_veg_final_sampdat$Annual.Precip + gut_r_veg_final_sampdat$Annual.Temp)
# Df SumOfSqs      F Pr(>F)    
# guts.jacc.geo.postsel.mat                                   3    2.024 1.7185 0.0001 ***
#   guts.jacc.env.final$Soil.Shann.Avg                          1    0.573 1.4609 0.0253 *  
#   guts.jacc.env.final$Ord_Cattle_PresenceRank                 2    1.545 1.9679 0.0001 ***
#   gut_r_veg_final_sampdat$Annual.Precip                       1    0.669 1.7031 0.0052 ** 
#   gut_r_veg_final_sampdat$Annual.Temp                         1    0.535 1.3628 0.0484 *  
#   gut_r_veg_final_sampdat[, 2]:gut_r_veg_final_sampdat[, 4]   1    0.857 2.1841 0.0002 ***
#   Residual                                                  187   73.402                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > summary(gut.jacc.fin.dbrda)


gut.jacc.fin.dbrda.results 
# Permutation test for dbrda under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 9999
# 
# Model: dbrda(formula = gut.jacc.dist ~ guts.jacc.geo.postsel.mat + guts.jacc.env.final$Soil.Shann.Avg + guts.jacc.env.final$Ord_Cattle_PresenceRank + gut_r_veg_final_sampdat[, 2] * gut_r_veg_final_sampdat[, 4] + gut_r_veg_final_sampdat$Annual.Precip + gut_r_veg_final_sampdat$Annual.Temp + guts.jacc.env.final$Temp.Variability)
# Df SumOfSqs      F Pr(>F)    
# guts.jacc.geo.postsel.mat                                   3    2.039 1.7352 0.0001 ***
#   guts.jacc.env.final$Soil.Shann.Avg                          1    0.630 1.6092 0.0074 ** 
#   guts.jacc.env.final$Ord_Cattle_PresenceRank                 2    1.529 1.9516 0.0001 ***
#   gut_r_veg_final_sampdat$Annual.Precip                       1    0.659 1.6814 0.0075 ** 
#   gut_r_veg_final_sampdat$Annual.Temp                         1    0.567 1.4468 0.0295 *  
#   guts.jacc.env.final$Temp.Variability                        1    0.541 1.3798 0.0463 *  
#   gut_r_veg_final_sampdat[, 2]:gut_r_veg_final_sampdat[, 4]   1    0.835 2.1312 0.0006 ***
#   Residual                                                  186   72.862                  


vif.gut.jacc.fin.dbrda <- vif.cca(gut.jacc.fin.dbrda)
vif.gut.jacc.fin.dbrda # Temp, precip, and X and Y have high VIFs. 
# guts.jacc.geo.postsel.matcentr_X                                          guts.jacc.geo.postsel.matcentr_Y 
# 9.614222                                                                 36.017668 
# guts.jacc.geo.postsel.matposMEM_5                                        guts.jacc.env.final$Soil.Shann.Avg 
# 3.305247                                                                  3.378674 
# guts.jacc.env.final$Ord_Cattle_PresenceRank.L                             guts.jacc.env.final$Ord_Cattle_PresenceRank.Q 
# 4.553174                                                                  7.801464 
# gut_r_veg_final_sampdat[, 2]P.vindex                                      gut_r_veg_final_sampdat[, 4]Sympatry 
# 12.045446                                                                  9.397523 
# gut_r_veg_final_sampdat$Annual.Precip                                       gut_r_veg_final_sampdat$Annual.Temp 
# 12.428185                                                                 54.262589 
# guts.jacc.env.final$Temp.Variability gut_r_veg_final_sampdat[, 2]P.vindex:gut_r_veg_final_sampdat[, 4]Sympatry 
# 4.307364                                                                 10.410232 

# Investigating a few likely correlations
precip_long.cor <- cor(guts.jacc.env.final$Annual.Precip, guts.jacc.geo.df$centr_X)
precip_long.cor #0.901! Very high
temp_lat.cor <- cor(guts.jacc.env.final$Annual.Temp, guts.jacc.geo.df$centr_Y)
temp_lat.cor #-0.8635219! Very high

# REMOVE TEMP AND PRECIP BECAUSE HIGH VIF AND CORR WITH LAT AND LONG
gut.jacc.fin.dbrda_POSTVIF <- dbrda(gut.jacc.dist ~ guts.jacc.geo.postsel.mat + guts.jacc.env.final$Soil.Shann.Avg + guts.jacc.env.final$Ord_Cattle_PresenceRank + gut_r_veg_final_sampdat[ ,2] * gut_r_veg_final_sampdat[ ,4] + + guts.jacc.env.final$Temp.Variability)
vif.gut.jacc.fin.dbrda_POSTVIF <- vif.cca(gut.jacc.fin.dbrda_POSTVIF)
vif.gut.jacc.fin.dbrda_POSTVIF
# Now, except for Phanaeus species and "patry", all VIFs are below 3.
# guts.jacc.geo.postsel.matcentr_X                                          guts.jacc.geo.postsel.matcentr_Y 
# 2.814403                                                                  3.004592 
# guts.jacc.geo.postsel.matposMEM_5                                      guts.jacc.env.final$Temp.Variability 
# 1.317170                                                                  1.510967 
# guts.jacc.env.final$Soil.Shann.Avg                             guts.jacc.env.final$Ord_Cattle_PresenceRank.L 
# 2.909812                                                                  1.225509 
# guts.jacc.env.final$Ord_Cattle_PresenceRank.Q                                      gut_r_veg_final_sampdat[, 2]P.vindex 
# 2.304351                                                                 10.973212 
# gut_r_veg_final_sampdat[, 4]Sympatry gut_r_veg_final_sampdat[, 2]P.vindex:gut_r_veg_final_sampdat[, 4]Sympatry 
# 7.533435                                                                  9.554460 

# Test significance of constraints
set.seed(93)                            
gut.jacc.fin.dbrda_POSTVIF.results <- anova.cca(gut.jacc.fin.dbrda, permutations = 9999, by = "margin")
gut.jacc.fin.dbrda_POSTVIF.results

# Permutation test for dbrda under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 9999
# 
# Model: dbrda(formula = gut.jacc.dist ~ guts.jacc.geo.postsel.mat + guts.jacc.env.final$Soil.Shann.Avg + guts.jacc.env.final$Ord_Cattle_PresenceRank + gut_r_veg_final_sampdat[, 2] * gut_r_veg_final_sampdat[, 4] + gut_r_veg_final_sampdat$Annual.Precip + gut_r_veg_final_sampdat$Annual.Temp + guts.jacc.env.final$Temp.Variability)
# Df SumOfSqs      F Pr(>F)    
# guts.jacc.geo.postsel.mat                                   3    2.039 1.7352 0.0001 ***
#   guts.jacc.env.final$Soil.Shann.Avg                          1    0.630 1.6092 0.0090 ** 
#   guts.jacc.env.final$Ord_Cattle_PresenceRank                 2    1.529 1.9516 0.0001 ***
#   gut_r_veg_final_sampdat$Annual.Precip                       1    0.659 1.6814 0.0042 ** 
#   gut_r_veg_final_sampdat$Annual.Temp                         1    0.567 1.4468 0.0286 *  
#   guts.jacc.env.final$Temp.Variability                        1    0.541 1.3798 0.0448 *  
#   gut_r_veg_final_sampdat[, 2]:gut_r_veg_final_sampdat[, 4]   1    0.835 2.1312 0.0007 ***
#   Residual                                                  186   72.862                                                      187   73.402           

# What if we took out the interaction effect?
gut.jacc.fin.noint.dbrda <- dbrda(gut.jacc.dist ~ guts.jacc.geo.postsel.mat + guts.jacc.env.final$Soil.Shann.Avg + guts.jacc.env.final$Ord_Cattle_PresenceRank + gut_r_veg_final_sampdat[ ,2] + gut_r_veg_final_sampdat[ ,4])
vif.gut.jacc.fin.noint.dbrda <- vif.cca(gut.jacc.fin.noint.dbrda) #now all 5 or lower. 

#### what if we did varpart on this model? #####
### July 1 #####
# final.jacc.guts.dbrda.varpart <- varpart(gut.jacc.dist, as.data.frame(guts.jacc.geo.postsel.mat), as.data.frame(guts.jacc.env.final[,c(2,4)]), as.data.frame(gut_r_veg_final_sampdat[ ,2]), as.data.frame(gut_r_veg_final_sampdat[ ,4]) )
# Nov. 28th-- needed  cattle, precip, both temps, and soil shannon
final.jacc.guts.dbrda.varpart <- varpart(gut.jacc.dist, as.data.frame(guts.jacc.geo.postsel.mat), as.data.frame(guts.jacc.env.final), as.data.frame(gut_r_veg_final_sampdat[ ,2]), as.data.frame(gut_r_veg_final_sampdat[ ,4]) )

par(mfrow = c(1,2))
#quartz()
plot(final.jacc.guts.dbrda.varpart, digits = 3, bg = c("yellow", "navyblue", "red", "lightblue"), Xnames = c("Geo.", "Env", "Species", "Patry"), id.size = 1.2) 
#### nothing too interesting here #####

plot(final.guts.dbrda.varpart, digits = 3, bg = c("yellow", "navyblue", "red", "lightblue"), Xnames = c("Geo.", "Env", "Species", "Patry"), id.size = 1.2) 
#### nothing too interesting here #####

# SAVE ALL OF THE RESULTS FOR ALL OF THE GUTS CONSIDERED TOGETHER
#objects saved below April 5, 2023
# save(guts.jacc.geo, guts.jacc.env.for.sel, guts.jacc.env.for.sel.results, guts.jacc.env.final, guts.wUF.geo, guts.wUF.env.for.sel, guts.wUF.env.for.sel.results, guts.wUF.env.final,
#      gut.jacc.varpart, gut.wUF.varpart, gut.jacc.mod.forsel, gut.jacc.mod.forsel_results, guts.jacc.geo.postsel.mat, gut.jacc.fin.dbrda, gut.jacc.fin.dbrda.results, vif.gut.jacc.fin.dbrda,
#      gut.jacc.fin.dbrda_POSTVIF, vif.gut.jacc.fin.dbrda_POSTVIF, gut.jacc.fin.dbrda_POSTVIF.results, guts.wUF.dbrda.full.vif, gut.wUF.mod.forsel, gut.wUF.mod.forsel_results, guts.wUF.geo.postsel.mat, 
#      gut.wUF.fin.dbrda, gut.wUF.fin.dbrda.results, final.guts.wUFdbrda.varpart, file= "dbRDAsVarPart_allGuts.RData")
# # SAVED NOV. 28, 2023
# save(guts.jacc.geo, guts.jacc.env.for.sel, guts.jacc.env.for.sel.results, guts.jacc.env.final,gut.jacc.varpart, gut.jacc.mod.forsel, gut.jacc.mod.forsel_results, guts.jacc.geo.postsel.mat, gut.jacc.fin.dbrda, gut.jacc.fin.dbrda.results, vif.gut.jacc.fin.dbrda,
#      gut.jacc.fin.dbrda_POSTVIF, vif.gut.jacc.fin.dbrda_POSTVIF, gut.jacc.fin.dbrda_POSTVIF.results, file= "dbRDAsVarPart_allGutsNov28_2023.RData")

### STOPPED HERE, NOV. 28, 2023 ####
################### SEPARATE MODELS FOR PV AND PD ###################
# BECAUSE WE FOUND AN INTERACTION EFFECT BETWEEN SPECIES AND PATRY, WE WANT TO RUN MODELS FOR PV AND PD SEPARATELY

#### SUBSET guts.jacc.geo BY PHANAEUS SPECIES AND REMOVE DBMEMS B/C THEY WERE CONSTRUCTED ON ALL COORDS. 
guts.jacc.geo.postsel.PV <- guts.jacc.geo.postsel.mat[c(1:14, 27:32, 53:59, 79:87, 99:109, 128:141, 160:181, 194:199), c(1:2)] #PV
class(guts.jacc.geo.postsel.PV)

guts.jacc.geo.postsel.PD <- guts.jacc.geo.postsel.mat[c(15:26, 33:52, 60:78, 88:98, 110:127, 142:159, 182:193), c(1:2)] #PD

#### SUBSET ENV FINAL BY PHANAEUS SPECIES 
dim(guts.jacc.env.final)
colnames(guts.jacc.env.final)
# REMOVE TEMP AND PRECIP FOR REASON ABOVE (high VIFs)
guts.jacc.env.final_PV <- guts.jacc.env.final[c(1:14, 27:32, 53:59, 79:87, 99:109, 128:141, 160:181, 194:199), c(3:4)] #PV

guts.jacc.env.final_PD <- guts.jacc.env.final[c(15:26, 33:52, 60:78, 88:98, 110:127, 142:159, 182:193),c(3:4)] #PD

#### SUBSET SAMPDAT BY PHANAEUS SPECIES
PV_sampdat <- gut_r_veg_final_sampdat[c(1:14, 27:32, 53:59, 79:87, 99:109, 128:141, 160:181, 194:199),] #PV
PV_sampdat[,4]

PD_sampdat <- gut_r_veg_final_sampdat[c(15:26, 33:52, 60:78, 88:98, 110:127, 142:159, 182:193),] #PD

####### MAKE PV MODEL #######
# (without dbMEMs)
set.seed(93)
PVsubset.jacc.final.dbrda <- dbrda(PV.jacc.dist ~ guts.jacc.geo.postsel.PV + guts.jacc.env.final_PV$Soil.Shann.Avg + guts.jacc.env.final_PV$Ord_Cattle_PresenceRank +PV_sampdat[,4])
PVsubset.jacc.final.dbrda.results <- anova.cca(PVsubset.jacc.final.dbrda, by = "margin", permutations =9999)
PVsubset.jacc.final.dbrda.results #geo, cattle, patry 

# #*                                               Df SumOfSqs      F Pr(>F)    
# guts.jacc.geo.postsel.PV                        2    1.425 1.7865 0.0003 ***
#   guts.jacc.env.final_PV$Soil.Shann.Avg           1    0.496 1.2430 0.1240    
# guts.jacc.env.final_PV$Ord_Cattle_PresenceRank  2    1.261 1.5818 0.0024 ** 
#   PV_sampdat[, 4]                                 1    0.699 1.7518 0.0041 ** 
#   Residual                                       82   32.696   

PVsubset.jacc.final.dbrda.vif <- vif.cca(PVsubset.jacc.final.dbrda) # all below 3 

###### MAKE PD MODEL #######
# (without dbMEMs)
set.seed(93)
PDsubset.jacc.final.dbrda <- dbrda(PD.jacc.dist ~ guts.jacc.geo.postsel.PD + guts.jacc.env.final_PD$Soil.Shann.Avg + guts.jacc.env.final_PD$Ord_Cattle_PresenceRank + PD_sampdat[,4])
PDsubset.jacc.final.dbrda.results <- anova.cca(PDsubset.jacc.final.dbrda, by = "margin", permutations =9999)
PDsubset.jacc.final.dbrda.results #patry, cattle, geo
# Df SumOfSqs      F Pr(>F)    
# guts.jacc.geo.postsel.PD                         2    1.438 1.8599 0.0003 ***
#   guts.jacc.env.final_PD$Soil.Shann.Avg            1    0.515 1.3309 0.0739 .  
# guts.jacc.env.final_PD$Ord_Cattle_PresenceRank   2    1.362 1.7611 0.0002 ***
#   PD_sampdat[, 4]                                  1    0.604 1.5615 0.0167 *  
#   Residual                                       103   39.828   


PDsubset.jacc.final.dbrda.vif <- vif.cca(PDsubset.jacc.final.dbrda) # all below 3
PDsubset.jacc.final.dbrda.vif # all below 6 (sympatry and soil shann avg 5.XX)

# last saved April 5, 2023
#save(PVsubset.jacc.final.dbrda, PVsubset.jacc.final.dbrda.results, PVsubset.jacc.final.dbrda.vif, PDsubset.jacc.final.dbrda,
    # PDsubset.jacc.final.dbrda.results, PDsubset.jacc.final.dbrda.vif, file = "PV_PDseparate_dbRDAs.RData")

############# ALL GUTS, ONLY IN SYMPATRY ##############
# JUNE 30, 2020
# BECAUSE WE FOUND AN INTERACTION EFFECT BETWEEN SPECIES AND PATRY, WE WANT TO RUN MODELS FOR JUST SYMPATRY

#### SUBSET guts.jacc.geo IN SYMPATRY VERSUS ALLOPATRY
guts.jacc.geo.postsel.SYMP <- guts.jacc.geo.postsel.mat[c(7:59, 67:109, 116:138), c(1:2)] 
dim(guts.jacc.geo.postsel.SYMP) #119,2 

guts.jacc.geo.postsel.ALLO <- guts.jacc.geo.postsel.mat[c(1:6, 60:66, 110:115, 139:199), c(1:2)] 
dim(guts.jacc.geo.postsel.ALLO) #80, 2

#### SUBSET ENV FINAL BY RANGE OVERLAP
dim(guts.jacc.env.final)
colnames(guts.jacc.env.final)
guts.jacc.env.final_SYMP <- guts.jacc.env.final[c(7:59, 67:109, 116:138), c(3:4)]

guts.jacc.env.final_ALLO <- guts.jacc.env.final[c(1:6, 60:66, 110:115, 139:199),c(3:4)] 

#### SUBSET SAMPDAT BY PHANAEUS SPECIES
SYMP_sampdat <- gut_r_veg_final_sampdat[c(7:59, 67:109, 116:138),]

ALLO_sampdat <- gut_r_veg_final_sampdat[c(1:6, 60:66, 110:115, 139:199),] 

#### DBRDAS ####

# IN SYMPATRY
set.seed(93)
SYMPsubset.jacc.final.dbrda <- dbrda(SYMP.jacc.dist ~ guts.jacc.geo.postsel.SYMP + guts.jacc.env.final_SYMP$Soil.Shann.Avg + guts.jacc.env.final_SYMP$Ord_Cattle_PresenceRank +SYMP_sampdat[,2])
SYMPsubset.jacc.final.dbrda.results <- anova.cca(SYMPsubset.jacc.final.dbrda, by = "margin", permutations =9999)
SYMPsubset.jacc.final.dbrda.results

# Df SumOfSqs      F Pr(>F)    
# guts.jacc.geo.postsel.SYMP                         2    1.671 2.1835 0.0001 ***
#   guts.jacc.env.final_SYMP$Soil.Shann.Avg            1    0.740 1.9324 0.0013 ** 
#   guts.jacc.env.final_SYMP$Ord_Cattle_PresenceRank   2    1.536 2.0059 0.0001 ***
#   SYMP_sampdat[, 2]                                  1    1.190 3.1100 0.0001 ***

SYMPsubset.jacc.vif <- vif.cca(SYMPsubset.jacc.final.dbrda)
SYMPsubset.jacc.vif #all below 5

###### varpart
SYMP_subset_vp <- varpart(SYMP.jacc.dist, as.data.frame(guts.jacc.geo.postsel.SYMP), as.data.frame(guts.jacc.env.final_SYMP), as.data.frame(SYMP_sampdat[,2]))

par(mfrow = c(1,2))
quartz()
SYMP_subset_vp.plot <-plot(SYMP_subset_vp, digits = 3, bg = c("yellow", "navyblue", "red"), Xnames = c("Geo.", "Envir.", "Species"), id.size = 1.2) 

# IN ALLOPATRY
set.seed(93)
ALLOsubset.jacc.final.dbrda <- dbrda(ALLO.jacc.dist ~ guts.jacc.geo.postsel.ALLO + guts.jacc.env.final_ALLO$Soil.Shann.Avg + guts.jacc.env.final_ALLO$Ord_Cattle_PresenceRank +ALLO_sampdat[,2])
ALLOsubset.jacc.final.dbrda.results <- anova.cca(ALLOsubset.jacc.final.dbrda, by = "margin", permutations =9999)
ALLOsubset.jacc.final.dbrda.results #shows that geo, cattle, and species are important

# guts.jacc.geo.postsel.ALLO                        2   1.4667 1.8418 0.0001 ***
#   guts.jacc.env.final_ALLO$Soil.Shann.Avg           1   0.4906 1.2322 0.1185    
# guts.jacc.env.final_ALLO$Ord_Cattle_PresenceRank  2   1.1356 1.4260 0.0089 ** 
#   ALLO_sampdat[, 2]                                 1   0.6599 1.6574 0.0076 ** 

ALLOsubset.jacc.vif <- vif.cca(ALLOsubset.jacc.final.dbrda)
ALLOsubset.jacc.vif # species high VIF, 
# guts.jacc.geo.postsel.ALLOcentr_X                  guts.jacc.geo.postsel.ALLOcentr_Y 
# 6.302495                                           4.262943 
# guts.jacc.env.final_ALLO$Soil.Shann.Avg guts.jacc.env.final_ALLO$Ord_Cattle_PresenceRank.L 
# 1.560815                                           2.913364 
# guts.jacc.env.final_ALLO$Ord_Cattle_PresenceRank.Q                          ALLO_sampdat[, 2]P.vindex 
# 4.590142                                           9.493895 

###### varpart?
ALLO_subset_vp <- varpart(ALLO.jacc.dist, as.data.frame(guts.jacc.geo.postsel.ALLO), as.data.frame(guts.jacc.env.final_ALLO), as.data.frame(ALLO_sampdat[,2]))

par(mfrow = c(1,2))
quartz()
ALLO_subset_vp.plot <-plot(ALLO_subset_vp, digits = 3, bg = c("yellow", "navyblue", "red"), Xnames = c("Geo.", "Envir.", "Species"), id.size = 1.2) 
#### nothing too interesting here #####

# Last saved April 5, 2023
# save(SYMPsubset.jacc.final.dbrda, SYMPsubset.jacc.final.dbrda.results, SYMPsubset.jacc.vif, SYMP_subset_vp,
#      ALLOsubset.jacc.final.dbrda, ALLOsubset.jacc.final.dbrda.results, ALLOsubset.jacc.vif, ALLO_subset_vp, file= "alloSymp_separateDbRDAs.RData")

