#################################################################################
#   Setting up phyloseq objects, rarefying, Shannon, and getting dissimilarities     
#               Phan Biogeo 2019 Project
# (THIS SCRIPT FOLLOWS the script called q2_to_phyloseq2.R )
################################################################################
## This script rarefies gut and soil samples, gets the Shannon diversity of the soil samples,
# and makes new phyloseq objects. It also finds quantitative Jaccard and weighted UniFrac
# distances for various subsets of gut samples (i.e. PV only, PD only, allopatry only,
# sympatry only)

# finds quantitative Jaccard and weighted UniFrac distances for the gut samples

##### NOTE: all_samp_physeq WAS RENAMED TO all_samp_physeq_raw HERE (SO DIFFERENT THAN  q2_to_phyloseq2.R)

library("phyloseq")
library("ape")
library("vegan")

#load(file = "ps_rare_diss_shann_.RData") #created in this file and can be used to run a lot of the code here too

#######################################################################################

otu.table.all = read.csv(file = "~/Desktop/2019_biogeo_phan_GM/_All_pops/phyloseq/nonrarefied_OTU_matrix.csv", sep=",", row.names=1)
otu.table.all = as.matrix(otu.table.all)
head(otu.table.all)

taxonomy.all = read.csv("~/Desktop/2019_biogeo_phan_GM/_All_pops/phyloseq/taxonomy_only.csv", sep=",", row.names = 1)
taxonomy.all = as.matrix(taxonomy.all)
head(taxonomy.all)

metadata.no.shan = read.csv("~/Desktop/2019_biogeo_phan_GM/_All_pops/phyloseq/all_samples/metadata_nosoilshannon_forR_csv.csv", row.names=1) #need row names = 1 so that OTU names are consistent across objects
#for this, removed q2 "categories" label
#has all samples and all metadata, EXCEPT soil Shannon 
head(metadata.no.shan)

phy_tree = read_tree("~/Desktop/2019_biogeo_phan_GM/_All_pops/phyloseq/tree.nwk") #Yay! It worked!

OTU = otu_table(otu.table.all, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy.all)
META = sample_data(metadata.no.shan)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	CHECKING THAT OTU NAMES CONSISTENT
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree) 

#	CHECK THAT SAMPLE NAMES ARE CONSISTENT
# Both use . not -.
sample_names(OTU)
sample_names(META)

#	MERGE INTO ONE PHYLOSEQ OBJECT
all_samp_physeq_raw	<-	phyloseq(OTU,	TAX,	META,	phy_tree)
dim(otu_table(all_samp_physeq_raw))
#still 23641   288,
######what phyloseq looks like!### numbers look correct

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 23641 taxa and 288 samples ]
#sample_data() Sample Data:       [ 287 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 23641 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 23641 tips and 23640 internal nodes ]

sample_data(all_samp_physeq_raw)
# END q2_to_phyloseq2.R

# Re-saved March 1, 2023
#save(otu.table.all, taxonomy.all, metadata.no.shan, phy_tree, OTU, TAX, META, all_samp_physeq_raw, file ="raw_phyloseq.dat.RData")

#################################################################################
#                 PHYLOSEQ OBJECTS BY SAMPLE TYPE
################################################################################

# CREATE SEPARATE PHYLOSEQ OBJECTS FOR GUT AND SOIL SAMPLES
# guts and controls and blanks
guts_names <- rownames(sample_data(all_samp_physeq_raw))[grep(".*PV|PD.*",rownames(sample_data(all_samp_physeq_raw)))]
gutsAndCont <- c(guts_names, "PCRblankA", "PCRblankB", "PCRblankC", "EXT.BLANK")
only.guts.phyloseq <- prune_samples(gutsAndCont, all_samp_physeq_raw) #Looks good. 254 samples!
dim(otu_table(only.guts.phyloseq)) #23641, 254
gut.otu.table.tp <- t(otu_table(only.guts.phyloseq))

# SOIL SAMPLES
only.soil.samples <- rownames(sample_data(all_samp_physeq_raw))[grep("*.S.*",rownames(sample_data(all_samp_physeq_raw)))]
only.soil.phyloseq <- prune_samples(only.soil.samples, all_samp_physeq_raw) #Looks good. 34 samples!
dim(otu_table(only.soil.phyloseq)) #23641, 34

# DROP GUT SAMPLES WITH LESS THAN 3507 READS AND SOIL SAMPLES WITH LESS THAN 26336 FOR RAREFYING
## GUT
gutsThresh <- 3507 #threshold

rare_gutsIndex <- which(colSums(otu_table(only.guts.phyloseq)) >= gutsThresh)
GutSampsToKeep <- colnames(otu_table(only.guts.phyloseq))[rare_gutsIndex]
GutSampsToKeep #199 samples

# Only keep samples above 3507 threshold
otu.gut.to.rarefy <- prune_samples(GutSampsToKeep,only.guts.phyloseq)

## SOIL
soilsThresh <- 26336
rare_soilsIndex <- which(colSums(otu_table(only.soil.phyloseq)) >= soilsThresh)
length(rare_soilsIndex) #33
soilSampsToKeep <- colnames(otu_table(only.soil.phyloseq))[rare_soilsIndex]
otu.soil.to.rarefy <- prune_samples(soilSampsToKeep, only.soil.phyloseq)

#################################################################################
#               RAREFYING  MULTIPLE TIMES, JACCARD DIST, AND SHANNON
################################################################################

######################################### GUTS #########################################

# FUNCTION TO GET OTU OF PHYLOSEQ FOR RAREFYING (OR REALLY WORKING IN VEGAN IN GENERAL)
###### (switches rows and columns for vegan!)

psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# GUT OTU TABLE OUT OF PHYLOSEQ
otu.gut.to.rarefy.veg <-psotu2veg(otu.gut.to.rarefy)
dim(otu.gut.to.rarefy.veg) #199 23641 (so columns are now ASVs)

# GETTING RID OF ASVs THAT ARE PRESENT IN NO GUT SAMPLES (to avoid memory problems with and to speed up rarefying)
otu.gut.to.rarefy.nozeros <- otu.gut.to.rarefy.veg[, colSums(otu.gut.to.rarefy.veg !=0) > 0]
dim(otu.gut.to.rarefy.nozeros) #199 rows and 1358 columns 
#check to make sure that you keep OTU names
colnames(otu.gut.to.rarefy.nozeros) 
# To check that this worked correctly:
length(which(colSums(otu.gut.to.rarefy.veg) == 0)) #this equals 22,283
(length(which(colSums(otu.gut.to.rarefy.veg) == 0)) + dim(otu.gut.to.rarefy.nozeros)[2]) == dim(otu.gut.to.rarefy.veg)[2] #TRUE

# another way of doing this!
otu.gut.to.rarefy.nozeros.sec.way <- otu.gut.to.rarefy.veg[,which(!apply(otu.gut.to.rarefy.veg, 2, FUN = function(x) {all(x == 0)}))]
dim(otu.gut.to.rarefy.nozeros.sec.way) #199 1358

set.seed(93)
# RAREFYING GUTS MULTIPLE TIMES AND THEN TAKING MEAN 
gut.r.3507.1k<-array(dim=c(199, 1358)) #199 rows, 1358 columns 
raw.rare<-list()
set.seed(93) 
for(j in 1:199){
  tempsamp<-array(dim=c(1000, 1358)) # rarefying 1000 times, columns match my number of columns 
  cat("\n",j,"of 199")
  for(i in 1:1000){ 
    tempsamp[i,]<-rrarefy(otu.gut.to.rarefy.nozeros[j,],3507) #3507 is where we want to rarefy to across gut samples
  }
  raw.rare[[j]]<-tempsamp
  gut.r.3507.1k[j,] <- apply(tempsamp,2,mean) # gets mean
}

dim(gut.r.3507.1k) #199 1358
sum(gut.r.3507.1k[,])/199 #equals 3507, so each sample was rarefied!
class(gut.r.3507.1k)

# RE-APPLY ROW AND COLUMN NAMES
rownames(gut.r.3507.1k) <- rownames(otu.gut.to.rarefy.nozeros) 
colnames(gut.r.3507.1k) <- colnames(otu.gut.to.rarefy.nozeros) 

# SEPARATE P. VINDEX AND P. DIFFORMIS TO GET DIFFERENT DISTANCES
# P.VINDEX
PV.r.data.allcol <- gut.r.3507.1k[c(1:14, 27:32, 53:59, 79:87, 99:109, 128:141, 160:181, 194:199),]
dim(PV.r.data.allcol)  #89 1358
#filter out empty columns 
PV.r.data <- PV.r.data.allcol[, colSums(PV.r.data.allcol !=0) > 0]
dim(PV.r.data) #89, 855

#P. DIFFORMIS
PD.r.data.allcol <- gut.r.3507.1k[c(15:26, 33:52, 60:78, 88:98, 110:127, 142:159, 182:193),]
dim(PD.r.data.allcol) # 110 1358
#filter out empty columns
PD.r.data <- PD.r.data.allcol[, colSums(PD.r.data.allcol !=0) > 0]
dim(PD.r.data) #110, 907

#### GET SYMPATRY ONLY #####
row.names(gut.r.3507.1k)[7]
SYMP.r.data.allcol <- gut.r.3507.1k[c(7:59, 67:109,116:138), ]
dim(SYMP.r.data.allcol) #119, 1358
#filter out empty columns 
SYMP.r.data <- SYMP.r.data.allcol[, colSums(SYMP.r.data.allcol !=0) > 0]
dim(SYMP.r.data) #119, 943

#### GET ALLOPATRY ONLY (JUNE 30, 2020)
ALLO.r.data.allcol <- gut.r.3507.1k[c(1:6, 60:66, 110:115, 139:199),]
dim(ALLO.r.data.allcol) #80, 1358
### filter out empty columns
ALLO.r.data <- ALLO.r.data.allcol[, colSums(ALLO.r.data.allcol !=0) > 0]
dim(ALLO.r.data) #80, 765

# GET JACCARD DISTANCE

# For all gut samples
gut.jacc.dist <- vegdist(gut.r.3507.1k, method = "jaccard")
str(gut.jacc.dist) #looks like we have sample names!

# For just P. vindex
PV.jacc.dist <- vegdist(PV.r.data, method = "jaccard")

# For just P. difformis
PD.jacc.dist <- vegdist(PD.r.data, method = "jaccard")
str(PD.jacc.dist)

# For just symp gut samples
SYMP.jacc.dist <- vegdist(SYMP.r.data, method = "jaccard")

# For just allo gut samples
ALLO.jacc.dist <- vegdist(ALLO.r.data, method = "jaccard")

######################################### SOIL #########################################

# SOIL OTU TABLE OUT OF PHYLOSEQ
otu.soil.to.rarefy.veg <-psotu2veg(otu.soil.to.rarefy)
dim(otu.soil.to.rarefy.veg) # 33 23641 (so columns are now ASVs)

# GETTING RID OF ASVs THAT ARE PRESENT IN NO SOIL SAMPLES (to avoid memory problems with and to speed up rarefying)
otu.soil.to.rarefy.nozeros <- otu.soil.to.rarefy.veg[, colSums(otu.soil.to.rarefy.veg !=0) > 0]
dim(otu.soil.to.rarefy.nozeros) #33 22365
#check to make sure that you keep OTU names
colnames(otu.soil.to.rarefy.nozeros) 

set.seed(93)

# RAREFYING SOIL MULTIPLE TIMES AND THEN TAKING MEAN 
soil.r.26336.1k<-array(dim=c(33, 22365)) #33 rows, 22365 columns in otu.soil.to.rarefy.nozeros
raw.rare<-list()
set.seed(93) 
for(j in 1:33){
  tempsamp<-array(dim=c(1000, 22365)) # rarefying 1000 times, columns match my number of columns 
  cat("\n",j,"of 33")
  for(i in 1:1000){ 
    tempsamp[i,]<-rrarefy(otu.soil.to.rarefy.nozeros[j,], 26336) #26336 is where we want to rarefy to. 
  }
  raw.rare[[j]]<-tempsamp
  soil.r.26336.1k[j,] <- apply(tempsamp,2,mean) #this step takes the mean. 
}

dim(soil.r.26336.1k) #33 22365
sum(soil.r.26336.1k[,])/33 #26336, means we rarefied correctly!

# RE-APPLY ROW AND COLUMN NAMES
rownames(soil.r.26336.1k) <- rownames(otu.soil.to.rarefy.nozeros) 
colnames(soil.r.26336.1k) <- colnames(otu.soil.to.rarefy.nozeros) 

# GET SHANNON DIVERSITY FROM RAREFIED SOIL
soil_shann <- diversity(soil.r.26336.1k, index = "shannon", MARGIN = 1)

#################################################################################
#                 CREATE NEW PHYLOSEQ OBJECTS WITH STUFF MADE ABOVE
################################################################################

###### 1. ALL GUT SAMPLES ######
# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
gut.r.3507.1k.tp <-t(gut.r.3507.1k)
dim(gut.r.3507.1k.tp) # 1358  199

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
gut_r_OTU = otu_table(gut.r.3507.1k.tp, taxa_are_rows = TRUE)

# READ IN NEW CSV FOR METADATA THAT HAS CORRECT SHANNON DATA COLLECTED ABOVE
metadata.shann = read.csv("~/Desktop/2019_biogeo_phan_GM/_All_pops/phyloseq/all_samples/metadata_total_forR_csv.CSV", row.names=1) #need row names = 1 so that OTU names are consistent across objects
head(metadata.shann) 
dim(metadata.shann) #288, 20

gut_r_META = sample_data(metadata.shann)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(gut_r_OTU)
taxa_names(phy_tree) 

#	MERGE INTO ONE PHYLOSEQ OBJECT
gut.r.phyloseq <-	phyloseq(gut_r_OTU,	TAX,	gut_r_META,	phy_tree)
gut.r.phyloseq
dim(otu_table(gut.r.phyloseq))

#####
###### 2. ONLY P.VINDEX SAMPLES ######

# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
PV.r.data.tp <-t(PV.r.data)
dim(PV.r.data.tp) # 855  89

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
PV_r_OTU = otu_table(PV.r.data.tp, taxa_are_rows = TRUE)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(PV_r_OTU)
taxa_names(phy_tree) 
setdiff(taxa_names(PV_r_OTU), taxa_names(phy_tree)) #this shows that all of the taxa names in PV_r_OTU are also in phytree!

#	MERGE INTO ONE PHYLOSEQ OBJECT
PV.r.phyloseq <-	phyloseq(PV_r_OTU,	TAX,	gut_r_META,	phy_tree)
PV.r.phyloseq #855 taxa and 89 samples. What we expect

#####
###### 3. ONLY P.DIFFORMIS SAMPLES ######

# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
PD.r.data.tp <-t(PD.r.data)
dim(PD.r.data.tp) # 907 110

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
PD_r_OTU = otu_table(PD.r.data.tp, taxa_are_rows = TRUE)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(PD_r_OTU)
taxa_names(phy_tree)

#	MERGE INTO ONE PHYLOSEQ OBJECT
PD.r.phyloseq <-	phyloseq(PD_r_OTU,	TAX,	gut_r_META,	phy_tree)
PD.r.phyloseq #907 taxa and 110 samples, WHAT WE expect

###### 4. ONLY SYMPATRIC SAMPLES ######
# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
SYMP.r.data.tp <-t(SYMP.r.data)
dim(SYMP.r.data.tp) # 943  119

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
SYMP_r_OTU = otu_table(SYMP.r.data.tp, taxa_are_rows = TRUE)

SYMP_r_META = sample_data(metadata.shann)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(SYMP_r_OTU)
taxa_names(phy_tree) #are these different?

#	MERGE INTO ONE PHYLOSEQ OBJECT
SYMP.r.phyloseq<-	phyloseq(SYMP_r_OTU,	TAX,	SYMP_r_META,	phy_tree)
SYMP.r.phyloseq
sample_data(SYMP.r.phyloseq) # looks good!

#################################################################################
#                      GET NEW, POST RAREFACTION UNIFRAC DISTANCES
################################################################################

# GUT UNIFRAC DISTANCES
guts.wUF.dist <- UniFrac(gut.r.phyloseq, weighted = TRUE)
guts.uwUF.dist <-UniFrac(gut.r.phyloseq, weighted = FALSE)

# P. VINDEX UNIFRAC
PV.wUF.dist <- UniFrac(PV.r.phyloseq, weighted = TRUE)
PV.uwUF.dist <-UniFrac(PV.r.phyloseq, weighted = FALSE)

# P. DIFFORMIS UNIFRAC
PD.wUF.dist <- UniFrac(PD.r.phyloseq, weighted = TRUE)
PD.uwUF.dist <-UniFrac(PD.r.phyloseq, weighted = FALSE)

# SYMPATRY ONLY UNIFRAC
SYMP.wUF.dist <- UniFrac(SYMP.r.phyloseq, weighted = TRUE)
SYMP.uwUF.dist <-UniFrac(SYMP.r.phyloseq, weighted = FALSE)

##### USEFUL FUNCTION:
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

gut_r_vg_sampdat_catord_shavg <- read.csv(file = "gut_r_vg_sampdat_catord_shavg.csv", row.names = 1, header = TRUE)

Ord_Cattle_PresenceRank <- factor(gut_r_vg_sampdat_catord_shavg$CattlePresenceRank, levels = c("low", "medium", "high"), labels = c("low", "medium", "high"), ordered = TRUE)
length(Ord_Cattle_PresenceRank)

gut_r_veg_final_sampdat <- cbind(gut_r_vg_sampdat_catord_shavg[,c(1:8, 10:21)], Ord_Cattle_PresenceRank)
colnames(gut_r_veg_final_sampdat)
gut_r_veg_final_sampdat$Ord_Cattle_PresenceRank

#################################################################################
#                      SAVE OBJECTS AS R DATA FILES
################################################################################

###### NOTE: gut_r_veg_final_sampdat has ordinal cattle variables 
##### re-saved March 27, 2023
#save(gut_r_veg_final_sampdat, all_samp_physeq_raw, psotu2veg, pssd2veg, gut.r.3507.1k, gut.r.3507.1k, gut.jacc.dist, soil.r.26336.1k, soil_shann, gut.r.phyloseq, guts.wUF.dist, guts.uwUF.dist, PV.r.data, PD.r.data, SYMP.r.data, ALLO.r.data, PV.jacc.dist, PD.jacc.dist, SYMP.jacc.dist, ALLO.jacc.dist, PV.r.phyloseq, PD.r.phyloseq, PV.wUF.dist, PV.uwUF.dist, PD.wUF.dist, PD.uwUF.dist, file = "savedObjectsMarch27_2023/ps_rare_diss_shann_.RData")

