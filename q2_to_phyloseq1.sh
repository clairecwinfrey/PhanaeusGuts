#!/bin/bash

####################################################################################################
# Preparing UNRAREFIED qiime2 files for import into phyloseq
# Table is already filtered to remove mito, chloroplasts, singletons, and reads without phylum.
# May 17, 2020
#following Robert Murdoch's instructions here: https://github.com/rwmurdoch/Manneheimia/blob/master/README.md
####################################################################################################

####################################################################################################
#                                  export table
####################################################################################################
qiime tools export \
--input-path filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro-nodoubl.qza \
--output-path phyloseq \ #note that this re-wrote former feature-table.biom that was rarefied.

#####convert to tsv table
biom convert \
-i phyloseq/feature-table.biom \
-o phyloseq/nonrarefied_OTU_table.txt \
--to-tsv \

####now, open up the otu_table.txt in a text editor and change #OTUID to OTUID (was #OTU ID, now is OTUID)


######## NOT NECESSARY TO RE-EXPORT TAXONOMY AND TREE BELOW; JUST KEPT SCRIPT FOR
#DOCUMENTATION PURPOSESE
####################################################################################################
#                                  export taxonomy
####################################################################################################

#####Now, export taxonomy table
qiime tools export \
--input-path seq_taxonomy.qza \
--output-path phyloseq \
#is called taxonomy.tsv
#has 33,488 features (i.e. all remaining after sepp. Qiime2 cannot filter out taxonomy)

####open this up in a text editor and change feature ID to OTUID

####################################################################################################
#                                  export tree
####################################################################################################

###this tree is rooted; Robert's was unrooted in tutorial
qiime tools export \
--input-path frag_insert_output/insertion-tree.qza \
--output-path phyloseq \

######now, we'll need to merge taxonomy and OTU tables in R and then output a
#merged file.

