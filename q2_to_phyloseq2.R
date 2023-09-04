# This script takes output from q2_to_phyloseq1.sh and continues to format it for exploration
# with phyloseq in R. The code used and steps taken here are adopted from R. Murdoch's code
# found at  https://github.com/rwmurdoch/Manneheimia/blob/master/qiime2_to_physeq1.R

getwd()
otu.all <- read.table(file = "phyloseq/nonrarefied_OTU_table.txt", header = TRUE)
head(otu.all)

tax <- read.csv(file="phyloseq/taxonomy_nospace.csv",header = TRUE)
#Robert originally had this as .tsv but I kept getting error: rror in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
#line 8458 did not have 3 elements. Was resolved when I opened in Excel, saved as csv, then axed the sep= value 
#use read.csv because this is identical to read.table except for csvs.
#ans this made it work!
head(tax)

merged_file.all <- merge(otu.all, tax, by.x = c("OTUID"), by.y = c("OTUID"))
head(merged_file.all)

write.table(merged_file.all, file = "combined_otu_tax_nonrarefied.txt", sep = '\t', col.names = TRUE, row.names = FALSE)

#now, open txt file that was just made in Excel and split into two files:
#1) for taxonomy only (only OTUID and taxonomic info columns)
#Now, use 
###For this file, is taxonomy_only.csv in phyloseq folder
#use data —> text-to-columns in Excel and separate on semicolon to get columns for kingdom, phylum, class, etc

######ADDED KINGDOM, Phylum, calss, etc. as headings 
#and 2)for the OTU matrix (containing only OTUID and abundances in each sample.) 
######new with non-rarefied: nonrarefied_OTU_matrix.csv