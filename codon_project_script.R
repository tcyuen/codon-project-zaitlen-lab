library('readr')
library('data.table')
library('dplyr')

#names of the files
variant_file <- 'variants.tsv.bgz'
height_file <- '50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'

#read the file and filter the rsid's of interest (synonymous)

rsid <- read_delim(variant_file, delim = "\t", col_names = TRUE, 
                   col_select = c("variant", "rsid", "consequence_category", 
                                  "ref", "alt"))

head(rsid, 5)
filtered_rsids <- filter(rsid, consequence_category == "synonymous")


head(filtered_rsids, 5)
#read in the other data and filter 
standing_height <- read_delim(height_file, delim = '\t', col_names = TRUE,
                              col_select = c("variant", "beta", "se", "tstat", 
                                             "pval"))

#filtering the standing height data
s_h_syn <- inner_join(filtered_rsids, standing_height, 
                           c('variant' = 'variant'))
head(s_h_syn, 5)

