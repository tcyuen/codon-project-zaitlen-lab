#Installing packages
#install.packages("R.oo")      #needed for R.utils package
#install.packages("R.utils")   #needed for gunzip function

#Install package to use read_tsv()
#install.packages("readr")  

#Load package
library(R.oo)
library(R.utils)
library(readr)

#Standing Height Phenotype
#standing_height_both_sexes <- gunzip("50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "standing_height_both_sexes.tsv")
standing_height_both_sexes <- read_tsv("standing_height_both_sexes.tsv")

#Variants 
#variants <- gunzip("variants.tsv.bgz", "variants1.tsv")
variants <- read_tsv("variants1.tsv")


#New data table that has rsid and consequence (syn v. non-syn) columns
rsid_cons <- variants[,-c(1,2,3,4,5,7,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
rsid_cons <-variants[,c("rsid","consequence","consequence_category")]


#Combining rsid/cons data frame with GWAS dataframe

blendedGWAS <-cbind(standing_height_both_sexes, rsid_cons)

