library('readr')
library('data.table')
library('dplyr')



#names of the files
variant_file <- 'variants.tsv.bgz'
height_file <- '50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'
chol_file <- '30690_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz'
bmi_file <- '21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz' #4121
asthma_file <- '22127.gwas.imputed_v3.both_sexes.tsv.bgz'
wbc_file <- '30000_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'
anxiety_file <- 'KRA_PSY_ANXIETY.gwas.imputed_v3.both_sexes.tsv.bgz' #10117
depress_file <- '20544_11.gwas.imputed_v3.both_sexes.v2.tsv.bgz' #11540
diabetes_file <- 'E4_DM1.gwas.imputed_v3.both_sexes.tsv.bgz' #9033
obesity_file <- 'E4_OBESITY.gwas.imputed_v3.both_sexes.tsv.bgz'#9057
ad_file <- 'AD.gwas.imputed_v3.both_sexes.tsv.bgz' #8450


files <- c(height_file, chol_file, bmi_file, asthma_file, wbc_file, 
           anxiety_file, depress_file, diabetes_file, obesity_file, ad_file)

#read the file and filter the rsid's of interest (synonymous)
#function used to sort the data
f <- function(x, pos) subset(x, consequence_category == "synonymous")
rsid <- read_delim_chunked(variant_file, DataFrameCallback$new(f),
                           delim = "\t", col_names = TRUE) 

#only keep the variant, rsid, and consequence category 
rsid <- rsid %>% select(c("variant", "rsid", "consequence_category"))


listofdata <- list()
#function to sort the data
synonymous <- function(x, pos) subset(x, variant %in% rsid$variant) 

#for loop to handle all the files
for (i in 1:length(files)){
  df <- read_delim_chunked(files[i], DataFrameCallback$new(synonymous),
                           delim = "\t", col_names = TRUE)
  #only keep certain categories that seem interesting
  df <- df %>% select(c("variant", "beta", "se", "tstat", "pval"))
  df <- inner_join(df, rsid, c('variant' = 'variant'))
  listofdata[[i]] <- df
  
}

