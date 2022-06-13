library('readr')
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

names <- c("height", "chol", "bmi", "asthma", "wbc", "anxiety", 
           "depress", "diabetes", "obesity", "alzheimers")

#read the file and filter the rsid's of interest (synonymous)
#function used to sort the data
f <- function(x, pos) subset(x, consequence_category == "synonymous")
rsid <- read_delim_chunked(variant_file, DataFrameCallback$new(f),
                           delim = "\t", col_names = TRUE) 

#only keep the variant, rsid, raf, alt, and consequence category 
#rsid <- rsid %>% dplyr::select(c("variant","ref","alt", "rsid", "consequence_category"))



#function to sort the data
synonymous <- function(x, pos) subset(x, variant %in% rsid$variant) 


#for loop to handle combine the files
for (i in 1:length(files)){
  df <- read_delim_chunked(files[i], DataFrameCallback$new(synonymous),
                           delim = "\t", col_names = TRUE)
  #only keep certain categories that seem interesting
  # df <- df %>% select(c("variant", "beta", "se", "tstat", "pval"))
  df <- inner_join(df, rsid, c('variant' = 'variant'))
  
  write.table(df,file = paste(c("combined_", names[i], "_data.tsv"),
                              collapse = ""), 
              sep = "\t", col.names = T, row.names = F)
}




# # Makes vector with the number of characters in the allele column test df
# ch_count_vec <-nchar(test$allele)
# 
# # Turns vector into data frame
# count_char <-data.frame(ch_count_vec)
# 
# # Combines data frames then filters rsid's with only two alleles
# test2 <-cbind(test, count_char)
# test_allele <-subset(test2, ch_count_vec == 3)
# 
# #removes data with no value for ancestral allele
# filtered_alleles <- test_allele[!(!is.na(test_allele$allele_1) & test_allele$allele_1==""), ]

#############################################################################################
# Testing only on the height_file
df <- read_delim_chunked(height_file, DataFrameCallback$new(synonymous),
                         delim = "\t", col_names = TRUE)

#only keep certain categories that seem interesting
#df <- df %>% dplyr::select(c("variant", "beta", "se", "tstat", "pval"))
df <- inner_join(df, rsid, c('variant' = 'variant'))

write.table(df, file = "test.tsv", sep = "\t", col.names = T, row.names = F)

checking <- read_delim("test.tsv", delim = "\t", col_names = T)

