library(readr)
library(dplyr)

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

# read the file and filter the rsid's of interest (synonymous)
# function used to sort the data
f <- function(x, pos) subset(x, consequence_category == "synonymous")
rsid <- read_delim_chunked(variant_file, DataFrameCallback$new(f),
                           delim = "\t", col_names = TRUE) 

delta_data <- read_delim("ancestral_and_derived_with_delta.tsv", delim = "\t", col_names = TRUE)

#remove any empty values for delta
delta_data <- delta_data[!is.na(delta_data$delta),]

table <- data.frame()
p_vals <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001)
#function to sort the data
synonymous <- function(x, pos) subset(x, variant %in% rsid$variant) 


#for loop to handle all the files
for (i in 1:length(files)){
  df <- read_delim_chunked(files[i], DataFrameCallback$new(synonymous),
                           delim = "\t", col_names = TRUE)
  
  df <- inner_join(df, rsid, c('variant' = 'variant'))
  combined_data <- inner_join(df, delta_data, c("rsid" = "input"))
  combined_data <- combined_data[!(!is.na(combined_data$beta) & combined_data$beta=="NaN"), ]
  combined_data$beta <- as.numeric(as.character(combined_data$beta))
  combined_data$pval <- as.numeric(as.character(combined_data$pval))
  #change the sign if the ref != ancestral
  combined_data$beta <- ifelse(combined_data$ref != combined_data$ancestral_allele,
                                      -combined_data$beta, combined_data$beta)
  
  
  
  output <- c()
  for (j in 0:length(p_vals)){
    if (j == 0){
      betas <- combined_data$beta
      
      delta_codon_freqs <- combined_data$delta
      
      lin_mod <- summary(lm(betas ~ delta_codon_freqs))
      
      
      stderr <- round(lin_mod$coef[2,2], 5)
      beta_coef <- round(lin_mod$coef[2,1], 5)
      tval <- round(lin_mod$coef[2,3], 5)
      pval <- round(lin_mod$coef[2,4], 5) 
      n <- length(betas)
      values <- c(beta_coef, stderr, tval, pval, n)
      
      cell <- paste(values, sep="", collapse = ", ")
      output <- append(output, cell, after = length(output))
    }
    else{
      betas <- combined_data[combined_data$pval < p_vals[j],]$beta
      delta_codon_freqs <- combined_data[combined_data$pval < p_vals[j],]$delta
      
      lin_mod <- summary(lm(betas ~ delta_codon_freqs))
      
      stderr <- round(lin_mod$coef[2,2], 5)
      beta_coef <- round(lin_mod$coef[2,1], 5)
      tval <- round(lin_mod$coef[2,3], 5)
      pval <- round(lin_mod$coef[2,4], 5) 
      n <- length(betas)
      values <- c(beta_coef, stderr, tval, pval, n)
      
      cell <- paste(values, sep="", collapse = ", ")
      output <- append(output, cell, after = length(output))
    }
  }
  table <- rbind(table, output)
  
}

colnames(table) <- c("all (Estimate, Std. Err, Tval, Pval, N)", p_vals)
rownames(table) <- c("Height", "Cholesterol", "BMI", "Asthma", 
                     "White Blood Cell", "Anxiety", "Depression", "Diabetes", 
                     "Obesity", "Alzheimer's Disease")

write.table(table, file = "lin_regress_data.tsv", sep = "\t", row.names = T, col.names = T)


