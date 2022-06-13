# shift in codon percentage
# (anc_codon - der_codon) / anc_codon
library(readr)
library(dplyr)

#names of the files
height_file <- 'combined_height_data.tsv'
chol_file <- 'combined_chol_data.tsv'
bmi_file <- 'combined_bmi_data.tsv' #4121
asthma_file <- 'combined_asthma_data.tsv'
wbc_file <- 'combined_wbc_data.tsv'
anxiety_file <- 'combined_anxiety_data.tsv' #10117
depress_file <- 'combined_depress_data.tsv' #11540
diabetes_file <- 'combined_diabetes_data.tsv' #9033
obesity_file <- 'combined_obesity_datatsv'#9057
ad_file <- 'combined_alzheimers_data.tsv' #8450


files <- c(height_file, chol_file, bmi_file, asthma_file, wbc_file, 
           anxiety_file, depress_file, diabetes_file, obesity_file, ad_file)

# read the file and filter the rsid's of interest (synonymous)
# function used to sort the data


delta_data <- read_delim("ancestral_and_derived_with_delta.tsv", delim = "\t", col_names = TRUE)

#remove any empty values for delta
delta_data <- delta_data[!is.na(delta_data$delta),]

table <- data.frame()
p_vals <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001)
#function to sort the data



#for loop to handle all the files
for (i in 1:length(files)){
  df <- read_delim(files[i], delim = "\t", col_names = TRUE)
  
  # combine the dataset with the delta codon freq
  combined_data <- inner_join(df, delta_data, c("rsid" = "input"))
  combined_data <- combined_data[!(!is.na(combined_data$beta) & combined_data$beta=="NaN"), ]
  combined_data$beta <- as.numeric(as.character(combined_data$beta))
  combined_data$pval <- as.numeric(as.character(combined_data$pval))
  
  #get the column index of the ancestral allele. 
  col_index <- grep("^anc_freq$", colnames(combined_data))
  
  
  change_percent <- apply(combined_data, MARGIN = 1, FUN = function(x){
    return((x[col_index] - x[(col_index + 1)])/x[col_index])
  })
  
  
  combined_data$change_percentage <- change_percent
  
  
  combined_data$pval <- as.numeric(as.character(combined_data$pval))
  #change the sign if the ref != ancestral
  combined_data$beta <- ifelse(combined_data$ref != combined_data$ancestral_allele,
                               -combined_data$beta, combined_data$beta)
  
  
  
  output <- c()
  for (j in 0:length(p_vals)){
    if (j == 0){
      betas <- abs(combined_data$beta)
      
      delta_percentage <- abs(combined_data$change_percentage)
      
      lin_mod <- summary(lm(betas ~ delta_percentage))
      
      
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
      betas <- abs(combined_data[combined_data$pval < p_vals[j],]$beta)
      delta_percentage <- abs(combined_data[combined_data$pval < p_vals[j],]$change_percentage)
      
      lin_mod <- summary(lm(betas ~ delta_percentage))
      
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
  abs_percentage_table <- rbind(abs_percentage_table, output)
  
}

colnames(abs_percentage_table) <- c("all (Estimate, Std. Err, Tval, Pval, N)", p_vals)

rownames(abs_percentage_table) <- c("Height", "Cholesterol", "BMI", "Asthma", 
                                "White Blood Cell", "Anxiety", "Depression", "Diabetes", 
                                "Obesity", "Alzheimer's Disease")

write.table(abs_percentage_table, file = "lin_regress_abs_change_percent.tsv", 
            sep = "\t", row.names = T, col.names = T)


