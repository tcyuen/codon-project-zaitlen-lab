# shift in codon percentage
# (anc_codon - der_codon) / anc_codon
library(readr)
library(dplyr)
library(ggplot2)

height_file <- '50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz'

f <- function(x, pos) subset(x, consequence_category == "synonymous" & consequence == "synonymous_variant")
rsid <- read_delim_chunked(variant_file, DataFrameCallback$new(f),
                           delim = "\t", col_names = TRUE)
# Testing only on the height_file
synonymous <- function(x, pos) subset(x, variant %in% rsid$variant) 
df <- read_delim_chunked(height_file, DataFrameCallback$new(synonymous),
                         delim = "\t", col_names = TRUE)

#only keep certain categories that seem interesting
#df <- df %>% dplyr::select(c("variant", "beta", "se", "tstat", "pval"))
df <- inner_join(df, rsid, c('variant' = 'variant'))

delta_data <- read_delim("ancestral_and_derived_with_delta.tsv", delim = "\t", col_names = TRUE)

#remove the NA for delta
delta_data <- delta_data[!is.na(delta_data$delta),]

combined_height_data <- inner_join(df, delta_data, c("rsid" = "input"))

#remove data without values
combined_height_data <- combined_height_data[!(!is.na(combined_height_data$beta) & combined_height_data$beta=="NaN"), ]

#change the type from char to numeric
combined_height_data$beta <- as.numeric(as.character(combined_height_data$beta))


anc_freq <- apply(combined_height_data, MARGIN = 1, FUN = function(x){
  anc_index <- which(codon_freq$codons == tolower(x[38]))
  return(round(codon_freq$freq[anc_index], 2))
})

combined_height_data$anc_freq <- anc_freq

der_freq <- apply(combined_height_data, MARGIN = 1, FUN = function(x){
  der_index <- which(codon_freq$codons == tolower(x[39]))
  return(round(codon_freq$freq[der_index], 2))
})


combined_height_data$der_freq <- der_freq


change_percent <- apply(combined_height_data, MARGIN = 1, FUN = function(x){
  anc_index <- which(codon_freq$codons == tolower(x[38]))
  der_index <- which(codon_freq$codons == tolower(x[39]))
  return((codon_freq$freq[anc_index] - codon_freq$freq[der_index])/codon_freq$freq[anc_index])
})


combined_height_data$change_percentage <- change_percent

write.table(combined_height_data, file = "height_data_percent_change.tsv", row.names = FALSE, 
            col.names = TRUE, sep = "\t")

#lin regression model

combined_height_data$beta <- ifelse(combined_height_data$ref != combined_height_data$ancestral_allele,
                                    -combined_height_data$beta, combined_height_data$beta)

betas <- combined_height_data$beta
mean(betas)
delta_percentage <- combined_height_data$change_percentage

lin_mod <- lm(betas ~ delta_percentage)
summary(lin_mod)

ggplotRegression(lin_mod)
#try partitioning the pval <0.05
combined_height_data$pval <- as.numeric(as.character(combined_height_data$pval))

betas <- combined_height_data[combined_height_data$pval < 0.05, ]$beta

delta_percentage <- combined_height_data[combined_height_data$pval < 0.05, ]$change_percentage

summary(lm(betas ~ delta_percentage))


lin_mod <- lm(betas ~ delta_percentage)

ggplotRegression(lin_mod)

#pval <0.01
betas <- combined_height_data[combined_height_data$pval < 0.01, ]$beta

delta_percentage <- combined_height_data[combined_height_data$pval < 0.01, ]$change_percentage


summary(lm(betas ~ delta_percentage))
lin_mod <- lm(betas ~ delta_percentage)

ggplotRegression(lin_mod)

#pval <0.001
betas <- combined_height_data[combined_height_data$pval < 0.001, ]$beta

delta_percentage <- combined_height_data[combined_height_data$pval < 0.001, ]$change_percentage


summary(lm(betas ~ delta_percentage))

lin_mod <- lm(betas ~ delta_percentage)

ggplotRegression(lin_mod)

#pval < 0.0001
betas <- combined_height_data[combined_height_data$pval < 0.0001, ]$beta

delta_percentage <- combined_height_data[combined_height_data$pval < 0.0001, ]$change_percentage

summary(lm(betas ~ delta_percentage))
lin_mod <- lm(betas ~ delta_percentage)

ggplotRegression(lin_mod)

#pval < 0.00001
betas <- combined_height_data[combined_height_data$pval < 0.00001, ]$beta

delta_percentage <- combined_height_data[combined_height_data$pval < 0.00001, ]$change_percentage

summary(lm(betas ~ delta_percentage))
lin_mod <- lm(betas ~ delta_percentage)

ggplotRegression(lin_mod)

#pval < 0.000001
betas <- combined_height_data[combined_height_data$pval < 0.000001, ]$beta

delta_percentage <- combined_height_data[combined_height_data$pval < 0.000001, ]$change_percentage

summary(lm(betas ~ delta_percentage))
lin_mod <- lm(betas ~ delta_percentage)

ggplotRegression(lin_mod)

#pval < 0.0000001
betas <- combined_height_data[combined_height_data$pval < 0.0000001, ]$beta

delta_percentage <- combined_height_data[combined_height_data$pval < 0.0000001, ]$change_percentage

summary(lm(betas ~ delta_percentage))
lin_mod <- lm(betas ~ delta_percentage)
ggplotRegression(lin_mod)
