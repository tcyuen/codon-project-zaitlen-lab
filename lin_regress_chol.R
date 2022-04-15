library('readr')
library('dplyr')

chol_file <- '30690_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz'

# Testing only on the height_file
synonymous <- function(x, pos) subset(x, variant %in% rsid$variant) 
df <- read_delim_chunked(chol_file, DataFrameCallback$new(synonymous),
                         delim = "\t", col_names = TRUE)

#only keep certain categories that seem interesting
df <- df %>% dplyr::select(c("variant", "beta", "se", "tstat", "pval"))
df <- inner_join(df, rsid, c('variant' = 'variant'))

#adding the ancestral allele to the data 
delta_data <- read_delim("ancestral_and_derived_with_delta.tsv", delim = "\t", col_names = TRUE)

combined_chol_data <- inner_join(df, delta_data, c("rsid" = "input"))

#get rid of empty values
combined_chol_data <- combined_chol_data[!(!is.na(combined_chol_data$beta) & combined_chol_data$beta=="NaN"), ]

#change the type from char to numeric
combined_chol_data$beta <- as.numeric(as.character(combined_chol_data$beta))

write.table(combined_chol_data, file = "final_combined_chol_data.tsv", row.names = FALSE, 
            col.names = TRUE, sep = "\t")



#change the sign if the ref != ancestral
combined_chol_data$beta <- ifelse(combined_chol_data$ref != combined_chol_data$ancestral_allele,
                                    -combined_chol_data$beta, combined_chol_data$beta)
betas <- combined_chol_data$beta

delta_codon_freqs <- combined_chol_data$delta

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

#run linear model on all beta and delta freqs
summary(lm(betas ~ delta_codon_freqs))

lin_mod <- lm(betas ~ delta_codon_freqs)

ggplotRegression(lin_mod)


#try partitioning the pval < 0.05
combined_chol_data$pval <- as.numeric(as.character(combined_chol_data$pval))

betas <- combined_chol_data[combined_chol_data$pval < 0.05, ]$beta

delta_codon_freqs <- combined_chol_data[combined_chol_data$pval < 0.05, ]$delta

summary(lm(betas ~ delta_codon_freqs))
lin_mod <- lm(betas ~ delta_codon_freqs)

ggplotRegression(lin_mod)

#pval <0.001
betas <- combined_chol_data[combined_chol_data$pval < 0.001, ]$beta

delta_codon_freqs <- combined_chol_data[combined_chol_data$pval < 0.001, ]$delta

summary(lm(betas ~ delta_codon_freqs))
lin_mod <- lm(betas ~ delta_codon_freqs)

ggplotRegression(lin_mod)

#pval < 0.0001
betas <- combined_chol_data[combined_chol_data$pval < 0.0001, ]$beta

delta_codon_freqs <- combined_chol_data[combined_chol_data$pval < 0.0001, ]$delta

summary(lm(betas ~ delta_codon_freqs))
lin_mod <- lm(betas ~ delta_codon_freqs)

ggplotRegression(lin_mod)

#pval < 0.00001
betas <- combined_chol_data[combined_chol_data$pval < 0.00001, ]$beta

delta_codon_freqs <- combined_chol_data[combined_chol_data$pval < 0.00001, ]$delta

summary(lm(betas ~ delta_codon_freqs))
lin_mod <- lm(betas ~ delta_codon_freqs)

ggplotRegression(lin_mod)

#pval < 0.000001
betas <- combined_chol_data[combined_chol_data$pval < 0.000001, ]$beta

delta_codon_freqs <- combined_chol_data[combined_chol_data$pval < 0.000001, ]$delta

summary(lm(betas ~ delta_codon_freqs))
lin_mod <- lm(betas ~ delta_codon_freqs)

ggplotRegression(lin_mod)

#pval < 0.0000001
betas <- combined_chol_data[combined_chol_data$pval < 0.0000001, ]$beta

delta_codon_freqs <- combined_chol_data[combined_chol_data$pval < 0.0000001, ]$delta

summary(lm(betas ~ delta_codon_freqs))
lin_mod <- lm(betas ~ delta_codon_freqs)

ggplotRegression(lin_mod)
