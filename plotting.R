library(gridExtra)
library(ggplot2)
library(gtools)
library(Biostrings)
library(readr)
library(dplyr)


#plot for codon usage table

freq_df <- as.data.frame(codon_freq[c(1,2,4)])
freq_df$codons <- toupper(freq_df$codons)
freq_df$aa[4:64] <- sapply(freq_df$aa[4:64], function(x) AMINO_ACID_CODE[[x]])

ggplot(freq_df[c(4:11, 14:30),], aes(factor(codons), freq, fill = aa)) + 
  geom_bar(stat="identity", position = "dodge") +
  labs(title="Codon Frequency by Amino Acid", x="Codon", y="Frequency", fill = "Amino Acid") + 
  theme(plot.title = element_text(size=22),
        axis.text = element_text(face="bold"))





#total change in codon frequency
hist(combined_data$delta, xlab = "Change in Codon Frequency", 
     main = "Distribution of Change in Codon Frequency", 
     col = "royalblue2")
abline(v = mean(combined_data$delta, na.rm = T), col = "red", lwd = 3)
text(x = mean(combined_data$delta, na.rm = T) * -8,                   # Add text for mean
     y = 11000,
     paste("Mean =", round(mean(combined_data$delta, na.rm = T), 2)),
     col = "red",
     cex = 2)





# statistically significant lm() height p< 0.001
#
combined_height_data <- read_tsv("final_combined_height_data.tsv")

#change the sign if the ref != ancestral
combined_height_data$beta <- ifelse(combined_height_data$ref != combined_height_data$ancestral_allele,
                                    -combined_height_data$beta, combined_height_data$beta)

betas <- combined_height_data[combined_height_data$pval < 0.001, ]$beta

delta_codon_freqs <- combined_height_data[combined_height_data$pval < 0.001, ]$delta

fit <- lm(betas ~ delta_codon_freqs)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = "Height (P < 0.001)",
       subtitle = paste("P =",signif(summary(fit)$coef[2,4], 5)),
       x = "Delta Codon Freq", y = "Effect Size of SNP") +
  theme(plot.title = element_text(size=22)) +
  annotate("text", x = 0, y = 250, 
           label = paste("Slope = ", signif(fit$coef[[2]], 5)),
           size = 10, col = "red")


# lm() chol p < 0.001
combined_chol_data <- read_tsv("final_combined_chol_data.tsv")

#change the sign if the ref != ancestral
combined_chol_data$beta <- ifelse(combined_chol_data$ref != combined_chol_data$ancestral_allele,
                                  -combined_chol_data$beta, combined_chol_data$beta)

betas <- combined_chol_data[combined_chol_data$pval < 0.001, ]$beta

delta_codon_freqs <- combined_chol_data[combined_chol_data$pval < 0.001, ]$delta

summary(lm(betas ~ delta_codon_freqs))
fit <- lm(betas ~ delta_codon_freqs)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(title = "Cholesterol (P < 0.001)",
       subtitle = paste(" P =",signif(summary(fit)$coef[2,4], 5)),
       x = "Delta Codon Freq", y = "Effect Size of SNP") +
  theme(plot.title = element_text(size=22)) +
  annotate("text", x = 0, y = 250, 
           label = paste("Slope = ", signif(fit$coef[[2]], 5)),
           size = 10, col = "red")




#adding the delta to the data 
delta_and_location <- read_delim("ancestral_and_derived_with_delta.tsv", delim = "\t", col_names = TRUE)

#remove any empty values for delta
delta_and_location <- delta_data[!is.na(delta_data$delta),]

#split with codons depending on location of snp
delta_and_location$first <- ifelse(delta_and_location$cds_start < 120, "Beginning", "End")
delta_and_location <- delta_and_location[!is.na(delta_and_location$delta),]

#plot the data split on location. 
ggplot(delta_and_location, aes(x = first, y = delta)) +
  geom_boxplot(fill = "green3", alpha = 0.3) +
  labs(title = "SNP Location and Codon Frequency", x = "Location", y = "Delta") +
  theme(plot.title = element_text(size=22)) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")
