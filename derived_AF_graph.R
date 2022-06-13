library(readr)
library(dplyr)
library(ggplot2)
library(Hmisc)

allele_freq <- read_tsv("combined_height_data.tsv", col_names = TRUE)

#get the derived AF 
derived_allele_freq <- ifelse(allele_freq$ref == allele_freq$ancestral_allele,
                              allele_freq$minor_AF.x, 1 - allele_freq$minor_AF.x)

allele_freq$derived_AF <- derived_allele_freq

write.table(allele_freq, file = "combined_height_data_derived_AF.tsv", 
            col.names = T,
            row.names = F,
            sep = "\t",)

cutoffs <- seq(from = 5, to = 100, by = 5)
percentages <- as.character(cutoffs)

allele_freq$bins <- cut(allele_freq$derived_AF, breaks =  
                          quantile(allele_freq$derived_AF, 
                                   probs = (0:20)/20, na.rm = T), 
                        include.lowest = T, 
                        labels = percentages)


#example of subsetting the data
subset(allele_freq, bins == "(-0.001,0.05]")
allele_freq[allele_freq$bins == "(0.05,0.1]",]

sum(is.na(allele_freq$bins))

#plot the data (boxplot)
ggplot(allele_freq, aes(x = bins, y = delta, fill = bins)) +
  geom_boxplot() +
  stat_summary(fun = "mean", col = "red") +
  theme(legend.position = "none", plot.title = element_text(size = 22)) + 
  labs(title = "Allele Frequency and Delta Codon Frequency",
       x = "Quantiles", y = "Delta Codon Frequency")


