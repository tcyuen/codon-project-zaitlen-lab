library(readr)
library(dplyr)
library(ggplot2)

allele_freq <- read_tsv("final_combined_height_data.tsv", col_names = TRUE)
#get the derived AF 
derived_allele_freq <- ifelse(allele_freq$ref == allele_freq$ancestral_allele,
                              allele_freq$minor_AF.x, 1 - allele_freq$minor_AF.x)

allele_freq$derived_AF <- derived_allele_freq

write.table(allele_freq, file = "combined_height_data_derived_AF.tsv", 
            col.names = T,
            row.names = F,
            sep = "\t",)

allele_freq$bins <- cut(allele_freq$derived_AF, breaks = 20)

#example of subsetting the data
subset(allele_freq, bins == "(-0.001,0.05]")
allele_freq[allele_freq$bins == "(0.05,0.1]",]

#plot the data
ggplot(allele_freq, aes(x = bins, y = delta)) +
  stat_summary(fun.y = "mean", geom = "bar")

#mean of derived AF delta
mean(allele_freq[allele_freq$bins == "(-0.001,0.05]",]$delta)

#things TO DO: 
# run permutation test on each bin cutoff
# change the x-axis and y-axis label 
# add a title to the plot
