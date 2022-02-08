# Makes vector with the number of characters in the allele column test df
ch_count_vec <-nchar(test$allele)

# Turns vector into data frame
count_char <-data.frame(ch_count_vec)

# Combines data frames then filters rsid's with only two alleles
test2 <-cbind(test, count_char)
test_allele <-subset(test2, ch_count_vec == 3)

typeof(test_allele$allele)