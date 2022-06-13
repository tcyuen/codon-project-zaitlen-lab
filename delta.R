library(readr)
library(tidyr)
library(dplyr)
library(stringr)

# read the file
data_from_file <- read_delim(file = "ancestral_and_codon_from_rsid.tsv", delim = "\t", col_names = TRUE)

#split the codons
split_codons <- separate(data_from_file, col = "codons", sep="/", into = c("codon_1", "codon_2"), remove=TRUE)

find_ancestral <- split_codons[,7:9]

#determine what codon is the ancestral
ancestral <- apply(find_ancestral, MARGIN = 1, FUN = function(x){
  if (str_detect(x[2], x[1])) {
    return(x[2])
  } else if (str_detect(x[3], x[1])) {
    return(x[3])
  } else if (str_detect(x[2], chartr("ATGC", "TACG", x[1]))) {
    return(x[2])
  } else{
    return(x[3])
  }
})

split_codons$ancestral_codon <- ancestral

#determine which codon is derived
get_derived <- split_codons[,c(8,9,12)]

derived <- apply(get_derived, MARGIN = 1, FUN = function(x){
  if (x[1] == x[3]){
    return(x[2])
  } else{
    return(x[1])
  }
})


split_codons$derived_codon <- derived

#find the change codon frequency
#delta = anc_freq - der_freq

codon_freq <- read_delim("codon_frequency_table.tsv", delim = "\t", col_names = TRUE)

#find anc and der freqs
anc_freq <- apply(split_codons, MARGIN = 1, FUN = function(x){
  anc_index <- which(codon_freq$codons == tolower(x[12]))
  return(round(codon_freq$freq[anc_index], 2))
})

split_codons$anc_freq <- anc_freq

der_freq <- apply(split_codons, MARGIN = 1, FUN = function(x){
  der_index <- which(codon_freq$codons == tolower(x[13]))
  return(round(codon_freq$freq[der_index], 2))
})

split_codons$der_freq <- der_freq

# calculate delta codon freq
delta <- apply(split_codons, MARGIN = 1, FUN = function(x){
  anc_index <- which(codon_freq$codons == tolower(x[12]))
  der_index <- which(codon_freq$codons == tolower(x[13]))
  return(round(codon_freq$freq[anc_index] - codon_freq$freq[der_index], 2))
})

split_codons$delta <- delta
split_codons$delta <- as.numeric(split_codons$delta)

#remove the unneded columns
split_codons <- select(split_codons, -c("id", "codon_1", "codon_2"))

write.table(split_codons, file = "ancestral_and_derived_with_delta.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



check_data <- read_delim("ancestral_and_derived_with_delta.tsv", delim = "\t", col_names = TRUE)



all_equal(check_data, split_codons)


