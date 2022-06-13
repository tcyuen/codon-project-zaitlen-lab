library(biomartr) # for cds
library(dplyr)
library(r2r) #for hash tables
library(data.table)
library(Biostrings) # for codon to amino acid

# check what kinds of databases and genomes are available
aval <- is.genome.available(db = "ensembl", organism = "Homo sapiens", details = TRUE)


human_cds_ensembl <- getCDS(
                            db = "ensembl",
                            organism = "Homo sapiens",
                            release = 75, #GRCh37 build
                            gunzip = FALSE,
                            path = file.path("_ensembl_downloads", "CDS")
                          )


# this is for data downloaded directly from ensembl
# ftp_download_cds <- "_ensembl_downloads/CDS/Homo_sapiens.GRCH37.cds.all.fa.gz"

human_cds <- read_cds(file = human_cds_ensembl, obj.type = "data.table")

# ftp_cds <- read_cds(file = ftp_download_cds, obj.type = "data.table")

# both of the files are the same seqs.
# all_equal(human_cds$seqs, ftp_cds$seqs)
# identical(human_cds$seqs, ftp_cds$seqs)


c1 <- c('a', 't', 'g', 'c')
c2 <- c('a', 't', 'g', 'c')
c3 <- c('a', 't', 'g', 'c')

#make all possible combinations of codons
codons <- apply(expand.grid(c1,c2,c3), 1, function(x) paste0(x, collapse=""))


codon_counts <- hashmap(default = 0)

#iterate through all of the cds_seqs
for (i in 1:length(human_cds$seqs)) {
  #split the seq into substrings of 3
  seq_split <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", human_cds$seqs[i])," ",T)
  for (codon in 1:length(codons)) {
    #iterate through the codons and count how many times they appear in the seq
    counts <- sum(unlist(seq_split) == codons[codon])
    #add the count to the dictionary
    codon_counts[[codons[codon]]] <- codon_counts[[codons[codon]]] + counts
  }
}




keys(codon_counts)
values(codon_counts)



#making the frequency table
#make a dataframe from the codons and their counts
codon_usage <- data.frame(codons = unlist(keys(codon_counts)), counts = unlist(values(codon_counts)))

#find the amino acid for each codon and append it to the dataframe. (uses Biostrings)
codon_usage$aa <- sapply(codon_usage$codons, function(x) GENETIC_CODE[[toupper(x)]])
#sort the dataframe (uses mixedorder() from gtools)
codon_usage <- codon_usage[mixedorder(codon_usage$aa),] 

#calc the relative frequency of the codon for each amino acid. (uses dplyr )
codon_usage <- codon_usage %>% 
  group_by(aa, codons) %>%
  summarize(counts = counts) %>%
  mutate(rel_freq = round(counts/sum(counts), 2))

#export the data
write.table(codon_usage, file = "codon_frequency_table.tsv", row.names = FALSE, 
            col.names = TRUE, sep = "\t")
