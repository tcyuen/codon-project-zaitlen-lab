# install.packages("biomartr", dependencies = TRUE)

library(biomartr)
library(dplyr)
library(r2r) #for hash tables
library(data.table)
# library(hash)
aval <- is.genome.available(db = "genbank", organism = "Homo sapiens", details = TRUE)

# counting_codons <- function(x){
#   data.table(x)[, .N, keyby = x]
# }
# 
# number <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", ftp_cds$seqs[1])," ",T)
# number2 <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", ftp_cds$seqs[2])," ",T)
# 
# testing <- counting_codons(unlist(number))
# testing2 <- counting_codons(unlist(number2))
# h <- hash(codons, c(0)*length(codons))
# h
# h[testing$x] <- h[testing$x] + testing$N
# h[testing2$x] <- h[testing2$x] + testing2$N
# class(c(h[testing$x]))
# listGenomes(db = "ensembl", type = "all", subset = NULL, details = FALSE)
# 
# human_cds_ensembl <- getCDS(
#                             db = "ensembl",
#                             organism = "Homo sapiens",
#                             release = 75,
#                             gunzip = FALSE,
#                             path = file.path("_ensembl_downloads", "CDS")
#                           )
# 
# another_cds_human <- getCDS(
#                             db = "ensembl",
#                             organism = "Homo sapiens",
#                             release = 100,
#                             gunzip = FALSE,
#                             path = file.path("_ensembl_downloads", "CDS")
#                           )

ftp_download_cds <- "_ensembl_downloads/CDS/Homo_sapiens.GRCH37.cds.all.fa.gz"

#human_cds <- read_cds(file = human_cds_ensembl, obj.type = "data.table")

ftp_cds <- read_cds(file = ftp_download_cds, obj.type = "data.table")

all_equal(human_cds$seqs, ftp_cds$seqs)
identical(human_cds$seqs, ftp_cds$seqs)


c1 <- c('a', 't', 'g', 'c')
c2 <- c('a', 't', 'g', 'c')
c3 <- c('a', 't', 'g', 'c')

codons <- apply(expand.grid(c1,c2,c3), 1, function(x) paste0(x, collapse=""))


codon_counts <- hashmap(default = 0)

#iterate through all of the cds_seqs
for (i in 1:length(ftp_cds$seqs)) {
  #split the seq into substrings of 3
  seq_split <- strsplit(gsub("([[:alnum:]]{3})", "\\1 ", ftp_cds$seqs[i])," ",T)
  for (codon in 1:length(codons)) {
    #iterate through the codons and count how many times they appear in the seq
    counts <- sum(unlist(seq_split) == codons[codon])
    codon_counts[[codons[codon]]] <- codon_counts[codons[codon]] + counts
  }
}




keys(codon_counts)
values(codon_counts)
