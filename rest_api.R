library(httr)
library(jsonlite)
library(xml2)
library(dplyr)
library(purrr)
library(magrittr)
library(tidyr)
library('readr')

# code used to practice how the things work. 

# truncated_data <- paste(shQuote(rsid$rsid[52901:53000], type = "cmd"), collapse=", ")
# 
# data_to_include <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
# 
# server <- "https://grch37.rest.ensembl.org"
# ext <- "/vep/human/id"
# cod <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_to_include)
# stop_for_status(cod)
# data <- fromJSON(toJSON(httr::content(cod)))
# data$id <- as.character(data$id)
# class(data$id)
# server <- "https://grch37.rest.ensembl.org"
# ext <- "/variation/homo_sapiens"
# aa <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_to_include)
# 
# stop_for_status(aa)
# 
# data_anc <-  fromJSON(toJSON(httr::content(aa)))
# class(ancestral_allele$input)




chunk_aa_query <- function(start,stop){
  truncated_data <- paste(shQuote(rsid$rsid[start:stop], type = "cmd"), collapse=", ")
  data_json <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
  server <- "https://grch37.rest.ensembl.org"
  ext <- "/variation/homo_sapiens"
  aa <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_json)
  
  stop_for_status(aa)
  data_anc <- data.frame(t(sapply(httr::content(aa),c)))

  
  data_anc <- t(data_anc)
  data_anc <- as.data.frame(data_anc)
  
  ancestral_allele <- data.frame()
  for (i in 1:length(data_anc$V1)) {
    print(i)
    
    if (length(data_anc[[1]][[i]]$mappings[[1]]$ancestral_allele) == 0) {
      
      name <- data.frame("input" = data_anc[[1]][[i]]$name, 
                         "ancestral_allele" = "")
      ancestral_allele <- rbind(ancestral_allele, name)
    }
    else{
      name <- data.frame("input" = data_anc[[1]][[i]]$name, 
                         "ancestral_allele" = data_anc[[1]][[i]]$mappings[[1]]$ancestral_allele)
      ancestral_allele <- rbind(ancestral_allele, name)
    }
  }
  return(ancestral_allele)
}

chunk_cod_query <- function(start,stop){
  truncated_data <- paste(shQuote(rsid$rsid[start:stop], type = "cmd"), collapse=", ")
  data_json <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
  server <- "https://grch37.rest.ensembl.org"
  ext <- "/vep/human/id"
  cod <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_json)
  stop_for_status(cod)

  query_data <- fromJSON(toJSON(httr::content(cod))) #puts into data.frame
  query_data <- query_data %>% dplyr::select(c("input", "transcript_consequences", "allele_string", "id"))
  query_data$id <- as.character(query_data$id)
  return(query_data)
}



# use this if you get a simple nested list back, otherwise inspect its structure
#head(data.frame(t(sapply(httr::content(r),c))))
# data <- fromJSON(toJSON(httr::content(r)))
# class(data)
# df <- data.frame(t(sapply(httr::content(r),c)))
# df$X1[[1]]$transcript_consequences[[1]]$codons
# df[[3]][[1]]$allele_string



#run the for loop to get the codons and anc_allele from ensembl
#might break at points, so just start the for loop at the point it stopped. 
#save it all to data.frame

all_data <- data.frame()

for (i in seq(from = 1, to = length(rsid$rsid), by = 200)){
  tryCatch({
  print(i)
  if (i+199 > length(rsid$rsid) ) {
    cod_data <- chunk_cod_query(i, length(rsid$rsid))
    aa_data <- chunk_aa_query(i, length(rsid$rsid))
    new_data <- inner_join(cod_data, aa_data, c("id"="input"))
  }
  else{
    cod_data <- chunk_cod_query(i, (i+199))
    aa_data <- chunk_aa_query(i, (i+199))
    new_data <- inner_join(cod_data, aa_data, c("id"="input"))
  }
  
  all_data <- rbind(all_data, new_data)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}



#get rid of snps with more than two alleles
all_data_filtered <- subset(all_data, nchar(allele_string) == 3)


#extract codons from list
codons <- data.frame()
for (i in 1:length(all_data_filtered$id)) {
  extracted <- map(all_data_filtered$transcript_consequences[i], "codons")
  unlisted <- unlist(extracted)
  if (length(unlisted) != 0 ) {
    parsed <- data.frame(id = c(all_data_filtered$id[i]),codons = unlisted)
    codons <- rbind(codons, parsed)
  }
}

#remove duplicate values
no_dupe_codons <- codons[!duplicated(codons[c("id","codons")]),]

#add to existing data
final_data <- inner_join(all_data_filtered, no_dupe_codons, c("id" = "id"))

#remove dupes again
final_data <- final_data[!duplicated(final_data[c("id","codons")]),]

#remove data with no ancestral allele
final_data <- final_data[!(!is.na(final_data$ancestral_allele) & final_data$ancestral_allele==""), ]

#format the data, remove the column transcript_consequences, change types
final_data <- final_data %>% select(!c(transcript_consequences))
final_data$input <- as.character(final_data$input)
final_data$allele_string <- as.character((final_data$allele_string))

#write out the data to save it
write.table(final_data, file = "ancestral_and_codon_from_rsid.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# read the file 
data_from_file <- read_delim(file = "ancestral_and_codon_from_rsid.tsv", delim = "\t", col_names = TRUE)

#checks if both data.frames are equal (from dplyr)
all_equal(data_from_file, final_data)
