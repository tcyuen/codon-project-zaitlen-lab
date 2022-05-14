library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(dplyr)

data_from_file <- read_tsv("ancestral_and_codon_from_rsid.tsv", col_names = TRUE)

pos_data <- data.frame()


for (i in seq(from = 1, to = length(rsid$rsid), by = 300)) {
  tryCatch({
    print(i)
    if (i+199 > length(rsid$rsid) ) {
      truncated_data <- paste(shQuote(rsid$rsid[i:length(rsid$rsid)], type = "cmd"), collapse=", ")
      data_json <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
      server <- "https://grch37.rest.ensembl.org"
      ext <- "/vep/human/id"
      cod <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_json)
      stop_for_status(cod)
      query_data <- fromJSON(toJSON(httr::content(cod))) #puts into data.frame
      query_data <- query_data %>% dplyr::select(c("input", "transcript_consequences", "id"))
      query_data$id <- as.character(query_data$id)
      
      
    }
    else{
      truncated_data <- paste(shQuote(rsid$rsid[i:(i+299)], type = "cmd"), collapse=", ")
      data_json <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
      server <- "https://grch37.rest.ensembl.org"
      ext <- "/vep/human/id"
      cod <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_json)
      stop_for_status(cod)
      query_data <- fromJSON(toJSON(httr::content(cod))) #puts into data.frame
      query_data <- query_data %>% dplyr::select(c("input", "transcript_consequences", "id"))
      query_data$id <- as.character(query_data$id)
      
      
    }
    pos_data <- rbind(pos_data, query_data)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}


