library(httr)
library(jsonlite)
library(xml2)

# truncated_data <- paste(shQuote(rsid$rsid[1:300], type = "cmd"), collapse=", ")
# 
# data_to_include <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
# 
# server <- "https://grch37.rest.ensembl.org"
# ext <- "/vep/human/id"
# r <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_to_include)
# stop_for_status(r)
# data <- fromJSON(toJSON(httr::content(r)))

chunk_query <- function(start,stop){
  truncated_data <- paste(shQuote(rsid$rsid[start:stop], type = "cmd"), collapse=", ")
  data_json <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
  server <- "https://grch37.rest.ensembl.org"
  ext <- "/vep/human/id"
  r <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = data_json)
  stop_for_status(r)
  query_data <- fromJSON(toJSON(httr::content(r)))
  return(query_data)
}

for (i in seq(from = 1, to = length(rsid$rsid), by = 200)){
  if (i == 1) {
    all_data <- chunk_query(i, (i+199))

    
  }
  
  if (i+200 > length(rsid$rsid) ) {
    
    new_data <- chunk_query(i, length(rsid$rsid))
    all_data <- rbind(all_data, new_data)
  }
  else{
    new_data <- chunk_query(i, (i+199))
    all_data <-rbind(all_data, new_data)
  }
  
  if (as.numeric(r[["headers"]][["x-ratelimit-reset"]]) < 100) {
    Sys.sleep(100)
  }
  
}

for (i in seq(from = 1, to = length(rsid$rsid), by = 200)){
  truncated_data <- paste(shQuote(rsid$rsid[i:(i+199)], type = "cmd"), collapse=", ")
  #print(truncated_data)
  
  }
# use this if you get a simple nested list back, otherwise inspect its structure
#head(data.frame(t(sapply(httr::content(r),c))))
data <- fromJSON(toJSON(httr::content(r)))
class(data)
df <- data.frame(t(sapply(httr::content(r),c)))
df$X1[[1]]$transcript_consequences[[1]]$codons
df[[3]][[1]]$allele_string
