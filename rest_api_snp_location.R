library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(jmuOutlier)
library(coin)

data_from_file <- read_tsv("ancestral_and_codon_from_rsid.tsv", col_names = TRUE)

pos_data <- data.frame()


for (i in seq(from = 1, to = length(rsid$rsid), by = 300)) {
  tryCatch({
    print(i)
    if (i+199 > length(rsid$rsid) ) {
      truncated_data <- paste(shQuote(rsid$rsid[i:length(rsid$rsid)], type = "cmd"), collapse=", ")
      query <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
      server <- "https://grch37.rest.ensembl.org"
      ext <- "/vep/human/id"
      cod <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = query)
      stop_for_status(cod)
      query_data <- fromJSON(toJSON(httr::content(cod))) #puts into data.frame
      query_data <- query_data %>% dplyr::select(c("input", "transcript_consequences", "id"))
      query_data$id <- as.character(query_data$id)
      
      
    }
    else{
      truncated_data <- paste(shQuote(rsid$rsid[i:(i+299)], type = "cmd"), collapse=", ")
      query <- paste('{ "ids" :[' , truncated_data, '] }',collapse = ", ", sep = "")
      server <- "https://grch37.rest.ensembl.org"
      ext <- "/vep/human/id"
      cod <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = query)
      stop_for_status(cod)
      query_data <- fromJSON(toJSON(httr::content(cod))) #puts into data.frame
      query_data <- query_data %>% dplyr::select(c("input", "transcript_consequences", "id"))
      query_data$id <- as.character(query_data$id)
      
      
    }
    pos_data <- rbind(pos_data, query_data)
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}


#unnest the json file
pos_df <- unnest(pos_data,transcript_consequences) %>% 
  subset(consequence_terms == "synonymous_variant") %>% filter(input %in% data_from_file$input)

#change the dataframe from lists to vectors
pos_df2 <- pos_df %>% select(input, cds_start, cds_end, transcript_id)

pos_df2 <- as.data.frame(lapply(pos_df2, unlist))

#get the smallest position of the rsid
smallest_pos <- pos_df2 %>% group_by(input) %>% slice_min(cds_start, n = 1, with_ties = F)

delta_and_location <- inner_join(delta_data, smallest_pos, c("input" = "input"))

delta_and_location$first <- ifelse(delta_and_location$cds_start < 120, TRUE, FALSE)
delta_and_location <- delta_and_location[!is.na(delta_and_location$delta),]

ggplot(delta_and_location, aes(x = first, y = delta)) +
  geom_point() + stat_smooth(method = "lm") +
  labs(title = "SNP Location and Codon", x = "First 40 codons", y = "Delta") +
  theme(plot.title = element_text(size=22))

length(delta_and_location$first[delta_and_location$first == TRUE])


mean(delta_and_location$delta[delta_and_location$first == F])


lm(delta_and_location$delta ~ delta_and_location$first)

perm.test(delta_and_location$delta[delta_and_location$first == T],
          delta_and_location$delta[delta_and_location$first == F], num.sim = 10000)


kruskal.test(delta_and_location$delta ~ delta_and_location$first)
wilcox.test(delta_and_location$delta ~ delta_and_location$first)
delta_and_location$first <- as.factor(delta_and_location$first)


wilcox_test(delta~first, data = delta_and_location)
oneway_test(delta~first, data = delta_and_location)
independence_test(delta~first, data = delta_and_location)
