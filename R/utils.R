source("R/prepare_data.R")
library(tidyverse)
library(preprocessCore)

# Generic aggregation function for columns of numeric or character type
#
agg_fun <- function(values){
  if(all(is.na(values))){
    return(NA)
  } else {
    if(class(values)=="character"){
      return(paste(unique(values),collapse=";"))
    } else {
      return(mean(values,na.rm = TRUE))
    }
  }
}

# Convenience functions for loading cached data
load_mq_annotated <- function(){
  if ( !file.exists("cache/mq_data_annotated.rds")){
    make_mq_annotated()
  }
  readRDS("cache/mq_data_annotated.rds")
}

load_interpro <- function(){
  if( !file.exists("cache/interpro_scan.rds") ){
    make_interpro()
  }
  readRDS("cache/interpro_scan.rds")
}

load_ceph_specific <- function(){
  if ( !file.exists("cache/cephalopod_specific.rds") ){
    make_ceph_specific()
  }
  readRDS("cache/cephalopod_specific.rds")
}

# For ease of use in downstream code we create the following function which transforms the LFQ raw data to log2 normalized
#
lognorm_lfq <- function(mq_data){
  lfq_mq_data <- mq_data %>% 
    select(contains('LFQ intensity ')) %>% 
    log2() %>% 
    as.matrix() 
  
  lognorm_mq_data <- normalize.quantiles(lfq_mq_data)
  colnames(lognorm_mq_data) <- colnames(lfq_mq_data)
  lognorm_mq_data <-cbind(prot_id = mq_data$prot_id, as.data.frame(lognorm_mq_data))  
  
  mq_data %>% select(-contains('LFQ intensity ')) %>% 
    left_join(lognorm_mq_data,by="prot_id")
}

# Clean sample identifiers so that they are reduced to codes in the sample table
# input data should have original sample names encoded in a variable "sample"
#
extract_sample_codes <- function(sample_names){
  str_match(sample_names,pattern = "([0-9][AP][0-9]|S[0-9][AB]|M[1-3])")[,2]
}


clean_desc <- function(data,pattern,replacement){
  data %>% mutate(Description = ifelse(grepl(Description,pattern=pattern,ignore.case = TRUE),replacement,Description))
}

load_sscr_proteins <- function(){
  mq_data <- lognorm_lfq(load_mq_annotated())
  ceph_specific <- load_ceph_specific()
  
  sscr_mq_data <- mq_data %>% 
    mutate(cysteine_density = cysteine_count/`Sequence length`) %>% 
    filter(Signal == "Y") %>% # 385 protein groups
    filter( (max_cysteines >= 5) | ( cysteine_density > 0.05 ) ) %>%
    filter(`Sequence length` < 200 ) %>% 
    mutate(is_ceph_specific = ifelse( prot_id %in% ceph_specific , TRUE, FALSE))
  
  
  
  # Write out sscr proteins for import to geneious
  #
  #write_tsv(sscr_mq_data,path = "raw_data/sscr_proteins/sscr_proteins.tsv")
  sscr_mq_data
}


serine_protease_treedata <- function(){
  species_abbreviations <- c(
    "Octopus bimaculoides"="Obi",
    "Doryteuthis opalescens"="Dop",
    "Sepia latimanus"="Slat",
    "Octopus kaurna"="Okau",
    "Hapalochlaena maculosa"="Hmac",
    "Heterololigo bleekeri"="Hble",
    "Euprymna scolopes"="Esco",
    "Abdopus aculeatus"="Aacu",
    "Octopus cyanea"="Ocya",
    "Pareledone turqueti"="Ptur",
    "Sepioteuthis australis"="Saus",
    "Loliolus noctiluca"="Lnoc",
    "Sepia pharaonis"="Spha",
    "Crassostrea gigas"="Cgig",
    "Octopus vulgaris"="Ovul",
    "Sepia officinalis"="Sepof",
    "Lottia gigantea"="Lotgi"
      )

  species_clades <- c(
    "Octopus bimaculoides"="Octopus",
    "Doryteuthis opalescens"="Squid",
    "Sepia latimanus"="Cuttlefish",
    "Octopus kaurna"="Octopus",
    "Hapalochlaena maculosa"="Octopus",
    "Heterololigo bleekeri"="Squid",
    "Euprymna scolopes"="Squid",
    "Abdopus aculeatus"="Octopus",
    "Octopus cyanea"="Octopus",
    "Pareledone turqueti"="Octopus",
    "Sepioteuthis australis"="Squid",
    "Loliolus noctiluca"="Squid",
    "Sepia pharaonis"="Cuttlefish",
    "Octopus vulgaris"="Octopus",
    "Sepia officinalis"="Cuttlefish",
    "Lottia gigantea"="Gastropod",
    "Crassostrea gigas"="Bivalve"
      )
  
  
    
  tip_data <- read_tsv("raw_data/serine_proteases/sp1_200_molluscs_idmap.tsv",col_names = c("ID","Species")) %>% 
    mutate(original_id=ID) %>% 
    mutate(ID = str_replace(ID,pattern = ">",replacement = "")) %>% 
    mutate(ID = ifelse( grepl(ID, pattern = "\\|.*\\|") , str_match(ID, pattern = "\\|(.*)\\|")[,2], ID)) %>%
    mutate(ID = str_replace_all(ID,pattern = "\\|",replacement = "_")) %>% 
    mutate(ID = str_replace_all(ID,pattern = "\\:",replacement = "_")) %>% 
    mutate(Order = species_clades[Species]) %>% 
    group_by(Species) %>% 
    mutate(simple_label = paste(species_abbreviations[Species],row_number(),sep="_")) %>% 
    ungroup()
  tip_data
}
