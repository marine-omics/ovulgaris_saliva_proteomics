# Combine data from the following sources into a master table for use by other scripts
#
# MaxQuant Searches
# InterProScan searches
# blast vs nr,swissprot,arachnoserver, uniprot_toxins
# Trinotate annotations
#

library(splitstackshape)
library(tidyverse)

## BLAST annotations. Based on 
# qaccver saccver ident evalue mismatch gapopen qstart qend sstart send stitle

read_blast <- function(file_name, columns = c("ID"=1,"Hit_ID"=2,"Evalue"=4,"Description"=11)){
  # place this behind the column argument to serve as "default if no other columns are given: c(1,2,5,12)
  blast_result <- read_tsv(file_name,col_names = FALSE)
  blast_result <- blast_result[,columns]
  names(blast_result) <- names(columns)
  blast_result <- filter(blast_result, Evalue < 1e-5)
  # There is sometimes more than one blast hit per protein so we need to pick the best E value for each
  blast_result <- blast_result %>% group_by(Hit_ID) %>% top_n(1,1/Evalue)
  blast_result
}

make_mq_annotated <- function(){
  blast_arachnoserver <- read_blast("raw_data/maxquant/all_proteins.arachnoserver.blastp") %>% 
    select(prot_id=ID,Hit_ID,Evalue,Description) %>% add_column(database="arachnoserver")

  blast_toxprot <- read_blast("raw_data/maxquant/all_proteins.toxprot.blastp") %>% 
    select(prot_id=ID,Hit_ID,Evalue,Description) %>% add_column(database="toxprot")

  blast_sp <- read_blast("raw_data/maxquant/all_proteins.sp.blastp", columns = c("ID"=1,"Hit_ID"=2,"Evalue"=6,"Description"=13)) %>% 
    select(prot_id=ID,Hit_ID,Evalue,Description) %>% add_column(database="sp")

  blast_nr <- read_blast("raw_data/maxquant/all_proteins.nr_besthit.blastp", columns = c("ID"=1,"Hit_ID"=2,"Evalue"=6,"Description"=13)) %>% 
    select(prot_id=ID,Hit_ID,Evalue,Description) %>% add_column(database="nr") 

  db_ranks <- c("sp" = 1,"arachnoserver"=2,"toxprot"=3,"nr"=4)
  
  all_blast <- rbind(blast_arachnoserver,blast_toxprot,blast_sp,blast_nr)
  best_blast <- all_blast %>% group_by(prot_id) %>% 
    top_n(1,1/Evalue) %>% 
    distinct() %>% 
    mutate(database_rank = db_ranks[database]) %>% 
    top_n(1,-database_rank)

  signalp_id_map <- read_tsv("raw_data/maxquant/all_proteins_ids.txt",col_names = c("prot_id","name"))
  signalp <- read_tsv("raw_data/maxquant/all_proteins.signalp",skip = 2,
                      col_names = c("name","D","Signal")) %>% left_join(signalp_id_map,by="name")
  
  trinotate <- read_tsv("raw_data/Transcriptome/hpc/trinity/trinotate_annotation_report.xls") %>% mutate(prot_id=paste("lcl|",prot_id,sep=""))

  cysteine_counts <- read_tsv("raw_data/maxquant/OVulgarisMQ_20172206_countaa.tsv",col_names = c("id","cysteine_count"))
  mwpi <- read_tsv("raw_data/maxquant/OVulgarisMQ_20172206_mwpi.tsv",col_names = c("id","mw","pi"))
  cysteine_density <- read_tsv("raw_data/maxquant/OVulgarisMQ_20172206_cysteine_density.tsv",col_names = c("id","max_cysteines","knot"))

  mq_data <- read_tsv("raw_data/maxquant/proteinGroupsMQ6All.txt",na = "0") %>% rownames_to_column("group_id")
  mq_data_long <- cSplit(mq_data,splitCols = c('Protein IDs'),sep = ";",direction = "long") %>% 
    rename(prot_id=`Protein IDs`)

  mq_data_long_annotated <- mq_data_long %>% 
    left_join(trinotate,by="prot_id") %>% 
    left_join(cysteine_counts,by=c("prot_id"="id")) %>% 
    left_join(cysteine_density,by=c("prot_id"="id")) %>% 
    left_join(mwpi,by=c("prot_id"="id")) %>% 
    left_join(best_blast) %>% 
    left_join(signalp) %>% 
    filter(!grepl(prot_id,pattern="CON_")) %>% 
    filter(!grepl(prot_id,pattern="REV_")) %>% 
    select(prot_id,
           group_id,
           gene_id,
           transcript_id,
           Hit_ID,
           Evalue,
           Description,
           database,
           Signal,
           Pfam,TmHMM,eggnog,Kegg,gene_ontology_blast,gene_ontology_pfam, # from Trinotate
           cysteine_count,max_cysteines,mw,pi,
          `Sequence length`,contains('iBAQ'),contains('LFQ')) %>% 
    distinct() # This final line is because there are unfortunately a few duplicates in the maxquant database
  

  saveRDS(mq_data_long_annotated,"cache/mq_data_annotated.rds")
}

make_interpro <- function(){

# We save interpro_scan separately from the main table. 
# This allows us to preserve the idea that there is one row per protein in the main table

  interpro_scan <- read_tsv("raw_data/maxquant/all_proteins.fasta.tsv",
                            col_names = c("prot_id","md5","length","ipr_analysis",
                                          "ipr_accession","ipr_description",
                                          "ipr_start","ipr_end","ipr_score",
                                          "ipr_status","ipr_date"
                            ))

# interpro_scan_agg <- interpro_scan %>% 
#   group_by(prot_id,ipr_analysis) %>% 
#   summarise(ipr_accession = agg_fun(ipr_accession), 
#             ipr_description = paste(unique(ipr_description),collapse=";")) %>% 
#   spread(ipr_analysis,ipr_description)

  saveRDS(interpro_scan,"cache/interpro_scan.rds")
}

make_ceph_specific <- function(){

  # Taxon Specific Proteins
  
  # To find taxon specific proteins we performed a BLAST search for every O. vulgaris protein identified by MaxQuant against the ncbi nr database. For these searches we allowed up to 500 matches to be returned for each query and retrieved ncbi taxon ids for each match. We then defined a taxon specific protein as one where there were no matches (E<1x10-3) or where all of the matches had taxon ids from within the specified taxon (octopodiformes or cephalopoda).
  # 
  # Read in the blast search results and taxon ids
  
  nr_blast <- read_tsv("raw_data/maxquant/all_proteins.nr.blastp",col_names = c("qaccver","saccver","staxids","sscinames","ident","evalue","mismatch","gapopen","qstart","qend","sstart","send","stitle"))
  
  octopodiformes_ids <- read_tsv("raw_data/maxquant/octopodiformes_taxids.txt",col_names = c("id")) %>% pull(id)
  cephalopoda_ids <- read_tsv("raw_data/maxquant/cephalopoda_taxids.txt",col_names = c("id")) %>% pull(id)
  
  all_ids <- read_tsv("raw_data/maxquant/all_proteins_ids.txt",col_names = FALSE) %>% pull(X1)
  all_ids_with_hits <- nr_blast %>% filter(grepl(qaccver,pattern="lcl")) %>% pull(qaccver) %>% unique()
  
  not_octopod_specific <- nr_blast %>% filter(!(staxids %in% octopodiformes_ids)) %>% filter(evalue<0.001) %>% pull(qaccver) %>% unique()
  not_ceph_specific <- nr_blast %>% filter(!(staxids %in% cephalopoda_ids)) %>% filter(evalue<0.001) %>% pull(qaccver) %>% unique()
  
  octopod_specific <- setdiff(all_ids,not_octopod_specific)
  octopod_specific_with_hits <- setdiff(all_ids_with_hits,not_octopod_specific)
  
  cephalopod_specific <- setdiff(all_ids,not_ceph_specific)

  saveRDS(cephalopod_specific,"cache/cephalopod_specific.rds")
  cephalopod_specific
}



