
library(tidyverse)
library(Rbowtie2)
library(Rsamtools)
library(RColorBrewer)
library(ggrepel)
library(Biostrings)
library(ggpubr)
library(BioCircos)
library(ShortRead)
library(R.utils)


  # reverse_complement_string
  rev_compl_string <- function(x) {
    as.character(reverseComplement(DNAStringSet(x)))
  }
   
  

  #Set working directory
  workingDIR <- "/home/dima/Desktop/BiocR/Genomed_1/372/372_Knorre_Nano300_Test-365509233/"
  setwd(workingDIR)

  #read file with Reference reads file generated from Puddu et al. 2019
  reference <- read_csv2("/home/dima/Desktop/BiocR/Genomed_1/Reference_Uptags.csv")
  
    
  u1 <- DNAString("GATGTCCACGAGGTCTCT")
  ?"DNAString"
  u1_rc <- reverseComplement(u1)
  
  u2 <- DNAString("CGTACGCTGCAGGTCGAC")
  ?"DNAString"
  u2_rc <- reverseComplement(u2)
  
   
  ##unpack files in cycle
  ####### skip this when done
  DIRlist <- paste(workingDIR, dir(), sep  = "")
  
  for (DIRs in DIRlist)
  {
   setwd(DIRs)
    targz_files <- dir() %>% str_subset(".gz") 
    for (tar_file in targz_files)
    {
      print(tar_file)
      gunzip(tar_file, skip = TRUE)
    }
  }
  
  ################
  
  
  # 
  # for (DIRs in DIRlist) {
  #   fastq_files <- dir() %>% str_subset(".gz", negate = TRUE) %>%
  #     str_subset(".fastq")
  #   for (fastq_file in fastq_files)
  #   {
  # strm <- FastqStreamer(fastq_file)
  strm <- FastqStreamer("Sov3_S6_L001_R1_001.fastq")
  
  out_table <- tibble(seq = character(), count = integer())
  
  repeat {
    fq <- yield(strm)
    if (length(fq) == 0) {break}
   t1 <- tables(sread(fq), n = 100000)  ## consider increase/decrease n here! 
   ## n limit the amplicons diversity
   print(c("min count = ", min(t1$top)))
   print(c("max count = ", max(t1$top)))
   
   attach_this_table <- tibble(seq_new = names(t1$top), count_new = t1$top) %>% 
     mutate(count_new = replace_na(count_new, 0))
   
   out_table <- out_table %>% 
     full_join(attach_this_table, by = c("seq" = "seq_new"), suffix = c("","new")) %>%
     mutate(count = replace_na(count, 0)) %>%
     mutate(count = count + count_new) %>% 
     select("seq","count")
     }
  # }
  # 
  
  raw_reads_count <- out_table
  read_count <- length(fq)
  max_coverage <- max(t1$top)

  write.csv(raw_reads_count, paste(fastq_file, 
                                   "raw_reads_count.csv", sep = "_"))

  
  filtered_reads <- raw_reads_count %>% 
    filter(str_detect(seq, as.character(u1)) | str_detect(seq, as.character(u1_rc))) %>%
  filter(str_detect(seq, as.character(u2)) | str_detect(seq, as.character(u2_rc))) 
  
  barcoded_reads <- filtered_reads %>% 
    mutate(barcode = str_match(seq, "GATGTCCACGAGGTCTCT\\s*(.*?)\\s*CGTACGCTGCAGGTCGAC")[,2]) %>%
    mutate(barcode = ifelse(is.na(barcode),
                            str_match(seq,"GTCGACCTGCAGCGTACG\\s*(.*?)\\s*AGAGACCTCGTGGACATC")[,2], 
                            barcode)) %>% 
    filter(!is.na(barcode)) %>% 
    mutate(barcode = ifelse(str_detect(seq, "GATGTCCACGAGGTCTCT"),
                            barcode, 
                            rev_compl_string(barcode)))
    
  
  final_result <- barcoded_reads %>% 
    group_by(barcode) %>%
    summarise(n = n(), count = count) %>%
    arrange(desc(count))
  
  
  reference_table <- unique(reference) %>%
      summarise(n = n(), ORF = Confirmed_deletion) %>%
      arrange(desc(n))
  
  unique(reference$UPTAG_seqs)
    
  output <- unique(reference) %>% 
    left_join(final_result, by = c("UPTAG_seqs" = "barcode")) %>% 
    filter(!is.na(count)) %>%
    arrange(desc(count)) 
  
  out_summary <- output %>% filter(!is.na(count)) 
  
  output %>% ggplot(aes(x = count)) +
    geom_histogram(binwidth = 10) +
    xlim(0,1000)
  
  write_csv2(output,"output_count.csv")
  
  total_count <- sum(output$count, na.rm = TRUE)
  raw_count <- sum(raw_reads_count$count.x, na.rm = TRUE) +
  sum(sum(raw_reads_count$count.y, na.rm = TRUE))
  hundred_and_more_ORFs <- sum(output$count > 10)
  hundred_and_more_ORFs_raw <- sum((raw_reads_count$count.x + raw_reads_count$count.y) > 10)
   
  as.character(u1_rc)
  as.character(u2_rc)
  

  raw_reads_count$seq %>% 
    # str_detect(as.character(u1_rc)) %>%
    str_detect(as.character(u1)) %>%
    sum()
  
################
## Reads Stats
################
summary <- list(reads = countFastq(fastq_files[1])[[1]]*2,
       repeated_barcodes = 2*sum(out1$count, na.rm = TRUE),
       primers_detectes = 2*sum(filtered_reads$count.x, na.rm = TRUE),
       unique_barcodes = length(final_result$barcode),
       ORFs = length(out_summary$ORF),
       final_mapped = sum(out_summary$count),
       median = median(out_summary$count))
         
summary
  
#############


dna2char(u1_rc)
  
  findU1 <- function (x) {
     vcountPattern("GATGTCCACGAGGTCTCT", sread(x)) 
  }
  
  u1filter <- srFilter(findU1)
  myfive_f1 <- fq[u1filter(fq)]
  sread(myfive_f1)
  
  findU2 <- function (x) {
     vcountPattern("CGTACGCTGCAGGTCGAC", sread(x)) 
     
  }
  
  u2filter <- srFilter(findU2)
  myfive_f2 <- myfive_f1[u2filter(myfive_f1)]
  sread(myfive_f2)
  

  
  
  
  t1 <- tables(sread(fq), n = 20000)
  min(t1$top)
  
  length(t1$top)
  
  ss1 <- tibble(t1$top, names(t1$top))
  
  
  
  sample <- FastqSampler("/home/dima/Desktop/NGS/Amplicons/D1n_S15/D1n_S15_R1_001.fastq")
  fq <- yield(sample)
  myfive <- fq[1:100]
  sread(myfive)
  quality(myfive)
  