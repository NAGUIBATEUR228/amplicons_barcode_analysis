library(tidyverse)
library("ggrepel")
library('patchwork')

folder_with_fastq<-'fin'
to<-read_tsv("../tables/tight_orfs.tsv")%>%
  transmute(sys_name=`Gene > Systematic Name`,std_name= `Gene > Standard Name`)

to_join<-function(x){
  y<-x%>%
    transmute(sys_name=Confirmed_deletion,count)%>%
    left_join(to,by=c("sys_name"))%>%
    group_by(sys_name,std_name)%>%
    summarise(count=sum(count))%>%
    ungroup()
  y$count[is.na(y$count)]<-0
  return(y)
}

"%+%" <- function(...){
  paste0(...)
}

workingDIR <- "C:/Users/zokmi/Desktop/study/coursework/mtx"
setwd(workingDIR)

seqdirs<-list.dirs(path='../'%+%folder_with_fastq,recursive = FALSE)%+%'/artem'
seqfiles<-unlist(map(seqdirs, ~paste(.x,dir(.x)%>%str_subset("output_count.csv$"),sep='/')))

tablenames<-str_remove_all(seqfiles,'_L00_R')%>%
  str_match_all("artem/\\s*(.*?)\\s*_output_count.csv")%>%
  map_chr(function(x){x[,2]%>%tolower()})


read_files<-function(type,abrev,join=FALSE){
  seqdirs<-get('seqdirs',envir = .GlobalEnv)
  tablenames<-get('tablenames',envir = .GlobalEnv)
  
  seqfiles<-unlist(map(seqdirs, ~paste(.x,dir(.x)%>%str_subset(type%+%".csv$"),sep='/')))
  
  chiffre<-tibble(tablename=tablenames,seqfile=seqfiles)
  
  if(join){
    map2(tablenames%+%abrev,seqfiles,
         function(x,y){
           assign(x,
                  read_csv(y,show_col_types = FALSE)%>%
                    to_join%>%
                    mutate(exp=x),
                  envir=.GlobalEnv)
         }
    )
  }
  else{map2(tablenames%+%abrev,seqfiles,
            function(x,y){
              assign(x,
                     read_csv(y,show_col_types = FALSE),
                     envir=.GlobalEnv)
            }
  )
  }
  return(chiffre)
}

chif<-read_files("output_count",'',TRUE)
read_files("output_dm_count",'dm',TRUE)
read_files("mixaled_dm",'mdm',TRUE)
read_files("blasted_dm",'bdm',TRUE)
read_files("mixaled",'m',TRUE)
read_files("blasted",'b',TRUE)
read_files('raw_reads_count','r')

tbs<-c(tablenames,tablenames%+%'dm',tablenames%+%'mdm',tablenames%+%'m')

counts<-purrr::reduce(map(tbs,get),full_join)%>%
  pivot_wider(names_from = exp,values_from = count)%>%
  mutate(across(where(is.numeric),
                ~ifelse(is.na(.x),0,.x)))%>%
  mutate(name=ifelse(is.na(std_name),sys_name,std_name))
colnames(counts)
  
jcounts<-counts%>%select(name, dck11:yg32m)

jcounts <- jcounts %>% pivot_longer(- name, 
                                   values_to = "count", 
                                   names_to = "Conditions") %>%
  mutate(Experiments = ifelse(str_detect(Conditions,'dck') | str_detect(Conditions,'sov'),str_sub(Conditions, 1, 4),str_sub(Conditions, 1, 3))) %>%
  group_by(name, Experiments) %>%
  summarise(count = sum(count)) %>% 
  ungroup()%>%
  pivot_wider(names_from = "Experiments", values_from= "count")
jcounts <- jcounts%>%dplyr::mutate(across(where(~is.numeric(.)), floor))
  
write_csv(jcounts,seqdirs[1]%+%'/../../jcounts_fin.csv')
