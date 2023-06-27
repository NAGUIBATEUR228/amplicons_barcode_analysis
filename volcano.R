library(tidyverse)
library("ggrepel")
library(RColorBrewer)
library(ggpubr)
library("DESeq2")
library('EnhancedVolcano')

### FUNCTIONS


#Set working directory
setwd("C:/Users/zokmi/Desktop/study/coursework/fin")
jcounts <- read.csv("jcounts_fin.csv")

jc <- as.data.frame(jcounts) 
jc<-jc%>%select(!starts_with(c('dck','sov','yd3')))

difexp<-function(jc,t,u){
rownames(jc) <- jc$name
cd <- jc %>% select(!name)%>%select(starts_with(t)|starts_with(u))
#cd<-jc%>%filter(rowSums(as.matrix(jc%>%select(starts_with('rd')|starts_with('pg'))))>=100)

colnames <- data.frame(experiment = colnames(cd)[1:ncol(cd)])
rownames(colnames) <- colnames[,1]
colnames <- colnames %>% mutate(condition = str_sub(experiment, 1,-2)) %>%
  mutate(rep = str_sub(experiment, -1, -1))  

### DDS

dds <- DESeqDataSetFromMatrix(countData = cd,
                              colData = colnames,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]

dds <- DESeq(dds)

dds$condition
res <- results(dds, contrast=c("condition",t,u))
return(res)
}
### res <- results(dds, contrast=c("condition","treated","untreated"))
res<-difexp(jc,'rd','yd')
plotMA(res, ylim=c(-5,5))
summary(res)
res


sg<-tibble(gene=character())
sg_tc<-character(0)
cut_tc<-character(0)

for (i in c('yg_yd','rd_yd','rg_rd','pg_pd','pg_rg')){
condition<-i
t=str_split_1(condition,'_')[1]
u=str_split_1(condition,'_')[2]
res<-difexp(jc,t,u)

res %>% as.data.frame() %>% 
  mutate(gene = rownames(.)) %>% 
  select(gene, log2FoldChange, pvalue) %>% 
  drop_na(pvalue) %>%as_tibble->res
res%>% mutate(pa=p.adjust(pvalue))%>%
  filter(pa < 10e-5)->cut

#write_tsv(cut%>%filter(log2FoldChange>0)%>%select(gene),'../GOen/'%+%condition%+%'_p.txt',col_names=F)
#write_tsv(cut%>%filter(log2FoldChange<0)%>%select(gene),'../GOen/'%+%condition%+%'_m.txt',col_names=F)
write_tsv(res%>%select(gene),'../GOen/'%+%condition%+%'_ref.txt',col_names=F)

if(nrow(cut)<20){
  tc<-cut$gene
}
if(nrow(cut)>=20){
  tc_p<-cut%>%filter(log2FoldChange>0)%>%arrange(desc(log2FoldChange))
  tc_p<-(tc_p%>%pull(gene))[1:min(nrow(tc_p),10)]
  tc_m<-cut%>%filter(log2FoldChange<0)%>%arrange(log2FoldChange)
  tc_m<-(tc_m%>%pull(gene))[1:min(nrow(tc_m),10)]
  tc<-c(tc_p,tc_m)%>%na.omit()%>%unique()
}
sg_add<-res%>%transmute(gene,!!condition :=log2FoldChange,!!(condition%+%'_pv') := -log(pvalue,base=10))
cut_tc<-c(cut_tc,cut$gene)%>%unique
sg_tc<-c(sg_tc,tc)%>%unique

sg<-full_join(sg,sg_add,by='gene')
}
hm<-sg%>%filter(gene%in%cut_tc)
sg<-sg%>%filter(gene%in%sg_tc)
#write_tsv(sg,'logFC_all.txt')
#write_tsv(sg%>%filter(gene%in%cut_tc),'logFC_cut.txt')
library("pheatmap")
ann <- sg %>% select(ends_with('_pv')) %>% as.data.frame()
rownames(ann) <- sg$gene
mat<-sg%>%select(!ends_with('_pv'))%>%select(!gene,)%>%as.matrix
rownames(mat)<-sg$gene
plot<-pheatmap(mat, annotation_row = ann, annotation_names_row=F, treeheight_col = 0, treeheight_row = 0,cluster_rows = T, cluster_cols = T,display_numbers = T,fontsize_number = 9, fontsize=10,number_color = 'black',number_format = "%.2f",border_color = 'black', angle_col='0',main='Топ LogFC гены')
ggsave("../top_LogFC_genes_fin.png",plot,scale=2.3)


ann <- hm %>% select(ends_with('_pv')) %>% as.data.frame()
rownames(ann) <- hm$gene
mat<-hm%>%select(!ends_with('_pv'))%>%select(!gene)%>%as.matrix
rownames(mat)<-hm$gene
plot<-pheatmap(mat, annotation_row = ann, annotation_names_row=F, treeheight_col = 1, treeheight_row = 50, cluster_rows = T, cluster_cols = T,display_numbers = F, fontsize=10, show_rownames=F, number_color = 'black', number_format = "none", border_color = 'black', cellwidth=50, angle_col='0',main='Топ LogFC гены')
ggsave("../top_LogFC_genes_all_fin.png",plot,scale=2)



sgd<-read_tsv('sgd_term.txt')
ab<-read_tsv('ab.txt',col_names = c('name'))
sgd%>%mutate(name=ifelse(is.na(std_name),sys_name,std_name))%>%select(name,DBID)%>%right_join(ab)%>%select(DBID)%>%na.omit%>%write_tsv('DBID.txt',col_names=F)
ab%>%filter(!((name %in% sgd$std_name)|(name %in% sgd$sys_name)))
sgd%>%mutate(name=ifelse(is.na(std_name),sys_name,std_name))%>%select(name,DBID)%>%right_join(read_tsv('goref.txt',col_names = 'name'))%>%select(DBID)%>%na.omit%>%write_tsv('DBID_ref.txt',col_names=F)

plo<-EnhancedVolcano(res,
                     #title=bquote(~"G\u00D7"*italic(rho)^'-'~vs~"D\u00D7"*italic(rho)^'-'),#bquote(~"D\u00D7"*italic(rho)^0~vs~D)
                     lab = rownames(res),
                     subtitle = element_blank(),
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote(~-Log[10] ~ italic(P)),
                     pCutoff = 10e-16,
                     FCcutoff = 2.0,
                     pointSize = 4.0,
                     labSize = 6.0,
                     legendPosition = 'none',
                     legendLabSize = 12,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     #xlim=c(-5,10),
                     widthConnectors = 0.75)
ggsave('volcano_pgpd.png',plo,width = 1280, height = 720, units = 'px', scale=3)


#dev.new(width = 1280, height = 720, units='px', scale=5)

res<-difexp(jc,'rd','yd')
