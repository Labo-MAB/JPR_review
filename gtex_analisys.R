#install.packages("stringr")
library(stringr)
dat.gene_tpm <- read.delim(file="./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", skip=2)
dat.annotation<- read.delim(file="./GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

newcolnames<-c()
for (n in colnames(dat.gene_tpm)[3:length(colnames(dat.gene_tpm))] ){
  n1<-gsub("\\.","-",n)
  newcolnames<-c(newcolnames,dat.annotation[which(dat.annotation$SAMPID==n1),"SMTSD"])
}
colnames(dat.gene_tpm)[3:length(colnames(dat.gene_tpm))]<-newcolnames

#POU5F1
pou5f1.fam<-dat.gene_tpm[grepl("POU5F1",dat.gene_tpm$Description),]
pou5f1.fam_t<-data.frame(t(pou5f1.fam[,3:length(colnames(pou5f1.fam))]))
colnames(pou5f1.fam_t)<-pou5f1.fam$Description
pou5f1.fam_t$tissue<-gsub("\\..*","",rownames(pou5f1.fam_t))
#get median of each tissue per gene to select the most expressed
gene_median_tissue<-list()
for (gene in pou5f1.fam$Description){
  median_list<-list()
  for (tis in unique(pou5f1.fam_t$tissue)){
    median_list[tis]<-median(pou5f1.fam_t[which(pou5f1.fam_t$tissue==tis),gene])
  }
  gene_median_tissue[gene]<-names(median_list[which.max(median_list)])
}
print(gene_median_tissue)
#$POU5F1P4
#[1] "Brain - Cerebellum"
#
#$POU5F1P7
#[1] "Adipose - Subcutaneous"
#
#$POU5F1P6
#[1] "Testis"
#
#$POU5F1
#[1] "Kidney - Medulla"
#
#$POU5F1P2
#[1] "Adipose - Subcutaneous"
#
#$POU5F1B
#[1] "Cervix - Endocervix"
#
#$POU5F1P5
#[1] "Skin - Not Sun Exposed (Suprapubic)"
#
#$POU5F1P3
#[1] "Nerve - Tibial"

#do boxplots for each gene-tissue
pou5f1.fam_t.tmp<-list()
for (i in 1:length(gene_median_tissue)){
  gene<-names(gene_median_tissue)[i]
  tis<-gene_median_tissue[gene]
  pou5f1.fam_t.tmp[[gene]]<-pou5f1.fam_t[which(pou5f1.fam_t$tissue==tis),gene]
}
boxplot(pou5f1.fam_t.tmp)
boxplot(pou5f1.fam_t.tmp,outline=FALSE)
#save
saveRDS(pou5f1.fam_t.tmp, file = "pou5f1_fam_t_tmp.rds")

#NANOG
nanog.fam<-dat.gene_tpm[grepl("NANOG",dat.gene_tpm$Description),]
nanog.fam_t<-data.frame(t(nanog.fam[,3:length(colnames(nanog.fam))]))
colnames(nanog.fam_t)<-nanog.fam$Description
nanog.fam_t$tissue<-gsub("\\..*","",rownames(nanog.fam_t))
#get median of each tissue per gene to select the most expressed
gene_median_tissue<-list()
for (gene in nanog.fam$Description){
  median_list<-list()
  for (tis in unique(nanog.fam_t$tissue)){
    median_list[tis]<-median(nanog.fam_t[which(nanog.fam_t$tissue==tis),gene])
  }
  gene_median_tissue[gene]<-names(median_list[which.max(median_list)])
}
print(gene_median_tissue)
#$NANOGNBP1
#[1] "Adipose - Subcutaneous"
#
#$NANOGP2
#[1] "Skin - Not Sun Exposed (Suprapubic)"
#
#$NANOGP3
#[1] "Adipose - Subcutaneous"
#
#$NANOGP11
#[1] "Brain - Cerebellum"
#
#$NANOGP4
#[1] "Brain - Cerebellum"
#
#$NANOGP5
#[1] "Brain - Cerebellum"
#
#$NANOGP6
#[1] "Liver"
#
#$NANOGNB
#[1] "Testis"
#
#$NANOG
#[1] "Testis"
#
#$NANOGP1
#[1] "Testis"
#
#$NANOGNBP2
#[1] "Testis"
#
#$NANOGP7
#[1] "Artery - Aorta"
#
#$NANOGP8
#[1] "Testis"
#
#$NANOGP10
#[1] "Adipose - Subcutaneous"
#
#$NANOGP9
#[1] "Skin - Sun Exposed (Lower leg)"
#
#$NANOGNBP3
#[1] "Brain - Nucleus accumbens (basal ganglia)"

#do boxplots for each gene-tissue
nanog.fam_t.tmp<-list()
for (i in 1:length(gene_median_tissue)){
  gene<-names(gene_median_tissue)[i]
  tis<-gene_median_tissue[gene]
  nanog.fam_t.tmp[[gene]]<-nanog.fam_t[which(nanog.fam_t$tissue==tis),gene]
}
boxplot(nanog.fam_t.tmp,las=2)
boxplot(nanog.fam_t.tmp,outline=FALSE,las=2)
#save
saveRDS(nanog.fam_t.tmp, file = "nanog_fam_t_tmp.rds")
