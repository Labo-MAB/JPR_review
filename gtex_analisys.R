library(stringr)
setwd("C:/Users/Lera/Documents/JPR review/POU5F1B_vis")
dat.gct <- read.delim(file="./gtex_median_tpm/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", skip=2)
#pou5f1
pou5f1.fam<-dat.gct[grepl("POU5F1",dat.gct$Description),]
pou5f1.fam_t<-data.frame(t(pou5f1.fam[,3:56]))
colnames(pou5f1.fam_t)<-pou5f1.fam$Description
# p2, b, 1, 3,4,7,6,5
#boxplot(pou5f1.fam_t$POU5F1P2, pou5f1.fam_t$POU5F1B,pou5f1.fam_t$POU5F1, pou5f1.fam_t$POU5F1P3,
#        pou5f1.fam_t$POU5F1P4,pou5f1.fam_t$POU5F1P7,pou5f1.fam_t$POU5F1P6, pou5f1.fam_t$POU5F1P5,
#        ) 
# pg, b, p3, p4, p6, p7, p2, p5
boxplot(pou5f1.fam_t$POU5F1,pou5f1.fam_t$POU5F1B,pou5f1.fam_t$POU5F1P3,pou5f1.fam_t$POU5F1P4,
        pou5f1.fam_t$POU5F1P6,pou5f1.fam_t$POU5F1P7,pou5f1.fam_t$POU5F1P2,pou5f1.fam_t$POU5F1P5, 
        col="white",outline=FALSE, las=2,names=c("POU5F1","POU5F1B","POU5F1P3","POU5F1P4","POU5F1P6","POU5F1P7","POU5F1P2","POU5F1P5"))

#nanog
nanog.fam<-dat.gct[grepl("NANOG",dat.gct$Description),]
nanog.fam_t<-data.frame(t(nanog.fam[,3:56]))
colnames(nanog.fam_t)<-nanog.fam$Description
nanog.fam_t_filt<-nanog.fam_t[,-which(str_detect(colnames(nanog.fam_t),"NANOGNB"))]
#boxplot(nanog.fam_t_filt$NANOGP4, nanog.fam_t_filt$NANOGP7,nanog.fam_t_filt$NANOGP2,nanog.fam_t_filt$NANOGP3,
#        nanog.fam_t_filt$NANOGP6,nanog.fam_t_filt$NANOG, nanog.fam_t_filt$NANOGP8,
#        nanog.fam_t_filt$NANOGP5,nanog.fam_t_filt$NANOGP10, nanog.fam_t_filt$NANOGP9,
#        nanog.fam_t_filt$NANOGP1,nanog.fam_t_filt$NANOGP11)
# pg, p8, p6, p5, p3, p10, p2, p9, p7, p1, p4, p11
boxplot(nanog.fam_t_filt$NANOG,nanog.fam_t_filt$NANOGP8,nanog.fam_t_filt$NANOGP6,nanog.fam_t_filt$NANOGP5,
        nanog.fam_t_filt$NANOGP3,nanog.fam_t_filt$NANOGP10,nanog.fam_t_filt$NANOGP2,nanog.fam_t_filt$NANOGP9,
        nanog.fam_t_filt$NANOGP7,nanog.fam_t_filt$NANOGP1,nanog.fam_t_filt$NANOGP4,nanog.fam_t_filt$NANOGP11,
        col="white",outline=FALSE, las=2,names=c("NANOG","NANOGP8","NANOGP6","NANOGP5","NANOGP3","NANOGP10","NANOGP2","NANOGP9","NANOGP7","NANOGP1","NANOGP4","NANOGP11"))

