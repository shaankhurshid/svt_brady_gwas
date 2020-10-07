library(data.table)
library(R.utils)

setwd("/medpop/afib/skhurshid/svt_brady_gwas/svt_brady_gwas/analysis2/") 
for (pheno in "Bradyarrhythmia_sinus_node_dysfunction"){
  file0<-NULL
  for (chr in c(1:22)){
    a<-fread(paste0("/medpop/afib/skhurshid/svt_brady_gwas/svt_brady_gwas/analysis2/temp/",pheno,".chr",chr,".UKBB.EUR.noRel.logistic"),header=T,data.table=F)
    a$ID3<-paste0(a$ID,a$NON_EFFECT_ALLELE,a$EFFECT_ALLELE,a$INFO)
    b<-read.table(paste0("/medpop/afib/lcweng/UKBB_all/SVT_Brady/UKbiobank_analysis/variantList_chr",chr,"_info04.txt.gz"),header=F,data.table=F)
    print(length(a$ID))
    print(length(unique(a$ID3)))
    print(length(b$V1))
    a<-a[!duplicated(a$ID3),-13]
    print(length(a$ID))
    file0<-rbind(file0,a)
    print(chr)
  }
  write.table(file0,file=paste0("/medpop/afib/skhurshid/svt_brady_gwas/svt_brady_gwas/analysis2/",pheno,".UKBB.EUR.noRel.final.logistic"),quote=F,row.names=F,sep="\t")
}


