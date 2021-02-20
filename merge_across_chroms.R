##################################
#
##combine results across each chromosome
#
##################################
library(R.utils)
library(data.table)

merge_across <- function(paths,phenos){
  for (iter in 1:length(paths)){
    
    path <- paths[iter]; pheno <- phenos[iter]
    setwd(path)
    file0<-NULL
    
  for (chr in c(1:22)){
    a<-fread(paste0(pheno,".chr",chr,".UKBB.EUR.noRel.logistic"),header=T,data.table=F)
    a$ID3<-paste0(a$ID,a$NON_EFFECT_ALLELE,a$EFFECT_ALLELE,a$INFO)
    b<-fread(paste0("/medpop/afib/lcweng/UKBB_all/SVT_Brady/UKbiobank_analysis/variantList_chr",chr,"_info04.txt.gz"),header=F,data.table=F)
    print(length(a$ID))
    print(length(unique(a$ID3)))
    print(length(b$V1))
    a<-a[!duplicated(a$ID3),-13]
    print(length(a$ID))
    file0<-rbind(file0,a)
    print(chr)
  }
  setwd('..')
  write.table(file0,file=paste0(pheno,".UKBB.EUR.noRel.final.logistic"),quote=F,row.names=F,sep="\t")
  }
}

paths <- c('/medpop/afib/skhurshid/svt_brady_gwas/svt_brady_gwas/analysis2/temp/',
           '/medpop/afib/skhurshid/svt_brady_gwas/svt_brady_gwas/analysis_snd_hard/temp/')
phenos <- c("Bradyarrhythmia_sinus_node_dysfunction",
            "Bradyarrhythmia_sinus_node_dysfunction_HARD_V2")

merge_across(paths=paths,phenos=phenos)


