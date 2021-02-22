##################################
#
##combine results within each chromosome
#
##################################
library(data.table)

merge_results <- function(paths,phenos){
  for (iter in 1:length(paths)){
    
    path <- paths[iter]; pheno <- phenos[iter]
    setwd(path)
    
    frq_file <- list.files(pattern="*.freq")
    frq_file1<-frq_file[grep(pheno,frq_file)]
    for (chr in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")){
      frq_file2<-frq_file1[grep(paste0(chr,"\\."),frq_file1)]
      info0<-fread(paste0("/broad/ukbb/imputed_v3/ukb_mfi_",chr,"_v3.txt"),header=F,data.table=F)
      number <- sapply(strsplit(as.character(frq_file2), split= "\\."), '[', 3)
      file0<-NULL
      for (num in number){
        log<-fread(paste0(pheno,".",chr,".",num,".assoc.",pheno,".glm.logistic"),header=T,data.table=F)
        frq<-fread(paste0(pheno,".",chr,".",num,".assoc.afreq"),header=T,data.table=F)
        
        info<-read.table(paste0("/medpop/afib/lcweng/UKBB_all/SVT_Brady/UKbiobank_analysis/chunk/variantList_",chr,"_info04.",num,".temp.txt"),header=F,as.is=T)
        info<-info0[info0$V2%in%info$V1,c(2,4,5,8)]
        names(info)<-c("ID","V4","V5","info")
        
        log2<-merge(log,frq,by=c("#CHROM","ID","REF","ALT"),suffixes=c(".log",".frq"))
        print(dim(log))
        print(dim(log2))
        print(dim(subset(log2,OBS_CT.log== OBS_CT.frq/2)))
        log2<-merge(log2,info,by.x=c("ID","REF","ALT"),by.y=c("ID","V5","V4"))
        print(dim(log2))
        #only select variants with good quality (info>=0.4)
        log2<-subset(log2,info>=0.4)
        
        ##reassign correct alt_allele
        log2$NON_EFFECT_ALLELE<-ifelse(log2$A1==log2$ALT,log2$REF,log2$ALT)
        log2$EFFECT_ALLELE<-ifelse(log2$A1==log2$ALT,log2$ALT,log2$REF)
        log2$EAF<-ifelse(log2$A1==log2$ALT,log2$ALT_FREQ,1-log2$ALT_FREQ)
        log2$BETA<-log(as.numeric(log2$OR))
        log2<-log2[,c(1,4,5,8,15,16,17,18,19,9,10,12)]
        names(log2)<-c("ID","CHR","POS","N","INFO","NON_EFFECT_ALLELE","EFFECT_ALLELE","EAF","BETA","OR","SE","P")
        file0<-rbind(file0,log2)
        print(chr)
        print(num)
      }
      write.table(file0,file=paste0(pheno,".",chr,".UKBB.EUR.noRel.logistic"),quote=F,row.names=F,sep="\t")
    }
  }
}

##################################
#
##combine results across each chromosome
#
##################################
library(R.utils)

merge_across <- function(paths,phenos){
  for (iter in 1:length(paths)){
    
    path <- paths[iter]; pheno <- phenos[iter]
    setwd(path)
    file0<-NULL
    
  for (chr in c(1:22)){
    a<-fread(paste0(pheno,".chr",chr,".UKBB.EUR.noRel.logistic"),header=T,data.table=F)
    a$ID3<-paste0(a$ID,a$NON_EFFECT_ALLELE,a$EFFECT_ALLELE,a$INFO)
    b<-read.table(paste0("/medpop/afib/lcweng/UKBB_all/SVT_Brady/UKbiobank_analysis/variantList_chr",chr,"_info04.txt.gz"),header=F)
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

#merge_results(paths=paths,phenos=phenos)
merge_across(paths=paths,phenos=phenos)


