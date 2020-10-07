# Script to create locus zoom plots on broad server

# Dependencies
library(data.table)

phenos <- c("DISTAL_INC","SND_REST")
for (i in phenos){
  data <- fread(paste0('/Volumes/medpop_afib/projects/SVT_BRADY/Preliminary/Meta_analysis/',i,'/METAANALYSIS1.TBL'))
  lz <- data[,c('MarkerName','P-value')]
  write.table(lz,file=paste0('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/svt_brady_gwas/',i,'_lz.txt'),
              row.names=F,quote=F,sep='\t')
}

data <- fread('/Volumes/medpop_afib/projects/SVT_BRADY/collaborator_GWAS_results/QC/AVNRT_ALL/CLEANED.BroadAF_AVNRT_ALL')
lz <- data[,c('SNP','PVAL')]
setnames(lz,'SNP','MarkerName'); setnames(lz,'PVAL','P-value')
write.table(lz,file=paste0('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/svt_brady_gwas/AVNRT_ALL_CLEANED_lz.txt'),
            row.names=F,quote=F,sep='\t')