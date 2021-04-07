# Script to create locus zoom plots on broad server

# Dependencies
library(data.table)

phenos <- c("SVT")
for (i in phenos){
  data <- fread(paste0('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/meta_outputs/',i,'.tsv'))
  lz <- data[,c('MarkerName','P-value')]
  write.table(lz,file=paste0('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/meta_outputs/',i,'_lz.txt'),
              row.names=F,quote=F,sep='\t')
}
