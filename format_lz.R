# Script to format meta-analysis results for locus zoom

# Depends
library(data.table)

# List paths
root <- '/medpop/afib/projects/SVT_BRADY/Preliminary/Meta_analysis/'
phenos <- c('SVT','WPW_EUR','AVNRT_EUR','SND_INC','SND_REST','DISTAL_INC','DISTAL_REST','PACER_EUR')
paths <- paste0(root,phenos,'/METAANALYSIS1.TBL')

# Loop and reformat
for (i in 1:length(paths)){
  tbl <- fread(paths[i])
  splits <- do.call(rbind,strsplit(tbl$MarkerName,':'))
  tbl[,':='(Chr = as.numeric(splits[,1]),
             Genpos = as.numeric(splits[,2]))]
  setkey(tbl,Chr,Genpos)
  write.table(tbl,file=paste0('/medpop/afib/skhurshid/svt_brady_gwas/meta_outputs/',phenos[i],'.tsv'),quote=F,sep='\t',row.names=F)
}