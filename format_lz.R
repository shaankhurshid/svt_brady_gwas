# Script to format meta-analysis results for locus zoom

# Depends
library(data.table)

# List paths
root <- '/medpop/afib/projects/SVT_BRADY/Preliminary/Meta_analysis/'
phenos <- c('SVT','WPW_EUR','AVNRT_EUR','SND_INC','SND_REST','DISTAL_INC','DISTAL_REST','PACER_EUR')
paths <- paste0(root,phenos,'/METAANALYSIS1.TBL')

# Loop and reformat
for (i in 1:1){
  file <- fread(paths[i])
  splits <- do.call(rbind,strsplit(file$MarkerName,':'))
  file[,':='(Chr = as.numeric(splits[,1]),
             Genpos = as.numeric(splits[,2]))]
  #write.table(file,'/medpop/afib/skhurshid/svt_brady_gwas/meta_outputs/',phenos[i],'.tsv',sep='\t',row.names=F)
}