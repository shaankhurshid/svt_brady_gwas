# Script to load phenotypes from broad server

# Dependencies
library(data.table)

# Load phenotype list
files_list <- list.files('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes')

# Load phenotypes
out <- list()
for (i in 1:length(files_list)){
  out[[i]] <- paste0('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/',files_list[i])
}
paths <- do.call(rbind,out)

data_long <- do.call(rbind,lapply(paths,fread,data.table=F))

# Widen
widen <- function(data){
  out <- list()
  n <- 0
  for (i in unique(data[,'sample_id'])){
    if((n-1) %% 1000 == 0){print(paste0('Just finished row ',n-1,' out of ',
                                        length(unique(data[,'sample_id'])),'!'))}
    n <- n + 1
    subset <- data[data[,'sample_id']==i,]
    out[[n]] <- data.frame(matrix(ncol=length(unique(subset[,'disease']))*3+1,nrow=0))
    out[[n]][1,1] <- i; names(out[[n]])[1] <- 'sample_id'
    m <- 2
    for (j in unique(subset[,'disease'])){
      out[[n]][1,m:(m+2)] <- c(subset[subset[,'disease']==j,'prevalent_disease'],
                               subset[subset[,'disease']==j,'incident_disease'],
                               subset[subset[,'disease']==j,'censor_date'])
      names(out[[n]])[m:(m+2)] <- c(paste0('prev_',j),paste0('incd_',j),paste0('censor_',j))
      m <- m+3
    }
  }
  return(do.call(rbind,out))
}

data_wide <- widen(data_long)

save(output,file='/mnt/ml4cvd/projects/skhurshid/svt_brady_gwas/phenotypes_wide_032520.RData')