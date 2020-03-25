# Script to data incident disease counts for dataset using BigQuery data (to be run on CGP)

# Dependencies
library(data.table)

### Step 1: Load and process phenotypes from big query
## Load data
data <- fread(file='/mnt/ml4cvd/projects/skhurshid/svt_brady_gwas/phenotypes_long.csv',data.table=F)

## Reduce columns to sample_id + phenotypes only
data[,-c('enroll_age')]

## Create function to convert to standard wide format
widen <- function(data){
  out <- list()
  n <- 0
  for (i in unique(data[,'sample_id'])){
    n <- n + 1
    subset <- data[data[,'sample_id']==i,]
    out[[n]] <- data.frame(matrix(ncol=length(unique(subset[,'disease']))*3+1,nrow=0))
    out[[n]][1,1] <- i; names(out[[n]])[1] <- 'sample_id'
    m <- 2
    for (j in unique(subset[,'disease'])){
      out[[n]][1,m:(m+2)] <- c(subset[subset[,'disease']==j,'prevalent_disease'],
                               subset[subset[,'disease']==j,'incident_disease'],
                               subset[subset[,'disease']==j,'survyears_disease'])
      names(out[[n]])[m:(m+2)] <- c(paste0('prev_',j),paste0('incd_',j),paste0('survyears_',j))
      m <- m+3
    }
  }
  return(do.call(rbind,out))
}

save(data,file='/mnt/ml4cvd/projects/skhurshid/svt_brady_gwas/phenotypes_wide_032520.RData')