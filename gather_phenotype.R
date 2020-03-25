# Script to data incident disease counts for dataset using BigQuery data (to be run on CGP)

# Dependencies
library(data.table)

### Step 1: Load and process phenotypes from big query
## Load data
data <- fread(file='/mnt/ml4cvd/projects/skhurshid/svt_brady_gwas/phenotypes_long.csv',data.table=F)

## Specify how many columns exist prior to "disease" (so they are preserved)
lead_columns <- 2

## Create function to convert to standard wide format
widen <- function(data){
  for (i in unique(data[,'disease'])){
    data[,paste0('incident',i)] <- rep(NA,nrow(data))
    data[,paste0('prevalent',i)] <- rep(NA,nrow(data))
    data[,paste0('survyears',i)] <- rep(NA,nrow(data))}
  for (i in 1:nrow(data)){
    data[,paste0('incident',data[,'disease'][i])][i] <- data[,'incident_disease'][i]
    data[,paste0('prevalent',data[,'disease'][i])][i] <- data[,'prevalent_disease'][i]
    data[,paste0('survyears',data[,'disease'][i])][i] <- data[,'survyears_disease'][i]
    if (i %% 1000 == 0){
    print(paste0('I am working on row # ',i,' out of ',nrow(data),' rows!'))}
  }
  return(data)
}

## Convert
data <- widen(data)

## Back to DT
setDT(data)

## Apply max function to get unique column values
## List of column names
out_list <- list()
i <- 1
while (i <= length(unique(data$disease))){
  out_list[[length(out_list)+1]] <- paste0('incident_',unique(data$disease)[i])
  out_list[[length(out_list)+1]] <- paste0('prevalent_',unique(data$disease)[i])
  out_list[[length(out_list)+1]] <- paste0('survyears_',unique(data$disease)[i])
  i <- i+1
}

out_list <- unlist(out_list)

c[,(out_list) := (lapply(.SD,max,na.rm=T)),
  .SDcols=names(data)[(lead_columns+1):length(names(data))],by=sample_id]

data <- unique(data,by='sample_id')

save(data,file='/mnt/ml4cvd/projects/skhurshid/svt_brady_gwas/phenotypes_wide_032520.RData')