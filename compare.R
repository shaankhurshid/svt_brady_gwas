# Depends
library(data.table)

# Create summarizing function
summarize <- function(path){
all_output <- data.table()
  for (i in list.files(path)[str_detect(list.files(path),'\\.gp')]){
    dataset <- fread(paste0(path,i))
    output <- data.table(phenotype=names(dataset)[3],events=sum(dataset[,3]==2),total=nrow(dataset))
    all_output <- rbind(all_output,output)
  }
return(all_output)
}

summarize_lcw <- function(path){
all_output <- data.table()
  for (i in list.files(path)[str_detect(list.files(path),'\\.gp') & !str_detect(list.files(path),'summary')]){
    dataset <- fread(paste0(path,i))
    output <- data.table(phenotype=str_remove(i,path),events=sum(dataset[,3]==2),total=nrow(dataset))
    all_output <- rbind(all_output,output)
  }
  return(all_output)
}

sk <- summarize(path='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/processed_phenotypes/')
lcw <- summarize_lcw(path='/Volumes/medpop_afib/lcweng/UKBB_all/SVT_Brady/UKbiobank_analysis/phenotype/')