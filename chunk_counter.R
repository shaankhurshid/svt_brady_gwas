# Script to count files
library(stringr)

# Function to count chunk files
# path = path to PLINK association output files
# phenotype = 
# chunk_list = numeric vector of how chunks are expected per chromosome (OPTIONAL if return_chunks is TRUE)
# return_chunks = set to TRUE if you want to see what chunks you have; otherwise output is just missing chunks

file_counter <- function(path,phenotype,chunk_list=NULL,return_chunks=FALSE){
  files_list <- list.files(path)
  missing_list <- list(); chunks_list <- list()
  for (i in 1:22){
    pattern1 <- paste0(phenotype,'\\.chr',i,'\\.[:digit:]+\\.assoc\\.',phenotype,'\\.glm\\.logistic$')
    pattern2 <- '\\.[:digit:]+\\.'
    model_files <- files_list[(str_detect(files_list,regex(pattern1)))]
    chunks <- str_extract(str_extract(model_files,regex(pattern2)),'[:digit:]+')
    chunks_list[[i]] <- as.numeric(chunks[order(as.numeric(chunks))])
    if (length(chunk_list) > 0){
      expected_chunks <- 1:chunk_list[i]
      missing_chunks <- as.numeric(expected_chunks[!(expected_chunks %in% chunks)])
    missing_list[[i]] <- missing_chunks[order(missing_chunks)]}
  }
if (return_chunks==TRUE){return(chunks_list)} else {return(missing_list)}
}

chunk_list <- c(55,60,50,49,45,43,40,39,31,34,34,33,25,23,21,23,20,20,16,16,10,10)

missing_snd_hard <- file_counter(path='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/svt_brady_gwas/analysis_snd_hard/temp',
                      phenotype='Bradyarrhythmia_sinus_node_dysfunction_HARD_V2',
                      chunk_list=chunk_list,return_chunks=FALSE)

missing_snd <- file_counter(path='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/svt_brady_gwas/analysis2/temp',
                      phenotype='Bradyarrhythmia_sinus_node_dysfunction',
                      chunk_list=chunk_list,return_chunks=FALSE)

