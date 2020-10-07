
############################
######## script to create processed phenotypes
############################
library(data.table)

# Function to create phenotypes
## TRAIT = single phenotype as string for which you want to build a dataset
## EXCLUDE ALL BOTH = vector of phenotype names as strings where you want to exclude both cases and controls with condition at any time
## EXCLUDE ALL CASES = vector of phenotype names as strings where you want to exclude cases with the condition at any time
## EXCLUDE ALL CONTROLS = vector of phenotype names as strings where you want to exclude controls with the condition at any time
## EXCLUDE INCIDENT CASES = vector of phenotype names as strings where you want to exclude cases with the condition BEFORE case diagnosis
## EXCLUDE FLEXIBLE = vector of phenotype names as strings where you want to exclude cases with the condition BEFORE case diagnosis AND controls with the condition at any time

create<-function(trait,exclude_all_both=NULL,exclude_all_cases=NULL,exclude_all_controls=NULL,
                 exclude_incident_cases=NULL,exclude_flexible=NULL,
                 pheno_path){
  ##phenotype file
  a<-fread(paste0(pheno_path,trait,'.csv'),header=T) #503629
  setnames(a,'has_disease',trait)
  a<-a[!is.na(get(trait))]
  print(paste0('Total N:',nrow(a)))
  
  ##sample qc file
  b<-fread("/Volumes/medpop_esp2/pradeep/UKBiobank/v2data/ukb_sqc_v2_7089.tsv",header=T) ##488377
  setkey(a,'sample_id'); setkey(b,'eid')
  ab <- a[b,nomatch=0]
  ab[,':='(array_UKBB = ifelse(genotyping_array=='UKBB',1,0))]
  print(paste0('N after merge with sample QC file:',nrow(ab)))
    
  ##white British, Irish, and other white
  white<-fread("/Volumes/medpop_afib/chaffin/ukbb_pca/UKB.BritishIrishOtherWhite.inclusion.list.v1.txt",header=T)
  ab[,':='(white = ifelse(sample_id %in% white$eid,1,0))]
  
  ##exclusion files
  drop<-fread("/Volumes/medpop_afib/skhurshid/phenotypes/withdrawals/w7089_20200820.csv",header=T)
  drop.1<-drop$sample_id
  ab[,':='(ex_drop = ifelse(sample_id %in% drop.1,1,0))]
  
  ##remove poor quality
  ab[,':='(ex_poor = ifelse(het_missing_outliers==1 | putative_sex_chromosome_aneuploidy==1,1,0),
           ex_sex = ifelse(Submitted_Gender==Inferred_Gender,0,1),
           ex_misKin = ifelse(ab$excluded_from_kinship_inference==1,1,0))]
  
  #high quality data
  ab <- ab[ab$ex_sex==0]
  print(paste0('N after removal of sex mismatch:',nrow(ab)))
  ab <- ab[ab$ex_poor==0]
  print(paste0('N after removal of poor:',nrow(ab)))
  ab <- ab[ab$ex_misKin==0]
  print(paste0('N after removal of missing kinship inference:',nrow(ab)))
  ab <- ab[ab$ex_drop==0]
  print(paste0('N after removal of withdrawn consent:',nrow(ab)))
  
  # Loop over "exclude all both" phenotypes - all individuals with exclusion phenotype at any time removed for both cases/controls
  if (length(exclude_all_both)!=0){
  for (i in exclude_all_both){
    exclude <- fread(paste0(pheno_path,i,'.csv'),header=T)
    setkey(ab,sample_id); setkey(exclude,sample_id)
    ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
    ab[,exclude := ifelse(c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd == 1)),1,0)]
       print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for cases/controls'))
       ab <- ab[exclude==0]
       ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
  }}
  
  # Loop over "exclude all cases" phenotypes - all individuals with exclusion phenotype at any time removed for cases
  if (length(exclude_all_cases)!=0){
  for (i in exclude_all_cases){
    exclude <- fread(paste0(pheno_path,i,'.csv'),header=T)
    setkey(ab,sample_id); setkey(exclude,sample_id)
    ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
    ab[,exclude := ifelse(c(c(!is.na(get(trait)) & get(trait)==1) & 
                              c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd==1))),1,0)]
    print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for cases only'))
    ab <- ab[exclude==0]
    ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
  }}
  
  # Loop over "exclude all controls" phenotypes - all individuals with exclusion phenotype at any time removed for controls
  if (length(exclude_all_controls)!=0){
   for (i in exclude_all_controls){
    exclude <- fread(paste0(pheno_path,i,'.csv'),header=T)
    setkey(ab,sample_id); setkey(exclude,sample_id)
    ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
    ab[,exclude := ifelse(c(c(get(trait)==0 | is.na(get(trait))) & 
                              c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd==1))),1,0)]
    print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for controls only'))
    ab <- ab[exclude==0]
    ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
  }}
  
  # Loop over "exclude incident" - only cases with exclusion phenotype before disease removed
  if (length(exclude_incident_cases)!=0){
    for (i in exclude_incident_cases){
      exclude <- fread(paste0(pheno_path,i,'.csv'),header=T)
      setkey(ab,sample_id); setkey(exclude,sample_id)
      ab[exclude,':='(exclude_disease = i.has_disease, exclude_prev = i.prevalent_dsease, exclude_censor = i.censor_date)]
      ab[,exclude := ifelse(c(c(!is.na(get(trait)) & get(trait)==1) & 
                                c(c(!is.na(exclude_disease) & (exclude_censor <= censor_date)) |
                                c(!is.na(exclude_disease) & exclude_prev==1))),1,0)]
      print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring before case diagnosis'))
      ab <- ab[exclude==0]
      ab <- ab[,!(c('exclude_prev','exclude_censor','exclude_disease','exclude'))]
    }}
  
  # Loop over "exclude flexible" - excludes any instance of exclusion phenotype among controls, and only exclusion phenotype prior to disease for cases
  if (length(exclude_flexible)!=0){
    for (i in exclude_flexible){
      exclude <- fread(paste0(pheno_path,i,'.csv'),header=T)
      setkey(ab,sample_id); setkey(exclude,sample_id)
      ab[exclude,':='(exclude_incd = i.incident_disease, exclude_disease = i.has_disease, exclude_prev = i.prevalent_disease, exclude_censor = i.censor_date)]
      ab[,exclude := ifelse(c(!is.na(get(trait)) & get(trait)==1),
                            ifelse(c(!is.na(exclude_prev) & exclude_prev==1),1,
                            ifelse(c(!is.na(exclude_incd) & (exclude_incd==1) & (exclude_censor <= censor_date)),1,0)),
                            ifelse(c(!is.na(exclude_disease) & exclude_disease==1),1,0))]
      print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring before case diagnosis or at any time for controls'))
      ab <- ab[exclude==0]
      ab <- ab[,!(c('exclude_disease','exclude_prev','exclude_incd','exclude_censor','exclude'))]
    }}
  
  #This file lists the pairs of individuals related up to the third degree in the data set. It is a plaintext file with space separated columns.
  rel<-fread("/Volumes/medpop_esp2/pradeep/UKBiobank/v2data/ukb708_rel_chr1_s488374.dat",header=T)
  pheno<-ab[,.SD,.SDcols=c('sample_id',trait)]
  names(pheno)<-c("ID","pheno")
  
  # Merge phenotypes with relative index
  setkey(rel,ID1); setkey(pheno,ID)
  rel <- rel[pheno,pheno.ID1:=i.pheno]
  setkey(rel,ID2)
  rel <- rel[pheno,pheno.ID2:=i.pheno]
  
  # Now create a list of people to remove for relatedness
  rel[,rel_id := ifelse(is.na(pheno.ID1) & pheno.ID2==0,ID1,
                        ifelse(is.na(pheno.ID2) & pheno.ID1==0,ID2,
                               ifelse(is.na(pheno.ID2),ID2,
                                      ifelse(c(pheno.ID2==1 & c(pheno.ID1==0 | is.na(pheno.ID1))),ID1,
                                             ifelse(pheno.ID1==0 & pheno.ID2==0,ID2,
                                                    ifelse(pheno.ID1==1,ID2,NA))))))]
  
  
  
  rel_id<-unique(rel$rel_id) #81838 ->81847 
  
  ab[,':='(male = ifelse(Inferred_Gender=='M',1,0),
           ex_rel = ifelse(sample_id %in% rel_id,1,0))]
  
  print(paste0('N related:',sum(ab$ex_rel)))
  
  #######
  ###test AF related PCs in each cleaned dataset
  #######
  
  #dim(subset(ab, used_in_pca_calculation==1 & ex_sex==0 & ex_poor==0 & ex_misKin==0 & white==1))
  #######
  ##all white, no relatives
  #######
  ab1<-ab[white==1 & ex_rel==0]
  form1<-formula(paste0(trait,"~enroll_age + ",paste0("PC",1:40,collapse="+"),"+ array_UKBB + male",collapse="+"))
  s1<-summary(glm(form1,family="binomial",data=ab1))$coefficients
  s1<-s1[substring(rownames(s1),1,2)=="PC",]
  ab1.1<-ab1[,.SD,.SDcols=c("enroll_age",trait)]
  caseN_1<-nrow(ab1.1[ab1.1[[trait]]==1 & !is.na(ab1.1[[trait]])])
  ctrlN_1<-nrow(ab1.1[ab1.1[[trait]]==0 & !is.na(ab1.1[[trait]])])
  allN_1<-nrow(ab1.1)
  male_1<-nrow(ab1[ab1[["male"]]==1 & !is.na(ab1.1[[trait]])])
  pcs1<-paste(rownames(s1[s1[,4]<0.05,]),collapse=",")
  
  #######
  ##create summary file
  #######
  t1<-c("case",caseN_1,"control",ctrlN_1,"all",allN_1,"mean_age",round(mean(ab1[!is.na(trait)]$enroll_age),2),"sd_age",round(sd(ab1[!is.na(trait)]$enroll_age),2),"male_N",male_1,"male%",round(mean(ab1[!is.na(trait)]$male)*100,2),"related-PCs",ifelse((pcs1!=""),pcs1,"None"))
  write.table(t1,file=paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/summary_",trait,".txt"),row.names=F,quote=F,sep="\t")
  
  
  #######
  ##create phenotype file
  #######
  ## Choose columns
  pheno<-ab1[,c("sample_id",trait,"enroll_age",rownames(s1)[1:5],"array_UKBB","male"),with=F]
  ## Format for PLINK
  setnames(pheno,"sample_id","FID")
  pheno[,':='(IID = FID)]
  pheno[,eval(trait) := ifelse(c(!is.na(get(trait)) & get(trait)==0),1,2)]
  setcolorder(pheno,c('FID','IID'))
  print(paste0('Final phenotype N: ',nrow(pheno)))
  write.table(pheno,file=paste0('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/processed_phenotypes/v2/',trait,".tsv"),sep="\t",col.names =T,row.names = F,quote = F)
}

#sqc<-ab4[,c(1,29,31:70,74,75,80)]#486553
#write.table(sqc,file="/medpop/afib/lcweng/UKBB_all/QC/sqc.tsv",sep="\t",col.names =T,row.names = F,quote = F)

#exclusion phenotypes
# SND - exclude valve disease, cardiac surgery, MI at or prior to SND (cases)
# SND - exclude valve disease, cardiac surgery, MI, SND-inclusive, DCD-inclusive, PM
# DCD - exclude valve disease, cardiac surgery, MI at or prior to DCD cases)
# DCD - exclude valve disease, cardiac surgery, MI, SND-inclusive, DCD-inclusive, PM (controls)
# PM - exclude valve disease, cardiac surgery, MI at or prior to PM (cases)
# PM - exclude valve disease, cardiac surgery, MI, SND-inclusive, DCD-inclusive, PM (controls)
# WPW - exclude HCM, Ebstein (cases)
# WPW - exclude HCM, Ebstein (controls)
# SVT - exclude HCM, Ebstein (cases)
# SVT - exclude HCM, Ebstein (controls)

create(trait="Bradyarrhythmia_AV_block_or_distal_conduction_disease",
       exclude_flexible=c("Cardiac_surgery","Myocardial_infarction","Valvular_disease_unspecified"),
       exclude_all_controls=c('Bradyarrhythmia_sinus_node_dysfunction','Bradyarrhythmia_Pacemaker_v2'),
       pheno_path = '/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/svt_brady/')

create(trait="Bradyarrhythmia_AV_block_or_distal_conduction_disease_HARD_V2",
       exclude_flexible=c("Cardiac_surgery","Myocardial_infarction","Valvular_disease_unspecified"),
       exclude_all_controls=c('Bradyarrhythmia_sinus_node_dysfunction','Bradyarrhythmia_Pacemaker_v2'),
       pheno_path = '/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/svt_brady/')

create(trait="Bradyarrhythmia_Pacemaker_v2",
       exclude_flexible=c("Cardiac_surgery","Myocardial_infarction","Valvular_disease_unspecified"),
       exclude_all_controls=c('Bradyarrhythmia_sinus_node_dysfunction','Bradyarrhythmia_AV_block_or_distal_conduction_disease'),
       pheno_path = '/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/svt_brady/')

create(trait="Bradyarrhythmia_sinus_node_dysfunction",
       exclude_flexible=c("Cardiac_surgery","Myocardial_infarction","Valvular_disease_unspecified"),
       exclude_all_controls=c('Bradyarrhythmia_AV_block_or_distal_conduction_disease','Bradyarrhythmia_Pacemaker_v2'),
       pheno_path = '/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/svt_brady/')

create(trait="Bradyarrhythmia_sinus_node_dysfunction_HARD_V2",
       exclude_flexible=c("Cardiac_surgery","Myocardial_infarction","Valvular_disease_unspecified"),
       exclude_all_controls=c('Bradyarrhythmia_AV_block_or_distal_conduction_disease','Bradyarrhythmia_Pacemaker_v2'),
       pheno_path = '/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/svt_brady/')

create(trait="Supraventricular_arrhythmia_WPW_v2",
       exclude_all_both=c("Hypertrophic_cardiomyopathy","Congenital_heart_disease_Ebstein_anomaly"),
       pheno_path = '/Volumes/mdpop_afib/skhurshid/phenotypes/2020_06/svt_brady/')

create(trait="Supraventricular_arrhythmia_SVT",
       exclude_all_both=c("Hypertrophic_cardiomyopathy","Congenital_heart_disease_Ebstein_anomaly"),
       pheno_path = '/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/svt_brady/')

