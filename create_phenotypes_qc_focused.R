############################
######## create phenotype file ----will need AF, age, sex, array, and pcs
############################

##batch, genetic sex, pc, and EUR are in ukb_sqc_v2.txt
##created by Mary
#sqc<-fread("/broad/ukbb/imputed/ukb_sqc_v2.txt",header=F)
#fam<-read.table('/path/to/ukb###_cal_chr1_v2_s488374.fam',sep=" ",header=F)
#sqc_7089<-as.data.frame(fam$V1)
#sqc_7089<-cbind.data.frame(fam$V1,sqc[,3:68])
#sqc_colnames<-read.table('/medpop/afib/lcweng/UKBB_all/QC/sqc_colnames.txt',header=F)
#colnames(sqc_7089)<-sqc_colnames$V1
#write.table(sqc_7089,'/medpop/esp2/pradeep/UKBiobank/v2data/ukb_sqc_v2_7089.tsv',sep="\t",col.names =T,row.names = F,quote = F)

############################
######## create phenotype file ----will need ID, phenotype, and baseline age
############################
library(data.table)

create<-function(trait){
  ##phenotype file
  a<-fread(paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/",trait,'.tab.tsv'),header=T,data.table=F) #503629
  names(a)[names(a)=='has_disease'] <- trait
  a<-a[!is.na(a[,3]),]
  print(paste0('Total N:',nrow(a)))
  
  ##sample qc file
  b<-fread("/Volumes/medpop_esp2/pradeep/UKBiobank/v2data/ukb_sqc_v2_7089.tsv",header=T,data.table=F) ##488377
  ab<-merge(a,b,by.x="sample_id",by.y="eid") #488374
  ab$array_UKBB<-ifelse(ab$genotyping_array=="UKBB",1,0)
  print(paste0('N after merge with sample QC file:',nrow(ab)))
  
  ##white British, Irish, and other white
  white<-read.table("/Volumes/medpop_afib/chaffin/ukbb_pca/UKB.BritishIrishOtherWhite.inclusion.list.v1.txt",header=T)
  ab$white<-ifelse(ab$sample_id%in%white$eid,1,0)
  
  ##exclusion files
  drop<-read.csv("/Volumes/medpop_afib/data/ukb9389/exclusion/w17488_20200204.csv",header=F)
  link<-read.csv("/Volumes/medpop_afib/data/ukb9389/ukb_app17488_app7089_link.csv",header=T)
  drop.1<-link[link$app17488%in%drop$V1,]$app7089
  ab$ex_drop<-ifelse(ab$sample_id%in%drop.1,1,0)
  
  ##remove poor quality
  ab$ex_sex<-ifelse(ab$Submitted_Gender==ab$Inferred_Gender,0,1) #378
  ab$ex_poor<-ifelse(ab$het_missing_outliers==1 | ab$putative_sex_chromosome_aneuploidy==1,1,0) #1620
  ab$ex_misKin<-ifelse(ab$excluded_from_kinship_inference==1,1,0) #977
  
  #high quality data
  ab <- ab[ab$ex_sex==0,]
  print(paste0('N after removal of sex mismatch:',nrow(ab)))
  ab <- ab[ab$ex_poor==0,]
  print(paste0('N after removal of poor:',nrow(ab)))
  ab <- ab[ab$ex_misKin==0,]
  print(paste0('N after removal of missing kinship inference:',nrow(ab)))
  ab <- ab[ab$ex_drop==0,]
  print(paste0('N after removal of withdrawn consent:',nrow(ab)))
  
  #This file lists the pairs of individuals related up to the third degree in the data set. It is a plaintext file with space separated columns.
  rel<-fread("/Volumes/medpop_esp2/pradeep/UKBiobank/v2data/ukb708_rel_chr1_s488374.dat",header=T,data.table=F)
  pheno<-ab[,c('sample_id',trait)]
  names(pheno)<-c("ID","pheno")
  rel<-merge(rel,pheno,by.x="ID1",by.y="ID",all.x=T)
  rel<-merge(rel,pheno,by.x="ID2",by.y="ID",suffixes=c(".ID1",".ID2"),all.x=T)
  rel$rel_ID<-NA
  rel$rel_ID[which(rel$pheno.ID1==1 )]<-rel$ID2[which(rel$pheno.ID1==1)]
  rel$rel_ID[which(rel$pheno.ID2==1 & (rel$pheno.ID1==0 | is.na(rel$pheno.ID1)))]<-rel$ID1[which(rel$pheno.ID2==1 & (rel$pheno.ID1==0 | is.na(rel$pheno.ID1)))]
  rel$rel_ID[which(rel$pheno.ID1==0 & rel$pheno.ID2==0)]<-rel$ID2[which(rel$pheno.ID1==0 & rel$pheno.ID2==0)]
  rel$rel_ID[which(rel$pheno.ID1==0 & is.na(rel$pheno.ID2))]<-rel$ID2[which(rel$pheno.ID1==0 & is.na(rel$pheno.ID2))]
  rel$rel_ID[which(rel$pheno.ID2==0 & is.na(rel$pheno.ID1))]<-rel$ID2[which(rel$pheno.ID2==0 & is.na(rel$pheno.ID1))]
  
  
  
  rel_id<-unique(rel$rel_ID) #81838 ->81847 
  ab$male<-ifelse(ab$Inferred_Gender=="M",1,0) 
  ab$ex_rel<-ifelse(ab$sample_id%in%rel_id,1,0) # 81837  ->81678
  print(paste0('N related:',sum(ab$ex_rel)))
  
  #######
  ###test AF related PCs in each cleaned dataset
  #######
  
  #dim(subset(ab, used_in_pca_calculation==1 & ex_sex==0 & ex_poor==0 & ex_misKin==0 & white==1))
  #######
  ##all white, no relatives
  #######
  ab1<-subset(ab, ex_rel==0 & white==1)
  form1<-formula(paste0(trait,"~enroll_age + ",paste0("PC",1:40,collapse="+"),"+ array_UKBB + male",collapse="+"))
  s1<-summary(glm(form1,family="binomial",data=ab1))$coefficients
  s1<-s1[substring(rownames(s1),1,2)=="PC",]
  ab1.1<-ab1[,colnames(ab1)%in%c("sample_id",trait)]
  caseN_1<-dim(ab1.1[ab1.1[,2]==1 & !is.na(ab1.1[,2]) ,])[1]
  ctrlN_1<-dim(ab1.1[ab1.1[,2]==0 & !is.na(ab1.1[,2]) ,])[1]
  allN_1<-dim(ab1.1[!is.na(ab1.1[,2]),])[1]
  male_1<-dim(ab1[ab1$male==1 & !is.na(ab1.1[,2]),])[1]
  pcs1<-paste(rownames(s1[s1[,4]<0.05,]),collapse=",")
  
  #######
  ##create summary file
  #######
  t1<-c("case",caseN_1,"control",ctrlN_1,"all",allN_1,"mean_age",round(mean(ab1[!is.na(ab1[,2]) ,]$enroll_age),2),"sd_age",round(sd(ab1[!is.na(ab1[,2]) ,]$enroll_age),2),"male_N",male_1,"male%",round(mean(ab1[!is.na(ab1[,2]) ,]$male)*100,2),"related-PCs",ifelse((pcs1!=""),pcs1,"None"))
  write.table(t1,file=paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/summary_",trait,".txt"),row.names=F,quote=F,sep="\t")
  
  
  #######
  ##create phenotype file
  #######
  pheno<-ab1[,c("sample_id",'trait',"male","enroll_age","array",rownames(s1[s1[,4]<0.05,]))]
  write.table(pheno,file=paste0('/Volumes/medpop_afib/skhurshid/svt_brady_gwas/',trait,".tsv"),sep="\t",col.names =T,row.names = F,quote = F)
}

#sqc<-ab4[,c(1,29,31:70,74,75,80)]#486553
#write.table(sqc,file="/medpop/afib/lcweng/UKBB_all/QC/sqc.tsv",sep="\t",col.names =T,row.names = F,quote = F)

pheno_list <- c('Bradyarrhythmia_AV_block_or_distal_conduction_disease',
                'Bradyarrhythmia_AV_block_or_distal_conduction_disease_HARD_V2',
                'Bradyarrhythmia_Pacemaker_v2',
                'Bradyarrhythmia_sinus_node_dysfunction',
                'Bradyarrhythmia_sinus_node_dysfunction_HARD_V2',
                'Supraventricular_arrhythmia_WPW_v2',
                'Supraventricular_arrhythmia_SVT')

for (i in pheno_list){
  create(i)
}
