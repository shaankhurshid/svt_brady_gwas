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
  a<-fread(paste0("/medpop/esp2/projects/UK_Biobank/Phenotype_Library/phenoV2/disease/",trait,'.tab.tsv'),header=T,data.table=F) #503629
  a<-a[!is.na(a[,2]),]
  print(paste0('Total N:',nrow(a)))
  
  ##sample qc file
  b<-fread("/medpop/esp2/pradeep/UKBiobank/v2data/ukb_sqc_v2_7089.tsv",header=T,data.table=F) ##488377
  ab<-merge(a,b,by.x="sample_id",by.y="eid") #488374
  ab$array_UKBB<-ifelse(ab$genotyping_array=="UKBB",1,0)
  print(paste0('N after merge with sample QC file:',nrow(ab)))
  
  ##white British, Irish, and other white
  white<-read.table("/medpop/afib/chaffin/ukbb_pca/UKB.BritishIrishOtherWhite.inclusion.list.v1.txt",header=T)
  ab$white<-ifelse(ab$sample_id%in%white$eid,1,0)
  
  ##exclusion files
  drop<-read.csv("/medpop/afib/data/ukb9389/exclusion/w1748_20170726.csv",header=F)
  link<-read.csv("/medpop/afib/data/ukb9389/ukb_app17488_app7089_link.csv",header=T)
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
  rel<-fread("/medpop/esp2/pradeep/UKBiobank/v2data/ukb708_rel_chr1_s488374.dat",header=T,data.table=F)
  pheno<-ab[,colnames(ab)%in%c("sample_id",trait)]
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
  ab$ex_rel<-ifelse(ab$sample_Id%in%rel_id,1,0) # 81837  ->81678
  print(paste0('N related:',sum(ab$ex_rel)))
  
  #######
  ###test AF related PCs in each cleaned dataset
  #######
  
  #dim(subset(ab, used_in_pca_calculation==1 & ex_sex==0 & ex_poor==0 & ex_misKin==0 & white==1))
  #######
  ##all white, no relatives
  #######
  ab1<-subset(ab, ex_rel==0 & white==1)
  form1<-formula(paste0(has_disease,"~enroll_age + ",paste0("PC",1:40,collapse="+"),"+ array_UKBB + male",collapse="+"))
  s1<-summary(glm(form1,family="binomial",data=ab1))$coefficients
  s1<-s1[substring(rownames(s1),1,2)=="PC",]
  ab1.1<-ab1[,colnames(ab1)%in%c("sample_id",trait)]
  caseN_1<-dim(ab1.1[ab1.1[,2]==1 & !is.na(ab1.1[,2]) ,])[1]
  ctrlN_1<-dim(ab1.1[ab1.1[,2]==0 & !is.na(ab1.1[,2]) ,])[1]
  allN_1<-dim(ab1.1[!is.na(ab1.1[,2]),])[1]
  male_1<-dim(ab1[ab1$male==1 & !is.na(ab1.1[,2]),])[1]
  pcs1<-paste(rownames(s1[s1[,4]<0.05,]),collapse=",")
  
  #######
  ##all white
  #######
  ab2<-subset(ab, white==1) 
  s2<-summary(glm(form1,family="binomial",data=ab2))$coefficients
  s2<-s2[substring(rownames(s2),1,2)=="PC",]
  ab2.1<-ab2[,colnames(ab2)%in%c("sample_id",trait)]
  caseN_2<-dim(ab2.1[ab2.1[,2]==1 & !is.na(ab2.1[,2]) ,])[1]
  ctrlN_2<-dim(ab2.1[ab2.1[,2]==0 & !is.na(ab2.1[,2]) ,])[1]
  allN_2<-dim(ab2.1[!is.na(ab1.1[,2]),])[1]
  male_2<-dim(ab2[ab2$male==1 & !is.na(ab2.1[,2]),])[1]
  pcs2<-paste(rownames(s2[s2[,4]<0.05,]),collapse=",")
  
  #######
  ##all individuals, no relatives
  #######
  ab3<-subset(ab, ex_rel==0 )
  s3<-summary(glm(form1,family="binomial",data=ab3))$coefficients
  s3<-s3[substring(rownames(s3),1,2)=="PC",]
  ab3.1<-ab3[,colnames(ab3)%in%c("sample_id",trait)]
  caseN_3<-dim(ab3.1[ab3.1[,2]==1 & !is.na(ab3.1[,2]) ,])[1]
  ctrlN_3<-dim(ab3.1[ab3.1[,2]==0 & !is.na(ab3.1[,2]) ,])[1]
  allN_3<-dim(ab3.1[ !is.na(ab3.1[,2]),])[1]
  male_3<-dim(ab3[ab3$male==1 &  !is.na(ab3.1[,2]),])[1]
  pcs3<-paste(rownames(s3[s3[,4]<0.05,]),collapse=",")
  
  #######
  ##all individuals
  #######
  ab4<-subset(ab, ex_sex==0 & ex_poor==0 & ex_misKin==0 )
  s4<-summary(glm(form1,family="binomial",data=ab4))$coefficients
  s4<-s4[substring(rownames(s4),1,2)=="PC",]
  ab4.1<-ab4[,colnames(ab4)%in%c("sample_id",trait)]
  caseN_4<-dim(ab4.1[ab4.1[,2]==1 & !is.na(ab4.1[,2]) ,])[1]
  ctrlN_4<-dim(ab4.1[ab4.1[,2]==0 & !is.na(ab4.1[,2]) ,])[1]
  allN_4<-dim(ab4.1[!is.na(ab4.1[,2]),])[1]
  male_4<-dim(ab4[ab4$male==1 & !is.na(ab4.1[,2]),])[1]
  pcs4<-paste(rownames(s4[s4[,4]<0.05,]),collapse=",")
  
  #######
  ##all British white, no relatives
  #######
  ab5<-subset(ab, ex_rel==0 & in_white_British_ancestry_subset==1)
  s5<-summary(glm(form1,family="binomial",data=ab5))$coefficients
  s5<-s5[substring(rownames(s5),1,2)=="PC",]
  ab5.1<-ab5[,colnames(ab5)%in%c("sample_id",trait)]
  caseN_5<-dim(ab5.1[ab5.1[,2]==1 & !is.na(ab5.1[,2]) ,])[1]
  ctrlN_5<-dim(ab5.1[ab5.1[,2]==0 & !is.na(ab5.1[,2]) ,])[1]
  allN_5<-dim(ab5.1[!is.na(ab5.1[,2]),])[1]
  male_5<-dim(ab5[ab5$male==1 & !is.na(ab5.1[,2]),])[1]
  pcs5<-paste(rownames(s5[s5[,4]<0.05,]),collapse=",")
  
  #######
  ##all British white
  #######
  ab6<-subset(ab, in_white_British_ancestry_subset==1) 
  s6<-summary(glm(form1,family="binomial",data=ab6))$coefficients
  s6<-s6[substring(rownames(s6),1,2)=="PC",]
  ab6.1<-ab6[,colnames(ab6)%in%c("sample_id",trait)]
  caseN_6<-dim(ab6.1[ab6.1[,2]==1 & !is.na(ab6.1[,2]) ,])[1]
  ctrlN_6<-dim(ab6.1[ab6.1[,2]==0 & !is.na(ab6.1[,2]) ,])[1]
  allN_6<-dim(ab6.1[!is.na(ab6.1[,2]) ,])[1]
  male_6<-dim(ab6[ab6$male==1 & !is.na(ab6.1[,2]) ,])[1]
  pcs6<-paste(rownames(s6[s6[,4]<0.05,]),collapse=",")
  
  
  #######
  ##create summary file
  #######
  t1<-c("case",caseN_1,"control",ctrlN_1,"all",allN_1,"mean_age",round(mean(ab1[!is.na(ab1[,2]) ,]$enroll_age),2),"sd_age",round(sd(ab1[!is.na(ab1[,2]) ,]$enroll_age),2),"male_N",male_1,"male%",round(mean(ab1[!is.na(ab1[,2]) ,]$male)*100,2),"related-PCs",ifelse((pcs1!=""),pcs1,"None"))
  t2<-c("case",caseN_2,"control",ctrlN_2,"all",allN_2,"mean_age",round(mean(ab2[!is.na(ab2[,2]) ,]$enroll_age),2),"sd_age",round(sd(ab2[!is.na(ab2[,2]) ,]$enroll_age),2),"male_N",male_2,"male%",round(mean(ab2[!is.na(ab2[,2]) ,]$male)*100,2),"related-PCs",ifelse((pcs2!=""),pcs2,"None"))
  t3<-c("case",caseN_3,"control",ctrlN_3,"all",allN_3,"mean_age",round(mean(ab3[!is.na(ab3[,2]) ,]$enroll_age),2),"sd_age",round(sd(ab3[!is.na(ab3[,2]) ,]$enroll_age),2),"male_N",male_3,"male%",round(mean(ab3[!is.na(ab3[,2]) ,]$male)*100,2),"related-PCs",ifelse((pcs3!=""),pcs3,"None"))
  t4<-c("case",caseN_4,"control",ctrlN_4,"all",allN_4,"mean_age",round(mean(ab4[!is.na(ab4[,2]) ,]$enroll_age),2),"sd_age",round(sd(ab4[!is.na(ab4[,2]) ,]$enroll_age),2),"male_N",male_4,"male%",round(mean(ab4[!is.na(ab4[,2]) ,]$male)*100,2),"related-PCs",ifelse((pcs4!=""),pcs4,"None"))
  t5<-c("case",caseN_5,"control",ctrlN_5,"all",allN_5,"mean_age",round(mean(ab5[!is.na(ab5[,2]) ,]$enroll_age),2),"sd_age",round(sd(ab5[!is.na(ab5[,2]) ,]$enroll_age),2),"male_N",male_5,"male%",round(mean(ab5[!is.na(ab5[,2]) ,]$male)*100,2),"related-PCs",ifelse((pcs5!=""),pcs5,"None"))
  t6<-c("case",caseN_6,"control",ctrlN_6,"all",allN_6,"mean_age",round(mean(ab6[!is.na(ab6[,2]) ,]$enroll_age),2),"sd_age",round(sd(ab6[!is.na(ab6[,2]) ,]$enroll_age),2),"male_N",male_6,"male%",round(mean(ab6[!is.na(ab6[,2]) ,]$male)*100,2),"related-PCs",ifelse((pcs6!=""),pcs6,"None"))
  
  
  t<-cbind(t4,t3,t6,t5,t2,t1)
  colnames(t)<-c("ALL","ALL-noRelatives","British_EUR","British_EUR-noRelatives","EUR","EUR-noRelatives")
  write.table(t,file=paste0("/medpop/afib/skhurshid/svt_brady_gwas/summary_",trait,".txt"),row.names=F,quote=F,sep="\t")
  
  
  #######
  ##create phenotype file
  #######
  pheno<-ab4[,colnames(ab4)%in%c("sample_id",trait,"enroll_age","ex_rel")]
  write.table(pheno,file=paste0('/medpop/afib/skhurshid/svt_brady_gwas/',trait,".tsv"),sep="\t",col.names =T,row.names = F,quote = F)
}

#sqc<-ab4[,c(1,29,31:70,74,75,80)]#486553
#write.table(sqc,file="/medpop/afib/lcweng/UKBB_all/QC/sqc.tsv",sep="\t",col.names =T,row.names = F,quote = F)


create("Bradyarrhythmia_AV_block_or_distal_conduction_disease")
