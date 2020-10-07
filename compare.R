# Script to compare LCW phenotypes with mine

# Dependencies
library(data.table)

# Load files
## Datasets
brady_lcw <- fread(file='/Volumes/medpop_afib/lcweng/UKBB_all/SVT_Brady/UKbiobank_analysis/phenotype/Bradyarrhythmia_Pacemaker_v2.tsv')
brady_sk <- fread(file='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/processed_phenotypes/Bradyarrhythmia_Pacemaker_v2.tsv')

## Disease dataset
brady <- fread(file='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/Bradyarrhythmia_Pacemaker_v2.tab.tsv')

## Exclusion datasets
valve <- fread(file='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/Valvular_disease_unspecified.tab.tsv')
mi <- fread(file='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/Myocardial_Infarction.tab.tsv')
ct <- fread(file='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/Cardiac_surgery.tab.tsv')
snd <- fread(file='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/Bradyarrhythmia_sinus_node_dysfunction.tab.tsv')
dcd <- fread(file='/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/Bradyarrhythmia_AV_block_or_distal_conduction_disease.tab.tsv')

# Define discrepant set
discrepancy <- brady_lcw[!(IID %in% brady_sk$IID)] # N=245, all cases
setkey(discrepancy,IID); setkey(brady,sample_id)
discrepancy[brady,disease_date := i.censor_date]

# Go through diseases
## Valve
has_valve <- valve[has_disease==1]
setkey(has_valve,sample_id); setkey(discrepancy,IID)
discrep_valve <- has_valve[discrepancy,ppm_date := i.disease_date]
valve_ppm <- discrep_valve[!is.na(ppm_date),c('prevalent_disease','censor_date','ppm_date')]

