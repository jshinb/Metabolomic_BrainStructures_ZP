wd = '/scratch/jshinb/run_gwas_regenie_20260326_zp/results/WM_metPC_loadingVars'
setwd(wd)

library(data.table)
library(dplyr)
library(stringr)

gwas_WM.PC1 = fread('gwas_WM_metPC_loadingVars_ADO1_WM.PC1.regenie.gz')
gwas_Met.PC3 = fread('gwas_WM_metPC_loadingVars_ADO1_Met.PC3.regenie.gz')

topSNP = 'rs174546'
assoc_wi_rs174546 = rbind(
gwas_WM.PC1 %>% filter(ID == topSNP),
gwas_Met.PC3 %>% filter(ID == topSNP))
print(assoc_wi_rs174546)  

#   CHROM   GENPOS       ID ALLELE0 ALLELE1   A1FREQ  INFO     N   TEST
#   <int>    <int>   <char>  <char>  <char>    <num> <num> <int> <char>
#1:    11 61569830 rs174546       T       C 0.663244     1   974    ADD
#2:    11 61569830 rs174546       T       C 0.661994     1   963    ADD
#
#        BETA       SE    CHISQ    LOG10P  EXTRA            P
#       <num>    <num>    <num>     <num> <lgcl>        <num>
#1: 0.0686652 0.059885  1.31473  0.599394     NA 2.515394e-01
#2: 1.9220000 0.224011 73.61500 17.022600     NA 9.492924e-18
