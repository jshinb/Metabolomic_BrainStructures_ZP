# rm(list=ls())
#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ----
#* [Put the full path of the input file]
project_specific_file = "~/Documents/scripts/Metabolomic_BrainStructures_ZP/project_specific_inputs.R"
## Excel file: Age_group_by_lipid_interaction_test_pvalues_2026-11-03.xlsx
## --Output: z_diff_and_interaction_assoc_test_results_2026-11-03.txt--
## Sent to Zdenka
## Get data from: ~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/get_data_for_interaction_test_2026-11-03.R
## Output_final: beta_diff_and_interaction_test_results_rank_within_age_group_2026-11-03.xlsx
#*****************************************************************************#

# 0. initialize ----
source(project_specific_file)
setwd(wd)
res_WM1 = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/TabS1_2026-03-17.txt")
head(res_WM1 )
res_WM1 = res_WM1 %>% filter(metaboliteid!='GlycA')

ind_ado = res_WM1$fdr.p2_WM.T1wSI_ados<0.05;print(sum(ind_ado))#650
ind_ado = ind_ado | res_WM1$fdr.p2_WM.MTR_ados <0.05;print(sum(ind_ado))#710# 
ind_ado = ind_ado | res_WM1$fdr.p2_WM.Vol_ados <0.05;print(sum(ind_ado))#711# after removing 'GlycA'

ind_ado1 = res_WM1$fdr.p2_WM.T1wSI_ados<0.05;print(sum(ind_ado1))#650
ind_ado1 = ind_ado1 & res_WM1$fdr.p2_WM.MTR_ados <0.05;print(sum(ind_ado1))#587# 
ind_ado1 = ind_ado1 & res_WM1$fdr.p2_WM.Vol_ados <0.05;print(sum(ind_ado1))#525# after removing 'GlycA'

ind_ado2 = res_WM1$fdr.p2_WM.T1wSI_ados<0.1;print(sum(ind_ado2))#672
ind_ado2 = ind_ado2 & res_WM1$fdr.p2_WM.MTR_ados <0.1;print(sum(ind_ado2))#617
ind_ado2 = ind_ado2 & res_WM1$fdr.p2_WM.Vol_ados <0.1;print(sum(ind_ado2))#564

ind_adult = res_WM1$fdr.p2_WM.T1wSI_adults<0.05;print(sum(ind_adult))#41
ind_adult = ind_adult | res_WM1$fdr.p2_WM.MTR_adults <0.05;print(sum(ind_adult))#63
ind_adult = ind_adult | res_WM1$fdr.p2_WM.Vol_adults <0.05;print(sum(ind_adult))#73

ind_adult2 = res_WM1$fdr.p2_WM.T1wSI_adults<0.1;print(sum(ind_adult2))#85
ind_adult2 = ind_adult2 & res_WM1$fdr.p2_WM.MTR_adults <0.1;print(sum(ind_adult2))#3
ind_adult2 = ind_adult2 & res_WM1$fdr.p2_WM.Vol_adults <0.1;print(sum(ind_adult2))#0

res_WM1 %>% filter(ind_ado) %>% select(beta_WM.T1wSI_ados,beta_WM.MTR_ados,beta_WM.Vol_ados) %>% cor(use='p')
res_WM1 %>% filter(ind_ado1) %>% select(beta_WM.T1wSI_ados,beta_WM.MTR_ados,beta_WM.Vol_ados) %>% cor(use='p')
res_WM1 %>% filter(ind_ado2) %>% select(beta_WM.T1wSI_ados,beta_WM.MTR_ados,beta_WM.Vol_ados) %>% cor(use='p')
res_WM1 %>% filter(ind_adult) %>% select(beta_WM.T1wSI_adults,beta_WM.MTR_adults,beta_WM.Vol_adults) %>% cor(use='p')
sig_lipids = res_WM1 %>% filter(ind_ado1) %>% pull(metaboliteid); sig_lipids%>% length#525

#
metabo_info = fread('data/metabodata_info_ados.txt')

# ----@
lipid.class = c('NeutL',
                'GPC',
                'SL',
                'NonL')#neutral lipids
nclass = length(lipid.class)*1
# cols = c('#E1AD01','#06D6A0','#118AB2','#8f00ff')#,'#EF476F',
cols <- c("#FFCE1B", "#06D6A0", "#118AB2", "#8f00ff")# barplot
class2.levels = c("Neutral lipids",
                  "Glycerophospholipids",
                  "Sphingolipids",
                  "Non-lipid metabolites")
metabo_info = fread('data/metabodata_info_ados.txt')
metabo_info = metabo_info %>% mutate(metabolite_id = str_remove(metabolite_id,"_nmol/mL"))
metabo_info = metabo_info%>% filter(metabolite_id %in% res_WM1$metaboliteid)#n=830
res_WM1 = res_WM1 %>% select(-`lipid class`)
res_WM1 = res_WM1 %>% 
  rename(metabolite_id = metaboliteid) %>%
  left_join(metabo_info)
res_WM1 = res_WM1 %>%
  mutate(class2 = case_when(
    class %in% c('TAG','DAG','MAG','CE') ~ class2.levels[1],
    class %in% c('PC','PE','PI','LPC','LPE') ~ class2.levels[2],
    class %in% c('SM','CER',"DCER","HCER","LCER") ~ class2.levels[3],
    class %in% c('Amino acids','Glycolysis related metabolites',
                 'Inflammation',
                 'Ketone bodies') ~ class2.levels[4]
  )) %>% mutate(class2 = factor(class2, class2.levels)) %>%
  rename(`lipid class` = class2)
class2.labels = paste(names(table(res_WM1 %>% filter(ind_ado1) %>% pull(`lipid class`)))," (n=",
                      table(res_WM1 %>% filter(ind_ado1) %>% pull(`lipid class`)),")",sep='')


# get the association results for the bmi-unadjusted analyses ----
## NOTE - modified the interaction tests - ranking within age-group to address 
## pipeline-differences

brain_phenos = c("Volume.WM.adjAgeSexICV", "T1wSI.WM.adjAgeSex", "MTR.WM.adjAgeSex")
run.test = function(){
  for (bi in brain_phenos){
    anal.date = '2026-11-03'
    load(file=paste0('data/data_for_interaction_tests_',anal.date,'.Rdata'),verbose = T)
    brain_phenos = c("Volume.WM.adjAgeSexICV", "T1wSI.WM.adjAgeSex", "MTR.WM.adjAgeSex")
    metabo.by.group.bmi <- c(); for (i in 1:nrow(res_WM1)){#i=1;mi = 'TAG54:7FA16:1'#res_WM1$metaboliteid[i]
      mi = res_WM1$metaboliteid[i]
      analdat = covdata_all %>% 
        left_join(brain_all %>% select(all_of(c('uniqueID',bi)))) %>%
        left_join(metabodata_all %>% select(all_of(c('uniqueID',mi))))
      analdat[['y']] <- analdat[[bi]]
      analdat[['x']] <- analdat[[mi]]
      analdat = analdat %>% data.frame
      rownames(analdat) = analdat$uniqueID
      
      #1. adjust metabolite for lipid medication status in adults
      fit1 = lmer(inormal(x) ~ on_lipid_low_med*Sex + (1|fam_id), data=analdat %>% filter(group=="adult"));summary(fit1)
      fit1 = lm(x ~ on_lipid_low_med*Sex, data=analdat %>% filter(group=="adult"))
      metabo.adjMed = data.frame(uniqueID = names(resid(fit1)),
                                 x.adjMed=resid(fit1)+ summary(fit1)$coef["(Intercept)","Estimate"])
      analdat = analdat %>% left_join(metabo.adjMed)
      analdat = analdat %>%
        mutate(x.adjMed = ifelse(group == "ado", x, x.adjMed)) %>%
        group_by(group) %>%
        mutate(
          rnt.y = inormal(y),
          rnt.x = inormal(x.adjMed)
        ) %>%
        ungroup()
      
      analdat %>% ggplot(aes(x=x, y=x.adjMed,color=group)) + geom_point(alpha=0.3)
      analdat %>% ggplot(aes(x=rnt.x, y=rnt.y,color=group)) + geom_point(alpha=0.3) +
        stat_cor() + geom_smooth(method="gam")
      
      ## ado fit ----@
      fit.ado = lmer(rnt.y~rnt.x + Sex + (1|fam_id),
                     data=analdat %>% filter(group=="ado"));
      summary(fit.ado)$coef["rnt.x",]
      
      ## adult fit ----@
      fit.adult = lm(rnt.y~rnt.x + Sex,
                     data=analdat %>% filter(group=="adult"));
      summary(fit.adult)$coef["rnt.x",]

      ## interaction test ----@
      fit.interact = lmer(rnt.y~rnt.x*group + Sex + (1|fam_id), 
                          data=analdat);summary(fit.interact)
      b.ado = summary(fit.interact)$coef["rnt.x","Estimate"]
      b.adult = b.ado + summary(fit.interact)$coef["rnt.x:groupadult","Estimate"]
      V = vcov(fit.interact)
      D = b.ado^2 - b.adult^2#>0 assoc. stronger in ado vs. adult
      g = c(2*b.ado, -2*b.adult)
      Vslopes = matrix(
        c(
          V["rnt.x", "rnt.x"],
          V["rnt.x", "rnt.x"] +
            V["rnt.x", "rnt.x:groupadult"],
          
          V["rnt.x", "rnt.x"] +
            V["rnt.x:groupadult", "rnt.x"],
          V["rnt.x", "rnt.x"] +
            V["rnt.x", "rnt.x:groupadult"] +
            V["rnt.x:groupadult", "rnt.x"] +
            V["rnt.x:groupadult", "rnt.x:groupadult"]
        ),
        nrow = 2,
        byrow = TRUE
      )
      SE_D = sqrt(t(g) %*% Vslopes %*% g)
      z.abs.diff.ado.adult = D/SE_D
      p_one_sided <- pnorm(z.abs.diff.ado.adult,lower.tail = F)
      
      fit.res.ado = summary(fit.ado)$coef["rnt.x",];
      fit.res.adult = summary(fit.adult)$coef["rnt.x",]
      fit.res.interact = summary(fit.interact)$coef["rnt.x:groupadult",]
      names(fit.res.ado) = paste(c('beta','se','df','t','pval'),'ado',sep='.')
      names(fit.res.adult) = paste(c('beta','se','t','pval'),'adult',sep='.')
      names(fit.res.interact) = paste(c('beta','se','df','t','pval'),'int',sep='.')
      
      metabo.by.group.bmi = 
        rbind(
          metabo.by.group.bmi,
          c(bname = bi,
            mname = mi,
            c(fit.res.ado,n.ado = nobs(fit.ado)), 
            c(fit.res.adult,n.adult = nobs(fit.adult)), 
            c(fit.res.interact,n.int = nobs(fit.interact)),
            c(abs.diff=D,abs.diff_se=SE_D,z.abs.diff=z.abs.diff.ado.adult,abs.diff.pval=p_one_sided)
          ))
    }
    metabo.by.group.bmi =  metabo.by.group.bmi %>% data.table
    metabo.by.group.bmi %>% head
    metabo.by.group.bmi = metabo.by.group.bmi %>% mutate(across(-c(1:2), as.numeric))
    metabo.by.group.bmi %>% arrange(pval.ado)
    metabo.by.group.bmi = metabo.by.group.bmi %>% dplyr::rename(metaboliteid=mname)
    metabo.by.group.bmi %>%
      write_tsv(paste0('results/assoc_lipid_',bi,'_stratified_interaction_unadjBMI_',Sys.Date(),".txt"))
  }
}
#run.test()
#"~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/get.assoc_unadjBMI_2026-10-01.R"
# anal.date = '2026-10-01'# rank across group
anal.date = '2026-11-03'# rank within group

int.files_unadjBMI=c(
  paste0('results/assoc_lipid_',brain_phenos[1],'_stratified_interaction_unadjBMI_',anal.date,".txt"),
  paste0('results/assoc_lipid_',brain_phenos[2],'_stratified_interaction_unadjBMI_',anal.date,".txt"),
  paste0('results/assoc_lipid_',brain_phenos[3],'_stratified_interaction_unadjBMI_',anal.date,".txt")
)

res.int.all_unadjBMI = rbind(
  fread(int.files_unadjBMI[1]),
  fread(int.files_unadjBMI[2]),
  fread(int.files_unadjBMI[3]))

res.int.all_unadjBMI = res.int.all_unadjBMI %>% 
  mutate(fdr.p.int = p.adjust(pval.int,method="fdr")) %>%
  mutate(fdr.p.abs.diff = p.adjust(abs.diff.pval,method="fdr")) %>%
  mutate(beta.abs.diff = abs(beta.ado)-abs(beta.adult)) # >0: if stronger in ados
res.int.all_unadjBMI %>% filter(fdr.p.int <0.05) %>% dim#566
res.int.all_unadjBMI %>% filter(fdr.p.abs.diff <0.05) %>% dim#291

res.int.all_unadjBMI %>% filter(fdr.p.int<0.05,beta.abs.diff<0) %>% dim #30 stronger in adults
res.int.all_unadjBMI %>% filter(fdr.p.int<0.05,beta.abs.diff>0) %>% dim #762 stronger in ados

# merge with res_WM1
head(res_WM1) 
dim(res_WM1)
dim(res.int.all_unadjBMI)/3
#[1] 830   9

# generate interaction test result table ----
# started: Mon, Nov 2 13:32
# long to wide:
res.int.all_unadjBMI = res.int.all_unadjBMI %>% 
  filter(metaboliteid != "GlycA") %>%
  mutate(fdr.p.int2 = p.adjust(pval.int)) %>%
  mutate(fdr.p.int2 = p.adjust( abs.diff))

interact_res_wide <- res.int.all_unadjBMI %>% reshape(
  direction = "wide", 
  idvar = "metaboliteid",         # Unique identifier for rows
  timevar = "bname",     # Column to split into headers
  v.names = setdiff(names(res.int.all_unadjBMI),c('bname','metaboliteid'))   # Column containing values
)
interact_res_wide %>% names
names(interact_res_wide )=
  str_replace(names(interact_res_wide ),"Volume.WM.adjAgeSexICV",'WM.Vol')
names(interact_res_wide )=
  str_replace(names(interact_res_wide ),"T1wSI.WM.adjAgeSex",'WM.T1wSI')
names(interact_res_wide )=
  str_replace(names(interact_res_wide ),"MTR.WM.adjAgeSex",'WM.MTR')

# same as res_WM1 ----
age_group_stratified_res = #same as res_WM1
  fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP_clean/manuscript/TabS1_2026-03-17.txt") %>% 
  filter(metaboliteid!="GlycA")

merged_assoc_res = age_group_stratified_res %>% 
  left_join(
    interact_res_wide, 
    join_by(metaboliteid)
  )

#*[NOTE: FDR correction for the interaction tests was applied within each phenotype]
selcol_names = unlist(
  str_split('metaboliteid
lipid class
beta_WM.Vol_ados
beta_WM.Vol_adults
pval.int.WM.Vol
fdr.p.int.WM.Vol
beta_WM.T1wSI_ados
beta_WM.T1wSI_adults
pval.int.WM.T1wSI
fdr.p.int.WM.T1wSI
beta_WM.MTR_ados
beta_WM.MTR_adults
pval.int.WM.MTR
fdr.p.int.WM.MTR',"\n"))

to_write = merged_assoc_res %>% 
  select(all_of(selcol_names))%>%
  filter(metaboliteid!='GlycA')

# group-specific tests of difference in age-group specific asso. estimates ----
## write output file 
to_write_sel_colnames=
  c("metaboliteid", "lipid class", "signif_ado",
    ## Volume
    "beta_WM.Vol_ados", "beta_WM.Vol_adults", 
    "pval_WM.Vol_diff", "fdr.p_WM.Vol_diff", 
    "pval.int.WM.Vol", "fdr.p.int.WM.Vol", 
    ## T1wSI
    "beta_WM.T1wSI_ados", "beta_WM.T1wSI_adults", 
    "pval_WM.T1wSI_diff", "fdr.p_WM.T1wSI_diff",
    "pval.int.WM.T1wSI", "fdr.p.int.WM.T1wSI",
    ## MTR
    "beta_WM.MTR_ados", "beta_WM.MTR_adults", 
    "pval_WM.MTR_diff", "fdr.p_WM.MTR_diff",
    "pval.int.WM.MTR", "fdr.p.int.WM.MTR"
  );to_write %>% 
  mutate(signif_ado = ifelse(metaboliteid%in%sig_lipids,TRUE,FALSE)) %>%
  relocate(metaboliteid,`lipid class`,signif_ado) %>% 
  left_join(
    age_group_stratified_res %>% 
      mutate(
        beta_WM.Vol_diff = beta_WM.Vol_ados - beta_WM.Vol_adults,
        se_WM.Vol_diff = sqrt(se_WM.Vol_ados^2 + se_WM.Vol_adults^2),
        z_WM.Vol_diff = (beta_WM.Vol_diff/se_WM.Vol_diff),
        pval_WM.Vol_diff = 2*pnorm(abs(z_WM.Vol_diff),lower.tail=F),
        ##
        beta_WM.T1wSI_diff = beta_WM.T1wSI_ados - beta_WM.T1wSI_adults,
        se_WM.T1wSI_diff = sqrt(se_WM.T1wSI_ados^2 + se_WM.T1wSI_adults^2),
        z_WM.T1wSI_diff = (beta_WM.T1wSI_diff/se_WM.T1wSI_diff),
        pval_WM.T1wSI_diff = 2*pnorm(abs(z_WM.T1wSI_diff),lower.tail=F),
        ##
        beta_WM.MTR_diff = beta_WM.MTR_ados - beta_WM.MTR_adults,
        se_WM.MTR_diff = sqrt(se_WM.MTR_ados^2 + se_WM.MTR_adults^2),
        z_WM.MTR_diff = (beta_WM.MTR_diff/se_WM.MTR_diff),
        pval_WM.MTR_diff = 2*pnorm(abs(z_WM.MTR_diff),lower.tail=F)
      ) %>%
      mutate(
        fdr.p_WM.Vol_diff = p.adjust(pval_WM.Vol_diff,method='fdr'),
        fdr.p_WM.T1wSI_diff = p.adjust(pval_WM.T1wSI_diff,method='fdr'),
        fdr.p_WM.MTR_diff = p.adjust(pval_WM.MTR_diff,method='fdr')
      ) %>% 
      select(metaboliteid,pval_WM.Vol_diff,fdr.p_WM.Vol_diff,
             pval_WM.T1wSI_diff,fdr.p_WM.T1wSI_diff,
             pval_WM.MTR_diff,fdr.p_WM.MTR_diff), join_by(metaboliteid) 
  ) %>% 
  select(all_of(to_write_sel_colnames)) %>%
  write_tsv(file=paste0('results/beta_diff_and_interaction_test_results_',Sys.Date(),'.txt'))

# compare interaction vs. heterogeneity test p-values ----
anal.date = "2026-11-03"; diff_int_res_fname = paste0('results/beta_diff_and_interaction_test_results_',anal.date,'.txt')
diff_int_res = fread(diff_int_res_fname)
p1=diff_int_res  %>% 
  ggplot(aes(x=-log10(pval_WM.Vol_diff),y=-log10(pval.int.WM.Vol))) + 
  geom_point() + 
  geom_abline(intercept = 0,slope = 1,col="darkblue")#+geom_smooth(method="lm")
p2=diff_int_res  %>% 
  ggplot(aes(x=-log10(pval_WM.T1wSI_diff),y=-log10(pval.int.WM.T1wSI))) + 
  geom_point() + 
  geom_abline(intercept = 0,slope = 1,col="darkblue")#+geom_smooth(method="lm")
p3=diff_int_res  %>% 
  ggplot(aes(x=-log10(pval_WM.MTR_diff),y=-log10(pval.int.WM.MTR))) + 
  geom_point() + 
  geom_abline(intercept = 0,slope = 1,col="darkblue")#+geom_smooth(method="lm")
plot.lims = c(0,max(-log10(diff_int_res$pval.int.WM.MTR)))
p1+p2+p3 +
  plot_layout(ncol=3,axes = 'collect') &
  coord_cartesian(xlim=plot.lims,ylim=plot.lims) &
  ggpubr::stat_cor(
    method = "spearman",
    label.x = plot.lims[1],
    label.y = plot.lims[2]*1.035,
    hjust = 0,
    vjust = 1
  )
## counts of signif/insignific lipids: (all) ----
rbind(
  diff_int_res %>% 
  mutate(
    diff.signif=ifelse(pval_WM.Vol_diff<0.05, 'p-diff<0.05',"p-diff>=0.05"),
    int.signif =ifelse(pval.int.WM.Vol<0.05,'p-int<0.05','p-int>=0.05')) %>%
  select(contains('.signif')) %>%
  table(),
  diff_int_res %>% 
    mutate(
      diff.signif=ifelse(pval_WM.T1wSI_diff<0.05, 'p-diff<0.05',"p-diff>=0.05"),
      int.signif =ifelse(pval.int.WM.T1wSI<0.05,'p-int<0.05','p-int>=0.05')) %>%
    select(contains('.signif')) %>%
    table(),
  diff_int_res %>% 
    mutate(
      diff.signif=ifelse(pval_WM.MTR_diff<0.05, 'p-diff<0.05',"p-diff>=0.05"),
      int.signif =ifelse(pval.int.WM.MTR<0.05,'p-int<0.05','p-int>=0.05')) %>%
    select(contains('.signif')) %>%
    table()
)
#              p-int<0.05 p-int>=0.05
# p-diff<0.05         136           4
# p-diff>=0.05         44         646
# p-diff<0.05         387          36
# p-diff>=0.05        107         300
# p-diff<0.05         464          16
# p-diff>=0.05         25         325

## counts of signif/insignific lipids: (525 lipids assoc. with all three in ados) ----
p1=diff_int_res %>% filter(metaboliteid%in% sig_lipids) %>%
  ggplot(aes(x=-log10(pval_WM.Vol_diff),y=-log10(pval.int.WM.Vol))) + 
  geom_point() + 
  geom_abline(intercept = 0,slope = 1,col="darkblue")#+geom_smooth(method="lm")
p2=diff_int_res %>% filter(metaboliteid%in% sig_lipids) %>%
  ggplot(aes(x=-log10(pval_WM.T1wSI_diff),y=-log10(pval.int.WM.T1wSI))) + 
  geom_point() + 
  geom_abline(intercept = 0,slope = 1,col="darkblue")#+geom_smooth(method="lm")
p3=diff_int_res %>% filter(metaboliteid%in% sig_lipids) %>%
  ggplot(aes(x=-log10(pval_WM.MTR_diff),y=-log10(pval.int.WM.MTR))) + 
  geom_point() + 
  geom_abline(intercept = 0,slope = 1,col="darkblue")#+geom_smooth(method="lm")
##
rbind(
  diff_int_res %>% filter(metaboliteid%in% sig_lipids) %>%
    mutate(
      diff.signif=ifelse(pval_WM.Vol_diff<0.05, 'p-diff<0.05',"p-diff>=0.05"),
      int.signif =ifelse(pval.int.WM.Vol<0.05,'p-int<0.05','p-int>=0.05')) %>%
    select(contains('.signif')) %>%
    table(),
  diff_int_res %>% filter(metaboliteid%in% sig_lipids) %>%
    mutate(
      diff.signif=ifelse(pval_WM.T1wSI_diff<0.05, 'p-diff<0.05',"p-diff>=0.05"),
      int.signif =ifelse(pval.int.WM.T1wSI<0.05,'p-int<0.05','p-int>=0.05')) %>%
    select(contains('.signif')) %>%
    table(),
  diff_int_res %>% filter(metaboliteid%in% sig_lipids) %>%
    mutate(
      diff.signif=ifelse(pval_WM.MTR_diff<0.05, 'p-diff<0.05',"p-diff>=0.05"),
      int.signif =ifelse(pval.int.WM.MTR<0.05,'p-int<0.05','p-int>=0.05')) %>%
    select(contains('.signif')) %>%
    table()
)
#              p-int<0.05 p-int>=0.05
# p-diff<0.05         110           1
# p-diff>=0.05         43         371
# p-diff<0.05         321          30
# p-diff>=0.05         87          87
# p-diff<0.05         391           8
# p-diff>=0.05         19         107
##
plot.lims = c(0,max(-log10(diff_int_res$pval.int.WM.MTR)))
p1+p2+p3 +
  plot_layout(ncol=3,axes = 'collect') &
  coord_cartesian(xlim=plot.lims,ylim=plot.lims) &
  ggpubr::stat_cor(
    method = "spearman",
    label.x = plot.lims[1],
    label.y = plot.lims[2]*1.035,
    hjust = 0,
    vjust = 1
  )
# numbers of lipids diff assoc with ados vs. adults ----
ind1 = to_write$fdr.p.int.WM.Vol <0.05;print(sum(ind1))#81
ind1 = ind1 | to_write$fdr.p.int.WM.T1wSI <0.05;print(sum(ind1))#327
ind1 = ind1 | to_write$fdr.p.int.WM.MTR <0.05;print(sum(ind1))#583

ind.sig = to_write$metaboliteid %in% sig_lipids;print(sum(ind.sig))
sum(ind.sig & ind1)#467

to_write %>% filter(ind.sig) %>% filter(fdr.p.int.WM.Vol<0.05) %>% dim#44
to_write %>% filter(ind.sig) %>% filter(fdr.p.int.WM.T1wSI<0.05) %>% dim#12
to_write %>% filter(ind.sig) %>% filter(fdr.p.int.WM.MTR<0.05) %>% dim#425
to_write%>% filter(ind1) %>% select(1:2)
table(to_write $fdr.p.int.WM.Vol<0.05,
      ind.sig,useNA = 'a')
  #       FALSE TRUE <NA>
  # FALSE   296  453    0
  # TRUE      9   72    0
  # <NA>      0    0    0
n1=res.int.all_unadjBMI %>% filter(bname=="Volume.WM.adjAgeSexICV") %>% filter(pval.int <0.05) %>% nrow#58
n2=res.int.all_unadjBMI %>% filter(bname=="T1wSI.WM.adjAgeSex") %>% filter(pval.int <0.05)%>% nrow#41
n3=res.int.all_unadjBMI %>% filter(bname=="MTR.WM.adjAgeSex") %>% filter(pval.int <0.05) %>% nrow#467
c(n1,n2,n3) %>% print;sum(n1,n2,n3) %>% print#141/122/564
# [1] 180 494 489
# [1] 1163

n1=res.int.all_unadjBMI %>% filter(bname=="Volume.WM.adjAgeSexICV") %>% 
  filter(fdr.p.int <0.05,metaboliteid%in% sig_lipids) %>% nrow#44
n2=res.int.all_unadjBMI %>% filter(bname=="T1wSI.WM.adjAgeSex") %>% 
  filter(fdr.p.int <0.05,metaboliteid%in% sig_lipids)%>% nrow#12
n3=res.int.all_unadjBMI %>% filter(bname=="MTR.WM.adjAgeSex") %>% 
  filter(fdr.p.int <0.05,metaboliteid%in% sig_lipids) %>% nrow#425
c(n1,n2,n3);sum(n1,n2,n3)
# [1]  72 248 350
# [1] 670

n1=res.int.all_unadjBMI %>% filter(bname=="Volume.WM.adjAgeSexICV") %>% 
  filter(fdr.p.int <0.05) %>% nrow#58
n2=res.int.all_unadjBMI %>% filter(bname=="T1wSI.WM.adjAgeSex") %>% 
  filter(fdr.p.int <0.05)%>% nrow#41
n3=res.int.all_unadjBMI %>% filter(bname=="MTR.WM.adjAgeSex") %>% 
  filter(fdr.p.int <0.05) %>% nrow#467
c(n1,n2,n3);sum(n1,n2,n3)
# [1]  81 306 405
# [1] 792

n1=res.int.all_unadjBMI %>% filter(bname=="Volume.WM.adjAgeSexICV") %>% 
  filter(fdr.p.int <0.05,metaboliteid%in% sig_lipids) %>% nrow#44
n2=res.int.all_unadjBMI %>% filter(bname=="T1wSI.WM.adjAgeSex") %>% 
  filter(fdr.p.int <0.05,metaboliteid%in% sig_lipids)%>% nrow#12
n3=res.int.all_unadjBMI %>% filter(bname=="MTR.WM.adjAgeSex") %>% 
  filter(fdr.p.int <0.05,metaboliteid%in% sig_lipids) %>% nrow#425
c(n1,n2,n3);sum(n1,n2,n3)#44+12+425
# [1]  72 248 350
# [1] 670

signif.int = res.int.all_unadjBMI %>% filter(fdr.p.int<0.05) %>% 
  pull(metaboliteid) %>% unique
signif.int %>% length#583
sum(sig_lipids %in% signif.int)#467


# scatter plot comparing estimates from ados vs. adults ---- 
anal.date = "2026-11-03"; diff_int_res_fname = paste0('results/beta_diff_and_interaction_test_results_',anal.date,'.txt')
diff_int_res = fread(diff_int_res_fname)

p1 = diff_int_res %>% #filter(metaboliteid %in% sig_lipids) %>% 
  mutate(signif = ifelse(pval_WM.Vol_diff<0.05, 'p.diff<0.05','p.diff>=0.05') ) %>%
  # mutate(signif = ifelse(pval.int.WM.MTR<0.05, 'p.int<0.05','p.int>=0.05') ) %>%
  ggplot(aes(x=beta_WM.Vol_adults,y=beta_WM.Vol_ados,fill = signif)) +
  geom_point(shape = 21) + geom_abline(intercept = 0, slope = 1) + 
  #geom_smooth(method='lm') + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0);p1

p2 = diff_int_res %>% #filter(metaboliteid %in% sig_lipids) %>% 
  mutate(signif = ifelse(pval_WM.T1wSI_diff<0.05, 'p.diff<0.05','p.diff>=0.05') ) %>%
  # mutate(signif = ifelse(pval.int.WM.MTR<0.05, 'p.int<0.05','p.int>=0.05') ) %>%
  ggplot(aes(x=beta_WM.T1wSI_adults,y=beta_WM.T1wSI_ados,fill = signif)) +
  geom_point(shape = 21) + geom_abline(intercept = 0, slope = 1) + 
  #geom_smooth(method='lm') + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0);p2

p3 = diff_int_res %>% #filter(metaboliteid %in% sig_lipids) %>% 
  mutate(signif = ifelse(pval_WM.MTR_diff<0.05, 'p.diff<0.05','p.diff>=0.05') ) %>%
  # mutate(signif = ifelse(pval.int.WM.MTR<0.05, 'p.int<0.05','p.int>=0.05') ) %>%
  ggplot(aes(x=beta_WM.MTR_adults,y=beta_WM.MTR_ados,fill = signif)) +
  geom_point(shape = 21) + geom_abline(intercept = 0, slope = 1) + 
  #geom_smooth(method='lm') + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0);p3

p1+p2+p3 + plot_layout(ncol=3,guides='collect')
