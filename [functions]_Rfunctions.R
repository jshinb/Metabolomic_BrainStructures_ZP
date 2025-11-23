# functions used for the projects

create_ggcorrplot = function(data, sel_colnames, new_colnames, ID_colname){
  
  df = data %>% select(all_of(sel_colnames))
  rownames.df = data[[ID_colname]]
  df = data.frame(df)
  rownames(df) = rownames.df 
  
  i=1; ind = !is.na(df[[1]])
  while(i<length(sel_colnames)){
    i <- i+1
    ind = ind | !is.na(df[[i]])
  }
  df = df[ind,]
  print(sum(ind))
  
  ## change the order of columns 
  names(df) = new_colnames
  cat(names(df),sep="\n")

  corr <- round(cor(df,use='p'), 2)
  p.mat <- cor_pmat(df)
  p.mat
  p_ggcorr = ggcorrplot(
    corr,
    # hc.order = TRUE,
    type = "lower",
    outline.color = "white",
    ggtheme = ggplot2::theme_bw,
    colors = c("#6D9EC1", "white", "#E46726"),
    # p.mat = p.mat,
    lab=TRUE,
  )
  ret = list(p_ggcorr=p_ggcorr, p.mat=p.mat, corr=corr)
  ret
}

inormal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

plot_heatmap = function(dat, axisNames)
{
  # define function to obtain lower triangle
  # get correlation matrix for both samples together
  cor = cor(dat, use = "p")
  cor = get_lower_tri(cor)
  # melt matrix
  melted = reshape2::melt(cor)
  # get rounded value
  melted$value_round = round(melted$value, digit = 2)
  melted$distance0 = abs(melted$value)
  
  # plot
  library(ggplot2)
  
  p = ggplot(data = melted)+
    geom_point(aes(x = Var1, y = Var2, shape = value, fill = value, size = distance0), 
               shape = 21, alpha = 0.7, colour = "white") +
    scale_fill_gradient2(low = "#82A0D8", high = "#8DDFCB", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab" ,name="Correlation", guide = "legend")+
    scale_size_continuous(range = c(1, 15), guide = "none")+
    geom_text(aes(Var1, Var2, label = value_round), color = "black", size = 4)+
    xlab("")+
    ylab("")+
    scale_x_discrete(labels = axisNames)+
    scale_y_discrete(labels = axisNames)+
    guides(fill = "none")+
    theme_bw()+
    theme(panel.border = element_blank(),
          axis.text.x = element_text(angle=-45,vjust = 0.5, hjust=0))
  
  return(p)
}

remove.outliers_grubbs <- function(dat, varnames){
  require(outliers)
  count.na0 <- count.na1 <- c()
  
  for(varname in varnames){
    x = dat[[varname]]
    count.na0 = c(count.na0,sum(is.na(x)))
    keep.going <- T
    while(keep.going){
      test <- grubbs.test(x,opposite = T)
      print(test)
      if(test$p.value<0.05){
        if(str_detect(test$alternative,"lowest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the lowest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------\n")
          x[which.min(x)] <- NA
        }else if(str_detect(test$alternative,"highest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the highest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.max(x)] <- NA
        }
      }
      test2 <- grubbs.test(x,opposite = F)
      print(test2)
      if(test2$p.value<0.05){
        if(str_detect(test2$alternative,"lowest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the lowest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.min(x)] <- NA
        }else if(str_detect(test2$alternative,"highest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the highest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.max(x)] <- NA
        }
      }
      
      keep.going=test$p.value<0.05 | test2$p.value<0.05
    }
    count.na1 = c(count.na1,sum(is.na(x)))
    dat[[varname]] <- x
  }
  counts.NA = data.frame(var=varnames,before=count.na0,after=count.na1)
  ret = list(counts.NA = counts.NA, clean_data=dat)
  ret
}

add_class2 = function(volcano_main_all,generation="ados"){
  if(generation=="ados"){
    metabodata_info = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/data/metabodata_info_ados.txt")
  }else{
    metabodata_info = fread("~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP/data/metabodata_info_adults.txt")
  }
  
  if(sum(names(volcano_main_all)=='class')==1){
    volcano_main_all = volcano_main_all %>% 
      dplyr::select(-class) 
  }
  volcano_main_all = volcano_main_all %>% 
    mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")) %>%
    left_join(
      metabodata_info %>%
        mutate(metabolite_id=str_remove(metabolite_id,"_nmol/mL")) %>% 
        dplyr::select(metabolite_id,class)
    )
  # class2.levels = c("CE", 
  #                   "TAG, DAG, or MAG", 
  #                   "PC, PE, or PI",
  #                   "LPC or LPE", 
  #                   "SM or ceramides",
  #                   "Nightingale")
  class2.levels = c("Cholesteryl esters",
                    "Acylglycerols",
                    "Phospholipids",
                    "Lysophospholipids",
                    "Sphingolipids",
                    "Inflammation/AA/metabolism")
  volcano_main_all[['class2']] <- class2.levels[1]
  volcano_main_all = volcano_main_all %>%
    mutate(class2 = ifelse(class %in% c('TAG','DAG','MAG'),class2.levels[2],class2)) %>%
    mutate(class2 = ifelse(class %in% c('PC','PE','PI'),class2.levels[3],class2)) %>%
    mutate(class2 = ifelse(class %in% c('LPC','LPE'),class2.levels[4],class2)) %>%
    mutate(class2 = ifelse(class %in% c('SM','CER',"DCER","HCER","LCER"),class2.levels[5],class2)) %>%
    mutate(class2 = ifelse(class %in% c('Amino acids','Glycolysis related metabolites','Inflammation','Ketone bodies'),
                           class2.levels[6],class2))
  table(volcano_main_all$class,volcano_main_all$class2,useNA = 'a')
  volcano_main_all = volcano_main_all %>%
    mutate(`lipid class` = factor(class2,levels=class2.levels)) 
  cols = c('#EF476F','#FFD166','#118AB2','#073B4C','#06D6A0','#8f00ff')
  
  class2.labels = paste(names(table(volcano_main_all$`lipid class`))," (n=",
                        table(volcano_main_all$`lipid class`),")",sep='')
  ret=(list(volcano_main_all=volcano_main_all,class2.labels=class2.labels,cols=cols))
  ret
}
