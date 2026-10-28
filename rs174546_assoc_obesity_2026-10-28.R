library(data.table)
library(dplyr)

gwascat = fread('~/Downloads/rs174546_associations_export.tsv')#accessed on Oct 28, 2026
#rs174546 - associated with lower odds of metabolomic syndrome

# riskAllele       "rs174546-?"             
# pValue           "5E-60"                  
# pValueAnnotation "-"                      
# riskFrequency    "NR"                     
# orValue          "-"                      
# beta             "0.0203186 unit decrease"
# ci               "[0.018-0.023]"          
# mappedGenes      "FADS2,FADS1"            
# traitName        "Metabolic syndrome"     
# efoTraits        "metabolic syndrome"     
# bgTraits         "-"                      
# accessionId      "GCST90444487"           
# locations        "11:61802358"            
# pubmedId         "39349817"               
# author           "Park S"         
