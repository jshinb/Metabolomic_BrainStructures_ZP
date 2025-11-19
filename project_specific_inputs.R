#************************** Project specific inputs **************************#
start_date = "2025-11-13"#*["YYYY-MM-DD"]
modify_date = "2025-11-13"#*["YYYY-MM-DD"]
analyst = "Jean Shin"#

# Progress
# 2025-11-13 
#* start cleaning up scripts (for better reproducibility)
#* addressing request of Zdenka:
#* correlation plots 
#** (sex-specific: correlation among adjisted and unadjusted brain phenotypes)

# set up working and sub-directories ----
# working directory(wd), file names for the brain and non-brain datasets, and missing value code
#*[directory where the R scripts and data files are stored]
#* The following directory is a 'cleaned' directory
wd = "~/Library/CloudStorage/OneDrive-Personal/Metabolomic_BrainStructures_ZP"
dir.create(wd)

# defile options ----
opt = list(
  data.dir = file.path(wd,'data'),
output.dir = file.path(wd,'outputs'),
script.dir = file.path(wd,'scripts'),
result.dir = file.path(wd,'results'),
group_name = 'SYS_ados'
)

# load packages ----
source("~/Documents/scripts/Metabolomic_BrainStructures_ZP/install_libs.R")
source("~/Documents/scripts/Metabolomic_BrainStructures_ZP/[functions]_Rfunctions.R")

