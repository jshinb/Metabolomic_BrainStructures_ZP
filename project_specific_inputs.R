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

# DATA INFORMATION ------------------------------------------------------------#
# working directory(wd), file names for the brain and non-brain datasets, and missing value code
#*[directory where the R scripts and data files are stored]
#* The following directory is a 'cleaned' directory
wd = "~/Library/CloudStorage/OneDrive-Personal/Lipidomic_BrainStructures_ZP"
dir.create(wd)

project.dir = c(data.dir = 'data',
output.dir = 'outputs',
script.dir = 'scripts',
result.dir = 'results')
names.project.dir = names(project.dir)

project.dir=file.path(wd,project.dir)
names(project.dir) = names.project.dir

lapply(project.dir,dir.create)
#------------------------------------------------------------------------------#

