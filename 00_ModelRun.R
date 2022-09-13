rm(list = ls())

# load libraries
require('rslurm')
require('sparta')
require('reshape2')

# Load data

visitData <- readRDS(file = "visitData_PlantBugs_210330.rds")

# regions dataframe


reg_data <- readRDS("../data/regions_non-crossed.rds")

# Set region aggregates
# region_aggs <- list(ENGLAND = names(reg_data)[!(names(reg_data) %in% c('SQ1_SQUARE','WALES','SCOTLAND','NORTHERN_IRELAND'))],
#                     GB = c(names(reg_data)[!(names(reg_data) %in% c('SQ1_SQUARE','WALES','SCOTLAND','NORTHERN_IRELAND'))], "WALES", "SCOTLAND"),
#                     UK = names(reg_data)[-1])

# if region aggregates already saved

region_aggs <- readRDS("../data/region_aggs_non-crossed.rds")

# Define function that loops through species
slurm_occDetFunc <- function(taxa_name){
  
  out <- occDetFunc(taxa_name = as.character(taxa_name),
                    occDetdata = visitData$occDetdata,
                    spp_vis = visitData$spp_vis,
                    write_results = TRUE,
                    n_chains = 3,
                    n_iterations = 32000,
                    burnin = 30000,
                    thinning = 6,
                    nyr = 2,
                    modeltype = c('ranwalk', 'halfcauchy', 'catlistlength'),
                    regional_codes = reg_data,
                    region_aggs = region_aggs,
                    return_data = FALSE,
                    seed = 123,
                    additional.parameters = "a",
                    allowSitesMultiRegions = TRUE,
                    rem_aggs_with_missing_regions=FALSE,
                    provenance = "2021_Francesca_Charlie_rerun")
  return(NULL)
}

# Create roster
pars <- data.frame(taxa_name = as.character(names(visitData[['spp_vis']])[-1]))

#####################################

# Create the job scipt and the R script needed to run the process on 
# lotus using slurm. Note: you can edit the templates used. These are
# found in the slurm folder in your R library (run '.Library' to find).
# You will need to add the command to load jaspy: module add jaspy
sjob <- slurm_apply(f = slurm_occDetFunc,
                    params = pars, 
                    jobname = 'PlantBugs',
                    nodes = nrow(pars), 
                    cpus_per_node = 1, 
                    submit = TRUE,
                    add_objects = c('visitData', 'reg_data', 'region_aggs'),
                    slurm_options = list(time = '167:59:00', 
                                         mem = 20000,
                                         partition = 'long-serial',
                                         error = '%a.err'))
