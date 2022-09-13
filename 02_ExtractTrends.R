################################################################
## Extract and plot trends in different arable land categories
## Authors: Rob Cooke & Francesca Mancini
## Date created: 2021-04-13
################################################################

# load the libraries
library(dplyr)
# Sys.setenv(JAGS_HOME="/usr/lib/JAGS")
# install.packages(c("rjags", "R2jags"))
library(rjags)
library(R2jags)
# remotes::install_github("https://github.com/03rcooke/wrappeR", ref = "main", force = TRUE)
library(wrappeR) # wrappeR: multi-species indicators
library(tidyr)
library(ggplot2)
library(HDInterval)
library(cowplot)
library(VennDiagram)
library(gridExtra)



# global parameters
startyear <- 1990
endyear <- 2019

# First we extract the occupancy data 
# for all taxonomic groups we are interested in ----

# Taxonomic group(s)
taxa <- c("Ladybirds", "Bees", "Carabids", "Spiders", 
          "PlantBugs", "Hoverflies"
          )

ntax <- length(taxa)

regions <- c("high", "no_agri", "low")

vers <- c("2021_Francesca", "2021_Francesca_bwars_rerun", "2021_Francesca_Charlie_rerun",
          "2021_Francesca_marlog_rerun", "2021_Francesca_Charlie_rerun", 
          "2021_Francesca"
          )

nregions <- length(regions)

# create roster - run for multiple groups & multiple regions
# ignore warning

# dir.create("/data/outputs/filtered")

roster <- createRoster(index = 1:(nregions * ntax),
                       modPath = "/data-s3/occmods/",
                       metaPath = "/data-s3/metadata/",
                       ver = rep(vers, each = nregions),
                       indicator = "all",
                       region = rep(regions, times = ntax),
                       nSamps = 999,
                       minObs = 1,
                       scaleObs = "region",
                       write = TRUE,
                       # create a folder called filtered in your directory
                       outPath = "/data/outputs/filtered/",
                       group = rep(taxa, each = nregions),
                       t0 = startyear,
                       tn = endyear)

# filtered data
filt_tax <- lapply(roster, wrappeR::applySamp, parallel = FALSE)



# model threshold quality ----
# This is the rules-of-thumb application 
# from Pocock et al - Rapid assessment of the 
# suitability of multi-species citizen science datasets 
# for occupancy trend analysis
# the function mod_qual also extracts convergence metrics

## observation metadata

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

allup <- function(x) {
  ot <- paste(toupper(stringr::word(x, 1)), sub(".*? ", "", x))
  ot
}

# load_rdata function
# loads an RData file, and assigns it to an object name
load_rdata <- function(fileName) {
  if(grepl("rdata$", fileName)){
  load(fileName)
  } 
  if(grepl("rds$", fileName)){
    out <- readRDS(fileName)
  }
  get(ls()[ls() != "fileName"])
}

# extract model quality metadata for each species
meta_high <- lapply(taxa, function(x) {
  
  load_rdata(paste0("/data/outputs/filtered/", x, "_all_high_samp.rdata")) %>% 
    .[[2]] %>% 
    dplyr::mutate(tax_grp = x)
  
}) %>% 
  dplyr::bind_rows(.) 

# ants, bees, carabids, wasps
meta_high <- meta_high  %>% 
  dplyr::mutate(species = firstup(Species_r_high)) %>% 
  dplyr::mutate(species = firstup(ifelse(tax_grp == "Spiders", 
                                 gsub("[[:punct:]]", "", meta_high$Species_r_high),
                                 species)))
#   # the ants are really inconsistent
#   dplyr::mutate(species = ifelse(tax_grp == "Ants", allup(Species_r_pa), species)) %>% 
#   dplyr::mutate(species = ifelse(tax_grp == "Ants" & Species_r_pa %in% c("tetramorium atratulum", "hypoponera punctatissima agg"), firstup(Species_r_pa), species))


meta_low <- lapply(taxa, function(x) {
  
  load_rdata(paste0("/data/outputs/filtered/", x, "_all_low_samp.rdata")) %>% 
    .[[2]] %>% 
    dplyr::mutate(tax_grp = x)
  
}) %>% 
  dplyr::bind_rows(.) 

# ants, bees, carabids, wasps
meta_low <- meta_low  %>% 
  dplyr::mutate(species = firstup(Species_r_low)) %>% 
  dplyr::mutate(species = firstup(ifelse(tax_grp == "Spiders", 
                                         gsub("[[:punct:]]", "", meta_low$Species_r_low),
                                         species)))


meta_no_agri <- lapply(taxa, function(x) {
  
  load_rdata(paste0("/data/outputs/filtered/", x, "_all_no_agri_samp.rdata")) %>% 
    .[[2]] %>% 
    dplyr::mutate(tax_grp = x)
  
}) %>% 
  dplyr::bind_rows(.) 

# ants, bees, carabids, wasps
meta_no_agri <- meta_no_agri  %>% 
  dplyr::mutate(species = firstup(Species_r_no_agri)) %>% 
  dplyr::mutate(species = firstup(ifelse(tax_grp == "Spiders", 
                                         gsub("[[:punct:]]", "", meta_no_agri$Species_r_no_agri),
                                         species)))


meta <- full_join(meta_high, meta_low, by = c("species", "tax_grp")) %>%
  full_join(meta_no_agri, by = c("species", "tax_grp"))

meta <- meta[-c(764,765),]

# meta_simple <- unique(meta[c("species", "tax_grp")])

## model quality

mod_qual <- function(roster, df) {
  
  # collapse roster into dataframe
  roster_df <- do.call(rbind, lapply(roster, as.data.frame))
  
  # extract species metrics for each species per group
  qual_df <- lapply(1:nrow(df), function(x) {
    
    # get metadata for species
    sp_meta <- df[x,]
    
    # subset roster for species
    sp_roster <- subset(roster_df, group == sp_meta$tax_grp)[1,]
    
    # identify model filename for species
    
    modFilePath <- file.path(sp_roster$modPath, sp_roster$group, 
                             "occmod_outputs", sp_roster$ver)
    modFiles_sp <- list.files(modFilePath, 
                              pattern = paste0("^", sp_meta$species,"(_|\\.)+"))
    
    modFiles_rdata <- modFiles_sp[grep(pattern = ".rdata$", modFiles_sp)]
    modFiles_rds <- modFiles_sp[grep(pattern = ".rds$", modFiles_sp)]
    
    if (length(modFiles_rdata) == 0 & length(modFiles_rds) == 
        0) 
      stop("Model files must be either .rds or .rdata")
    else {
      if (length(modFiles_rdata) == 0) {
        filetype <- "rds"
        modFiles <- modFiles_rds
      }
      else {
        filetype <- "rdata"
        modFiles <- modFiles_rdata
      }
    }
    
    if(sp_roster$group %in% c("Ladybirds", "Hoverflies")) {
      list_of_file_names <- gsub(pattern = paste0("\\.", filetype), 
                        repl = "", modFiles)
      
      list_of_file_names <- gsub("_[[:digit:]]{1}$", "", 
                                 list_of_file_names)
      iterations <- regmatches(list_of_file_names, regexpr("[[:digit:]]+$", 
                                                           list_of_file_names))
      first_iter <- min(as.numeric(iterations))
      last_iter <- max(as.numeric(iterations))
      
      
      modFiles_first <- modFiles[which(grepl(paste0("_", first_iter, "_")
                                             , modFiles)==TRUE)]
      modFiles <- modFiles[which(grepl(paste0("_", last_iter, "_"), 
                                       modFiles)==TRUE)]
      
      mod_first <- load_rdata(paste0(modFilePath, "/",
                               modFiles_first))
      first_year_index <- sp_roster$t0 - mod_first$min_year + 1
      last_year_index <- length(mod_first$min_year:mod_first$max_year)
      
      mod <- load_rdata(paste0(modFilePath, "/",
                               modFiles))$out
    
      } else {
      mod <- load_rdata(paste0(modFilePath, "/",
                               modFiles))
      
      # extract species metrics from model
      first_year_index <- sp_roster$t0 - mod$min_year + 1
      last_year_index <- length(mod$min_year:mod$max_year)
      }
    
    # extract species metrics from model
    first_last_year_converged <- sum(mod$BUGSoutput$summary[c(paste0("psi.fs.r_high[", 
                                                                     first_year_index,
                                                                     "]"), 
                                                          paste0("psi.fs.r_low[", 
                                                                 first_year_index,
                                                                 "]"), 
                                                          paste0("psi.fs.r_no_agri[", 
                                                                 first_year_index,
                                                                 "]"),
                                                          paste0("psi.fs.r_high[",
                                                                 last_year_index, "]"),
                                                          paste0("psi.fs.r_low[", 
                                                                 last_year_index, "]"),
                                                          paste0("psi.fs.r_no_agri[",
                                                                 last_year_index, "]")),
                                                        "Rhat"] > 1.1) ==0
    
    years_conv_high <-
      length(which(mod$BUGSoutput$summary[
        paste0("psi.fs.r_high[", first_year_index:last_year_index, "]"),
        "Rhat"] < 1.1)) 
    
    years_conv_low <-
      length(which(mod$BUGSoutput$summary[ 
        paste0("psi.fs.r_low[", first_year_index:last_year_index, "]"),
        "Rhat"] < 1.1))
    
    years_conv_no_agri <-
      length(which(mod$BUGSoutput$summary[
        paste0("psi.fs.r_no_agri[", first_year_index:last_year_index, "]"),
        "Rhat"] < 1.1))
    
    # prop_year_conv <-
    #   length(which(mod$BUGSoutput$summary[grep(paste(
    #     c("psi.fs.r_high", "psi.fs.r_low", "psi.fs.r_no_agri"),
    #     collapse = "|"
    #   ),
    #   rownames(mod$BUGSoutput$summary),
    #   value = TRUE),
    #   "Rhat"] <= 1.1)) / length(mod$BUGSoutput$summary[grep(paste(
    #     c("psi.fs.r_high", "psi.fs.r_low", "psi.fs.r_no_agri"),
    #     collapse = "|"
    #   ),
    #   rownames(mod$BUGSoutput$summary),
    #   value = TRUE), "Rhat"])
    # 
    mean_occ_prec <- mean(1/(mod$BUGSoutput$summary[
      c(paste0("psi.fs.r_high[", first_year_index:last_year_index, "]"), 
        paste0("psi.fs.r_low[", first_year_index:last_year_index, "]"),
        paste0("psi.fs.r_no_agri[", first_year_index:last_year_index, "]")),
    "sd"]))

    if(sp_roster$group %in% c("Ladybirds", "Hoverflies")) {
      out <- as.data.frame(attr(mod_first, "metadata")$analysis$spp_Metrics) %>%
        dplyr::mutate(species = sp_meta$species, 
                      tax_grp = sp_meta$tax_grp,
                      first_last_converged = first_last_year_converged,
                      years_converged_high = years_conv_high,
                      years_converged_low = years_conv_low,
                      years_converged_no_agri = years_conv_no_agri,
                      precision = mean_occ_prec)
    }
    else{
    out <- as.data.frame(attr(mod, "metadata")$analysis$spp_Metrics) %>%
      dplyr::mutate(species = sp_meta$species, 
                    tax_grp = sp_meta$tax_grp,
                    first_last_converged = first_last_year_converged,
                    years_converged_high = years_conv_high,
                    years_converged_low = years_conv_low,
                    years_converged_no_agri = years_conv_no_agri,
                    precision = mean_occ_prec)
    }
    
    print(paste("Complete:", sp_meta$species))
    
    return(out)
    
  }) %>% 
    # bind across species
    dplyr::bind_rows(.)
  
  return(qual_df)
  
}

mod_met <- mod_qual(roster = roster, df = meta)

saveRDS(mod_met, "/data/outputs/mod_met.rds")

# plot convergence stats

ggplot(data = mod_met %>%
         select(tax_grp, species, 
                years_converged_high,
                years_converged_low,
                years_converged_no_agri) %>%
         gather(key = "region", value = "years_converged",
                years_converged_high, years_converged_low, years_converged_no_agri)) +
  geom_histogram(aes(x = years_converged, fill = region),
                 binwidth=1, position="dodge") +
  facet_wrap(~tax_grp)


first_last_year_conv_summary <- mod_met %>%
  group_by(tax_grp, first_last_converged) %>%
  tally()

ggplot(data = mod_met) +
  geom_histogram(aes(x = log(precision)),
                 binwidth = 1) +
  # xlim(0,500) +
  facet_wrap(~tax_grp) +
  geom_vline(aes(xintercept = log(87)), linetype = "dashed")



pass_keep <- mod_met %>% 
  # EqualWt
  dplyr::mutate(pass = ifelse(prop_abs >= 0.990, 
                              P90 >= 3.1, 
                              P90 >= 6.7)) 
# %>%
# dplyr::mutate(pass = ifelse(converged >= 0.7,
#                             TRUE, FALSE))

pass_keep_summary <- pass_keep %>%
  group_by(tax_grp) %>%
  count(pass#, converged
        )

# only those that converged in first and last year are kept
pass_first_last_converged_keep <- mod_met %>% 
  # EqualWt
  dplyr::mutate(pass = ifelse(prop_abs >= 0.990, 
                              P90 >= 3.1, 
                              P90 >= 6.7)) %>%
  dplyr::mutate(converged = ifelse(first_last_converged == TRUE,
                              TRUE, FALSE))

pass__first_last_converged_keep_summary <- pass_first_last_converged_keep %>%
  group_by(tax_grp) %>%
  count(pass, converged)

# only those that have precision more than 10 are kept
pass_precision_keep <- mod_met %>% 
  # EqualWt
  dplyr::mutate(pass = ifelse(prop_abs >= 0.990, 
                              P90 >= 3.1, 
                              P90 >= 6.7)) %>%
  dplyr::mutate(converged = ifelse(precision >= 87,
                              TRUE, FALSE)) 

pass_precision_keep_summary <- pass_precision_keep %>%
  group_by(tax_grp) %>%
  count(pass, converged)

saveRDS(pass_keep, "/data/outputs/pass_keep.rds")
saveRDS(pass_first_last_converged_keep, "/data/outputs/pass_first_last_converged_keep.rds")
saveRDS(pass_precision_keep, "/data/outputs/pass_precision_keep.rds")

pass_keep <- readRDS("/data/outputs/pass_keep.rds") %>%
  dplyr::filter(pass == TRUE)

pass_first_last_converged_keep <- readRDS("/data/outputs/pass_first_last_converged_keep.rds") %>%
  dplyr::filter(pass == TRUE & converged == TRUE)

# pass_precision_keep <- readRDS("/data/outputs/pass_precision_keep.rds") %>%
#   dplyr::filter(pass == TRUE & converged == TRUE)


pass_keep <- tolower(pass_keep$species)
pass_first_last_converged_keep <- tolower(pass_first_last_converged_keep$species)
# pass_precision_keep <- tolower(pass_precision_keep$species)

## occupancy estimates ----
# Here we read in prepared occupancy data

# read in filtered data
occ_read <- function(x, reg) {
  
  load_rdata(paste0("/data/outputs/filtered/", x, "_all_", reg, "_samp.rdata")) %>% 
    .[[1]] %>% 
    dplyr::mutate(tax_grp = x)
  
}

# focal years 1990 to 2019
foc_yrs <- paste0("year_", startyear:endyear)

#   ## high coverage of arable land
#   
  occ_high <- lapply(taxa, function(x) occ_read(x = x, reg = "high")) %>%
    dplyr::bind_rows(.) %>% 
    # trim to focal years 1990 to 2019
    dplyr::select(foc_yrs, iteration, species, tax_grp) %>% 
    # add region aggregate name
    dplyr::mutate(agri = "high")
  
  spp_high <- data.frame(spp = unique(occ_high$species))
  
  ## low coverage of arable land
  
  occ_low <- lapply(taxa, function(x) occ_read(x = x, reg = "low")) %>% 
    dplyr::bind_rows(.) %>% 
    # trim to focal years 1990 to 2019
    dplyr::select(foc_yrs, iteration, species, tax_grp)  %>% 
    # add region aggregate name
    dplyr::mutate(agri = "low")
  
  spp_low <- data.frame(spp = unique(occ_low$species))
  
  ## no arable land
  occ_no_agri <- lapply(taxa, function(x) occ_read(x = x, reg = "no_agri")) %>% 
    dplyr::bind_rows(.) %>% 
    # trim to focal years 1990 to 2019
    dplyr::select(foc_yrs, iteration, species, tax_grp)  %>% 
    # add region aggregate name
    dplyr::mutate(agri = "no_agri")
  
  spp_no_agri <- data.frame(spp = unique(occ_no_agri$species))
  

  # remove species that didn't pass model quality threshold
  occ_high_pass <- dplyr::filter(occ_high, species %in% pass_keep & species != "harmonia axyridis")
  occ_high_pass_converged <- dplyr::filter(occ_high, species %in% pass_first_last_converged_keep)
  # occ_high_pass_precision <- dplyr::filter(occ_high, species %in% pass_precision_keep)
  occ_low_pass <- dplyr::filter(occ_low, species %in% pass_keep & species != "harmonia axyridis")
  occ_low_pass_converged <- dplyr::filter(occ_low, species %in% pass_first_last_converged_keep)
  # occ_low_pass_precision <- dplyr::filter(occ_low, species %in% pass_precision_keep)
  occ_no_agri_pass <- dplyr::filter(occ_no_agri, species %in% pass_keep & species != "harmonia axyridis")
  occ_no_agri_pass_converged <- dplyr::filter(occ_no_agri, species %in% pass_first_last_converged_keep)
  # occ_no_agri_pass_precision <- dplyr::filter(occ_no_agri, species %in% pass_precision_keep)

  spp_high_pass <- data.frame(spp = unique(occ_high_pass$species))
  spp_low_pass <- data.frame(spp = unique(occ_low_pass$species))
  spp_no_agri_pass <- data.frame(spp = unique(occ_no_agri_pass$species))
  
  

occ_out_pass <- list(occ_high_pass, occ_low_pass, occ_no_agri_pass)
occ_out_pass_converged <- list(occ_high_pass_converged,
                               occ_low_pass_converged,
                               occ_no_agri_pass_converged)
# occ_out_pass_precision <- list(occ_high_pass_precision, 
#                                occ_low_pass_precision, 
#                                occ_no_agri_pass_precision)



## trend lines by taxa and coverage of arable land

# function to calculate geometric mean
geomean <- function(x) exp(mean(log(x)))

# collapse occ_out list into a dataframe 

occ_out_pass_df <- do.call(rbind, occ_out_pass)
occ_out_pass_converged_df <- do.call(rbind, occ_out_pass_converged)
# occ_out_pass_precision_df <- do.call(rbind, occ_out_pass_precision)


# remove old objects to gain some memory

rm(occ_high, occ_high_pass, occ_high_pass_converged, occ_high_pass_precision,
   occ_low, occ_low_pass, occ_low_pass_converged, occ_low_pass_precision,
   occ_no_agri, occ_no_agri_pass, occ_no_agri_pass_converged, occ_no_agri_pass_precision,
   occ_out_pass, occ_out_pass_converged, occ_out_pass_precision)

gc() 

# summarise_occ takes a subset of the occ_out dataframe
# by taxon group, summarises occupancy across species by iteration
# transforms the dataframe into a long format
# then summarises (mean, median and geometric mean)
# across iterations and calculates CI by taxa and arable land cover

# get number of sites for each taxa

names(roster) <- rep(taxa, each = 3)


summarise_occ <- function(df, taxa, roster) {
  
  roster_tmp <- roster[taxa][[1]]
  
  filepath <- paste0(roster_tmp$modPath, taxa, 
                     "/occmod_outputs/", roster_tmp$ver)
  

  if(taxa %in% c("Ladybirds", "Hoverflies")) {
    
    modFiles <- list.files(filepath)
    
    list_of_file_names <- gsub(pattern = "\\.rdata", 
                               repl = "", modFiles)
    
    list_of_file_names <- gsub("_[[:digit:]]{1}$", "", 
                               list_of_file_names)
    iterations <- regmatches(list_of_file_names, regexpr("[[:digit:]]+$", 
                                                         list_of_file_names))
    first_iter <- min(as.numeric(iterations))
    last_iter <- max(as.numeric(iterations))
    
    
    modFiles_first <- modFiles[which(grepl(paste0("_", first_iter, "_")
                                           , modFiles)==TRUE)][1]
    # modFiles <- modFiles[which(grepl(paste0("_", last_iter, "_"), 
    #                                  modFiles)==TRUE)]
    
    out <- load_rdata(paste0(filepath, "/",modFiles_first))
    nsites <- out$nsites
    # first_year_index <- sp_roster$t0 - mod_first$min_year + 1
    # last_year_index <- length(mod_first$min_year:mod_first$max_year)
    # 
    # mod <- load_rdata(paste0(modFilePath, "/",
    #                          modFiles))$out
    
  } else {
  
  out <- load_rdata(
    list.files(paste0(roster_tmp$modPath, taxa, 
                      "/occmod_outputs/", roster_tmp$ver),
               full.names = TRUE)[1]
  )
  nsites <- out$nsites
  }
  

  df %>%
    filter(tax_grp == taxa) %>%
    mutate_at(
      vars(starts_with("year_")), 
           list(~case_when(
             . == 0 ~ 1/(nsites*5),
             TRUE ~ .))) %>%
    # either max from Nick's code or change 0s to 1/(n_sites * 5)
    dplyr::group_by(agri, iteration) %>%
    dplyr::summarise_at(
      vars(starts_with("year_")),
      geomean) %>%    # geomean
    ungroup() %>%
    gather(year_1990:year_2019, key = "year", value = "occupancy") %>%
    mutate(year = as.numeric(gsub("year_", "", year)))%>%
    group_by(agri, year) %>%
    summarise(
      mean = mean(occupancy),   # this one
      median = median(occupancy),
      # gmean = geomean(occupancy),
      low_ci = HDInterval::hdi(occupancy,
                              credMass = 0.95)[[1]],
      upp_ci = HDInterval::hdi(occupancy,
                              credMass = 0.95)[[2]]) %>%
    mutate(tax_grp = taxa)
}


# apply the function to all taxa and bind together

occ_pass_summary <- lapply(unique(occ_out_pass_df$tax_grp), 
                      summarise_occ, 
                      df = occ_out_pass_df,
                      roster = roster) %>%
  dplyr::bind_rows(.)

# saveRDS(occ_pass_summary, "/data/outputs/occ_pass_summary.rds")
saveRDS(occ_pass_summary, "/data/outputs/occ_pass_summary_no_harlequin.rds")


occ_pass_converged_summary <- lapply(unique(occ_out_pass_converged_df$tax_grp),
                           summarise_occ,
                           df = occ_out_pass_converged_df,
                           roster = roster) %>%
  dplyr::bind_rows(.)
saveRDS(occ_pass_converged_summary, "/data/outputs/occ_pass_converged_summary.rds")

# occ_pass_precision_summary <- lapply(unique(occ_out_pass_precision_df$tax_grp), 
#                            summarise_occ, 
#                            df = occ_out_pass_precision_df,
#                            roster = roster) %>%
#   dplyr::bind_rows(.)


# then plot

# n_species_annotation <- data.frame(n_species = c(bee_spp$species,
#                                                  carabids_spp$species,
#                                                  ladybirds_spp$species,
#                                                  plantbug_spp$species),
#                                    tax_grp = c("Bees", "Carabids",
#                                                "Ladybirds", "PlantBugs"))

occ_pass_by_agri <- ggplot(
  data = occ_pass_summary) +
  geom_line(aes(x = as.numeric(year), y = mean,
                colour = agri)) +
  geom_ribbon(aes(
    x = as.numeric(year),
    ymin = low_ci,
    ymax = upp_ci,
    fill = agri
  ),
  alpha = 0.2,
  show.legend = FALSE) +
  facet_wrap( ~ tax_grp,
              ncol = 2,
              scales = "free_x") +
  # colour
  scale_colour_manual(
    values = c("firebrick1", "goldenrod1", "dodgerblue"),
    labels = c("High", "Low", "No arable"),
    name = "Coverage of\narable land"
  ) +
  scale_fill_manual(
    values = c("firebrick1", "goldenrod1", "dodgerblue"),
    labels = c("High", "Low", "No arable")
  ) +
  # force lower limit to 0
  expand_limits(y = 0) +
  # axis labels
  labs(x = "Year", y = "Average occupancy") +
  # geom_text(data = n_species_annotation,
  #           aes(x = 2005, y = 0.7,
  #           label = paste0("n species = ", n_species))) +
  theme_light() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)#,
        #legend.position = c(0.75, 0.1) # c(0,0) bottom left, c(1,1) top-right
        ) +
  ggtitle("Pass Rule-of-thumb")

# occ_pass_converged_by_agri <- ggplot(
#   data = occ_pass_converged_summary) +
#   geom_line(aes(x = as.numeric(year), y = mean,
#                 colour = agri)) +
#   geom_ribbon(aes(
#     x = as.numeric(year),
#     ymin = low_ci,
#     ymax = upp_ci,
#     fill = agri
#   ),
#   alpha = 0.2,
#   show.legend = FALSE) +
#   facet_wrap( ~ tax_grp,
#               ncol = 2,
#               scales = "free_x") +
#   # colour
#   scale_colour_manual(
#     values = c("firebrick1", "goldenrod1", "dodgerblue"),
#     labels = c("High", "Low", "No arable"),
#     name = "Coverage of\narable land"
#   ) +
#   scale_fill_manual(
#     values = c("firebrick1", "goldenrod1", "dodgerblue"),
#     labels = c("High", "Low", "No arable")
#   ) +
#   # force lower limit to 0
#   expand_limits(y = 0) +
#   # axis labels
#   labs(x = "Year", y = "Average occupancy") +
#   # geom_text(data = n_species_annotation,
#   #           aes(x = 2005, y = 0.7,
#   #           label = paste0("n species = ", n_species))) +
#   theme_light() +
#   theme(text = element_text(size = 11),
#         axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)#,
#         #legend.position = c(0.75, 0.1) # c(0,0) bottom left, c(1,1) top-right
#   )+
#   ggtitle("Pass Rule-of-thumb and first/last year converged")
# 
# occ_pass_precision_by_agri <- ggplot(
#   data = occ_pass_precision_summary) +
#   geom_line(aes(x = as.numeric(year), y = mean,
#                 colour = agri)) +
#   geom_ribbon(aes(
#     x = as.numeric(year),
#     ymin = low_ci,
#     ymax = upp_ci,
#     fill = agri
#   ),
#   alpha = 0.2,
#   show.legend = FALSE) +
#   facet_wrap( ~ tax_grp,
#               ncol = 2,
#               scales = "free_x") +
#   # colour
#   scale_colour_manual(
#     values = c("firebrick1", "goldenrod1", "dodgerblue"),
#     labels = c("High", "Low", "No arable"),
#     name = "Coverage of\narable land"
#   ) +
#   scale_fill_manual(
#     values = c("firebrick1", "goldenrod1", "dodgerblue"),
#     labels = c("High", "Low", "No arable")
#   ) +
#   # force lower limit to 0
#   expand_limits(y = 0) +
#   # axis labels
#   labs(x = "Year", y = "Average occupancy") +
#   # geom_text(data = n_species_annotation,
#   #           aes(x = 2005, y = 0.7,
#   #           label = paste0("n species = ", n_species))) +
#   theme_light() +
#   theme(text = element_text(size = 11),
#         axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)#,
#         #legend.position = c(0.75, 0.1) # c(0,0) bottom left, c(1,1) top-right
#   )+
#   ggtitle("Pass Rule-of-thumb and precision threshold")

# # save: trend plot
png("/data/plots/all_taxa_pass_trends.png", 
    height = 16, width = 14, units = "cm", res = 320)

occ_pass_by_agri

dev.off()

# png("/data/plots/all_taxa_pass_converged_trends.png", 
#     height = 16, width = 14, units = "cm", res = 320)
# 
# occ_pass_converged_by_agri
# 
# dev.off()
# 
# png("/data/plots/all_taxa_pass_precision_trends.png", 
#     height = 16, width = 14, units = "cm", res = 320)
# 
# occ_pass_precision_by_agri
# 
# dev.off()

# png("/data/plots/all_taxa_pass_precision87_trends.png", 
#     height = 16, width = 14, units = "cm", res = 320)
# 
# occ_pass_precision_by_agri
# 
# dev.off()

## Trend change ----
# Here we calculate annual growth rate between the first and last year

# function to calculate annual growth rate per species
trend_change <- function(occ_df, taxa) {
  
  out <- lapply(taxa, function(x) {
    
    roster_tmp <- roster[taxa][[1]]
    
    mod <- load_rdata(
      list.files(paste0(roster_tmp$modPath, taxa, 
                        "/occmod_outputs/", roster_tmp$ver),
                 full.names = TRUE)[1]
    )
    nsites <- mod$nsites
    
    
    spp_df <- occ_df %>% 
      dplyr::filter(tax_grp == x) %>% 
      # only non-na columns
      dplyr::select(which(colMeans(is.na(.)) < 1)) %>% 
      # select first and last year
      dplyr::select(1, (ncol(.) - 4), iteration, species, tax_grp, agri) %>%
      # record first year
      dplyr::mutate(first_year = colnames(.)[1]) %>% 
      dplyr::mutate(first_year = as.numeric(gsub("year_", "", first_year))) %>% 
      # record last year
      dplyr::mutate(last_year = colnames(.)[2]) %>% 
      dplyr::mutate(last_year = as.numeric(gsub("year_", "", last_year))) %>% 
      # set column names
      setNames(c("occ_first", "occ_last", "iteration", "species", "tax_grp", "agri", "first_year", "last_year")) %>% 
      # number of years
      dplyr::mutate(nyr = last_year - first_year) %>% 
      # replace zero with very small number - can't divide by zero
      dplyr::mutate(occ_first = ifelse(occ_first == 0, 1/(nsites*5), occ_first)) %>%
      # 1/(n_sites * 5)
      # annual growth rate - see occurrenceChange()
      dplyr::mutate(change = (((occ_last / occ_first) ^ (1 / nyr)) - 1) * 100)
    
  }) %>% 
    # bind across taxonomic groups
    dplyr::bind_rows(.)
  
}

# trend_out <- lapply(1:length(occ_out), function(p) {
#   
#   occ_out_p <- occ_out[[p]]
  
  trend_pass_high <- trend_change(occ_df = filter(occ_out_pass_df, agri == "high"), 
                                  taxa = taxa)
  # trend_pass_converged_high <- trend_change(occ_df = filter(occ_out_pass_converged_df, 
  #                                                           agri == "high"), 
  #                                           taxa = taxa)
  # trend_pass_precision_high <- trend_change(occ_df = filter(occ_out_pass_precision_df, 
  #                                                           agri == "high"), 
  #                                           taxa = taxa)
  trend_pass_low <- trend_change(occ_df = filter(occ_out_pass_df, agri == "low"), 
                                 taxa = taxa)
  # trend_pass_converged_low <- trend_change(occ_df = filter(occ_out_pass_converged_df, 
  #                                                          agri == "low"), 
  #                                          taxa = taxa)
  # trend_pass_precision_low <- trend_change(occ_df = filter(occ_out_pass_precision_df, 
  #                                                          agri == "low"), 
  #                                          taxa = taxa)
  trend_pass_no_agri <- trend_change(occ_df = filter(occ_out_pass_df, agri == "no_agri"), 
                                     taxa = taxa)
  # trend_pass_converged_no_agri <- trend_change(occ_df = filter(occ_out_pass_converged_df,
  #                                                              agri == "no_agri"),
  #                                              taxa = taxa)
  # trend_pass_precision_no_agri <- trend_change(occ_df = filter(occ_out_pass_precision_df, 
  #                                                              agri == "no_agri"), 
  #                                              taxa = taxa)
  
#   return(list(trend_pa,
#               trend_unp))
#   
# })

  trend_pass_out <- list(trend_pass_high, trend_pass_low, trend_pass_no_agri)
  # trend_pass_converged_out <- list(trend_pass_converged_high, 
  #                                  trend_pass_converged_low, 
  #                                  trend_pass_converged_no_agri)
  # trend_pass_precision_out <- list(trend_pass_precision_high, 
  #                                  trend_pass_precision_low, 
  #                                  trend_pass_precision_no_agri)
  
saveRDS(trend_pass_out, "/data/outputs/trend_pass_out.rds")
  
  
# function to calculate annual growth rate per taxonomic group
tax_change <- function(trend_df) {
  
  trend_overall <- trend_df %>% 
    dplyr::mutate(tax_grp = "Overall")
  
  trend_df <- dplyr::bind_rows(trend_df, trend_overall)
  
  chng_tax <- trend_df %>% 
    dplyr::group_by(tax_grp, iteration) %>% 
    dplyr::summarise(change = mean(change, na.rm = TRUE))
}

# chng_out <- lapply(1:length(trend_out), function(p) {
#   
#   trend_out_p <- trend_out[[p]]
  
  chng_pass_high <- tax_change(trend_df = trend_pass_out[[1]])
  # chng_pass_converged_high <- tax_change(trend_df = trend_pass_converged_out[[1]])
  # chng_pass_precision_high <- tax_change(trend_df = trend_pass_precision_out[[1]])
  chng_pass_low <- tax_change(trend_df = trend_pass_out[[2]])
  # chng_pass_converged_low <- tax_change(trend_df = trend_pass_converged_out[[2]])
  # chng_pass_precision_low <- tax_change(trend_df = trend_pass_precision_out[[2]])
  chng_pass_no_agri <- tax_change(trend_df = trend_pass_out[[3]])
  # chng_pass_converged_no_agri <- tax_change(trend_df = trend_pass_converged_out[[3]])
  # chng_pass_precision_no_agri <- tax_change(trend_df = trend_pass_precision_out[[3]])
  
#   return(list(chng_pa,
#               chng_unp))
# })

  chng_pass_out <- list(chng_pass_high, chng_pass_low, chng_pass_no_agri)
  # chng_pass_converged_out <- list(chng_pass_converged_high, 
  #                                 chng_pass_converged_low, 
  #                                 chng_pass_converged_no_agri)
  # chng_pass_precision_out <- list(chng_pass_precision_high, 
  #                                 chng_pass_precision_low, 
  #                                 chng_pass_precision_no_agri)
  
  
# plot trends ----
  
# put the three agri categories together

  chng_pass_high$agri <- "high"
  # chng_pass_converged_high$agri <- "high"
  # chng_pass_precision_high$agri <- "high"
  
chng_pass_low$agri <- "low"
# chng_pass_converged_low$agri <- "low"
# chng_pass_precision_low$agri <- "low"

chng_pass_no_agri$agri <- "no_agri"
# chng_pass_converged_no_agri$agri <- "no_agri"
# chng_pass_precision_no_agri$agri <- "no_agri"



chng_pass_comb <- dplyr::bind_rows(chng_pass_high, chng_pass_low, chng_pass_no_agri) %>% 
  # order levels
  dplyr::mutate(agri = factor(agri, levels = c("no_agri", "low", "high")))

# chng_pass_converged_comb <- dplyr::bind_rows(chng_pass_converged_high, 
#                                              chng_pass_converged_low, 
#                                              chng_pass_converged_no_agri) %>% 
#   # order levels
#   dplyr::mutate(agri = factor(agri, levels = c("no_agri", "low", "high")))
# 
# chng_pass_precision_comb <- dplyr::bind_rows(chng_pass_precision_high, 
#                                              chng_pass_precision_low, 
#                                              chng_pass_precision_no_agri) %>% 
#   # order levels
#   dplyr::mutate(agri = factor(agri, levels = c("no_agri", "low", "high")))


chng_pass_comb_avg <- chng_pass_comb %>% 
  dplyr::group_by(agri, tax_grp) %>% 
  dplyr::summarise(mean = mean(change, na.rm = TRUE),
                   median = median(change, na.rm = TRUE),
                   low_ci = HDInterval::hdi(change,
                                            credMass = 0.95)[[1]],
                   upp_ci = HDInterval::hdi(change,
                                            credMass = 0.95)[[2]])

# chng_pass_converged_comb_avg <- chng_pass_converged_comb %>% 
#   dplyr::group_by(agri, tax_grp) %>% 
#   dplyr::summarise(mean = mean(change, na.rm = TRUE),
#                    median = median(change, na.rm = TRUE),
#                    low_ci = HDInterval::hdi(change,
#                                             credMass = 0.95)[[1]],
#                    upp_ci = HDInterval::hdi(change,
#                                             credMass = 0.95)[[2]])
# 
# chng_pass_precision_comb_avg <- chng_pass_precision_comb %>% 
#   dplyr::group_by(agri, tax_grp) %>% 
#   dplyr::summarise(mean = mean(change, na.rm = TRUE),
#                    median = median(change, na.rm = TRUE),
#                    low_ci = HDInterval::hdi(change,
#                                             credMass = 0.95)[[1]],
#                    upp_ci = HDInterval::hdi(change,
#                                             credMass = 0.95)[[2]])
# 

chng_pass_comb_avg$tax_grp <- factor(chng_pass_comb_avg$tax_grp,      # Reordering group factor levels
                            levels = c(taxa, "Overall"))

# chng_pass_converged_comb_avg$tax_grp <- factor(chng_pass_converged_comb_avg$tax_grp,      # Reordering group factor levels
#                                      levels = c(taxa, "Overall"))
# 
# chng_pass_precision_comb_avg$tax_grp <- factor(chng_pass_precision_comb_avg$tax_grp,      # Reordering group factor levels
#                                      levels = c(taxa, "Overall"))
# 

chng_pass_comb$tax_grp <- factor(chng_pass_comb$tax_grp,      # Reordering group factor levels
                            levels = c(taxa, "Overall"))
# chng_pass_converged_comb$tax_grp <- factor(chng_pass_converged_comb$tax_grp,      # Reordering group factor levels
#                                  levels = c(taxa, "Overall"))
# chng_pass_precision_comb$tax_grp <- factor(chng_pass_precision_comb$tax_grp,      # Reordering group factor levels
#                                  levels = c(taxa, "Overall"))
# 

# save(chng_pass_comb, chng_pass_comb_avg, file = "/data/outputs/trends.rds")
save(chng_pass_comb, chng_pass_comb_avg, file = "/data/outputs/trends_no_harlequin.rds")
# save(chng_pass_converged_comb, chng_pass_converged_comb_avg, file = "/data/outputs/trends_converged.rds")

chng_pass_plot <- ggplot(chng_pass_comb, aes(x = agri, y = change, colour = agri)) +
  # zero line
  geom_hline(yintercept = 0, colour = "grey", lty = 2, lwd = 1) +
  # plot iterations
  geom_jitter(alpha = 0.2, width = 0.3) +
  # plot mean and ci
  geom_pointrange(data = chng_pass_comb_avg, 
                  aes(ymax = upp_ci, ymin = low_ci, y = mean), 
                  size = 1.5, shape = 20, 
                  colour = rep(c("dodgerblue3", "goldenrod", "firebrick3"), 
                               each = (ntax + 1))) +
  # colour
  scale_colour_manual(values = c("dodgerblue", "goldenrod1", "firebrick1")) +
  # x labels
  scale_x_discrete(labels = c("No arable", "Low", "High")) +
  # axis labels
  labs(x = "", y = "Annual growth rate") +
  facet_wrap(~tax_grp) +
  theme_light() +
  theme(legend.position = "none",
        text = element_text(size = 11))+
  ggtitle("Pass Rule-of-thumb")


# chng_pass_converged_plot <- ggplot(chng_pass_converged_comb, 
#                                    aes(x = agri, y = change, colour = agri)) +
#   # zero line
#   geom_hline(yintercept = 0, colour = "grey", lty = 2, lwd = 1) +
#   # plot iterations
#   geom_jitter(alpha = 0.2, width = 0.3) +
#   # plot mean and ci
#   geom_pointrange(data = chng_pass_converged_comb_avg, 
#                   aes(ymax = upp_ci, ymin = low_ci, y = mean), 
#                   size = 1.5, shape = 20, 
#                   colour = rep(c("dodgerblue3", "goldenrod", "firebrick3"), 
#                                each = (ntax + 1))) +
#   # colour
#   scale_colour_manual(values = c("dodgerblue", "goldenrod1", "firebrick1")) +
#   # x labels
#   scale_x_discrete(labels = c("No arable", "Low", "High")) +
#   # axis labels
#   labs(x = "", y = "Annual growth rate") +
#   facet_wrap(~tax_grp) +
#   theme_light() +
#   theme(legend.position = "none",
#         text = element_text(size = 11))+
#   ggtitle("Pass Rule-of-thumb and first/last year converged")
# 
# chng_pass_precision_plot <- ggplot(chng_pass_precision_comb, 
#                                    aes(x = agri, y = change, colour = agri)) +
#   # zero line
#   geom_hline(yintercept = 0, colour = "grey", lty = 2, lwd = 1) +
#   # plot iterations
#   geom_jitter(alpha = 0.2, width = 0.3) +
#   # plot mean and ci
#   geom_pointrange(data = chng_pass_precision_comb_avg, 
#                   aes(ymax = upp_ci, ymin = low_ci, y = mean), 
#                   size = 1.5, shape = 20, 
#                   colour = rep(c("dodgerblue3", "goldenrod", "firebrick3"), 
#                                each = (ntax + 1))) +
#   # colour
#   scale_colour_manual(values = c("dodgerblue", "goldenrod1", "firebrick1")) +
#   # x labels
#   scale_x_discrete(labels = c("No arable", "Low", "High")) +
#   # axis labels
#   labs(x = "", y = "Annual growth rate") +
#   facet_wrap(~tax_grp) +
#   theme_light() +
#   theme(legend.position = "none",
#         text = element_text(size = 11))+
#   ggtitle("Pass Rule-of-thumb and precision threshold")
# 

# # save: growth rate plot
png("/data/plots/all_taxa_pass_growth_rate.png", 
    height = 13, width = 15, units = "cm", res = 320)

chng_pass_plot

dev.off()

# png("/data/plots/all_taxa_pass_converged_growth_rate.png", 
#     height = 13, width = 15, units = "cm", res = 320)
# 
# chng_pass_converged_plot
# 
# dev.off()
# 
# png("/data/plots/all_taxa_pass_precision_growth_rate.png", 
#     height = 13, width = 15, units = "cm", res = 320)
# 
# chng_pass_precision_plot
# 
# dev.off()

# png("/data/plots/all_taxa_pass_precision87_growth_rate.png", 
#     height = 13, width = 15, units = "cm", res = 320)
# 
# chng_pass_precision_plot
# 
# dev.off()

## Difference in trends between regions ----

diff_pass_summary <- chng_pass_comb %>%
  spread(agri, change) %>%
  rowwise() %>%
  mutate(high_low = high - low,
         high_no_agri = high - no_agri,
         low_no_agri = low - no_agri) %>%
  select(-c(high, low, no_agri)) %>%
  group_by(tax_grp) %>%
  summarise(mean_highlow = mean(high_low),
            lowCI_highlow = HDInterval::hdi(high_low,
                                              credMass = 0.95)[[1]],
            uppCI_highlow = HDInterval::hdi(high_low,
                                              credMass = 0.95)[[2]],
            mean_highnoagri = mean(high_no_agri),
            lowCI_highnoagri = HDInterval::hdi(high_no_agri,
                                              credMass = 0.95)[[1]],
            uppCI_highnoagri = HDInterval::hdi(high_no_agri,
                                              credMass = 0.95)[[2]],
            mean_lownoagri = mean(low_no_agri),
            lowCI_lownoagri = HDInterval::hdi(low_no_agri,
                                                  credMass = 0.95)[[1]],
            uppCI_lownoagri = HDInterval::hdi(low_no_agri,
                                                  credMass = 0.95)[[2]]) %>% 
    group_by(tax_grp) %>% 
    gather(key = diff, value = value, mean_highlow:uppCI_lownoagri) %>%
    mutate(metric = sub("_.*", "", diff),
           regions = sub(".*_", "", diff)) %>%
    select(-diff) %>% 
    spread(key = metric, value = value) %>%
    ungroup()


diff_pass <- chng_pass_comb %>%
  spread(agri, change) %>%
  rowwise() %>%
  mutate(high_low = high - low,
         high_noagri = high - no_agri,
         low_noagri = low - no_agri)%>%
  select(-c(high, low, no_agri)) %>% 
  group_by(tax_grp, iteration) %>% 
  gather(key = diff, value = value, high_low:low_noagri) %>%
  mutate(regions = sub("_", "", diff)) %>%
  select(-diff) %>%
  ungroup() %>%
  rename(diff = value)


# saveRDS(diff_pass, "/data/outputs/trends_diff.rds")
saveRDS(diff_pass_summary, "/data/outputs/trends_diff_summary_no_harlequin.rds")
saveRDS(diff_pass, "/data/outputs/trends_diff_no_harlequin.rds")


## which species in each subset? ----

pass_keep <- readRDS("/data/outputs/pass_keep.rds") %>%
  dplyr::filter(pass == TRUE)

pass_first_last_converged_keep <- readRDS("/data/outputs/pass_first_last_converged_keep.rds") %>%
  dplyr::filter(pass == TRUE & converged == TRUE)

# pass_precision_keep <- readRDS("/data/outputs/pass_precision_keep.rds") %>%
#   dplyr::filter(pass == TRUE & converged == TRUE)


pass_keep_hov <- pass_keep %>%
  filter(tax_grp == "Hoverflies") %>%
  group_by(species) %>%
  summarise() %>%
  select(species=species)

pass_first_last_converged_keep_hov <- pass_first_last_converged_keep %>%
  filter(tax_grp == "Hoverflies") %>%
  group_by(species) %>%
  summarise() %>%
  select(species=species)

# pass_precision_keep_bees <- pass_precision_keep %>%
#   filter(tax_grp == "Bees") %>%
#   group_by(species) %>%
#   summarise() %>%
#   select(species=species)



# Generate plot
v <- venn.diagram(list(rt=pass_keep_hov$species, 
                       conv=pass_first_last_converged_keep_hov$species
                       # , prec = pass_precision_keep_bees$species
                       ),
                  fill = c("orange", "blue"
                           # , "yellow"
                           ),
                  alpha = c(0.5, 0.5
                            # , 0.5
                            ), cat.cex = 1.5, cex=1,
                  filename=NULL)

# have a look at the default plot
grid.newpage()
grid.draw(v)

# get species names in the graph
overlaps <- calculate.overlap(list(rt=pass_keep_hov$species, 
                                   conv=pass_first_last_converged_keep_hov$species
                                   # ,
                                   # prec = pass_precision_keep_bees$species
                                   )) 

overlaps <- overlaps[ sort(names(overlaps)) ] 

spp_in_rt <- overlaps$a1
spp_rt_conv <- overlaps$a2
# spp_rt_prec <- overlaps$a4
spp_all <- overlaps$a5

# what is different in the trends?


occ_out_pass_summary <- readRDS("/data/outputs/occ_pass_summary_no_harlequin.rds")

occ_out_pass_converged_summary <- readRDS("/data/outputs/occ_pass_converged_summary.rds")

# then plot


occ_pass_converged_hov_by_agri <- ggplot(
  data = occ_out_pass_converged_summary %>% 
    filter(tax_grp == "Hoverflies")) +
  geom_line(aes(x = as.numeric(year), y = mean,
                colour = agri)) +
  geom_ribbon(aes(
    x = as.numeric(year),
    ymin = low_ci,
    ymax = upp_ci,
    fill = agri
  ),
  alpha = 0.2,
  show.legend = FALSE) +
  # colour
  scale_colour_manual(
    values = c("firebrick1", "goldenrod1", "dodgerblue"),
    labels = c("High", "Low", "No cropland"),
    name = "Cropland\ncover"
  ) +
  scale_fill_manual(
    values = c("firebrick1", "goldenrod1", "dodgerblue"),
    labels = c("High", "Low", "No cropland")
  ) +
  # force lower limit to 0
  expand_limits(y = c(0, 0.3)) +
  # axis labels
  labs(x = "Year", y = "Average occupancy") +
  theme_light() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)#,
        #legend.position = c(0.75, 0.1) # c(0,0) bottom left, c(1,1) top-right
  )+
  ggtitle("Occupancy for species in converged subset")




occ_pass_hov_by_agri <- ggplot(
  data = occ_out_pass_summary %>% 
    filter(tax_grp == "Hoverflies")) +
  geom_line(aes(x = as.numeric(year), y = mean,
                colour = agri)) +
  geom_ribbon(aes(
    x = as.numeric(year),
    ymin = low_ci,
    ymax = upp_ci,
    fill = agri
  ),
  alpha = 0.2,
  show.legend = FALSE) +
  # colour
  scale_colour_manual(
    values = c("firebrick1", "goldenrod1", "dodgerblue"),
    labels = c("High", "Low", "No cropland"),
    name = "Cropland\n cover"
  ) +
  scale_fill_manual(
    values = c("firebrick1", "goldenrod1", "dodgerblue"),
    labels = c("High", "Low", "No cropland")
  ) +
  # force lower limit to 0
  expand_limits(y = c(0, 0.3)) +
  # axis labels
  labs(x = "Year", y = "Average occupancy") +
  theme_light() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)#,
        #legend.position = c(0.75, 0.1) # c(0,0) bottom left, c(1,1) top-right
  )+
  ggtitle("Occupancy for species in rules-of-thumb subset")


png("/data/plots/occ_rt_vs_converged.png", 
    height = 12, width = 15, units = "cm", res = 320)

plot_grid(occ_pass_hov_by_agri,
          occ_pass_converged_hov_by_agri,
          # ,
          # occ_pass_precision_bees_by_agri,
          ncol = 1, nrow = 2)

dev.off()


