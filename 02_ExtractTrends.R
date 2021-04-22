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
# remotes::install_github("https://github.com/03rcooke/wrappeR", 
#                         ref = "main", force = TRUE)
library(wrappeR) # wrappeR: multi-species indicators
library(tidyr)
library(ggplot2)
library(HDInterval)


# global parameters
startyear <- 1990
endyear <- 2019

# First we extract the occupancy data 
# for all taxonomic groups we are interested in ----

# Taxonomic group(s)
taxa <- c("Bees", "Carabids", "Spiders", "PlantBugs", "Ladybirds")

ntax <- length(taxa)

regions <- c("high", "no_agri", "low")

vers <- c("2021_Francesca_bwars_rerun", "2021_Francesca_Charlie_rerun",
          "2021_Francesca_marlog_rerun", "2021_Francesca_Charlie_rerun",
          "2021_Francesca")

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
  load(fileName)
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
meta_high <- meta_high %>% 
  dplyr::mutate(species = firstup(Species_r_high)) 
# %>% 
#   # the ants are really inconsistent
#   dplyr::mutate(species = ifelse(tax_grp == "Ants", allup(Species_r_pa), species)) %>% 
#   dplyr::mutate(species = ifelse(tax_grp == "Ants" & Species_r_pa %in% c("tetramorium atratulum", "hypoponera punctatissima agg"), firstup(Species_r_pa), species))

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
    
    # read model for species
    if(file.exists(paste0(sp_roster$modPath, sp_roster$group,
                          "/occmod_outputs/", sp_roster$ver, "/",
                          sp_meta$species, ".rdata"))){
    mod <- load_rdata(paste0(sp_roster$modPath, sp_roster$group, 
                             "/occmod_outputs/", sp_roster$ver, "/", 
                             sp_meta$species, ".rdata"))
    
    # extract species metrics from model
    out <- as.data.frame(attr(mod, "metadata")$analysis$spp_Metrics) %>%
      dplyr::mutate(species = sp_meta$species, 
                    tax_grp = sp_meta$tax_grp)
    
    print(paste("Complete:", sp_meta$species))
    }
    if(file.exists(paste0(sp_roster$modPath, sp_roster$group,
                                   "/occmod_outputs/", sp_roster$ver, "/",
                                   sp_meta$species, ".rds"))){
      mod <- readRDS(paste0(sp_roster$modPath, sp_roster$group, 
                               "/occmod_outputs/", sp_roster$ver, "/", 
                               sp_meta$species, ".rds"))
      
      # extract species metrics from model
      out <- as.data.frame(attr(mod, "metadata")$analysis$spp_Metrics) %>%
        dplyr::mutate(species = sp_meta$species, 
                      tax_grp = sp_meta$tax_grp)
      
      print(paste("Complete:", sp_meta$species))
    }
    return(out)
    
  }) %>% 
    # bind across species
    dplyr::bind_rows(.)
  
  return(qual_df)
  
}

mod_met <- mod_qual(roster = roster, df = meta_high)

pass_keep <- mod_met %>% 
  # EqualWt
  dplyr::mutate(pass = ifelse(prop_abs >= 0.990, P90 >= 3.1, P90 >= 6.7))

saveRDS(pass_keep, "/data/outputs/pass_keep.rds")

pass_keep <- readRDS("/data/outputs/pass_keep.rds") %>% 
  dplyr::filter(pass == TRUE)

pass_keep <- tolower(pass_keep$species)

## occupancy estimates ----
# Here we read in prepared occupancy data

# read in filtered data
occ_read <- function(x, reg) {
  
  load_rdata(paste0("/data/outputs/filtered/", x, "_all_", reg, "_samp.rdata")) %>% 
    .[[1]] %>% 
    dplyr::mutate(tax_grp = x)
  
}

# focal years 1990 to 2019
foc_yrs <- paste0("year_", startyear:2019)

# reg_pairs <- list(GB = regions[1:2],
#                   ENG = regions[3:4],
#                   SCO = regions[5:6],
#                   WAL = regions[7:8])

# occ_out <- lapply(1:length(reg_pairs), function(p) {
#   
#   reg_p <- reg_pairs[[p]]
#   
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
  
  # # unify GB protected and GB unprotected species lists
  # all_spp <- list(spp_high, spp_low, spp_no_agri) %>%
  #   Reduce(inner_join, .)
  # # all_spp <- dplyr::inner_join(spp_high, spp_low) %>%
  # #   dplyr::inner_join(spp_no_agri)
  # 
  # # update to match species lists
  # occ_high <- dplyr::filter(occ_high, species %in% all_spp$spp)
  # occ_low <- dplyr::filter(occ_low, species %in% all_spp$spp)
  # occ_no_agri <- dplyr::filter(occ_no_agri, species %in% all_spp$spp)
  # 
  # remove species that didn't pass model quality threshold
  occ_high <- dplyr::filter(occ_high, species %in% pass_keep)
  occ_low <- dplyr::filter(occ_low, species %in% pass_keep)
  occ_no_agri <- dplyr::filter(occ_no_agri, species %in% pass_keep)
  
#   return(list(occ_high, occ_low, occ_no_agri))
#   
# })
  
occ_out <- list(occ_high, occ_low, occ_no_agri)

# GB_spp <- occ_out[[1]] %>% 
#   dplyr::distinct(species, .keep_all = TRUE) %>% 
#   dplyr::count(tax_grp)

# ENG_spp <- occ_out[[2]][[1]] %>% 
#   dplyr::distinct(species, .keep_all = TRUE) %>% 
#   dplyr::count(tax_grp)
# 
# SCOT_spp <- occ_out[[3]][[1]] %>% 
#   dplyr::distinct(species, .keep_all = TRUE) %>% 
#   dplyr::count(tax_grp)
# 
# WAL_spp <- occ_out[[4]][[1]] %>% 
#   dplyr::distinct(species, .keep_all = TRUE) %>% 
#   dplyr::count(tax_grp)

## trend lines by taxa and coverage of arable land

# function to calculate geometric mean
geomean <- function(x) exp(mean(log(x)))

# collapse occ_out list into a dataframe 

occ_out_df <- do.call(rbind, occ_out)

# summarise_occ takes a subset of the occ_out dataframe
# by taxon group, summarises occupancy across species by iteration
# transforms the dataframe into a long format
# then summarises (mean, median and geometric mean)
# across iterations and calculates CI by taxa and arable land cover

summarise_occ <- function(df, taxa) {

  df %>%
    filter(tax_grp == taxa) %>%
    dplyr::group_by(agri, iteration) %>%
    dplyr::summarise_at(
      vars(starts_with("year_")),
      mean, na.rm = TRUE) %>%
    ungroup() %>%
    gather(year_1990:year_2019, key = "year", value = "occupancy") %>%
    mutate(year = as.numeric(gsub("year_", "", year)))%>%
    group_by(agri, year) %>%
    summarise(
      mean = mean(occupancy),
      median = median(occupancy),
      gmean = geomean(occupancy),
      low_ci = HDInterval::hdi(occupancy,
                              credMass = 0.95)[[1]],
      upp_ci = HDInterval::hdi(occupancy,
                              credMass = 0.95)[[2]]) %>%
    mutate(tax_grp = taxa)
}


# apply the function to all taxa and bind together

occ_summary <- lapply(unique(occ_out_df$tax_grp), 
                      summarise_occ, 
                      df = occ_out_df) %>%
  dplyr::bind_rows(.)

 
# then plot

occ_by_agri <- ggplot(
  data = occ_summary) +
  geom_line(aes(x = as.numeric(year), y = gmean,
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
  theme_light() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
        legend.position = c(0.75, 0.1) # c(0,0) bottom left, c(1,1) top-right
        )


# # save: trend plot
png("/data/plots/all_taxa_trends.png", 
    height = 16, width = 14, units = "cm", res = 320)

occ_by_agri

dev.off()


## Trend change ----
# Here we calculate annual growth rate between the first and last year

# function to calculate annual growth rate per species
trend_change <- function(occ_df, taxa) {
  
  out <- lapply(taxa, function(x) {
    
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
      dplyr::mutate(occ_first = ifelse(occ_first == 0, 1e-07, occ_first)) %>% 
      # annual growth rate - see occurrenceChange()
      dplyr::mutate(change = (((occ_last / occ_first) ^ (1 / nyr)) - 1) * 100)
    
  }) %>% 
    # bind across taxonomic groups
    dplyr::bind_rows(.)
  
}

# trend_out <- lapply(1:length(occ_out), function(p) {
#   
#   occ_out_p <- occ_out[[p]]
  
  trend_high <- trend_change(occ_df = occ_out[[1]], taxa = taxa)
  trend_low <- trend_change(occ_df = occ_out[[2]], taxa = taxa)
  trend_no_agri <- trend_change(occ_df = occ_out[[3]], taxa = taxa)
  
#   return(list(trend_pa,
#               trend_unp))
#   
# })

  trend_out <- list(trend_high, trend_low, trend_no_agri)
  
# function to calculate annual growth rate per taxonomic group
tax_change <- function(trend_df) {
  
  trend_overall <- trend_df %>% 
    dplyr::mutate(tax_grp = "Overall")
  
  trend_df <- dplyr::bind_rows(trend_df, trend_overall)
  
  chng_tax <- trend_df %>% 
    dplyr::group_by(tax_grp, iteration) %>% 
    dplyr::summarise(change = mean(change))
}

# chng_out <- lapply(1:length(trend_out), function(p) {
#   
#   trend_out_p <- trend_out[[p]]
  
  chng_high <- tax_change(trend_df = trend_out[[1]])
  chng_low <- tax_change(trend_df = trend_out[[2]])
  chng_no_agri <- tax_change(trend_df = trend_out[[3]])
  
#   return(list(chng_pa,
#               chng_unp))
# })

  chng_out <- list(chng_high, chng_low, chng_no_agri)
  
# plot trends ----
  
# put the three agri categories together

chng_high$agri <- "high"
chng_low$agri <- "low"
chng_no_agri$agri <- "no_agri"

chng_comb <- dplyr::bind_rows(chng_high, chng_low, chng_no_agri) %>% 
  # order levels
  dplyr::mutate(agri = factor(agri, levels = c("no_agri", "low", "high")))

chng_comb_avg <- chng_comb %>% 
  dplyr::group_by(agri, tax_grp) %>% 
  dplyr::summarise(mean = mean(change),
                   median = median(change),
                   low_ci = HDInterval::hdi(change,
                                            credMass = 0.95)[[1]],
                   upp_ci = HDInterval::hdi(change,
                                            credMass = 0.95)[[2]])

chng_comb_avg$tax_grp <- factor(chng_comb_avg$tax_grp,      # Reordering group factor levels
                            levels = c("Bees", "Carabids", "Ladybirds", "PlantBugs",
                                       "Spiders", "Overall"))

  
chng_comb$tax_grp <- factor(chng_comb$tax_grp,      # Reordering group factor levels
                            levels = c("Bees", "Carabids", "Ladybirds", "PlantBugs",
                                       "Spiders", "Overall"))

chng_plot <- ggplot(chng_comb, aes(x = agri, y = change, colour = agri)) +
  # zero line
  geom_hline(yintercept = 0, colour = "grey", lty = 2, lwd = 1) +
  # plot iterations
  geom_jitter(alpha = 0.2, width = 0.3) +
  # plot mean and ci
  geom_pointrange(data = chng_comb_avg, 
                  aes(ymax = upp_ci, ymin = low_ci, y = mean), 
                  size = 1.5, shape = 20, 
                  colour = rep(c("dodgerblue3", "goldenrod", "firebrick3"), each = 6)) +
  # colour
  scale_colour_manual(values = c("dodgerblue", "goldenrod1", "firebrick1")) +
  # x labels
  scale_x_discrete(labels = c("No arable", "Low", "High")) +
  # axis labels
  labs(x = "", y = "Annual growth rate") +
  facet_wrap(~tax_grp) +
  theme_light() +
  theme(legend.position = "none",
        text = element_text(size = 11))

# # save: growth rate plot
png("/data/plots/all_taxa_growth_rate.png", 
    height = 13, width = 15, units = "cm", res = 320)

chng_plot

dev.off()


##  Species richness -----
# Here we calculate species richness across and per year

# function to calculate species richness per year and across years
sprich <- function(occ_df) {
  
  occ_overall <- occ_df %>% 
    dplyr::mutate(tax_grp = "Overall")
  
  occ_df <- dplyr::bind_rows(occ_df, occ_overall)
  
  spr_df <- occ_df %>% 
    dplyr::group_by(tax_grp, iteration, agri) %>% 
    # species richness per year
    dplyr::summarise_at(vars(dplyr::starts_with("year_")), sum) %>% 
    dplyr::ungroup() %>% 
    # species richness across years
    dplyr::mutate(spr = rowMeans(dplyr::select(., dplyr::starts_with("year_")), 
                                 na.rm = TRUE))
  
}

# spr_out <- lapply(1:length(occ_out), function(p) {
#   
#   occ_out_p <- occ_out[[p]]
  
  spr_high <- sprich(occ_df = occ_out[[1]])
  spr_low <- sprich(occ_df = occ_out[[2]])
  spr_no_agri <- sprich(occ_df = occ_out[[3]])
  
  spr_out <- list(spr_high, spr_low, spr_no_agri)
  
#   return(list(spr_pa,
#               spr_unp))
#   
# })

  
  spr_comb <- dplyr::bind_rows(spr_high, spr_low, spr_no_agri) %>% 
    # order levels
    dplyr::mutate(agri = factor(agri, levels = c("no_agri", "low", "high")))
  
  spr_comb_avg <- spr_comb %>% 
    dplyr::group_by(agri, tax_grp) %>% 
    dplyr::summarise(mean = mean(spr),
                     median = median(spr),
                     low_ci = HDInterval::hdi(spr,
                                              credMass = 0.95)[[1]],
                     upp_ci = HDInterval::hdi(spr,
                                              credMass = 0.95)[[2]])
  
  spr_comb$tax_grp <- factor(spr_comb$tax_grp,      # Reordering group factor levels
                             levels = c("Bees", "Carabids", "Ladybirds", "PlantBugs",
                                         "Spiders", "Overall"))
  
  
  spr_comb_avg$tax_grp <- factor(spr_comb_avg$tax_grp,      # Reordering group factor levels
                              levels = c("Bees", "Carabids", "Ladybirds", "PlantBugs",
                                         "Spiders", "Overall"))
  
  
  
  spr_plot <- ggplot(spr_comb, aes(x = agri, y = spr, colour = agri)) +
    # plot iterations
    geom_jitter(alpha = 0.2, width = 0.3) +
    # plot mean and ci
    geom_pointrange(data = spr_comb_avg, 
                    aes(ymax = upp_ci, ymin = low_ci, y = mean), 
                    size = 1.5, shape = 20, 
                    colour = rep(c("dodgerblue3", "goldenrod", "firebrick3"), each = 6)) +
    # colour
    scale_colour_manual(values = c("dodgerblue", "goldenrod1", "firebrick1")) +
    # x labels
    scale_x_discrete(labels = c("No arable", "Low", "High")) +
    # axis labels
    labs(x = "", y = "Species richness") +
    facet_wrap(~tax_grp) +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(size = 11))
  
  
  # # save: species richness plot
  png("/data/plots/all_taxa_sp_rich.png", 
      height = 13, width = 15, units = "cm", res = 320)
  
  spr_plot
  
  dev.off()
  
  