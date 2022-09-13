########################################
## identify agricultural sites
## Author: Francesca Mancini
## Date created: 2020-10-06
## Date modified:
########################################

library(raster)
library(rgdal)
library(dplyr)
# to install BRCmap from zip file
# remotes::install_local("BRCmap_0.10.3.5.zip")
library(BRCmap)
library(reshape2)

# read the LCM 1990-2015 change 

LCM1990 <- raster("/data/data/LCMchange_1990-2015/data/07b6e5e9-b766-48e5-a28c-5b3e35abecc0/LCC_GB_1990_to_2015.tif",
                  band = 1)

LCM2015 <- raster("/data/data/LCMchange_1990-2015/data/07b6e5e9-b766-48e5-a28c-5b3e35abecc0/LCC_GB_1990_to_2015.tif",
                  band = 2)

# crop raster to England
# 
# england <- readOGR("/data/data/England_shp/England_AL4-AL4.shp")
# 
# england <- spTransform(england, CRS = LCM1990@crs)
# 
# 
# plot(LCM2015)
# lines(england)
# 
# LCM1990_Eng <- mask(crop(LCM1990, extent(england)), england)  # this took a while!
# LCM2015_Eng <- mask(crop(LCM2015, extent(england)), england)  # this took a while!
# 
# plot(LCM2015_Eng)
# lines(england)
# 
# rm(LCM2015)
# rm(LCM1990)
# 
# 
# writeRaster(LCM1990_Eng, filename = "/data/data/LCM1990_eng.tif")
# writeRaster(LCM2015_Eng, filename = "/data/data/LCM2015_eng.tif")
# 


# reclassify raster so that arable gets 1, and everything else gets 0
LCM1990 <- reclassify(LCM1990, matrix(c(0,NA,1,0,2,1,3,0,4,0,5,0,6,0), ncol = 2, byrow = TRUE))
LCM2015 <- reclassify(LCM2015, matrix(c(0,NA,1,0,2,1,3,0,4,0,5,0,6,0), ncol = 2, byrow = TRUE))

# aggregate rasters to 1Km
LCM1990_1km <- aggregate(LCM1990, fact = 1000/25, fun = mean)
LCM2015_1km <-  aggregate(LCM2015, fact = 1000/25, fun = mean)


# subtract the two rasters to calculate change in agriculture cover
LCM90_15_change <- LCM1990_1km - LCM2015_1km

# all the cells that have changed by more than 10% are excluded
LCM_nochange <- LCM2015_1km

LCM_nochange[LCM90_15_change<=-0.1 | LCM90_15_change>=0.1] <- NA

writeRaster(LCM_nochange, "/data/data/LCM_agri_nochange.tif")


## calculate the % coverage for categories 4, 5 and 6
# reclassify raster so that water, built-up areas and other get 1 and the rest gets 0
LCM1990 <- raster("/data/data/LCMchange_1990-2015/data/07b6e5e9-b766-48e5-a28c-5b3e35abecc0/LCC_GB_1990_to_2015.tif",
                  band = 1)

LCM2015 <- raster("/data/data/LCMchange_1990-2015/data/07b6e5e9-b766-48e5-a28c-5b3e35abecc0/LCC_GB_1990_to_2015.tif",
                  band = 2)

LCM1990_other <- reclassify(LCM1990, matrix(c(0,NA,1,0,2,0,3,0,4,1,5,1,6,1), ncol = 2, byrow = TRUE))
LCM2015_other <- reclassify(LCM2015, matrix(c(0,NA,1,0,2,0,3,0,4,1,5,1,6,1), ncol = 2, byrow = TRUE))

# aggregate rasters to 1Km
LCM1990_other_1km <- aggregate(LCM1990_other, fact = 1000/25, fun = mean)
LCM2015_other_1km <-  aggregate(LCM2015_other, fact = 1000/25, fun = mean)


# subtract the two rasters to calculate change in agriculture cover
LCM90_15_other_change <- LCM1990_other_1km - LCM2015_other_1km

# all the cells that have changed by more than 10% are excluded
LCM_other_nochange <- LCM2015_other_1km

LCM_other_nochange[LCM90_15_other_change<=-0.1 | LCM90_15_other_change>=0.1] <- NA

writeRaster(LCM_other_nochange, "/data/data/LCM_other_nochange.tif")


# convert ratser to point and exract dataframe
LCM_agri <- raster("/data/data/LCM_agri_nochange.tif")

agri_df <- as.data.frame(rasterToPoints(LCM_agri))

LCM_other <- raster("/data/data/LCM_other_nochange.tif")

other_df <- as.data.frame(rasterToPoints(LCM_other))

# merge agri and other

sites_df <- left_join(agri_df, other_df, by = c("x", "y"))
sites_df <- na.exclude(sites_df)

# create categories

# use distribution of coverage across sites to define categories
# look at distribution of coverage
hist(sites_df$LCM_agri_nochange)
summary(sites_df$LCM_agri_nochange)
quantile(sites_df$LCM_agri_nochange, 
         probs = c(0.33, 0.66))

plot(ecdf(sites_df$LCM_agri_nochange), 
     main = "ECDF of arable cover (all sites)")

# given the ECDF a sensible categorisation is:
# "no agriculture" = 0, "low agriculture" = 0-50, "high agriculture" = 50-100
# the ECDF of coverage of arable land in sites with
# surveys was also checked and it was very similar to that for all sites

sites_df <- sites_df %>%
  mutate(agri_cover = case_when(LCM_agri_nochange > 0.5 ~ "high",
                                LCM_agri_nochange <= 0.5 & LCM_agri_nochange > 0 ~ "low",
                                LCM_agri_nochange == 0 ~ "no_agri"),
         other_cover = case_when(LCM_other_nochange >= 0.75 ~ "high",
                                 LCM_other_nochange < 0.75 & LCM_other_nochange > 0.25 ~ "medium",
                                 LCM_other_nochange <= 0.25 ~ "low"))

# exclude all those sites with high to medium coverage of wetland, built-up areas and other

sites_df<- filter(sites_df, LCM_other_nochange <= 0.25)


# use BRCmap to convert x and y to gridref

sites_df$gr <- gr_num2let(sites_df$x, sites_df$y, keep_precision = FALSE)
# save

write.csv(sites_df, "/data/data/sites.csv", row.names = FALSE)


## make csv file with regions for occ mod updates

agri <- read.csv("/data/data/sites.csv", stringsAsFactors = FALSE) %>%
  dplyr::rename(easting = x, northing = y, grid_ref = gr) %>% 
  # get 4-figure grid reference from 6-figure grid reference
  dplyr::mutate(grid_ref = BRCmap::reformat_gr(.$grid_ref, 
                                               prec_out = 1000)) %>%
  select(grid_ref, agri_cover)

pa <- read.csv("/data/data/pro_stat_grid_ref.csv", stringsAsFactors = FALSE)%>%
  select(grid_ref, prot)

grid_ref <- read.csv("/data/data/gr_ref.csv", stringsAsFactors = FALSE) %>%
  mutate(region = case_when(region == "ENG" ~ "ENGLAND",
                            region == "SCO" ~ "SCOTLAND",
                            region == "WAL" ~ "WALES",
                            region == "NIR" ~ "NORTHERN_IRELAND")) %>%
  select(grid_ref, region)


# inital master sheet did not have all sites
# correct and add the extra sites

sites <- readRDS("/data/data/sq1km_UK_regions.rds") %>%
  tidyr::gather(key = "region", value = "number", 
         ENGLAND, WALES, SCOTLAND, NORTHERN_IRELAND) %>%
  filter(number == 1) #%>%
  # mutate(region = case_when(region == "NORTHERN_IRELAND" ~ "NORTHERN.IRELAND",
  #                           TRUE ~ region))

sites <- sites[, -3]

extra_sites <- subset(sites, !(sites$SQ1_SQUARE %in% grid_ref$grid_ref)) %>%
  dplyr::rename(grid_ref = SQ1_SQUARE)

grid_ref <- dplyr::bind_rows(grid_ref, extra_sites)

# reformat all as wide

grid_ref_wide <- grid_ref %>%
  mutate(value = 1) %>%
  tidyr::spread(key = region, value = value, fill = 0)

agri_wide <- agri %>%
  mutate(value = 1) %>%
  tidyr::spread(key = agri_cover, value = value, fill = 0)

pa_wide <- pa %>%
  mutate(value = 1) %>%
  tidyr::spread(key = prot, value = value, fill = 0)


# test every site has value of 1 for one column only

max(rowSums(grid_ref_wide[,2:5]))

max(rowSums(agri_wide[,2:4]))

max(rowSums(pa_wide[,2:3]))

# put them all together

regions_df <- full_join(grid_ref_wide, agri_wide,
                        by = "grid_ref") %>%
  replace(is.na(.), 0) %>%
  full_join(pa_wide, by = "grid_ref") %>%
  replace(is.na(.), 0)


# regions_df has more rows than grid_ref_wide

which(!(regions_df$grid_ref) %in% grid_ref_wide$grid_ref)
# there are 196 sites in the agri data frame that 
# are not in the country data frame 
# I am keeping them in as I don't want to loose sites


# pa_agri_regions <- grid_ref %>%
#   left_join(agri, by = c("grid_ref" = "gr")) 
# pa_agri_regions <- pa_agri_regions[, c(1,4,9)]
# 
# pa_agri_regions <- pa_agri_regions %>%
#   left_join(pa, by = c("grid_ref", "region"))
# 
# pa_agri_regions <- pa_agri_regions[, c(1,2,3,7)]
# 
# # write.csv(pa_agri_regions, "/data/data/agri_pa_regions.csv", row.names = FALSE)
# 
# pa_agri_regions <- read.csv( "/data/data/agri_pa_regions.csv", 
#                              stringsAsFactors = FALSE)
# 
# 
# pa_agri_regions_complete <- bind_rows(pa_agri_regions, extra_sites)
# # write.csv(pa_agri_regions_complete,
# #           "/data/data/agri_pa_regions.csv", row.names = FALSE)
# 
# 
# 
# # Load up the regions data frame
# agri_reg <- read.csv('/data/data/agri_pa_regions.csv', stringsAsFactors = FALSE)
# #Replace NAs with none
# agri_reg <- agri_reg %>%
#   mutate(agri_cover = tidyr::replace_na(agri_cover, "agriNA"))%>%
#   mutate(prot = tidyr::replace_na(prot, "protNA"))
# # # Create a region which is a concatenation of all regions
# # agri_reg$concat <- paste(agri_reg$region,
# #                          agri_reg$agri_cover,
# #                          agri_reg$prot,sep = '_')
# # Create a value column for when we cast the table
# agri_reg$value <- 1
# 
# # Now we build the regions dataframe
# # Convert long table to wide
# regions_df <- reshape2::dcast(data = agri_reg, 
#                               formula = grid_ref ~ region ~ agri_cover ~ prot, 
#                               value.var = 'value')
# # And now replace NAs with 0. Not strictly necessary for sparta, but it is neater
# regions_df[is.na(regions_df)] <- 0

# Now we have our regions dataframe, we need to create the aggregates
# # First create a long master lookup table of all regions, with the same headers
# regions_lookup <-
#   merge(data.frame(region = unique(agri_reg$region)),
#         data.frame(agri_cover = unique(agri_reg$agri_cover))) %>%
#   merge(data.frame(prot = unique(agri_reg$prot)))
# regions_lookup$concat <- paste(regions_lookup$region,regions_lookup$agri_cover,
#                                regions_lookup$prot,sep = '_')
# regions_lookup[] <- lapply(regions_lookup, as.character)
# regions_lookup <- reshape2::melt(regions_lookup,id.vars = c('concat'))
# # Add GB to the table
# GB <- regions_lookup[regions_lookup$value %in% c('SCOTLAND',
#                                                  'ENGLAND',
#                                                  'WALES'),]
# GB$value <- 'GB'
# regions_lookup <- bind_rows(regions_lookup, GB)
# 
# # add UK
# 
# UK <- regions_lookup[regions_lookup$value %in% c('SCOTLAND',
#                                                  'ENGLAND',
#                                                  'WALES',
#                                                  'NORTHERN_IRELAND'),]
# UK$value <- 'UK'
# regions_lookup <- bind_rows(regions_lookup, UK)
# 
# # add country-protected for Rob
# country_prot <- regions_lookup[regions_lookup$value %in% c('SCOTLAND',
#                                                            'ENGLAND',
#                                                            'WALES',
#                                                            'NORTHERN_IRELAND',
#                                                            'unp', 'pa'),]
# country_prot <- country_prot %>%
#   group_by(concat) %>%
#   dplyr::summarise(variable = "region_pa", 
#             value = paste(value, collapse="_")) %>%
#   filter(!(value %in% c('SCOTLAND','ENGLAND','WALES', 'NORTHERN_IRELAND')))
# 
# regions_lookup <- bind_rows(regions_lookup, country_prot)
# 
# # Just in case one of the possible permutations was missing (e.g. SCO_no_agri_unp),
# # let's remove it from this lookup table
# perms_in_data <- regions_lookup$concat[regions_lookup$concat %in% names(regions_df)[-1]]
# regions_lookup <- regions_lookup[regions_lookup$concat %in% perms_in_data,]
# 
# # Now, work out all the regions we might be looking at e.g. SCO, no_agri, or GB
# regions <- unique(regions_lookup$value)
# regions <- regions[regions!='none']
# # Generate a list of all concatenated sub-regions which lie in each broad region
# region_aggs <- lapply(regions, FUN = function(region){
#   c(regions_lookup$concat[regions_lookup$value == region])
# })
# # And give them the right names
# names(region_aggs) <- regions


region_aggs <- list(UK = as.character(c("ENGLAND","SCOTLAND", 
                           "WALES", "NORTHERN_IRELAND")),
                    GB = as.character(c("ENGLAND","SCOTLAND", "WALES")))

# Now we're done. Save the outputs
# saveRDS(object = regions_df, file = '/data/data/regions_df.rds')
# saveRDS(object = region_aggs, file = '/data/data/region_aggs.rds')
# saveRDS(object = region_aggs, file = '/data/data/region_aggs_country_pa.rds')
saveRDS(object = regions_df, file = '/data/data/regions_non-crossed.rds')
saveRDS(object = grid_ref_wide, file = '/data/data/regions_country.rds')
saveRDS(object = region_aggs, file = '/data/data/region_aggs_non-crossed.rds')
