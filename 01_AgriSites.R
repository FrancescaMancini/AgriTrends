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
  rename(easting = x, northing = y) %>% 
  # get 4-figure grid reference from 6-figure grid reference
  dplyr::mutate(gr = BRCmap::reformat_gr(.$gr, prec_out = 1000))

pa <- read.csv("/data/data/pro_stat_grid_ref.csv", stringsAsFactors = FALSE)

grid_ref <- read.csv("/data/data/gr_ref.csv", stringsAsFactors = FALSE)

pa_agri_regions <- grid_ref %>%
  left_join(agri, by = c("grid_ref" = "gr")) 
pa_agri_regions <- pa_agri_regions[, c(1,4,9)]

pa_agri_regions <- pa_agri_regions %>%
  left_join(pa, by = c("grid_ref", "region"))

pa_agri_regions <- pa_agri_regions[, c(1,2,3,7)]

write.csv(pa_agri_regions, "/data/data/agri_pa_regions.csv", row.names = FALSE)



# inital master sheet did not have all sites
# correct and add the extra sites

sites <- readRDS("/data/data/sq1km_UK_regions.rds") %>%
  gather(key = "region", value = "number", 
         ENGLAND, WALES, SCOTLAND, NORTHERN_IRELAND) %>%
  filter(number == 1) %>%
  mutate(region = case_when(region == "ENGLAND" ~ "ENG",
                            region == "SCOTLAND" ~ "SCO",
                            region == "WALES" ~ "WAL",
                            region == "NORTHERN_IRELAND" ~ "NIR"))

sites <- sites[, -3]
extra_sites <- subset(sites, !(sites$SQ1_SQUARE %in% pa_agri_regions$grid_ref)) %>%
  rename(grid_ref = SQ1_SQUARE)


pa_agri_regions_complete <- bind_rows(pa_agri_regions, extra_sites)
write.csv(pa_agri_regions_complete, 
          "/data/data/agri_pa_regions.csv", row.names = FALSE)

# Load up the regions data frame
agri_reg <- read.csv('/data/data/agri_pa_regions.csv', stringsAsFactors = FALSE)
#Replace NAs with none
agri_reg[is.na(agri_reg)] <- 'none'
# Create a region which is a concatenation of all regions
agri_reg$concat <- paste(agri_reg$region,agri_reg$agri_cover,agri_reg$prot,sep = '_')
# Create a value column for when we cast the table
agri_reg$value <- 1

# Now we build the regions dataframe
# Convert long table to wide
regions_df <- reshape2::dcast(data = agri_reg, formula = grid_ref ~ concat, value.var = 'value')
# And now replace NAs with 0. Not strictly necessary for sparta, but it is neater
regions_df[is.na(regions_df)] <- 0

# Now we have our regions dataframe, we need to create the aggregates
# First create a long master lookup table of all regions, with the same headers
regions_lookup <-
  merge(data.frame(region = unique(agri_reg$region)),
        data.frame(agri_cover = unique(agri_reg$agri_cover))) %>%
  merge(data.frame(prot = unique(agri_reg$prot)))
regions_lookup$concat <- paste(regions_lookup$region,regions_lookup$agri_cover,
                               regions_lookup$prot,sep = '_')
regions_lookup[] <- lapply(regions_lookup, as.character)
regions_lookup <- melt(regions_lookup,id.vars = c('concat'))
# Add GB to the table
GB <- regions_lookup[regions_lookup$value %in% c('SCO','ENG','WAL'),]
GB$value <- 'GB'
regions_lookup <- bind_rows(regions_lookup, GB)
# add country-protected for Rob
country_prot <- regions_lookup[regions_lookup$value %in% c('SCO','ENG','WAL',
                                                           'unp', 'pa'),]
country_prot <- country_prot %>%
  group_by(concat) %>%
  dplyr::summarise(variable = "region_pa", 
            value = paste(value, collapse="_")) %>%
  filter(!(value %in% c('SCO','ENG','WAL')))

regions_lookup <- bind_rows(regions_lookup, country_prot)

# Just in case one of the possible permutations was missing (e.g. SCO_no_agri_unp),
# let's remove it from this lookup table
perms_in_data <- regions_lookup$concat[regions_lookup$concat %in% names(regions_df)[-1]]
regions_lookup <- regions_lookup[regions_lookup$concat %in% perms_in_data,]

# Now, work out all the regions we might be looking at e.g. SCO, no_agri, or GB
regions <- unique(regions_lookup$value)
regions <- regions[regions!='none']
# Generate a list of all concatenated sub-regions which lie in each broad region
region_aggs <- lapply(regions, FUN = function(region){
  c(regions_lookup$concat[regions_lookup$value == region])
})
# And give them the right names
names(region_aggs) <- regions

# Now we're done. Save the outputs
saveRDS(object = regions_df, file = '/data/data/regions_df.rds')
saveRDS(object = region_aggs, file = '/data/data/region_aggs.rds')
saveRDS(object = region_aggs, file = '/data/data/region_aggs_country_pa.rds')
