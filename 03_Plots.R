####################################
## Generate figures
####################################

library(ggplot2)
library(ggridges)
library(dplyr)
library(wesanderson)
library(raster)
library(cowplot)
library(tidyr)
library(grid)
library(gridExtra)


# map

sites_df <- read.csv("/data/data/sites.csv")

## plot map of no_agri, low and high sites

arable_cover <- ggplot(data = sites_df) +
  geom_point(aes(x = x, y = y,
                 colour = agri_cover),
             shape = 15, size = 0.1) +
  scale_colour_manual(
    values = wes_palette("FantasticFox1"),
    labels = c("High", "Low", "No cropland"),
    name = "Cropland\ncover"
  ) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  coord_quickmap() +
  theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

png("/data/plots/arable_cover.png", 
    height = 16, width = 14, units = "cm", res = 320)

arable_cover

dev.off()

## latitude by region

lat_by_class <- ggplot(data = sites_df) +
  geom_boxplot(aes(x = agri_cover, y = y,
                   fill = agri_cover), show.legend = FALSE) +
  scale_fill_manual(
    values = wes_palette("FantasticFox1")) +
  scale_x_discrete(labels = c("High", "Low", "No cropland")) +
  xlab("Cropland cover region") +
  ylab("Latitude") +
  theme_classic() +
  theme(text = element_text(size = 11))

png("/data/plots/latitude_by_region.png", 
    height = 9, width = 10, units = "cm", res = 320)

lat_by_class

dev.off()


# elevation by region

elev <- raster("/data/data/Elevation_UK_1km.tif")

sites_df$elev <- raster::extract(elev, sites_df[c("x","y")])


elev_by_class <- ggplot(data = sites_df) +
  geom_boxplot(aes(x = agri_cover, y = elev,
                   fill = agri_cover), show.legend = FALSE) +
  scale_fill_manual(
    values = wes_palette("FantasticFox1")) +
  scale_x_discrete(labels = c("High", "Low", "No cropland")) +
  xlab("Cropland cover region") +
  ylab("Elevation") +
  theme_classic() +
  theme(text = element_text(size = 11))

png("/data/plots/elevation_by_region.png", 
    height = 9, width = 10, units = "cm", res = 320)

elev_by_class

dev.off()


png("/data/plots/elevation_and_latitude_by_region.png", 
    height = 9, width = 14, units = "cm", res = 320)

plot_grid(lat_by_class,
          elev_by_class,
          ncol = 2, nrow = 1)

dev.off()


# occupancy by region

occ_pass_summary <- readRDS("/data/outputs/occ_pass_summary_no_harlequin.rds")

# New facet label names for taxa
taxa_grp <- c("Bees", "Carabids", "Hoverflies", "Ladybirds", "Plant Bugs", "Spiders")
names(taxa_grp) <- c("Bees", "Carabids", "Hoverflies", "Ladybirds", "PlantBugs", "Spiders")

occ_by_region <- ggplot(
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
              scales = "free_x", 
              labeller = labeller(tax_grp = taxa_grp)) +
  # colour
  scale_colour_manual(
    values = wes_palette("FantasticFox1"),
    labels = c("High", "Low", "No cropland"),
    name = "Cropland cover"
  ) +
  scale_fill_manual(
    values = wes_palette("FantasticFox1"),
    labels = c("High", "Low", "No cropland")
  ) +
  scale_x_continuous(breaks = seq(1990, 2020, 5)) +
  # force lower limit to 0
  expand_limits(y = 0) +
  # axis labels
  labs(x = "Year", y = "Average occupancy") +
  # geom_text(data = n_species_annotation,
  #           aes(x = 2005, y = 0.7,
  #           label = paste0("n species = ", n_species))) +
  theme_classic() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
        legend.position = "top") 


png("/data/plots/all_taxa_occ_by_region.png", 
    height = 16, width = 14, units = "cm", res = 320)

occ_by_region

dev.off()


# trends by region
taxa_grp2 <- c("Bees", "Carabids", "Hoverflies", "Ladybirds", "Plant Bugs", "Spiders", "Overall")
names(taxa_grp2) <- c("Bees", "Carabids", "Hoverflies", "Ladybirds", "PlantBugs", "Spiders", "Overall")

load("/data/outputs/trends_no_harlequin.rds")

pal <- wes_palette("FantasticFox1", n = 3)

trends_by_region <- ggplot(chng_pass_comb, aes(x = agri, y = change, colour = agri)) +
  # zero line
  geom_hline(yintercept = 0, colour = "grey", lty = 2, lwd = 1) +
  # plot iterations
  geom_jitter(alpha = 0.2, width = 0.3) +
  # plot mean and ci
  geom_pointrange(data = chng_pass_comb_avg, 
                  aes(ymax = upp_ci, ymin = low_ci, y = mean), 
                  size = 1, shape = 20, 
                  colour = "grey") +
  # colour
  scale_colour_manual(values = c(pal[3], pal[2], pal[1])) +
  # x labels
  scale_x_discrete(labels = c("No cropland", "Low", "High")) +
  # axis labels
  labs(x = "", y = "Annual growth rate") +
  facet_wrap(~tax_grp, nrow = 4, ncol = 2,
             labeller = labeller(tax_grp = taxa_grp2)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 11))


png("/data/plots/all_taxa_trends_by_region.png", 
    height = 16, width = 14, units = "cm", res = 320)

trends_by_region

dev.off()


# trends by region converged only
taxa_grp2 <- c("Bees", "Carabids", "Hoverflies", "Ladybirds", "Plant Bugs", "Spiders", "Overall")
names(taxa_grp2) <- c("Bees", "Carabids", "Hoverflies", "Ladybirds", "PlantBugs", "Spiders", "Overall")

load("/data/outputs/trends_converged.rds")

pal <- wes_palette("FantasticFox1", n = 3)

trends_by_region <- ggplot(chng_pass_converged_comb, aes(x = agri, y = change, colour = agri)) +
  # zero line
  geom_hline(yintercept = 0, colour = "grey", lty = 2, lwd = 1) +
  # plot iterations
  geom_jitter(alpha = 0.2, width = 0.3) +
  # plot mean and ci
  geom_pointrange(data = chng_pass_converged_comb_avg, 
                  aes(ymax = upp_ci, ymin = low_ci, y = mean), 
                  size = 1, shape = 20, 
                  colour = "grey") +
  # colour
  scale_colour_manual(values = c(pal[3], pal[2], pal[1])) +
  # x labels
  scale_x_discrete(labels = c("No cropland", "Low", "High")) +
  # axis labels
  labs(x = "", y = "Annual growth rate") +
  facet_wrap(~tax_grp, nrow = 4, ncol = 2,
             labeller = labeller(tax_grp = taxa_grp2)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 11))


png("/data/plots/all_taxa_trends_by_region_converged.png", 
    height = 16, width = 14, units = "cm", res = 320)

trends_by_region

dev.off()


# difference in trends

diff_pass <- readRDS("/data/outputs/trends_diff_summary_no_harlequin.rds") %>%
  mutate(magnitude = case_when(mean > 0 & lowCI > 0 & uppCI > 0 ~ "Positive",
                               mean < 0 & lowCI < 0 & uppCI < 0 ~ "Negative",
                               TRUE ~ "No difference"))

diff_plot <- ggplot(diff_pass,
                    aes(x = regions, y = mean,
                        colour = magnitude)) +
  geom_hline(yintercept = 0, colour = "grey", lty = 2, lwd = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = lowCI,
                    ymax = uppCI,
                    colour = magnitude), width = .2) +
  facet_wrap(~tax_grp, ncol = 2) +
  # colour
  scale_colour_manual(values = c("darkcyan", "cornsilk3", "coral2"),
                      name = "") +
  # x labels
  scale_x_discrete(labels = c("High - Low", "High - No cropland", "Low - No cropland")) +
  # axis labels
  labs(x = "", y = "Difference in annual growth rates") +
  theme_classic() +
  theme(text = element_text(size = 10))
  
png("/data/plots/diff_trend_by_region.png", 
    height = 17, width = 16, units = "cm", res = 320)

diff_plot

dev.off()


# new plot with posterior distribution

diff_pass <- readRDS("/data/outputs/trends_diff_no_harlequin.rds") 

diff_pass_summary <- diff_pass %>%
  filter(regions == "highlow") %>%
  group_by(tax_grp) %>%
  summarise(mean = mean(diff),
            lowCI = HDInterval::hdi(diff,
                                    credMass = 0.95)[[1]],
            uppCI = HDInterval::hdi(diff,
                                    credMass = 0.95)[[2]]
            )

diff_plot <- ggplot(diff_pass %>%
                      filter(regions == "highlow"),
                    aes(x = diff, y = tax_grp,
                        fill = stat(x))) +
  geom_density_ridges_gradient(color = "transparent", show.legend = FALSE,
                               jittered_points = TRUE, position = "raincloud",
                               alpha = 0.7, scale = 0.9, point_color = "lightgrey",
                               rel_min_height = 0.001) +
  geom_point(inherit.aes = FALSE, data = diff_pass_summary, 
             aes(x= mean, y = tax_grp),show.legend = FALSE, color = "black",
             position = position_nudge(y = -0.25))+
  geom_errorbarh(inherit.aes = FALSE, data = diff_pass_summary, 
                 aes(xmin = lowCI, xmax = uppCI, y = tax_grp), 
                 height = 0.1, show.legend = FALSE, color = "black",
                 position = position_nudge(y = -0.25)) +
  geom_vline(xintercept = 0, colour = "black", lty = 2, lwd = 1) +
  scale_fill_viridis_c(option = "C") +
  # axis labels
  labs(y = "", x = "Difference in annual growth rates between high and low cropland cover") +
  # scale_y_discrete(labels = labeller(tax_grp = taxa_grp2)) +
  theme_classic() +
  theme(text = element_text(size = 10))



png("/data/plots/diff_trends_high_low.png", 
    height = 10, width = 15, units = "cm", res = 320)

diff_plot

dev.off()


# for converged species only

load("/data/outputs/trends_converged.rds")

diff_pass_converged <- chng_pass_converged_comb %>%
  spread(agri, change) %>%
  rowwise() %>%
  mutate(high_low = high - low,
         high_noagri = high - no_agri,
         low_noagri = low - no_agri)%>%
  dplyr::select(-c(high, low, no_agri)) %>% 
  group_by(tax_grp, iteration) %>% 
  gather(key = diff, value = value, high_low:low_noagri) %>%
  mutate(regions = sub("_", "", diff)) %>%
  dplyr::select(-diff) %>%
  ungroup() %>%
  rename(diff = value)


diff_pass$subset <- "Rules-of-thumb"
diff_pass_converged$subset <- "Converged"

diff <- rbind(diff_pass, diff_pass_converged)

diff_summary <- diff %>%
  filter(regions == "highlow") %>%
  group_by(subset, tax_grp) %>%
  summarise(mean = mean(diff),
            lowCI = HDInterval::hdi(diff,
                                    credMass = 0.95)[[1]],
            uppCI = HDInterval::hdi(diff,
                                    credMass = 0.95)[[2]]
  )


diff_plot <- ggplot(diff_summary,
                    aes(x = tax_grp, y = mean,
                        group = subset, colour = subset)) +
  geom_point(position=position_dodge(width=0.8)) +
  geom_errorbar(aes(ymin = lowCI,
                    ymax = uppCI,
                    colour = subset),
                position=position_dodge(width=0.8), width = .2) +
  # facet_wrap(~tax_grp, ncol = 4) +
  scale_colour_manual(name = "Species subset",
                      values = c("cadetblue4", "chocolate 3")) +
  geom_hline(yintercept = 0, linetype = 4) +
  xlab("") +
  ylab("Mean difference in annual growth rates\nbetween high and low cropland cover") +
  theme_classic() +
  theme(text = element_text(size = 10))

png("/data/plots/diff_trends_converged_high_low.png", 
    height = 11, width = 15, units = "cm", res = 320)

diff_plot

dev.off()



## Difference in trends between regions by species ----

trend_pass_out <- readRDS("/data/outputs/trend_pass_out.rds")

trend_pass_out_df <- do.call(rbind, trend_pass_out) %>%
  dplyr::select(iteration, species, tax_grp, agri, change)

spp_diff_pass_summary <- trend_pass_out_df %>%
  spread(agri, change) %>%
  rowwise() %>%
  mutate(high_low = high - low,
         high_no_agri = high - no_agri,
         low_no_agri = low - no_agri) %>%
  dplyr::select(-c(high, low, no_agri)) %>%
  group_by(tax_grp, species) %>%
  summarise(mean_highlow = mean(high_low, na.rm = TRUE),
            lowCI_highlow = HDInterval::hdi(high_low,
                                            credMass = 0.95)[[1]],
            uppCI_highlow = HDInterval::hdi(high_low,
                                            credMass = 0.95)[[2]],
            mean_highnoagri = mean(high_no_agri, na.rm = TRUE),
            lowCI_highnoagri = HDInterval::hdi(high_no_agri,
                                               credMass = 0.95)[[1]],
            uppCI_highnoagri = HDInterval::hdi(high_no_agri,
                                               credMass = 0.95)[[2]],
            mean_lownoagri = mean(low_no_agri, na.rm = TRUE),
            lowCI_lownoagri = HDInterval::hdi(low_no_agri,
                                              credMass = 0.95)[[1]],
            uppCI_lownoagri = HDInterval::hdi(low_no_agri,
                                              credMass = 0.95)[[2]]) %>% 
  group_by(tax_grp, species) %>% 
  gather(key = diff, value = value, mean_highlow:uppCI_lownoagri) %>%
  mutate(metric = sub("_.*", "", diff),
         regions = sub(".*_", "", diff)) %>%
  dplyr::select(-diff) %>% 
  spread(key = metric, value = value) %>%
  ungroup()



spp_diff_plot <- ggplot(spp_diff_pass_summary %>%
                          filter(regions == "highlow") %>%
                          na.omit() %>%
                          arrange(tax_grp, mean) %>%
                          mutate(order = row_number())) +
  geom_point(aes(x = order, y = mean, color = mean), size = .5) +
  geom_errorbar(aes(x = order, ymin = lowCI,
                    ymax = uppCI,
                    color = mean),
                alpha = .4) +
  geom_hline(yintercept = 0, linetype = 4) +
  scale_color_viridis_c(name = "Effect size", option = "C") +
  facet_wrap(~tax_grp, scales = "free",
             nrow = 3, ncol = 2) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size = 10))

#create common x and y labels

y_grob <- textGrob("Difference in annual growth rates between high and low cropland", 
                   gp=gpar(fontface="bold", col="black", fontsize=11), rot=90)

x_grob <- textGrob("Species", 
                   gp=gpar(fontface="bold", col="black", fontsize=11))

#add to plot and save

png("/data/plots/species_diff_plot.png", 
    height = 14, width = 16, units = "cm", res = 320)

grid.arrange(arrangeGrob(spp_diff_plot, 
                         left = y_grob, bottom = x_grob))

dev.off()


