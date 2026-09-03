library(terra)
library(sf)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# --- Landcover legend ---
legend <- c(
  "0"   = "Unclassified",
  "20"  = "Water",
  "31"  = "Snow_Ice",
  "32"  = "Rock_Rubble",
  "33"  = "Exposed_Barren_Land",
  "40"  = "Bryoids",
  "50"  = "Shrubs",
  "80"  = "Wetland",
  "81"  = "Wetland-Treed",
  "100" = "Herbs",
  "210" = "Coniferous",
  "220" = "Broadleaf",
  "230" = "Mixedwood"
)

# Colors roughly matching typical VLCE2 land cover palettes
legend_colors <- c(
  "Unclassified"        = "grey80",
  "Water"                = "#1f78b4",
  "Snow_Ice"             = "#ffffff",
  "Rock_Rubble"          = "#b2b2b2",
  "Exposed_Barren_Land"  = "#8C510A",
  "Bryoids"              = "#c2b280",
  "Shrubs"               = "#F16913",
  "Wetland"              = "#66c2a5",
  "Wetland-Treed"        = "#238b45",
  "Herbs"                = "#a6d854",
  "Coniferous"           = "#014f2a",
  "Broadleaf"            = "#7fbf7b",
  "Mixedwood"            = "#41ab5d"
)

# --- File paths for each year ---
years <- 2007:2022

tif_paths <- setNames(
  paste0("./geospatial_assets/CA_For_", years, "/CA_forest_VLCE2_", years, ".tif"),
  as.character(years)
)

# --- Points of interest (lat/lon, WGS84) ---
pts_df <- janp_main |>
  select(location, ecoregion, latitude, longitude) |> distinct() |>
  rename(site = location, lat = latitude, lon = longitude) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# --- Function: buffer + mask + return raster as data.frame for one site/year ---
get_buffer_df <- function(tif_path, pt_geom, buffer_m, year_label, site_label) {

  r <- rast(tif_path)
  pt_proj <- st_transform(pt_geom, crs(r))
  pt_buff <- st_buffer(pt_proj, dist = buffer_m)
  pt_vect <- vect(pt_buff)
  r_crop <- crop(r, pt_vect)
  r_mask <- mask(r_crop, pt_vect)
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "code"
  df$code <- as.character(df$code)
  df$class <- recode(df$code, !!!legend)
  df$year <- year_label
  df$site <- site_label

  list(df = df, buff_sf = st_sf(site = site_label, year = year_label, geometry = pt_buff))
}

# --- Loop over every site x year combination ---
combos <- expand_grid(site = pts_df$site, year = names(tif_paths))

all_results <- pmap(combos, function(site, year) {
  pt_geom <- pts_df$geometry[pts_df$site == site]
  get_buffer_df(
    tif_path   = tif_paths[[year]],
    pt_geom    = pt_geom,
    buffer_m   = 150,
    year_label = year,
    site_label = site
  )
})

# Combine raster cell data and buffer outlines
plot_df   <- map_dfr(all_results, "df")
buff_all  <- map_dfr(all_results, "buff_sf")
native_crs <- crs(rast(tif_paths[1]))

pct_df <- plot_df %>%
  count(site, year, class) %>%
  group_by(site, year) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

prop_df <- plot_df |>
  count(site, year, class) |>
  group_by(site, year) |>
  mutate(proportion = n / sum(n)) |>
  ungroup() |>
  left_join(janp_main |> select(location, ecoregion) |> distinct(), by = c("site" = "location"))

prop_df_year <- prop_df |>
  group_by(year, class, ecoregion) |>
  summarise(n = sum(n), .groups = "drop") |>
  group_by(year) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()

all_classes <- unique(plot_df$class)

prop_df_complete <- prop_df |>
  complete(site, year, class = all_classes, fill = list(n = 0, proportion = 0))

pdf <- prop_df_complete |>
  left_join(janp_main |> select(location) |> distinct() |> drop_na(), by = c("site" = "location")) |>
  group_by(site, year) |>
  mutate(total = sum(n)) |>
  ungroup() |>
  mutate(year = factor(year)) |>
  filter(!n == 0) |>
  select(-n, -total, -ecoregion)

Methods

To assess whether landcover composition changed between 2007 and 2022, we conducted a permutational multivariate analysis of variance (PERMANOVA) using Bray–Curtis dissimilarities calculated from the proportional cover of each landcover class at each site. Elevation was included as a continuous predictor together with sampling year and their interaction (year × elevation). Because each site was sampled in both years, site was treated as a blocking factor (strata = site) to constrain permutations within sites, thereby accounting for the paired sampling design. Due to the small number of paired observations, the complete set of 1,023 possible permutations was evaluated. Community differences were visualized using principal coordinates analysis (PCoA) based on the Bray–Curtis dissimilarity matrix, with trajectories connecting repeated observations from the same site.

Results

Landcover composition differed significantly as a function of year and elevation (PERMANOVA: F = 10.70, R² = 0.667, P = 0.002), with the model explaining 66.7% of the total variation in community composition. Sequential tests indicated that elevation explained the largest proportion of variation in landcover composition (R² = 0.616, F = 29.62, P = 0.002), while year explained a smaller but significant component of variation (R² = 0.053, F = 2.57, P = 0.002). The interaction between year and elevation was not significant (F = -0.093, P = 1.00), indicating that the magnitude and direction of compositional change between 2007 and 2022 did not vary systematically along the elevational gradient. The PCoA ordination supported these results. Sites were primarily separated along the first ordination axis according to differences in landcover composition associated with elevation, while paired observations from 2007 and 2022 exhibited relatively short, consistent shifts in ordination space. Although many sites changed in composition over time, there was no clear evidence that high- and low-elevation sites followed different compositional trajectories, consistent with the non-significant year × elevation interaction.

Discussion

Landcover composition changed significantly between 2007 and 2022, demonstrating that measurable landscape change occurred across the study area over the 15-year period. As expected, elevation was the dominant driver of landcover composition, reflecting the strong environmental gradients that structure vegetation communities in Jasper National Park. However, the absence of a significant year × elevation interaction suggests that compositional changes were not systematically greater at either lower or higher elevations. In other words, while landcover changed through time, the overall pattern of change was relatively consistent across the elevational gradient.

These results provide important context for subsequent analyses examining temporal changes in bird communities. A primary objective of this study is to determine whether vegetation changes, such as conifer forest expansion at lower elevations or shrub expansion near treeline, have contributed to observed changes in avian occupancy and abundance. The PERMANOVA indicates that vegetation composition has shifted over time but does not identify which landcover classes were responsible for these changes. Consequently, these results alone do not provide evidence for either forest encroachment or alpine shrubbification.

To address these ecological hypotheses directly, subsequent analyses will quantify temporal changes in individual landcover classes using generalized linear mixed models. Specifically, changes in conifer cover will be evaluated to test for forest expansion at lower elevations, while changes in shrub cover will be examined as evidence of alpine shrubbification. These models will determine whether changes in specific vegetation types vary with elevation and will provide more ecologically interpretable predictors for modelling changes in bird distributions.

Taken together, the PERMANOVA establishes that landcover composition has changed over time, while the class-specific analyses will determine whether those changes correspond to the vegetation transitions hypothesized to influence bird communities. This hierarchical approach distinguishes broad community-level change from changes in individual habitat components, providing a stronger ecological basis for linking vegetation dynamics to subsequent analyses of avian community change.

