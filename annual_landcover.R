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
