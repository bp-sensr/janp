elev_trend <- function(species) {

  # 1. Abundance-weighted mean elevation per year
  wme <- plot_dat |>
    left_join(janp_main |> select(location, elevation) |> distinct(), by = "location") |>
    group_by(location) |>
    mutate(elevation = first(na.omit(elevation))) |>
    ungroup() |>
    filter(species_code == species,
           !is.na(elevation),
           lambda_hat > 0) |>
    group_by(year) |>
    summarise(
      wt_mean_elev = sum(lambda_hat * elevation) / sum(lambda_hat),
      .groups = "drop"
    ) |>
    mutate(year_c = year - mean(year))

  # 2. Linear regression
  fit <- lm(wt_mean_elev ~ year_c, data = wme)

  # 3. Plot
  p <- ggplot(wme, aes(x = year, y = wt_mean_elev)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title    = paste0(species, " — abundance-weighted mean elevation over time"),
      subtitle = paste0("Slope = ",
                        round(coef(fit)["year_c"], 2),
                        " m yr⁻¹  |  p = ",
                        round(tidy(fit)$p.value[2], 3)),
      x = "Year", y = "Weighted mean elevation (m)"
    ) +
    theme_bw()

  list(data = wme, fit = fit, plot = p)
}

# Usage
result <- elev_trend("DEJU")
result$plot
tidy(result$fit)
glance(result$fit)
