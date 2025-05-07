library(glmmTMB)
library(readxl)
library(tidyverse)
read_xlsx("data/cmon_voorraden_20240724.xlsx", sheet = "volledig") |>
  select(
    location = "Bodemlocaties_Naam", date = "Bodemlocaties_Datum",
    landuse = "Landgebruikscategorie Cmon", stock = "OrgC voorraad_0_30"
  ) |>
  filter(!is.na(.data$stock)) |>
  mutate(
    date = as.Date(.data$date),
    landuse2 = ifelse(
      str_detect(.data$landuse, "Ruimtebeslag"), "Ruimtebeslag", .data$landuse
    )
  ) -> organic_stock_0_30
ggplot(organic_stock_0_30, aes(x = landuse, y = stock)) +
  geom_boxplot() +
  coord_flip()
ggplot(organic_stock_0_30, aes(x = landuse2, y = stock)) +
  geom_boxplot() +
  coord_flip()

ggplot(organic_stock_0_30, aes(x = landuse2, y = stock)) +
  geom_violin() +
  geom_jitter(alpha = 0.2, position = position_jitter(width = 0.1)) +
  coord_flip()

summary(organic_stock_0_30)
organic_stock_0_30 |>
  count(landuse)

glmmTMB(
  stock ~ 0 + landuse2 + diag(0 + landuse2 | location),
  data = organic_stock_0_30, family = lognormal(link = "log")
) |>
  summary()
