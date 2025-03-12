# Plotting soil water retention curves

library(tidyverse)
library(udunits2)

# Read in csv
wrc <- read_csv("data_clean/rwc-swp.csv")
str(wrc)

# Make quick plot in raw units
ggplot(wrc) +
  geom_point(aes(x = rwc, y = swp,
                 color = population))

# Convert to Van Genuchten units
wrc <- wrc |> 
  mutate(psi_hm = ud.convert(swp, "MPa", "cmH2O") |> 
           ud.convert("cm", "m"))

# Plot along Van Genuchten axes
ggplot(wrc) +
  geom_point(aes(x = -psi_hm, y = rwc,
                 color = population)) +
  scale_x_log10(name = "Pressure Head (-m)",
                breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0.01, 0.1, 1, 10, 100, 1000))
