# ============================================================
# STEP SELECTION FUNCTION ANALYSIS — APPENDIX SCRIPT
# Roe Deer SSF Methodology
Gráinne Margaret Donohue
# ------------------------------------------------------------
# Study design:
#   Two contrasting 3-day periods in July 2024
#   Heat event:  July 20-22
#   Cool period: July 4-6
#
# Model structure (fit_issf via amt package):
#   case_ ~ temp_scaled * canopy_scaled
#          + temp_scaled * hli_scaled
#          + temp_scaled * water_scaled
#          + temp_scaled * soil_scaled
#          + strata(step_id_)
#
# Covariates:
#   - Weather station temperature (Klotten A, joined by timestamp)
#   - Canopy density (Copernicus 10m)
#   - Heat Load Index (derived from DEM via spatialEco)
#   - Distance to water (OSM water bodies and waterways)
#   - Soil moisture (SLU dataset)
#
# Required input files (place in working directory subfolders):
#   Data/
#     	raw GPS fixes, hot period with buffer
#    	raw GPS fixes, cool period with buffer
#     	hourly weather station data
#     	logger GPS coordinates
#     	individual logger CSV files
#   Rasters/
#     	Copernicus canopy raster (ETRS LAEA)
#   	Digital elevation model
#    	SLU soil moisture raster
#   	OSM shapefiles for water
#
# Output files are saved to Outputs/ and Rasters/ subfolders
# ============================================================


# ============================================================
# SECTION 1 — LIBRARIES AND SETUP
# ============================================================

library(amt)
library(terra)
library(tidyverse)
library(survival)
library(lubridate)
library(sf)
library(data.table)
library(spatialEco)   # for hli()
library(Hmisc)        # for rcorr()
library(corrplot)     # for correlogram plots
library(car)          # for vif()
library(patchwork)    # for combining plots
library(adehabitatHR) # for KDE home ranges
library(automap)      # for autofitVariogram in temperature regression
library(gstat)        # for krige() in temperature regression

# Set working directory — update this path to your own machine
## setwd(add your path here)

# Create output folders if they don't exist
dir.create("Outputs", showWarnings = FALSE)
dir.create("Rasters", showWarnings = FALSE)
dir.create("Data",    showWarnings = FALSE)


# ============================================================
# SECTION 2 — GPS DATA CLEANING
# Source: V1 script
# Raw GPS files contain two mixed timestamp formats and
# coordinates in WGS84 that must be reprojected to SWEREF99
# ============================================================

# ------------------------------------------------------------
# 2a — CLEAN HOT PERIOD GPS DATA
# ------------------------------------------------------------

gps_hot_raw <- read.csv("Data/gps_data_2024_hot.csv",
                        header           = TRUE,
                        stringsAsFactors = FALSE)

# The raw file has an extra empty column between timestamp and x
colnames(gps_hot_raw) <- c("id", "timestamp", "time_extra", "x", "y")

# Remove rows with missing coordinates
gps_hot_raw <- gps_hot_raw %>%
  filter(!is.na(x), x != "", !is.na(y), y != "")

# Two timestamp formats exist in the file:
# Format 1 — DD/MM/YYYY with time in a separate column
# Format 2 — ISO 8601: 2024-07-22T12:01:21.000Z
format2_rows <- grepl("T", gps_hot_raw$timestamp)

gps_hot_f1 <- gps_hot_raw[!format2_rows, ] %>%
  mutate(timestamp = as.POSIXct(
    paste(timestamp, time_extra),
    format = "%d/%m/%Y %H:%M:%S",
    tz = "UTC"))

gps_hot_f2 <- gps_hot_raw[format2_rows, ] %>%
  mutate(timestamp = as.POSIXct(
    timestamp,
    format = "%Y-%m-%dT%H:%M:%OS",
    tz = "UTC"))

gps_hot_combined <- bind_rows(
  gps_hot_f1 %>% select(id, timestamp, x, y),
  gps_hot_f2 %>% select(id, timestamp, x, y)
) %>% arrange(id, timestamp)

# Convert coordinates to numeric (fix decimal separator if needed)
gps_hot_combined$x <- as.numeric(gsub(",", ".", gps_hot_combined$x))
gps_hot_combined$y <- as.numeric(gsub(",", ".", gps_hot_combined$y))

# The raw coordinates are stored as latitude/longitude (WGS84)
# with x and y columns swapped — swap them back before reprojecting
gps_hot_combined <- gps_hot_combined %>%
  rename(latitude = x, longitude = y)

gps_hot_combined <- gps_hot_combined %>%
  select(id, timestamp, x, y) %>%
  filter(!is.na(timestamp)) %>%
  distinct(id, timestamp, .keep_all = TRUE) %>%
  arrange(id, timestamp)

write.csv(gps_hot_combined,
          "Data/gps_data_2024_hot_clean.csv",
          row.names = FALSE)

cat("Hot period GPS cleaned:", nrow(gps_hot_combined), "fixes,",
    n_distinct(gps_hot_combined$id), "individuals\n")

# ------------------------------------------------------------
# 2b — CLEAN COOL PERIOD GPS DATA
# Same procedure as hot period
# ------------------------------------------------------------

gps_cool_raw <- read.csv("Data/gps_data_2024_cool.csv",
                         header           = TRUE,
                         stringsAsFactors = FALSE)

colnames(gps_cool_raw) <- c("id", "timestamp", "time_extra", "x", "y")

gps_cool_raw <- gps_cool_raw %>%
  filter(!is.na(x), x != "", !is.na(y), y != "")

format2_rows <- grepl("T", gps_cool_raw$timestamp)

gps_cool_f1 <- gps_cool_raw[!format2_rows, ] %>%
  mutate(timestamp = as.POSIXct(
    paste(timestamp, time_extra),
    format = "%d/%m/%Y %H:%M:%S",
    tz = "UTC"))

gps_cool_f2 <- gps_cool_raw[format2_rows, ] %>%
  mutate(timestamp = as.POSIXct(
    timestamp,
    format = "%Y-%m-%dT%H:%M:%OS",
    tz = "UTC"))

gps_cool_combined <- bind_rows(
  gps_cool_f1 %>% select(id, timestamp, x, y),
  gps_cool_f2 %>% select(id, timestamp, x, y)
) %>% arrange(id, timestamp)

gps_cool_combined$x <- as.numeric(gsub(",", ".", gps_cool_combined$x))
gps_cool_combined$y <- as.numeric(gsub(",", ".", gps_cool_combined$y))

gps_cool_combined <- gps_cool_combined %>%
  rename(latitude = x, longitude = y)

gps_cool_combined <- gps_cool_combined %>%
  select(id, timestamp, x, y) %>%
  filter(!is.na(timestamp)) %>%
  distinct(id, timestamp, .keep_all = TRUE) %>%
  arrange(id, timestamp)

write.csv(gps_cool_combined,
          "Data/gps_data_2024_cool_clean.csv",
          row.names = FALSE)

cat("Cool period GPS cleaned:", nrow(gps_cool_combined), "fixes,",
    n_distinct(gps_cool_combined$id), "individuals\n")


# ============================================================
# SECTION 3 — LOGGER DATA EXTRACTION AND COORDINATE JOINING
# 39 microclimate loggers stored as individual semicolon-
# separated CSV files. Logger ID extracted from filename.
# Coordinates stored in a separate file and joined by ID.
# ============================================================

# ------------------------------------------------------------
# 3a — READ ALL LOGGER FILES
# ------------------------------------------------------------

logger_folder <- "Data/logger_data"
logger_files  <- list.files(path     = logger_folder,
                            pattern  = "*.csv",
                            full.names = TRUE)

cat("Logger files found:", length(logger_files), "\n")

# Extract logger ID from filename
extract_logger_id <- function(filepath) {
  parts     <- strsplit(basename(filepath), "_")[[1]]
  clean_pos <- which(tolower(parts) == "clean")
  if (length(clean_pos) == 0) return(NA)
  return(parts[clean_pos + 1])
}

# Read all logger files
read_logger_file <- function(filepath) {
  logger_id <- extract_logger_id(filepath)
  df <- tryCatch(
    read.csv(filepath, sep = ",", header = TRUE,
             stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)
  df$logger_id <- logger_id
  return(df)
}

logger_raw <- map_dfr(logger_files, read_logger_file)

cat("Total logger rows loaded:", nrow(logger_raw), "\n")
cat("Loggers found:", n_distinct(logger_raw$logger_id), "\n")

# ------------------------------------------------------------
# 3b — CLEAN AND AGGREGATE TO HOURLY MEANS
# Column mapping (from raw logger files):
#   T1 = soil temperature, T2 = 2cm, T3 = 15cm,
#   T4 = air temperature at 170cm (used in V1 interpolation),
#   VWC = volumetric water content (soil moisture)
# ------------------------------------------------------------

logger_raw <- logger_raw %>%
  rename(
    timestamp     = DT,
    soil_temp     = T1,
    temp_2cm      = T2,
    temp_15cm     = T3,
    air_temp      = T4,
    soil_moisture = VWC
  ) %>%
  mutate(
    timestamp     = as.POSIXct(timestamp,
                               format = "%Y-%m-%d %H:%M:%S",
                               tz = "UTC"),
    soil_temp     = as.numeric(soil_temp),
    air_temp      = as.numeric(air_temp),
    soil_moisture = as.numeric(soil_moisture)
  )

# Aggregate from 15-minute to hourly — covers both study periods with one-day buffer on each side
logger_hourly <- logger_raw %>%
  filter(timestamp >= as.POSIXct("2024-07-03 00:00:00", tz = "UTC"),
         timestamp <= as.POSIXct("2024-07-23 23:59:59", tz = "UTC")) %>%
  mutate(hour_timestamp = floor_date(timestamp, unit = "hour")) %>%
  group_by(logger_id, hour_timestamp) %>%
  summarise(
    air_temp      = mean(air_temp,      na.rm = TRUE),
    soil_temp     = mean(soil_temp,     na.rm = TRUE),
    soil_moisture = mean(soil_moisture, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(timestamp = hour_timestamp)

cat("Hourly logger rows:", nrow(logger_hourly), "\n")
cat("Loggers with data:", n_distinct(logger_hourly$logger_id), "\n")

# ------------------------------------------------------------
# 3c — JOIN LOGGER COORDINATES
# Logger coordinates are stored separately in WGS84 and must be reprojected to SWEREF99
# ------------------------------------------------------------

logger_coords <- read.csv("Data/Grimso_Loggers.csv", sep = ";") %>%
  select(Name, Longitude, Latitude) %>%
  rename(logger_id = Name,
         longitude  = Longitude,
         latitude   = Latitude)

# Remove leading zeros from logger IDs to match file-extracted IDs
logger_coords$logger_id <- sub("^0+", "",
                               as.character(logger_coords$logger_id))

logger_coords_sf <- st_as_sf(logger_coords,
                             coords = c("longitude", "latitude"),
                             crs    = 4326)
logger_coords_sf <- st_transform(logger_coords_sf, crs = 3006)

logger_coords$x <- st_coordinates(logger_coords_sf)[, 1]
logger_coords$y <- st_coordinates(logger_coords_sf)[, 2]

# Join coordinates to hourly logger data
logger_hourly <- logger_hourly %>%
  left_join(logger_coords %>% select(logger_id, x, y),
            by = "logger_id")

write.csv(logger_hourly,
          "Data/logger_data_july2024_with_coords.csv",
          row.names = FALSE)

cat("Logger data with coordinates saved!\n")
cat("Loggers with coordinates:",
    sum(!is.na(logger_hourly$x)) / nrow(logger_hourly) * 100,
    "% matched\n")
# ============================================================
# SECTION 4 — MICROCLIMATE TEMPERATURE REGRESSION
# Response variable: maximum logger air temperature (15cm above ground) across July 2024.
# Predictors: canopy density, HLI, elevation, soil moisture — extracted from raster layers at each logger location.
# Coefficients applied across rasters to produce a spatially continuous predicted temperature surface for the study area.

# LIMITATION: elevation produced unrealistic predictions (up to 38.5°C, exceeding the observed logger maximum of 30.5°C) in high-elevation, low-canopy pixels. 
# DEPENDENCY: this section uses canopy_r, hli_r, soil_r and dem_crop built in Section 8. Run Section 8 first, or load the rasters manually using the commented lines in Section 4b.
# ============================================================
# ------------------------------------------------------------
# 4a — BUILD LOGGER SPATIAL OBJECT
# Maximum air temperature per logger across July 2024
# aggregated from logger_hourly (built in Section 3)
# ------------------------------------------------------------

logger_means <- logger_hourly %>%
  group_by(logger_id, x, y) %>%
  summarise(
    max_air_temp       = max(air_temp,       na.rm = TRUE),
    mean_soil_moisture = mean(soil_moisture, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(x), !is.na(y),
         !is.na(max_air_temp))

logger_sf_interp <- st_as_sf(logger_means,
                              coords = c("x", "y"),
                              crs    = 3006)

cat("Logger locations for regression:", nrow(logger_sf_interp), "\n")
cat("Maximum air temperature range:",
    round(range(logger_sf_interp$max_air_temp), 2), "\n")

saveRDS(logger_sf_interp, "Data/logger_sf_interp.rds")

# ------------------------------------------------------------
# 4b — EXTRACT COVARIATE VALUES AT LOGGER LOCATIONS
# If running this section in isolation, load rasters first:
# canopy_r <- rast("Rasters/canopy_density.tif")
# hli_r    <- rast("Rasters/heat_load_index.tif")
# soil_r   <- rast("Rasters/soil_moisture_clean.tif")
# dem_crop <- rast("Rasters/dem_crop.tif")
# ------------------------------------------------------------

# Extract canopy, HLI and soil moisture from covariates stack
logger_covs <- terra::extract(covariates,
                               vect(logger_sf_interp),
                               bind = FALSE)

# Extract elevation separately — dem_crop is not in the
# covariates stack as elevation is not an SSF covariate
logger_elev <- terra::extract(dem_crop,
                               vect(logger_sf_interp),
                               bind = FALSE)

cat("Covariate values extracted at", nrow(logger_covs),
    "logger locations\n")
cat("Missing values per column:\n")
print(colSums(is.na(logger_covs)))

# Build regression dataframe
logger_coords_xy <- st_coordinates(logger_sf_interp)

reg_data <- data.frame(
  logger_id       = logger_sf_interp$logger_id,
  max_air_temp    = logger_sf_interp$max_air_temp,
  x               = logger_coords_xy[, 1],
  y               = logger_coords_xy[, 2],
  canopy_density  = logger_covs$canopy_density,
  heat_load_index = logger_covs$heat_load_index,
  soil_moisture   = logger_covs$soil_moisture,
  elevation       = logger_elev[, 2]
) %>% drop_na()

cat("Logger locations with complete data:", nrow(reg_data), "\n")
print(summary(reg_data))

# ------------------------------------------------------------
# 4c — FIT LINEAR REGRESSION
# Response:   maximum air temperature (°C) per logger
# Predictors: canopy density, HLI, elevation, soil moisture
# Replicates approach from prior Grimsö thesis —
# canopy density substituted for basal area
# ------------------------------------------------------------

cat("Fitting linear regression...\n")

temp_model <- lm(max_air_temp ~
                   canopy_density  +
                   heat_load_index +
                   elevation       +
                   soil_moisture,
                 data = reg_data)

summary(temp_model)

sink("Outputs/temperature_regression_results.txt")
cat("Temperature Regression — Grimso July 2024\n")
cat("Response: maximum air temperature at logger locations\n")
cat("Predictors: canopy density, HLI, elevation, soil moisture\n\n")
summary(temp_model)
sink()

cat("Regression results saved\n")

# ------------------------------------------------------------
# 4d — PRODUCE PREDICTED TEMPERATURE MAP
# Regression coefficients applied pixel-by-pixel across
# the raster stack to produce a landscape-wide surface.
# See LIMITATION note in section header regarding elevation.
# ------------------------------------------------------------

cat("Producing predicted temperature map...\n")

intercept   <- coef(temp_model)["(Intercept)"]
coef_canopy <- coef(temp_model)["canopy_density"]
coef_hli    <- coef(temp_model)["heat_load_index"]
coef_elev   <- coef(temp_model)["elevation"]
coef_soil   <- coef(temp_model)["soil_moisture"]

cat("Model coefficients:\n")
cat("Intercept:      ", round(intercept,   4), "\n")
cat("Canopy density: ", round(coef_canopy, 4), "\n")
cat("HLI:            ", round(coef_hli,    4), "\n")
cat("Elevation:      ", round(coef_elev,   4), "\n")
cat("Soil moisture:  ", round(coef_soil,   4), "\n")

predicted_temp <- intercept          +
  (coef_canopy * canopy_r)  +
  (coef_hli    * hli_r)     +
  (coef_elev   * dem_crop)  +
  (coef_soil   * soil_r)

names(predicted_temp) <- "predicted_temperature"

print(predicted_temp)
cat("Predicted temperature range:",
    round(range(values(predicted_temp), na.rm = TRUE), 2), "\n")

writeRaster(predicted_temp,
            "Rasters/predicted_temperature.tif",
            overwrite = TRUE)

# ------------------------------------------------------------
# 4e — PLOT PREDICTED TEMPERATURE SURFACE
# Replicates Figure 4 layout from prior Grimso thesis:
# four predictor maps alongside the predicted temperature map
# ------------------------------------------------------------

pred_temp_df <- as.data.frame(predicted_temp, xy = TRUE)
colnames(pred_temp_df) <- c("x", "y", "temperature")

logger_coords_df <- data.frame(
  x    = st_coordinates(logger_sf_interp)[, 1],
  y    = st_coordinates(logger_sf_interp)[, 2],
  temp = logger_sf_interp$max_air_temp
)

ggplot() +
  geom_raster(data = pred_temp_df,
              aes(x = x, y = y, fill = temperature)) +
  scale_fill_gradientn(
    colours = c("#313695", "#74add1", "#ffffbf",
                "#fdae61", "#d73027"),
    name    = "Temp (°C)") +
  geom_point(data = logger_coords_df,
             aes(x = x, y = y),
             size = 2, shape = 21,
             fill = "white", colour = "black", stroke = 1.2) +
  coord_equal() +
  labs(title    = "Predicted Maximum Air Temperature — Grimsö Study Area",
       subtitle = "Linear regression of maximum logger temperature against landscape covariates, July 2024",
       x        = "Easting (SWEREF99)",
       y        = "Northing (SWEREF99)") +
  theme_classic(base_size = 12, base_family = "Times New Roman")

ggsave("Outputs/predicted_temperature_map.png",
       width = 10, height = 12, dpi = 300)

# ------------------------------------------------------------
# 4f — VALIDATE — OBSERVED VS PREDICTED AT LOGGER LOCATIONS
# ------------------------------------------------------------

reg_data$predicted_temp <- predict(temp_model, newdata = reg_data)

ggplot(reg_data, aes(x = max_air_temp, y = predicted_temp)) +
  geom_point(colour = "firebrick", size = 3) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "gray50") +
  geom_smooth(method = "lm", colour = "steelblue",
              fill   = "steelblue", alpha = 0.2) +
  labs(title    = "Observed vs Predicted Temperature at Logger Locations",
       subtitle = "Dashed line shows perfect prediction (1:1)",
       x        = "Observed maximum air temperature (°C)",
       y        = "Predicted maximum air temperature (°C)") +
  theme_classic(base_size = 12, base_family = "Times New Roman")

ggsave("Outputs/observed_vs_predicted_temperature.png",
       width = 8, height = 6, dpi = 300)

cat("Predicted temperature map and validation plot saved\n")

# ============================================================
# SECTION 5 — WEATHER STATION DATA
# Klotten A weather station provides hourly air temperature, joined to GPS and SSF steps by timestamp
# ============================================================

weather_2024 <- read.csv("Data/Klotten_A_July2024.csv",
                         sep              = ";",
                         skip             = 9,
                         header           = TRUE,
                         stringsAsFactors = FALSE)

weather_2024 <- weather_2024 %>%
  select(Datum, Tid..UTC., Lufttemperatur) %>%
  rename(date        = Datum,
         time        = Tid..UTC.,
         temperature = Lufttemperatur)

weather_2024$timestamp <- as.POSIXct(
  paste(weather_2024$date, weather_2024$time),
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC")

weather_2024 <- weather_2024 %>%
  select(timestamp, temperature) %>%
  filter(!is.na(temperature), !is.na(timestamp))

# Filter to study periods with one-day buffer on each side
weather_hot <- weather_2024 %>%
  filter(timestamp >= as.POSIXct("2024-07-19 00:00:00", tz = "UTC"),
         timestamp <= as.POSIXct("2024-07-23 23:59:59", tz = "UTC"))

weather_cool <- weather_2024 %>%
  filter(timestamp >= as.POSIXct("2024-07-03 00:00:00", tz = "UTC"),
         timestamp <= as.POSIXct("2024-07-07 23:59:59", tz = "UTC"))

saveRDS(weather_2024, "Data/weather_2024.rds")

cat("Weather station data loaded\n")
cat("Hot period temperature range:",
    range(weather_hot$temperature, na.rm = TRUE), "\n")
cat("Cool period temperature range:",
    range(weather_cool$temperature, na.rm = TRUE), "\n")


# ============================================================
# SECTION 6 — PERIOD SELECTION PLOTS
# Daily temperature summaries from Klotten A and logger network, used to justify choice of heat event and cool period dates
# ============================================================

# Daily summary from weather station
daily_2024 <- weather_2024 %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(date) %>%
  summarise(
    max_temp  = max(temperature,  na.rm = TRUE),
    mean_temp = mean(temperature, na.rm = TRUE),
    min_temp  = min(temperature,  na.rm = TRUE)
  )

# Daily mean and max from logger network
logger_daily <- logger_hourly %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(date, logger_id) %>%
  summarise(daily_max  = max(air_temp,  na.rm = TRUE),
            daily_mean = mean(air_temp, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(date) %>%
  summarise(logger_max  = mean(daily_max,  na.rm = TRUE),
            logger_mean = mean(daily_mean, na.rm = TRUE),
            .groups = "drop")

daily_2024 <- daily_2024 %>%
  left_join(logger_daily, by = "date")

# Plot 1 — Mean temperature
p1 <- ggplot(daily_2024, aes(x = date)) +
  annotate("rect",
           xmin = as.Date("2024-07-20"), xmax = as.Date("2024-07-22"),
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "firebrick") +
  annotate("rect",
           xmin = as.Date("2024-07-04"), xmax = as.Date("2024-07-06"),
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "steelblue") +
  geom_hline(aes(yintercept = 15, linetype = "15°C threshold"),
             colour = "gray40") +
  geom_line(aes(y = mean_temp, colour = "Klotten A"),
            linewidth = 0.8) +
  geom_line(aes(y = logger_mean, colour = "Logger network"),
            linewidth = 0.8) +
  annotate("text", x = as.Date("2024-07-21"), y = 21,
           label = "Heat event", colour = "firebrick", size = 3) +
  annotate("text", x = as.Date("2024-07-05"), y = 21,
           label = "Cool period", colour = "steelblue", size = 3) +
  scale_colour_manual(name   = "Temperature record",
                      values = c("Klotten A"       = "purple",
                                 "Logger network"  = "darkgreen")) +
  scale_linetype_manual(name   = "Threshold",
                        values = c("15°C threshold" = "dashed")) +
  labs(title = "Daily Mean Temperature — Klotten A & Logger Network 2024",
       x = "Date", y = "Mean temperature (°C)") +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(legend.position = "right")

# Plot 2 — Max temperature
p2 <- ggplot(daily_2024, aes(x = date)) +
  annotate("rect",
           xmin = as.Date("2024-07-20"), xmax = as.Date("2024-07-22"),
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "firebrick") +
  annotate("rect",
           xmin = as.Date("2024-07-04"), xmax = as.Date("2024-07-06"),
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "steelblue") +
  geom_hline(aes(yintercept = 15, linetype = "15°C threshold"),
             colour = "gray40") +
  geom_line(aes(y = max_temp,   colour = "Klotten A"),
            linewidth = 0.8) +
  geom_line(aes(y = logger_max, colour = "Logger network"),
            linewidth = 0.8) +
  scale_colour_manual(name   = "Temperature record",
                      values = c("Klotten A"       = "purple",
                                 "Logger network"  = "darkgreen")) +
  scale_linetype_manual(name   = "Threshold",
                        values = c("15°C threshold" = "dashed")) +
  labs(title = "Daily Maximum Temperature — Klotten A & Logger Network 2024",
       x = "Date", y = "Maximum temperature (°C)") +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(legend.position = "right")

ggsave("Outputs/temperature_candidates_combined.png",
       plot  = p1 / p2,
       width = 10, height = 8, dpi = 300)

cat("Period selection plots saved\n")


# ============================================================
# SECTION 7 — RASTER PREPARATION
# All rasters projected to SWEREF99, cropped to combined GPS extent, and resampled to 10m resolution
# ============================================================

# ------------------------------------------------------------
# HEAT LOAD INDEX (derived from DEM)
# hli() from spatialEco takes DEM directly, force.hemisphere = "northern" is required for Sweden
# ------------------------------------------------------------

cat("Processing Heat Load Index...\n")
dem       <- rast("Rasters/dem.tif")
dem_sw    <- project(dem, "EPSG:3006", method = "bilinear")
dem_crop  <- crop(dem_sw, combined_extent)

hli_terra <- hli(dem_crop,
                 check            = TRUE,
                 force.hemisphere = "northern")

hli_r     <- resample(hli_terra, template_10m, method = "bilinear")
names(hli_r) <- "heat_load_index"
writeRaster(hli_r, "Rasters/heat_load_index.tif", overwrite = TRUE)
cat("HLI saved\n")


# ------------------------------------------------------------
# CONFIRM ALL RASTERS ALIGN
# ------------------------------------------------------------

covariates <- c(canopy_r, water_r, hli_r, soil_r)
names(covariates) <- c("canopy_density", "distance_to_water",
                       "heat_load_index", "soil_moisture")

cat("\nRaster stack check:\n")
print(covariates)


# ============================================================
# SECTION 8 — GPS PERIOD FILTERING
# Filter cleaned GPS data to the core 3-day periods, roedeer_05 excluded, not present in either period
# ============================================================

gps_hot$timestamp  <- as.POSIXct(gps_hot$timestamp,
                                  format = "%Y-%m-%d %H:%M:%S",
                                  tz = "UTC")
gps_cool$timestamp <- as.POSIXct(gps_cool$timestamp,
                                  format = "%Y-%m-%d %H:%M:%S",
                                  tz = "UTC")

# Core period filtering
gps_hot_final <- gps_hot %>%
  filter(timestamp >= as.POSIXct("2024-07-20 00:00:00", tz = "UTC"),
         timestamp <= as.POSIXct("2024-07-22 23:59:59", tz = "UTC"),
         id != "roedeer_05")

gps_cool_final <- gps_cool %>%
  filter(timestamp >= as.POSIXct("2024-07-04 00:00:00", tz = "UTC"),
         timestamp <= as.POSIXct("2024-07-06 23:59:59", tz = "UTC"),
         id != "roedeer_05")

cat("HOT PERIOD — core 3 days (July 20-22):\n")
cat("Fixes:", nrow(gps_hot_final),
    "| Individuals:", n_distinct(gps_hot_final$id), "\n")

cat("COOL PERIOD — core 3 days (July 4-6):\n")
cat("Fixes:", nrow(gps_cool_final),
    "| Individuals:", n_distinct(gps_cool_final$id), "\n")

# Verify no fixes fall outside the raster extent
raster_ext <- ext(covariates)

check_outside <- function(gps_df, period) {
  n_out <- sum(gps_df$x < raster_ext[1] | gps_df$x > raster_ext[2] |
                 gps_df$y < raster_ext[3] | gps_df$y > raster_ext[4])
  cat(period, "— fixes outside raster extent:", n_out, "\n")
}

check_outside(gps_hot_final,  "Hot period")
check_outside(gps_cool_final, "Cool period")


# ============================================================
# SECTION 9 — JOIN WEATHER STATION TEMPERATURE TO GPS DATA
# GPS timestamps rounded to nearest hour to match hourly weather station readings before joining
# ============================================================

join_temp <- function(gps_df) {
  gps_df %>%
    mutate(timestamp_hour = round_date(timestamp, unit = "hour")) %>%
    left_join(weather_2024 %>%
                rename(timestamp_hour = timestamp,
                       weather_temp   = temperature),
              by = "timestamp_hour")
}

gps_hot_final  <- join_temp(gps_hot_final)
gps_cool_final <- join_temp(gps_cool_final)

cat("Temperature join — hot period missing:",
    sum(is.na(gps_hot_final$weather_temp)), "\n")
cat("Temperature join — cool period missing:",
    sum(is.na(gps_cool_final$weather_temp)), "\n")


# ============================================================
# SECTION 10 — MOVEMENT TRACK AND STEP GENERATION
# Tracks resampled to 1-hour fix rate (tolerance ±10 min)
# Bursts with fewer than 3 fixes removed
# Individuals with fewer than 10 fixes after resampling removed
# 20 random steps generated per real step
# ============================================================

# Helper to build track, resample, and generate steps
build_ssf <- function(gps_df, period_label) {

  cat("\nBuilding", period_label, "track...\n")

  track <- mk_track(gps_df,
                    .x  = x,
                    .y  = y,
                    .t  = timestamp,
                    id  = id,
                    crs = 3006)

  track_resampled <- track %>%
    nest(data = -id) %>%
    mutate(
      data = map(data, ~ track_resample(
        .x,
        rate      = hours(1),
        tolerance = minutes(10)
      )),
      data = map(data, filter_min_n_burst, min_n = 3)
    )

  # Report fixes per individual after resampling
  track_resampled %>%
    mutate(n_fixes = map_int(data, nrow)) %>%
    select(id, n_fixes) %>%
    print(n = Inf)

  # Remove individuals with too few fixes to fit a model
  track_resampled <- track_resampled %>%
    mutate(n_fixes = map_int(data, nrow)) %>%
    filter(n_fixes > 10) %>%
    select(-n_fixes)

  cat("Individuals retained:", nrow(track_resampled), "\n")

  # Generate real steps and 20 random alternatives per step
  ssf <- track_resampled %>%
    mutate(steps = map(data, ~ steps_by_burst(.x) %>%
                         random_steps(n_control = 20))) %>%
    unnest(cols = steps)

  cat("Steps generated:", nrow(ssf),
      "| Real:", sum(ssf$case_),
      "| Random:", sum(!ssf$case_), "\n")

  return(ssf)
}

ssf_hot  <- build_ssf(gps_hot_final,  "hot period")
ssf_cool <- build_ssf(gps_cool_final, "cool period")


# ============================================================
# SECTION 11 — COVARIATE EXTRACTION AND SCALING
# Spatial covariates extracted at step endpoints
# Weather station temperature joined by step end time
# ============================================================

# Extract spatial covariates at step endpoints
cat("Extracting spatial covariates — hot period...\n")
ssf_hot <- ssf_hot %>%
  extract_covariates(covariates) %>%
  drop_na(canopy_density, distance_to_water,
          heat_load_index, soil_moisture)

cat("Extracting spatial covariates — cool period...\n")
ssf_cool <- ssf_cool %>%
  extract_covariates(covariates) %>%
  drop_na(canopy_density, distance_to_water,
          heat_load_index, soil_moisture)

cat("Hot period rows after extraction:", nrow(ssf_hot), "\n")
cat("Cool period rows after extraction:", nrow(ssf_cool), "\n")

# Join weather station temperature to steps
# Use step end time (t2_) rounded to nearest hour
join_temp_to_steps <- function(ssf_df) {
  ssf_df %>%
    mutate(timestamp_hour = round_date(t2_, unit = "hour")) %>%
    left_join(weather_2024 %>%
                rename(timestamp_hour = timestamp),
              by = "timestamp_hour") %>%
    select(-timestamp_hour)
}

ssf_hot  <- join_temp_to_steps(ssf_hot)
ssf_cool <- join_temp_to_steps(ssf_cool)

cat("Temperature join to steps — hot missing:",
    sum(is.na(ssf_hot$temperature)), "\n")
cat("Temperature join to steps — cool missing:",
    sum(is.na(ssf_cool$temperature)), "\n")

# Scale all covariates
# Scaling parameters are saved as attributes on each column
# and retrieved later for the habitat suitability maps
ssf_hot <- ssf_hot %>%
  mutate(
    temp_scaled   = scale(temperature),
    canopy_scaled = scale(canopy_density),
    water_scaled  = scale(distance_to_water),
    hli_scaled    = scale(heat_load_index),
    soil_scaled   = scale(soil_moisture)
  )

ssf_cool <- ssf_cool %>%
  mutate(
    temp_scaled   = scale(temperature),
    canopy_scaled = scale(canopy_density),
    water_scaled  = scale(distance_to_water),
    hli_scaled    = scale(heat_load_index),
    soil_scaled   = scale(soil_moisture)
  )

cat("Covariates scaled for both periods\n")


# ============================================================
# SECTION 12 — MULTICOLLINEARITY CHECKS
# Pearson correlation matrix and VIF both run on real steps
# Correlations above 0.7 flagged: none found
# ============================================================

# ------------------------------------------------------------
# PEARSON CORRELATION MATRIX
# ------------------------------------------------------------

run_correlation <- function(ssf_data, period_name) {

  cat("\n============================\n")
  cat("CORRELATION —", period_name, "\n")
  cat("============================\n")

  cov_matrix <- ssf_data %>%
    filter(case_ == TRUE) %>%
    select(temperature, canopy_density,
           heat_load_index, distance_to_water,
           soil_moisture) %>%
    as.matrix()

  cor_results <- rcorr(cov_matrix, type = "pearson")
  cor_r        <- round(cor_results$r, 2)
  cor_p        <- round(cor_results$P, 3)

  cat("\nPearson r matrix:\n"); print(cor_r)
  cat("\nP-value matrix:\n");   print(cor_p)

  high_cor <- which(abs(cor_r) > 0.7 & cor_r != 1, arr.ind = TRUE)
  if (nrow(high_cor) == 0) {
    cat("\nNo correlations above |0.7| — no problematic collinearity\n")
  } else {
    cat("\nCorrelations above |0.7|:\n")
    for (i in 1:nrow(high_cor)) {
      cat(rownames(cor_r)[high_cor[i,1]], "vs",
          colnames(cor_r)[high_cor[i,2]], ":",
          cor_r[high_cor[i,1], high_cor[i,2]], "\n")
    }
  }
  return(list(r = cor_r, p = cor_p))
}

cor_hot  <- run_correlation(ssf_hot,  "HOT PERIOD July 20-22")
cor_cool <- run_correlation(ssf_cool, "COOL PERIOD July 4-6")

write.csv(cor_hot$r,  "Outputs/correlation_r_hot.csv")
write.csv(cor_hot$p,  "Outputs/correlation_p_hot.csv")
write.csv(cor_cool$r, "Outputs/correlation_r_cool.csv")
write.csv(cor_cool$p, "Outputs/correlation_p_cool.csv")

# Correlograms
png("Outputs/correlogram_combined.png",
    width = 2800, height = 1400, res = 300)
par(mfrow = c(1, 2))
corrplot(cor_hot$r,
         method = "color", type = "upper", addCoef.col = "black",
         tl.col = "black", tl.srt = 45,
         title  = "Hot Period (Jul 20-22)",
         mar    = c(0, 0, 2, 0))
corrplot(cor_cool$r,
         method = "color", type = "upper", addCoef.col = "black",
         tl.col = "black", tl.srt = 45,
         title  = "Cool Period (Jul 4-6)",
         mar    = c(0, 0, 2, 0))
par(mfrow = c(1, 1))
dev.off()

cat("Correlograms saved\n")

# ------------------------------------------------------------
# VARIANCE INFLATION FACTORS
# ------------------------------------------------------------

vif_hot <- vif(glm(case_ ~
                     canopy_scaled + hli_scaled +
                     water_scaled  + soil_scaled,
                   data   = ssf_hot %>% filter(case_ == TRUE),
                   family = binomial))

vif_cool <- vif(glm(case_ ~
                      canopy_scaled + hli_scaled +
                      water_scaled  + soil_scaled,
                    data   = ssf_cool %>% filter(case_ == TRUE),
                    family = binomial))

sink("Outputs/vif_results.txt")
cat("VIF Results — Hot Period:\n");  print(round(vif_hot,  3))
cat("\nVIF Results — Cool Period:\n"); print(round(vif_cool, 3))
cat("\nAll values < 5 indicate acceptable collinearity\n")
sink()

cat("VIF check complete\n")


# ============================================================
# SECTION 13 — MODEL FITTING
# Conditional logistic regression via fit_issf (amt)
# One model per period, strata defined by step_id_
# Interaction terms allow temperature to modify each covariate
# ============================================================

cat("Fitting hot period SSF model...\n")
model_hot <- ssf_hot %>%
  fit_issf(
    case_ ~
      temp_scaled * canopy_scaled +
      temp_scaled * hli_scaled    +
      temp_scaled * water_scaled  +
      temp_scaled * soil_scaled   +
      strata(step_id_),
    model = TRUE
  )

cat("Fitting cool period SSF model...\n")
model_cool <- ssf_cool %>%
  fit_issf(
    case_ ~
      temp_scaled * canopy_scaled +
      temp_scaled * hli_scaled    +
      temp_scaled * water_scaled  +
      temp_scaled * soil_scaled   +
      strata(step_id_),
    model = TRUE
  )

summary(model_hot)
summary(model_cool)

# Save models and data for downstream use
saveRDS(model_hot,  "Data/model_hot.rds")
saveRDS(model_cool, "Data/model_cool.rds")
saveRDS(ssf_hot,    "Data/ssf_hot.rds")
saveRDS(ssf_cool,   "Data/ssf_cool.rds")

sink("Outputs/model_hot_results.txt")
cat("SSF Model — Hot Period July 20-22 2024\n\n")
summary(model_hot)
sink()

sink("Outputs/model_cool_results.txt")
cat("SSF Model — Cool Period July 4-6 2024\n\n")
summary(model_cool)
sink()

cat("Both models fitted and saved\n")


# ============================================================
# SECTION 14 — COEFFICIENT COMPARISON PLOTS
# Coefficients and 95% CI extracted for hot and cool periods
# Plotted as forest plot — filled points = significant
# ============================================================

extract_coefs <- function(model, period) {
  data.frame(
    covariate = names(coef(model$model)),
    estimate  = coef(model$model),
    ci_low    = confint(model$model)[, 1],
    ci_high   = confint(model$model)[, 2],
    period    = period
  ) %>% filter(!is.na(estimate))
}

coefs_combined <- bind_rows(
  extract_coefs(model_hot,  "Hot period"),
  extract_coefs(model_cool, "Cool period")
) %>%
  mutate(
    significant = ifelse(ci_low > 0 | ci_high < 0,
                         "Significant", "Not significant"),
    covariate_clean = dplyr::recode(
      covariate,
      "temp_scaled"               = "Temperature",
      "canopy_scaled"             = "Canopy density",
      "hli_scaled"                = "Heat load index",
      "water_scaled"              = "Distance to water",
      "soil_scaled"               = "Soil moisture",
      "temp_scaled:canopy_scaled" = "Temp x Canopy",
      "temp_scaled:hli_scaled"    = "Temp x HLI",
      "temp_scaled:water_scaled"  = "Temp x Water",
      "temp_scaled:soil_scaled"   = "Temp x Soil moisture"
    )
  )

# All coefficients
ggplot(coefs_combined,
       aes(x      = reorder(covariate_clean, estimate),
           y      = estimate,
           colour = period,
           shape  = significant)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.3, linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
  scale_colour_manual(values = c("Hot period"  = "firebrick",
                                 "Cool period" = "steelblue"),
                      name = "Period") +
  scale_shape_manual(values = c("Significant"     = 16,
                                "Not significant" = 1),
                     name = "Significance") +
  coord_flip() +
  labs(title = "SSF Coefficients — Hot vs Cool Period",
       x     = "Covariate",
       y     = "Coefficient estimate") +
  theme_classic(base_size = 12, base_family = "Times New Roman")

ggsave("Outputs/coefficients_hot_vs_cool.png",
       width = 10, height = 7, dpi = 300)

# Main effects only
ggplot(coefs_combined %>% filter(!grepl(":", covariate)),
       aes(x      = reorder(covariate_clean, estimate),
           y      = estimate,
           colour = period,
           shape  = significant)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.3, linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
  scale_colour_manual(values = c("Hot period"  = "firebrick",
                                 "Cool period" = "steelblue"),
                      name = "Period") +
  scale_shape_manual(values = c("Significant"     = 16,
                                "Not significant" = 1),
                     name = "Significance") +
  coord_flip() +
  labs(title    = "SSF Main Effects — Hot vs Cool Period",
       subtitle = "Filled points = significant at p < 0.05",
       x        = "Covariate",
       y        = "Coefficient estimate") +
  theme_classic(base_size = 12, base_family = "Times New Roman")

ggsave("Outputs/coefficients_main_effects.png",
       width = 9, height = 5, dpi = 300)

# Interaction terms only
ggplot(coefs_combined %>% filter(grepl(":", covariate)),
       aes(x      = reorder(covariate_clean, estimate),
           y      = estimate,
           colour = period,
           shape  = significant)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.3, linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
  scale_colour_manual(values = c("Hot period"  = "firebrick",
                                 "Cool period" = "steelblue"),
                      name = "Period") +
  scale_shape_manual(values = c("Significant"     = 16,
                                "Not significant" = 1),
                     name = "Significance") +
  coord_flip() +
  labs(title = "SSF Interaction Terms — Hot vs Cool Period",
       x     = "Interaction term",
       y     = "Coefficient estimate") +
  theme_classic(base_size = 12, base_family = "Times New Roman")

ggsave("Outputs/coefficients_interactions.png",
       width = 9, height = 5, dpi = 300)

cat("Coefficient plots saved\n")

# Summary table
summary_table <- coefs_combined %>%
  select(covariate_clean, period, estimate, ci_low, ci_high) %>%
  mutate(
    estimate = round(estimate, 3),
    ci_low   = round(ci_low,   3),
    ci_high  = round(ci_high,  3),
    ci       = paste0("(", ci_low, ", ", ci_high, ")")
  ) %>%
  select(covariate_clean, period, estimate, ci) %>%
  pivot_wider(names_from  = period,
              values_from = c(estimate, ci))

write.csv(summary_table, "Outputs/model_summary_table.csv",
          row.names = FALSE)

cat("Summary table saved\n")


# ============================================================
# SECTION 15 — RSS PLOTS
# Log-RSS plotted across the observed range of each covariate
# All other covariates held at their mean (scaled = 0)
# Reference location x2 set to mean of all covariates
# ============================================================

# Helper — build x1 dataframe for one covariate across its range
make_x1 <- function(var_name, var_values) {
  n  <- length(var_values)
  df <- data.frame(
    temp_scaled   = rep(0, n),
    canopy_scaled = rep(0, n),
    hli_scaled    = rep(0, n),
    water_scaled  = rep(0, n),
    soil_scaled   = rep(0, n)
  )
  df[[var_name]] <- var_values
  return(df)
}

x2_base <- data.frame(temp_scaled = 0, canopy_scaled = 0,
                      hli_scaled  = 0, water_scaled  = 0,
                      soil_scaled = 0)

# Function to produce and save one RSS plot
rss_plot <- function(model, ssf_data, var_name, var_label,
                     colour, period_label) {
  x1  <- make_x1(var_name,
                 seq(min(ssf_data[[var_name]], na.rm = TRUE),
                     max(ssf_data[[var_name]], na.rm = TRUE),
                     length.out = 100))
  rss <- log_rss(model, x1, x2_base, ci = "se")
  x_col <- paste0(var_name, "_x1")

  ggplot(rss$df, aes(x = .data[[x_col]], y = log_rss)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                alpha = 0.2, fill = colour) +
    geom_line(colour = colour, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "gray50") +
    labs(title    = paste("Log-RSS —", var_label),
         subtitle = period_label,
         x        = paste(var_label, "(scaled)"),
         y        = "Log RSS") +
    theme_classic(base_size = 12, base_family = "Times New Roman")
}

covs <- list(
  list(var = "temp_scaled",   label = "Temperature",       colour = "firebrick"),
  list(var = "canopy_scaled", label = "Canopy density",    colour = "darkgreen"),
  list(var = "hli_scaled",    label = "Heat load index",   colour = "orange"),
  list(var = "water_scaled",  label = "Distance to water", colour = "steelblue"),
  list(var = "soil_scaled",   label = "Soil moisture",     colour = "brown")
)

for (cov in covs) {
  p_hot  <- rss_plot(model_hot,  ssf_hot,
                     cov$var, cov$label, cov$colour, "Hot period July 20-22")
  p_cool <- rss_plot(model_cool, ssf_cool,
                     cov$var, cov$label, cov$colour, "Cool period July 4-6")
  combined <- p_hot + p_cool
  fname <- paste0("Outputs/rss_", gsub("_scaled", "", cov$var),
                  "_combined.png")
  ggsave(fname, combined, width = 14, height = 5, dpi = 300)
}

cat("RSS plots saved\n")


# ============================================================
# SECTION 16 — RSS PLOTS (canopy x temperature)
# Following Alston et al. — three RSS curves per period,
# each at a different temperature percentile (10th, 50th, 90th)
# Shows how canopy selection changes across the temperature range
# Reference x2 set to mean of all covariates
# ============================================================

alston_rss <- function(model, ssf_data, period_label) {

  # Temperature percentiles (10th, 50th, 90th)
  temp_pcts <- quantile(ssf_data$temp_scaled,
                        c(0.10, 0.50, 0.90), na.rm = TRUE)

  # Convert back to Celsius for plot labels
  temp_mean <- attr(ssf_data$temp_scaled, "scaled:center")
  temp_sd   <- attr(ssf_data$temp_scaled, "scaled:scale")
  temp_c    <- round(temp_pcts * temp_sd + temp_mean, 1)

  cat(period_label, "temperature percentiles (°C):\n")
  cat("10th:", temp_c[1], " 50th:", temp_c[2], " 90th:", temp_c[3], "\n")

  canopy_range <- seq(min(ssf_data$canopy_scaled, na.rm = TRUE),
                      max(ssf_data$canopy_scaled, na.rm = TRUE),
                      length.out = 100)

  x2 <- data.frame(canopy_scaled = 0, temp_scaled = 0,
                   hli_scaled = 0, water_scaled = 0, soil_scaled = 0)

  make_x1_alston <- function(temp_val) {
    data.frame(canopy_scaled = canopy_range,
               temp_scaled   = temp_val,
               hli_scaled    = 0,
               water_scaled  = 0,
               soil_scaled   = 0)
  }

  labels <- c(
    paste0("10th pct (", temp_c[1], "°C)"),
    paste0("50th pct (", temp_c[2], "°C)"),
    paste0("90th pct (", temp_c[3], "°C)")
  )

  rss_df <- map2_dfr(as.list(temp_pcts), labels, function(tp, lbl) {
    rss <- log_rss(model, make_x1_alston(tp), x2, ci = "se")
    rss$df %>%
      mutate(
        temp_level = lbl,
        rss        = exp(log_rss),
        rss_lwr    = exp(lwr),
        rss_upr    = exp(upr)
      )
  }) %>%
    mutate(temp_level = factor(temp_level, levels = labels))

  ggplot(rss_df,
         aes(x = canopy_scaled_x1, y = rss,
             colour = temp_level, fill = temp_level)) +
    geom_ribbon(aes(ymin = rss_lwr, ymax = rss_upr),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dotted",
               colour = "gray50") +
    scale_colour_manual(values = c("steelblue", "gray50", "firebrick"),
                        name = "Temperature") +
    scale_fill_manual(values = c("steelblue", "gray50", "firebrick"),
                      name = "Temperature") +
    labs(title    = paste("Canopy Selection —", period_label),
         subtitle = "RSS curves at 10th, 50th and 90th temperature percentiles",
         x        = "Canopy density (scaled)",
         y        = "Relative Selection Strength (RSS)") +
    theme_classic(base_size = 12, base_family = "Times New Roman")
}

p_alston_hot  <- alston_rss(model_hot,  ssf_hot,  "Hot period July 20-22")
p_alston_cool <- alston_rss(model_cool, ssf_cool, "Cool period July 4-6")

ggsave("Outputs/rss_alston_canopy_combined.png",
       p_alston_hot + p_alston_cool +
         plot_annotation(
           title    = "Canopy Selection Across Temperature Gradient",
           subtitle = "Roe deer — Grimsö July 2024") +
         plot_layout(guides = "collect"),
       width = 14, height = 6, dpi = 300)

cat("Alston-style RSS plots saved\n")


# ============================================================
# SECTION 17 — STEP LENGTH ANALYSIS
# Compares step lengths between hot and cool periods and tests for a relationship between step length and temperature
# ============================================================

real_steps_hot  <- ssf_hot  %>% filter(case_ == TRUE) %>%
  mutate(period = "Hot (July 20-22)")
real_steps_cool <- ssf_cool %>% filter(case_ == TRUE) %>%
  mutate(period = "Cool (July 4-6)")
real_steps_combined <- bind_rows(real_steps_hot, real_steps_cool)

cat("Hot period step lengths:\n");  print(summary(real_steps_hot$sl_))
cat("Cool period step lengths:\n"); print(summary(real_steps_cool$sl_))

step_ttest <- t.test(real_steps_hot$sl_, real_steps_cool$sl_)
cat("\nT-test hot vs cool step lengths:\n"); print(step_ttest)

# Step length distribution
ggplot(real_steps_combined,
       aes(x = sl_, fill = period, colour = period)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("Hot (July 20-22)"  = "firebrick",
                               "Cool (July 4-6)" = "steelblue"),
                    name = "Period") +
  scale_colour_manual(values = c("Hot (July 20-22)"  = "firebrick",
                                 "Cool (July 4-6)" = "steelblue"),
                      name = "Period") +
  labs(title = "Step Length Distribution — Hot vs Cool Period",
       x     = "Step length (metres)",
       y     = "Density") +
  theme_classic(base_size = 12, base_family = "Times New Roman")

ggsave("Outputs/steplength_distribution.png",
       width = 8, height = 5, dpi = 300)

# Step length vs temperature scatter
ggplot(real_steps_combined,
       aes(x = temperature, y = sl_, colour = period)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_colour_manual(values = c("Hot (July 20-22)"  = "firebrick",
                                 "Cool (July 4-6)" = "steelblue"),
                      name = "Period") +
  labs(title = "Step Length vs Temperature",
       x     = "Temperature (°C)",
       y     = "Step length (metres)") +
  theme_classic(base_size = 12, base_family = "Times New Roman")

ggsave("Outputs/steplength_vs_temperature.png",
       width = 8, height = 5, dpi = 300)

cat("Step length analysis complete\n")


# ============================================================
# SECTION 18 — INDIVIDUAL-LEVEL MODELS
# Separate SSF fitted per individual to assess consistency of selection direction across the population
# Temperature excluded from individual models (too few steps per individual for interaction terms to be estimable)
# ============================================================

run_individual_models <- function(ssf_data, period_name) {

  cat("\n", period_name, "individual models:\n")

  models <- ssf_data %>%
    group_by(id) %>%
    nest() %>%
    mutate(model = map(data, ~ tryCatch(
      fit_issf(
        case_ ~
          canopy_scaled +
          hli_scaled    +
          water_scaled  +
          soil_scaled   +
          strata(step_id_),
        data  = .x,
        model = TRUE
      ),
      error = function(e) NULL
    )))

  cat("Models fitted:",
      sum(!sapply(models$model, is.null)), "of", nrow(models), "\n")

  coefs <- models %>%
    filter(!sapply(model, is.null)) %>%
    mutate(
      canopy_coef = map_dbl(model, ~ coef(.x$model)["canopy_scaled"]),
      hli_coef    = map_dbl(model, ~ coef(.x$model)["hli_scaled"]),
      water_coef  = map_dbl(model, ~ coef(.x$model)["water_scaled"]),
      soil_coef   = map_dbl(model, ~ coef(.x$model)["soil_scaled"])
    ) %>%
    select(id, canopy_coef, hli_coef, water_coef, soil_coef)

  print(coefs)

  cat("Selecting for canopy:",
      round(mean(coefs$canopy_coef > 0) * 100, 1), "%\n")

  return(coefs)
}

ind_coefs_hot  <- run_individual_models(ssf_hot,  "HOT PERIOD")
ind_coefs_cool <- run_individual_models(ssf_cool, "COOL PERIOD")

# Plot individual canopy coefficients
ind_coefs_hot$period  <- "Hot (July 20-22)"
ind_coefs_cool$period <- "Cool (July 4-6)"
ind_coefs_combined    <- bind_rows(ind_coefs_hot, ind_coefs_cool)

ind_coefs_long <- ind_coefs_combined %>%
  pivot_longer(cols      = c(canopy_coef, hli_coef, water_coef, soil_coef),
               names_to  = "covariate",
               values_to = "estimate") %>%
  mutate(covariate = dplyr::recode(
    covariate,
    "canopy_coef" = "Canopy density",
    "hli_coef"    = "Heat load index",
    "water_coef"  = "Distance to water",
    "soil_coef"   = "Soil moisture"))

ggplot(ind_coefs_long,
       aes(x = id, y = estimate, fill = period, colour = period)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
  scale_fill_manual(values = c("Hot (July 20-22)"  = "firebrick",
                               "Cool (July 4-6)" = "steelblue"),
                    name = "Period") +
  scale_colour_manual(values = c("Hot (July 20-22)"  = "firebrick",
                                 "Cool (July 4-6)" = "steelblue"),
                      name = "Period") +
  facet_wrap(~ covariate, scales = "free_y", ncol = 2) +
  coord_flip() +
  labs(title = "Individual Coefficients — Hot vs Cool Period",
       x     = "Individual",
       y     = "Coefficient estimate") +
  theme_classic() +
  theme(strip.text = element_text(size = 9))

ggsave("Outputs/individual_all_covariates.png",
       width = 12, height = 8, dpi = 300)

cat("Individual models complete\n")


# ============================================================
# SECTION 19 — HABITAT SUITABILITY MAP 
# Model coefficients applied to scaled raster layers to produce a continuous habitat suitability surface.
# Scaling parameters extracted from the fitted SSF data to ensure the rasters are on the same scale as the model.
# Main effects only — temperature set to period mean (scaled = 0)
# ============================================================

# Extract period-specific scaling parameters
get_scale <- function(ssf_data, var) {
  list(
    mean = attr(ssf_data[[var]], "scaled:center"),
    sd   = attr(ssf_data[[var]], "scaled:scale")
  )
}

scale_raster <- function(r, sc) (r - sc$mean) / sc$sd

# Hot period
canopy_s_hot <- scale_raster(canopy_r, get_scale(ssf_hot, "canopy_scaled"))
hli_s_hot    <- scale_raster(hli_r,    get_scale(ssf_hot, "hli_scaled"))
water_s_hot  <- scale_raster(water_r,  get_scale(ssf_hot, "water_scaled"))
soil_s_hot   <- scale_raster(soil_r,   get_scale(ssf_hot, "soil_scaled"))

# Cool period
canopy_s_cool <- scale_raster(canopy_r, get_scale(ssf_cool, "canopy_scaled"))
hli_s_cool    <- scale_raster(hli_r,    get_scale(ssf_cool, "hli_scaled"))
water_s_cool  <- scale_raster(water_r,  get_scale(ssf_cool, "water_scaled"))
soil_s_cool   <- scale_raster(soil_r,   get_scale(ssf_cool, "soil_scaled"))

# Extract main effect coefficients
get_coef <- function(model, name) coef(model$model)[name]

hab_hot <- (get_coef(model_hot, "canopy_scaled") * canopy_s_hot) +
  (get_coef(model_hot, "hli_scaled")    * hli_s_hot)    +
  (get_coef(model_hot, "water_scaled")  * water_s_hot)  +
  (get_coef(model_hot, "soil_scaled")   * soil_s_hot)
names(hab_hot) <- "suitability_hot"

hab_cool <- (get_coef(model_cool, "canopy_scaled") * canopy_s_cool) +
  (get_coef(model_cool, "hli_scaled")    * hli_s_cool)    +
  (get_coef(model_cool, "water_scaled")  * water_s_cool)  +
  (get_coef(model_cool, "soil_scaled")   * soil_s_cool)
names(hab_cool) <- "suitability_cool"

writeRaster(hab_hot,  "Rasters/habitat_suitability_hot.tif",  overwrite = TRUE)
writeRaster(hab_cool, "Rasters/habitat_suitability_cool.tif", overwrite = TRUE)

cat("Habitat suitability rasters saved for ArcGIS\n")

# Plot both maps side by side
hab_hot_df  <- as.data.frame(hab_hot,  xy = TRUE) %>%
  drop_na() %>% rename(suitability = suitability_hot)
hab_cool_df <- as.data.frame(hab_cool, xy = TRUE) %>%
  drop_na() %>% rename(suitability = suitability_cool)

shared_min <- min(hab_hot_df$suitability, hab_cool_df$suitability)
shared_max <- max(hab_hot_df$suitability, hab_cool_df$suitability)

hab_map <- function(df, period_label) {
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = suitability)) +
    scale_fill_gradient2(low = "steelblue", mid = "lightyellow",
                         high = "firebrick", midpoint = 0,
                         limits = c(shared_min, shared_max),
                         name   = "Suitability") +
    coord_equal() +
    labs(title = period_label,
         x = "Easting (SWEREF99)", y = "Northing (SWEREF99)") +
    theme_classic()
}

ggsave("Outputs/habitat_suitability_hot_vs_cool.png",
       hab_map(hab_hot_df,  "Hot period — July 20-22") +
         hab_map(hab_cool_df, "Cool period — July 4-6") +
         plot_annotation(
           title = "Predicted Habitat Suitability — Hot vs Cool Period",
           subtitle = "Roe deer — Grimsö July 2024"),
       width = 14, height = 12, dpi = 300)

cat("Habitat suitability map saved\n")


# ============================================================
# SECTION 20 — KDE HOME RANGES
# 95% kernel density estimate per individual per period using adehabitatHR::kernelUD with href bandwidth
# Plotted over habitat suitability surfaces
# Shapefiles exported
# ============================================================

gps_hot_sf  <- st_as_sf(
  gps_hot %>%
    mutate(timestamp = as.POSIXct(timestamp,
                                  format = "%Y-%m-%d %H:%M:%S",
                                  tz = "UTC")) %>%
    filter(timestamp >= as.POSIXct("2024-07-20 00:00:00", tz = "UTC"),
           timestamp <= as.POSIXct("2024-07-22 23:59:59", tz = "UTC")),
  coords = c("x", "y"), crs = 3006)

gps_cool_sf <- st_as_sf(
  gps_cool %>%
    mutate(timestamp = as.POSIXct(timestamp,
                                  format = "%Y-%m-%d %H:%M:%S",
                                  tz = "UTC")) %>%
    filter(timestamp >= as.POSIXct("2024-07-04 00:00:00", tz = "UTC"),
           timestamp <= as.POSIXct("2024-07-06 23:59:59", tz = "UTC")),
  coords = c("x", "y"), crs = 3006)

# KDE function — returns sf polygon per individual at 95% contour
calc_kde <- function(gps_sf, period_name) {

  cat("\nCalculating KDE for", period_name, "...\n")

  coords <- st_coordinates(gps_sf)
  sp_pts <- SpatialPointsDataFrame(
    coords      = coords,
    data        = data.frame(id = gps_sf$id),
    proj4string = CRS("+init=epsg:3006")
  )

  kde    <- kernelUD(sp_pts[, "id"], h = "href", extent = 2)
  kde_95 <- getverticeshr(kde, percent = 95)

  # Use st_transform to correctly set CRS (not st_set_crs)
  kde_sf <- st_as_sf(kde_95) %>%
    st_transform(3006) %>%
    rename(individual = id) %>%
    mutate(period  = period_name,
           area_ha = round(as.numeric(st_area(.)) / 10000, 1))

  cat("Individuals:", nrow(kde_sf), "\n")
  print(kde_sf %>% select(individual, area_ha) %>% st_drop_geometry())

  return(kde_sf)
}

kde_hot  <- calc_kde(gps_hot_sf,  "Hot period")
kde_cool <- calc_kde(gps_cool_sf, "Cool period")

# Colour palette for individuals
individual_colours <- c(
  "roedeer_02" = "#E41A1C",
  "roedeer_03" = "#377EB8",
  "roedeer_04" = "#4DAF4A",
  "roedeer_06" = "#984EA3",
  "roedeer_07" = "#FF7F00",
  "roedeer_08" = "#A65628",
  "roedeer_09" = "#F781BF"
)

# Habitat suitability + individual KDE polygons
kde_hab_map <- function(hab_df, kde_sf, period_label) {
  ggplot() +
    geom_raster(data = hab_df, aes(x = x, y = y, fill = suitability)) +
    scale_fill_gradient2(low = "steelblue", mid = "lightyellow",
                         high = "firebrick", midpoint = 0,
                         limits = c(shared_min, shared_max),
                         name   = "Suitability") +
    geom_sf(data = kde_sf, aes(colour = individual),
            fill = NA, linewidth = 0.8, inherit.aes = FALSE) +
    scale_colour_manual(values = individual_colours, name = "Individual") +
    coord_sf() +
    labs(title = period_label,
         subtitle = "95% KDE home ranges per individual",
         x = "Easting (SWEREF99)", y = "Northing (SWEREF99)") +
    theme_classic()
}

ggsave("Outputs/habitat_suitability_kde_individual.png",
       kde_hab_map(hab_hot_df,  kde_hot,  "Hot period — July 20-22") +
         kde_hab_map(hab_cool_df, kde_cool, "Cool period — July 4-6") +
         plot_layout(guides = "collect"),
       width = 16, height = 12, dpi = 300)

# Export KDE polygons as shapefiles for ArcGIS
st_write(kde_hot,  "Outputs/kde_individual_hot.shp",  delete_dsn = TRUE)
st_write(kde_cool, "Outputs/kde_individual_cool.shp", delete_dsn = TRUE)
st_write(st_union(kde_hot)  %>% st_sf() %>% mutate(period = "Hot"),
         "Outputs/kde_merged_hot.shp",  delete_dsn = TRUE)
st_write(st_union(kde_cool) %>% st_sf() %>% mutate(period = "Cool"),
         "Outputs/kde_merged_cool.shp", delete_dsn = TRUE)

cat("KDE maps and shapefiles saved\n")

# ============================================================
# SECTION 21 — ZONAL STATISTICS
# Home range overlap with predicted habitat suitability
# Calculates % of each KDE home range in positive, neutral, and negative suitability terrain for both periods
# ============================================================

# Load habitat suitability rasters
suit_hot  <- rast("Rasters/habitat_suitability_hot.tif")
suit_cool <- rast("Rasters/habitat_suitability_cool.tif")

# Ensure CRS matches
hr_hot_sf  <- st_transform(hr_hot_sf,  crs = 3006)
hr_cool_sf <- st_transform(hr_cool_sf, crs = 3006)

# Function to calculate suitability overlap statistics
# for each individual home range
calc_suitability_overlap <- function(hr_sf, suit_rast, period_name) {
  results <- map_dfr(1:nrow(hr_sf), function(i) {
    hr_single <- hr_sf[i, ]
    vals <- terra::extract(suit_rast, vect(hr_single))[[2]]
    vals <- vals[!is.na(vals)]
    total <- length(vals)
    data.frame(
      period       = period_name,
      individual   = hr_sf$id[i],
      n_pixels     = total,
      pct_positive = round(sum(vals > 0)   / total * 100, 1),
      pct_neutral  = round(sum(vals >= -0.2 & vals <= 0.2) / total * 100, 1),
      pct_negative = round(sum(vals < 0)   / total * 100, 1),
      mean_suit    = round(mean(vals), 3)
    )
  })
  return(results)
}

# Run for both periods
overlap_hot  <- calc_suitability_overlap(hr_hot_sf,  suit_hot,  "Hot period")
overlap_cool <- calc_suitability_overlap(hr_cool_sf, suit_cool, "Cool period")

# Combine results
overlap_all <- bind_rows(overlap_hot, overlap_cool)

# Rename columns for clean table output
overlap_table <- overlap_all %>%
  rename(
    Period             = period,
    Individual         = individual,
    `Pixels (n)`       = n_pixels,
    `% Positive`       = pct_positive,
    `% Neutral`        = pct_neutral,
    `% Negative`       = pct_negative,
    `Mean suitability` = mean_suit
  )

# Save outputs
write.csv(overlap_table,
          "Outputs/home_range_suitability_overlap_formatted.csv",
          row.names = FALSE)

cat("Zonal statistics saved!\n")
cat("Outputs/home_range_suitability_overlap_formatted.csv\n")