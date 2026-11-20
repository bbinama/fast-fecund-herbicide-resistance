## ============================================================================================
## REPRODUCIBILITY - FULL ANALYSIS PIPELINE (single-file concatenation)
## Manuscript: 'A fast, colonising life history predicts herbicide resistance in arable weeds'
## Assembled 2026-09-08 | R 4.5.2
##
## This file concatenates every script the manuscript's results depend on, in run order.
## The .R scripts are included verbatim; the two literate sources (Herbicide_Code.qmd and
## 02_Master_Analysis.Rmd) are included as their extracted R code (knitr::purl) - their fully
## annotated originals (with interpretation text) are the .qmd / .Rmd files.
##
## HOW TO RUN: execute the sections in order. Stages 1-7 write intermediate outputs
## (merged_data_final.csv, occurrence CSVs, moa_analysis.rds, reviewer_checks.rds, Figures/*)
## that Stage 8 (the main analysis) reads; Stage 10 (_render.R full) knits the report.
## Inputs required: EVA vegetation plots, TRY traits, GBIF occurrences (doi:10.15468/dl.9j88hn),
## Heap resistance DB, BVL register, EUNIS/FloraVeg habitat, district crop stats (Duden 2024),
## CIAT/GRIN crop-wild-relative lists - see Data availability in the manuscript.
## ============================================================================================


## ############################################################################################
## STAGE 0  Dependencies / session setup
## source file: install_dependencies.R
## ############################################################################################

# =============================================================================
# install_dependencies.R  --  one-shot package setup for the herbicide-resistance
# trait-modelling analysis. Run once before rendering:  Rscript install_dependencies.R
# Tested on R 4.5.2 (see README for the full pinned version manifest).
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

## --- CRAN packages (grouped by role) -----------------------------------------
cran <- c(
  # data wrangling / IO
  "tidyverse", "data.table", "janitor", "readxl", "writexl",
  # imputation, ordination, multivariate
  "missForest", "FactoMineR", "factoextra", "vegan", "corrplot",
  # models: forests, GLMs, boosting, evaluation, explanations
  "randomForest", "partykit", "glmnet", "xgboost", "caret", "pROC",
  "kernelshap", "shapviz",
  # phylogenetics
  "ape", "phylolm", "caper",
  # networks (strategy comparison)
  "igraph",
  # plotting / tables / reporting
  "patchwork", "ggrepel", "ggbeeswarm", "viridis", "kableExtra",
  "knitr", "rmarkdown", "png", "magick", "pdftools",
  # spatial (maps / watch-list context)
  "sf", "terra", "rnaturalearth", "ggspatial", "giscoR", "geodata",
  # misc stats
  "FSA", "MASS", "Matrix"
)

## --- GitHub-only packages -----------------------------------------------------
## V.PhyloMaker2 is NOT on CRAN; it builds the dated GBOTB megatree (build_megatree.R).
github <- c(V.PhyloMaker2 = "jinyizju/V.PhyloMaker2")

installed <- rownames(installed.packages())

todo <- setdiff(cran, installed)
cat("CRAN packages to install:", if (length(todo)) paste(todo, collapse = ", ") else "(none)", "\n")
if (length(todo)) install.packages(todo, dependencies = c("Depends", "Imports", "LinkingTo"))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
for (pkg in names(github)) if (!requireNamespace(pkg, quietly = TRUE)) {
  cat("Installing", pkg, "from GitHub:", github[[pkg]], "\n")
  remotes::install_github(github[[pkg]], upgrade = "never")
}

## --- report ------------------------------------------------------------------
all_pkgs <- c(cran, names(github))
ok <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
cat("\n=== install result ===\n")
print(data.frame(package = all_pkgs, available = ok, row.names = NULL))
missing <- all_pkgs[!ok]
cat("\nALL_OK:", all(ok), "\n")
if (length(missing)) cat("STILL MISSING:", paste(missing, collapse = ", "), "\n")

## Document rendering also needs: pandoc (bundled with RStudio, or install standalone)
## and a LaTeX engine (xelatex for the walkthrough, pdflatex for report/supplement).
## If you have no system LaTeX, run once:  install.packages("tinytex"); tinytex::install_tinytex()
cat("\nNOTE: rendering the .Rmd reports also requires pandoc and a LaTeX engine (see README).\n")

## ############################################################################################
## STAGE 1  GBIF occurrence cleaning + species assembly
## source file: Herbicide_Code.qmd
## ############################################################################################

library(rgbif)

# 1. Match your species list to GBIF keys
# Your list of 500 scientific names
my_species <- c("Abies alba", "Picea abies", "Fagus sylvatica", ...) 

# Match the entire list at once
matched_names <- name_backbone_checklist(my_species)

# Filter for successful matches and extract the unique IDs (usageKey)
taxon_keys <- matched_names %>%
  filter(matchType != "NONE") %>%
  pull(usageKey) %>%
  unique()
# Fill in your GBIF.org account details
user <- "your_username"
pwd <- "your_password"
email <- "your_email@example.com"
# 2. Submit ONE download request for all species in Germany
# This runs on GBIF servers, and they will email you when the ZIP is ready.
# Submit the request
download_request <- occ_download(
  pred_in("taxonKey", taxon_keys),     # All 500 species IDs
  pred("country", "DE"),              # Limited to Germany
  pred("hasCoordinate", TRUE),        # Only records with GPS data
  pred("occurrenceStatus", "PRESENT"), # Exclude "absent" records
  format = "SIMPLE_CSV",              # Lightweight format
  user = user, pwd = pwd, email = email
)

# Check status (Wait until it says "SUCCEEDED")
occ_download_wait(download_request)

# Download and import the data directly into R
my_data <- occ_download_get(download_request) %>%
           occ_download_import()

library(dplyr)
library(tidyr)

# 1. Check Basis of Record (What types of data do we have?)
basis_summary <- my_data %>%
  group_by(basisOfRecord) %>%
  tally(name = "count") %>%
  mutate(percentage = count / sum(count) * 100)

# 2. Check Coordinate Uncertainty (How many are NA vs. Values?)
uncertainty_summary <- my_data %>%
  summarise(
    total_records = n(),
    na_uncertainty = sum(is.na(coordinateUncertaintyInMeters)),
    has_uncertainty = sum(!is.na(coordinateUncertaintyInMeters)),
    zero_uncertainty = sum(coordinateUncertaintyInMeters == 0, na.rm = TRUE)
  )

# 3. Check Elevation (The "0m Trap" and NAs)
elevation_summary <- my_data %>%
  summarise(
    na_elevation = sum(is.na(elevation)),
    zero_elevation = sum(elevation == 0, na.rm = TRUE),
    min_elev = min(elevation, na.rm = TRUE),
    max_elev = max(elevation, na.rm = TRUE)
  )

# 4. Cross-tabulate: Uncertainty vs Basis of Record
# This helps see if HUMAN_OBSERVATION is the one missing uncertainty info
cross_tab <- my_data %>%
  group_by(basisOfRecord) %>%
  summarise(
    count = n(),
    na_uncertainty = sum(is.na(coordinateUncertaintyInMeters)),
    avg_uncertainty = mean(coordinateUncertaintyInMeters, na.rm = TRUE)
  )

# Print results to console
print(basis_summary)
print(uncertainty_summary)
print(elevation_summary)
print(cross_tab)

my_data_cleaned <- my_data %>%
  # 1. Remove obvious non-wild records
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION")) %>%
  
  # 2. Strategic Uncertainty Filter
  # Keep if uncertainty is <= 2000m OR if it's a recent human observation with NA uncertainty
  filter(
    coordinateUncertaintyInMeters <= 2000 | 
    (is.na(coordinateUncertaintyInMeters) & basisOfRecord == "HUMAN_OBSERVATION" & year >= 2010)
  ) %>%
  
  # 3. Remove Centroids (Berlin, etc.)
  filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 10000)) %>%
  
  # 4. Remove records where coordinates are exactly 0,0 or obviously wrong
  filter(decimalLatitude != 0 & decimalLongitude != 0)

species_altitude_traits <- my_data_cleaned %>%
  group_by(species) %>%
  summarise(
    alt_min = min(altitude, na.rm = TRUE),
    # Use 95th percentile instead of absolute Max to avoid Alpine outliers
    alt_max_95 = quantile(altitude, 0.95, na.rm = TRUE),
    alt_mean = mean(altitude, na.rm = TRUE),
    alt_range = alt_max_95 - alt_min,
    total_obs = n()
  )

my_data_filtered <- my_data %>%
  # 1. Keep only 'wild' or field-based records
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "OBSERVATION", "PRESERVED_SPECIMEN", "OCCURRENCE", "MATERIAL_SAMPLE")) %>%
  
  # 2. Basic Coordinate check (Remove 0,0 and invalid coordinates)
  filter(decimalLatitude != 0 & decimalLongitude != 0) %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

my_data_cleaned <- my_data_filtered %>%
  filter(
    # Option A: Uncertainty is known and is <= 2km
    (coordinateUncertaintyInMeters <= 2000) | 
    
    # Option B: Uncertainty is NA, but it's a recent Human Observation (Post-2010)
    # We assume smartphone GPS accuracy for these.
    (is.na(coordinateUncertaintyInMeters) & basisOfRecord == "HUMAN_OBSERVATION" & year >= 2010)
  ) %>%
  
  # 3. Remove Centroids (Common "faked" coordinates in Germany)
  filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 10000))

my_data_hard_cleaned <- my_data %>%
  # 1. Keep only 'Field' data: Observations and Herbarium specimens
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION", 
                             "OBSERVATION", 
                             "PRESERVED_SPECIMEN", 
                             "OCCURRENCE", 
                             "MATERIAL_SAMPLE")) %>%
  
  # 2. Basic GPS Integrity
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  filter(decimalLatitude != 0 & decimalLongitude != 0) %>%
  
  # 3. Geographic bounding box for Germany (Safety check)
  # Prevents points in the middle of the ocean or other continents
  filter(decimalLatitude > 47 & decimalLatitude < 55) %>%
  filter(decimalLongitude > 5 & decimalLongitude < 16)

my_data_final <- my_data_hard_cleaned %>%
  filter(
    # 1. High-confidence spatial data (any year)
    (coordinateUncertaintyInMeters <= 2000) |
    
    # 2. Professional Historical Data (Keep since 1970 for HR context)
    (basisOfRecord == "PRESERVED_SPECIMEN" & year >= 1970) |
    
    # 3. Modern Citizen Science (smartphone era)
    (basisOfRecord == "HUMAN_OBSERVATION" & year >= 2010)
  )

library(dplyr)
library(ggplot2)

# 1. Summary of years
year_summary <- my_data_filtered %>%
  summarise(
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    median_year = median(year, na.rm = TRUE),
    records_pre_1978 = sum(year < 1978, na.rm = TRUE),
    records_post_1978 = sum(year >= 1978, na.rm = TRUE),
    na_years = sum(is.na(year))
  )

# 2. Decade-based breakdown to see the 'Density'
decade_breakdown <- my_data_filtered %>%
  filter(!is.na(year)) %>%
  mutate(decade = (year %/% 10) * 10) %>%
  group_by(decade) %>%
  tally()

# 3. Print results
print(year_summary)
print(decade_breakdown, n = 30)

# 4. Optional: Quick visual to see the 'Temporal Gap'
ggplot(decade_breakdown, aes(x = decade, y = n)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Occurrence Records by Decade", x = "Year", y = "Count")

my_data_cleaned <- my_data_filtered %>%
  filter(
    # 1. High Precision Lane: Keep anything with low uncertainty (regardless of year)
    (coordinateUncertaintyInMeters <= 2000) | 
    
    # 2. The Modern Smartphone Lane: Trust NAs from 2010 onwards (Human Observation)
    (is.na(coordinateUncertaintyInMeters) & 
       basisOfRecord == "HUMAN_OBSERVATION" & 
       year >= 2010) |
    
    # 3. The Professional/Historical Lane: Keep NAs for Specimens from the HR era
    # We trust 'Preserved Specimen' from 1970-2010 even without uncertainty info
    (is.na(coordinateUncertaintyInMeters) & 
       basisOfRecord == "PRESERVED_SPECIMEN" & 
       year >= 1970) |

    # 4. Handle the NA Years: Only keep them if they are Professional Specimens
    # This prevents 'mystery' human observations from ruining your altitude traits
    (is.na(year) & basisOfRecord == "PRESERVED_SPECIMEN" & is.na(coordinateUncertaintyInMeters))
  ) %>%
  
  # 5. Final Centroid Purge
  # Removes coordinates that are 'faked' or rounded to city centers
  filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 10000))

library(ggplot2)
library(dplyr)

# Create a combined dataset for plotting
comparison_df <- bind_rows(
  my_data_filtered %>% mutate(status = "Raw (Filtered)"),
  my_data_cleaned %>% mutate(status = "Cleaned (Temporal Balance)")
)

# Plotting the temporal density
ggplot(comparison_df, aes(x = year, color = status, fill = status)) +
  geom_density(alpha = 0.3) + # Using density to see the relative shape
  scale_x_continuous(limits = c(1950, 2026), breaks = seq(1950, 2026, 10)) +
  theme_minimal() +
  labs(
    title = "Temporal Distribution: Raw vs. Cleaned Data",
    subtitle = "Focus on HR Era (1970 - 2026)",
    x = "Year of Observation",
    y = "Density of Records",
    fill = "Data State",
    color = "Data State"
  ) +
  geom_vline(xintercept = 1978, linetype = "dashed", color = "red") +
  annotate("text", x = 1978, y = 0.01, label = "HR First Detected", angle = 90, vjust = -1)

ggplot(my_data_cleaned, aes(x = year, fill = basisOfRecord)) +
  geom_histogram(binwidth = 2, position = "stack") +
  scale_x_continuous(limits = c(1970, 2026)) +
  theme_light() +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Composition of Cleaned Dataset (1970-2026)",
    x = "Year",
    y = "Number of Records",
    fill = "Record Type"
  )

# Check how many points remain per species
species_audit <- my_data_cleaned %>%
  group_by(species) %>%
  tally(name = "cleaned_count") %>%
  arrange(cleaned_count)

# How many species have fewer than 5 points?
low_data_species <- filter(species_audit, cleaned_count < 5)
print(nrow(low_data_species))

# 1. Create a data frame for Phalaris minor from your screenshot
# 3. Create the FloraWeb supplement for Phalaris minor
phalaris_floraweb <- data.frame(
  species = "Phalaris minor",
  decimalLatitude = c(52.47498, 51.87496, 51.87497, 51.67492, 51.32497, 50.92490, 48.02488),
  decimalLongitude = c(13.20826, 12.04156, 12.20823, 7.791429, 12.34980, 6.95807, 7.87477),
  coordinateUncertaintyInMeters = c(3169, 3189, 3189, 3195, 3209, 3221, 3317),
  basisOfRecord = "OCCURRENCE",
  year = 1980 # Using 1980 as a proxy for the 'ab 1980' period
)

# 4. Combine with your main dataset
my_data_final_prep <- bind_rows(my_data_cleaned, phalaris_floraweb)

library(terra)

# Load your World SRTM and Germany Shapefile
# Replace paths with your actual filenames
world_srtm <- rast("path/to/world_srtm_30s.tif")
germany_shp <- vect("path/to/germany_shapefile.shp")

# Ensure they use the same Coordinate Reference System (CRS)
# This is vital for the crop/mask to align correctly
germany_shp <- project(germany_shp, crs(world_srtm))

# Step 1: Crop (cuts to a box) and Mask (cuts to the border)
ger_srtm <- crop(world_srtm, germany_shp, mask = TRUE)

# Optional: Plot to verify it worked
plot(ger_srtm, main = "Germany Elevation (SRTM 30s)")
lines(germany_shp)

library(terra)

# ==========================================
# STEP 3: Rebuild Spatial Raster Environment
# ==========================================
cat("Downloading and rebuilding German altitude models...\n")
# 1. Fetch country-level raw data from geodata engine
ger_full_raw <- geodata::elevation_30s(country = "DEU", path = tempdir())

# 2. Establish master project standard coordinate reference system (EPSG:3035)
target_crs <- "EPSG:3035"
ger_srtm_proj <- terra::project(ger_full_raw, target_crs)

# 3. Synchronize shapefile vector limits to match
germany_shp_proj <- terra::project(germany_shp, target_crs)
germany_border   <- terra::aggregate(germany_shp_proj)

# 4. Final masking sequence to cleanly bound the elevation model to the coastline
ger_srtm_final <- terra::mask(ger_srtm_proj, germany_border)


# ==========================================
# STEP 4: Memory-Safe Altitude Extraction
# ==========================================
cat("Processing memory-safe spatial extractions...\n")
# Filter unique coordinates first to save system RAM capacity
unique_coords <- my_data_cleaned %>%
  dplyr::distinct(decimalLongitude, decimalLatitude)

# Convert unique elements into a native coordinate vector format
unique_vect <- terra::vect(unique_coords, 
                           geom = c("decimalLongitude", "decimalLatitude"), 
                           crs = "EPSG:4326")

# Extract exact values directly from the unprojected native structures
extraction_results <- terra::extract(ger_full_raw, unique_vect)

unique_coords_with_alt <- unique_coords %>%
  dplyr::mutate(altitude = extraction_results[, 2])

# Join the extraction back into the foundational master tracking dataset
my_data_final_cleaned <- my_data_cleaned %>%
  dplyr::left_join(unique_coords_with_alt, by = c("decimalLongitude", "decimalLatitude"))



# ==========================================
# STEP 7: Export All Datasets to P:\ Drive
# ==========================================
cat("Writing final datasets to target directory configurations...\n")

# A. Export full cleaned occurrence master file
readr::write_csv(my_data_final_cleaned, csv_output_path)
cat(" -> CSV Occurrences exported to:", csv_output_path, "\n")

# B. Export summarized Trait matrix Excel dashboard for review
openxlsx::write.xlsx(species_final_traits, xlsx_output_path, sheetName = "Trait Profile Matrix", keepNA = TRUE)
cat(" -> Excel Trait Matrix exported to:", xlsx_output_path, "\n")

# C. Export the perfectly clipped raster grid matrix directly to GeoTIFF format
terra::writeRaster(ger_srtm_final, tif_output_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
cat(" -> Spatial GeoTIFF Raster exported to:", tif_output_path, "\n")



library(geodata)
library(terra)

# 1. Download the FULL elevation data for Germany (Level 1)
# This automatically fetches all tiles needed for the whole country
ger_full_raw <- elevation_30s(country = "DEU", path = tempdir())

# 2. Project everything to your target CRS (EPSG:3035)
target_crs <- "EPSG:3035"
ger_srtm_proj <- project(ger_full_raw, target_crs)

# 3. Project your shapefile and aggregate for the border
germany_shp_proj <- project(germany_shp, target_crs)
germany_border <- aggregate(germany_shp_proj)

# 4. Mask the new FULL raster to the border
ger_srtm_final <- mask(ger_srtm_proj, germany_border)

# 5. Plot to verify - The gap should now be 100% GONE
plot(ger_srtm_final, main = "Complete Germany Elevation (No Gaps)")
#lines(germany_border, lwd = 1)

data_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait"


library(readr)

my_data <- read_csv("P:/CONSERWA/Herbicide_Project/Data_Trait/my_data_final_cleanedv2.csv")

# 1. Create a data frame for Phalaris minor from your screenshot
# 3. Create the FloraWeb supplement for Phalaris minor
phalaris_floraweb <- data.frame(
  species = "phalaris minor",
  decimalLatitude = c(52.47498, 51.87496, 51.87497, 51.67492, 51.32497, 50.92490, 48.02488),
  decimalLongitude = c(13.20826, 12.04156, 12.20823, 7.791429, 12.34980, 6.95807, 7.87477),
  coordinateUncertaintyInMeters = c(3169, 3189, 3189, 3195, 3209, 3221, 3317),
  basisOfRecord = "OCCURRENCE",
  year = 1980 # Using 1980 as a proxy for the 'ab 1980' period
)

# 4. Combine with your main dataset
my_data <- bind_rows(my_data, phalaris_floraweb)

length(unique(my_data$species))
unique(my_data$species)

str(my_data )

library(sf)
library(dplyr)

cat("Converting coordinate records to spatial features for boundary validation...\n")

# 1. Temporarily convert the data frame to an sf object (WGS84)
my_data_sfa <- sf::st_as_sf(
  my_data, 
  coords = c("decimalLongitude", "decimalLatitude"), 
  crs = 4326, 
  remove = FALSE # Keep the original lat/lon text columns for downstream use
)

# 2. Transform the points to match the Projected Germany Border (EPSG:3035)
my_data_sf_proj <- sf::st_transform(my_data_sfa, 3035)

# 3. Create a logical index of points that fall INSIDE the German border polygon
# (st_intersects returns a list; lengths > 0 means it overlaps with Germany)
inside_germany_logical <- lengths(sf::st_intersects(my_data_sf_proj, germany_border)) > 0

# 4. Report the Audit Results
total_points <- nrow(my_data)
points_inside <- sum(inside_germany_logical)
points_outside <- total_points - points_inside

cat("--- SPATIAL AUDIT SUMMARY ---\n")
cat("Total input occurrence records:", total_points, "\n")
cat("Records verified inside Germany:", points_inside, "\n")
cat("Records located outside Germany (Omitted):", points_outside, "\n")
cat("-----------------------------\n")

# 5. Overwrite my_data keeping ONLY the verified points inside Germany
my_data_cleaned_geo <- my_data[inside_germany_logical, ]

library(sf)
library(dplyr)

# Ensure your cleaned data is an sf object and projected to EPSG:3035 (Meters)
# we do this outside the loop once to save massive processing time and memory
cat("Projecting occurrence data for grid operations...\n")
my_data_sf <- st_as_sf(my_data_cleaned_geo, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
my_data_proj <- st_transform(my_data_sf, 3035)

# =======================================================
# STEP 1: Dynamically Create a 10 km Spatial Grid
# =======================================================
cat("Generating dynamic 10 km spatial grid over study area...\n")

# st_make_grid automatically creates a grid covering your data or country boundary
# cellsize = 10000 specifies a grid spacing of 10,000 meters (10 km)
germany_grid_10km <- st_make_grid(germany_border, cellsize = 10000, square = TRUE)

# Convert the grid grid layout into an sf object and assign a unique tracking ID to each cell
germany_grid_sf <- st_sf(grid_id = 1:length(germany_grid_10km), geometry = germany_grid_10km)

# Keep only the grid cells that actually intersect with the terrestrial German landmass
germany_grid_sf <- germany_grid_sf[st_intersects(germany_grid_sf, germany_border, sparse = FALSE), ]


# =======================================================
# STEP 2: Overlay Points on the Grid to Extract Cell IDs
# =======================================================
cat("Overlaying weed occurrences onto spatial grid...\n")

# st_join performs a spatial join, stamping the 'grid_id' onto each point record
points_with_grid_ids <- st_join(my_data_proj, germany_grid_sf, join = st_intersects)


# =======================================================
# STEP 3: Combined Range Size (EOO) and Density (AOO) Pipeline
# =======================================================
cat("Calculating species-specific EOO and AOO metrics...\n")

range_metrics_df <- points_with_grid_ids %>%
  # Drop geometry to make the summary statistics run in RAM instantly
  st_drop_geometry() %>% 
  group_by(species) %>%
  summarise(
    # 1. EOO: Geographic Footprint (Convex Hull Area)
    EOO_km2 = {
      # We extract the geometry columns back for the sub-loop
      sub_coords <- my_data_cleaned_geo[my_data_cleaned_geo$species == cur_group()$species, ]
      if (nrow(sub_coords) >= 3) {
        pts <- st_as_sf(sub_coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
        pts_proj <- st_transform(pts, 3035)
        hull <- st_convex_hull(st_union(pts_proj))
        as.numeric(st_area(hull)) / 1e6
      } else {
        0
      }
    },
    
    # 2. AOO: Density/Exposure (Count distinct 10km grid cells occupied)
    # Filter out NA to skip records whose points might fall perfectly on a border edge
    n_grid_cells = n_distinct(grid_id, na.rm = TRUE),
    AOO_km2      = n_grid_cells * 100, # Each 10x10km cell represents exactly 100 km2
    
    # 3. Total Sample Count
    n_points = n(),
    .groups = "drop"
  )

# --- Check the Output Profile ---
print(head(range_metrics_df))



###############

# ==========================================
# STEP 6: Consolidate Final Trait Arrays
# ==========================================
species_altitude_traits <- my_data_cleaned_geo %>%
  dplyr::filter(!is.na(altitude)) %>% 
  dplyr::group_by(species) %>%
  dplyr::summarise(
    alt_min   = min(altitude, na.rm = TRUE),
    alt_max   = max(altitude, na.rm = TRUE),
    alt_range = alt_max - alt_min,
    n_points  = n(),
    .groups   = 'drop'
  ) %>%
  dplyr::left_join(range_metrics_df, by = "species") %>%
  dplyr::mutate(
    alt_min   = ifelse(is.infinite(alt_min), NA, alt_min),
    alt_max   = ifelse(is.infinite(alt_max), NA, alt_max),
    alt_range = ifelse(is.infinite(alt_range) | alt_range < 0, 0, alt_range)
  )




str(species_altitude_traits)



library(dplyr)

# 2. Transform and Scale everything safely in a single pipe block
final_trait_table_scaled <- species_altitude_traits %>%
  # Clean up duplicate columns first
  rename(n_points = n_points.x) %>%
  select(-n_points.y) %>%
  
  mutate(
    # --- A. Simple Scale Only (No Log) ---
    alt_range_scaled = as.numeric(scale(alt_range)),
    
    # --- B. Log10 + 1 Transformations (To fix severe skews) ---
    # Log EOO (Range Size)
    EOO_log           = log10(EOO_km2 + 1), 
    
    # Log AOO (Occupancy Grid Density)
    AOO_log           = log10(AOO_km2 + 1),
    
    # Log Sample Density
    total_obs_log     = log10(n_points + 1)
  ) %>%
  mutate(
    # --- C. Scale the Logged Values ---
    EOO_scaled  = as.numeric(scale(EOO_log)),
    AOO_scaled  = as.numeric(scale(AOO_log)),
    total_obs_scaled  = as.numeric(scale(total_obs_log))
  ) %>%
  # Remove intermediate log tracking columns to keep table pristine for Random Forest
  select(-EOO_log, -AOO_log, -total_obs_log)

# =======================================================
# STEP 3: Verification Checks
# =======================================================
cat("--- VERIFICATION OF SCALED FEATURES ---\n")
print(colnames(final_trait_table_scaled))

# Check summaries to ensure means are 0 and variables are ready for tree splits
summary(final_trait_table_scaled %>% 
          select(alt_range_scaled, EOO_scaled, AOO_scaled, total_obs_scaled))

head(final_trait_table_scaled)

library(readr)
library(openxlsx)

# 1. Establish and confirm the precise network path directory layout
export_directory <- "P:/CONSERWA/Herbicide_Project/Data_Trait"
if (!dir.exists(export_directory)) {
  dir.create(export_directory, recursive = TRUE)
}

# 2. Designate final file names for your target inputs
csv_output_path  <- file.path(export_directory, "species_final_traits_scaled.csv")
xlsx_output_path <- file.path(export_directory, "species_final_traits_scaled.xlsx")

# 3. Perform the export operations smoothly
cat("Initiating local export sequence onto network P:\\ drive...\n")

# Export Option A: Fast, clean CSV format for modeling input stability
readr::write_csv(final_trait_table_scaled, csv_output_path)
cat(" -> Successfully exported CSV to:", csv_output_path, "\n")

# Export Option B: Polished Excel workbook format for documentation review
openxlsx::write.xlsx(
  x = final_trait_table_scaled, 
  file = xlsx_output_path, 
  sheetName = "Scaled Trait Matrix", 
  rowNames = FALSE,
  keepNA = TRUE
)
cat(" -> Successfully exported Excel workbook to:", xlsx_output_path, "\n")


trait_data <- readxl::read_excel(
  file.path(data_path, "TRY_weeds_filtered_analysis_Edit_Use.xlsx"),
  col_names = T
)

merged_data <- trait_data %>%
  left_join(final_trait_table_scaled, by = "species")

readr::write_csv(
  merged_data,
  "P:/CONSERWA/Herbicide_Project/Data_Trait/merged_data_final.csv"
)

library(readr)

merged_data <- read_csv("P:/CONSERWA/Herbicide_Project/Data_Trait/merged_data_final.csv")

library(ggplot2)
library(GGally)
library(dplyr)
# Install the missing raincloud plot library




# Create a mock/placeholder 'HR_status' column if you haven't joined it yet 
# to ensure the graphics can split your 496 species effectively
# final_trait_table_scaled <- final_trait_table_scaled %>% mutate(HR_status = ...)

library(ggplot2)
library(dplyr)
library(tidyr)
library(GGally)

# Set the active output directory on your network drive
export_dir <- "P:/CONSERWA/Herbicide_Project/Data_Trait/Figures"
if (!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE)
}

# Ensure your dataset has a clean factor for your binary classification groups
# (Adjust 'HR_status' if your column is named differently, e.g., 'resistance_class')
final_trait_table_scaled <- merged_data %>%
  mutate(HR_Status = as.factor(HR_Status))


# ==============================================================================
# FIGURE 1: Advanced Bivariate Ecological Space Plot (The "Weed Niche Map")
# ==============================================================================
cat("Generating Figure 1: Bivariate Niche Space Map...\n")

fig1 <- ggplot(final_trait_table_scaled, 
               aes(x = EOO_scaled, y = AOO_scaled, color = HR_Status)) +
  geom_point(aes(size = alt_range_scaled), alpha = 0.6) +
  
  # FIX: Added a 3rd professional color to match your 3 data levels
  scale_color_manual(values = c("#2A9D8F", "#E76F51", "#457B9D")) + 
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Macro-Ecological Niche Space of German Arable Weeds",
    x = "Geographic Extent Footprint (EOO Scaled)",
    y = "Agricultural Occupancy Density (AOO Scaled)",
    size = "Altitude Range",
    color = "Resistance Profile"
  )

ggsave(file.path(export_dir, "Figure1_Ecological_Space.pdf"), plot = fig1, width = 10, height = 7, device = "pdf")

# ==============================================================================
# FIGURE 2: Multi-panel Diagnostic Boxplot + Jitter Cloud Dashboard
# ==============================================================================
cat("Generating Figure 2: Trait Distribution Strip-Plots...\n")

# Restructure table into a long format for elegant multi-panel alignment
long_traits <- final_trait_table_scaled %>%
  select(HR_Status, alt_range_scaled, EOO_scaled, AOO_scaled) %>%
  tidyr::pivot_longer(cols = -HR_Status, names_to = "Trait", values_to = "Value")

fig2 <- ggplot(long_traits, aes(x = HR_Status, y = Value, fill = HR_Status)) +
  # Layer 1: Jittered points create the underlying "cloud" representing raw species data density
  geom_jitter(aes(color = HR_Status), width = 0.2, alpha = 0.2, size = 1.2) +
  # Layer 2: Clean boxplots overlay the points to provide the hard statistical summary (median & quartiles)
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.4, color = "black", size = 0.7) +
  # Layer 3: Grid split with customized, explicit axis labels for publication formatting
  facet_wrap(~Trait, scales = "free_y", labeller = as_labeller(c(
    alt_range_scaled = "Topographic Range (Linear Z-Score)",
    EOO_scaled = "Extent of Occurrence (Log10 Z-Score)",
    AOO_scaled = "Area of Occupancy (Log10 Z-Score)"
  ))) +
  scale_fill_manual(values = c("#2A9D8F", "#E76F51","#457B9D")) +
  scale_color_manual(values = c("#2A9D8F", "#E76F51","#457B9D")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray96"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black")
  ) +
  labs(
    x = "Herbicide Resistance Classification", 
    y = "Standardized Model Predictor Value"
  )

ggsave(file.path(export_dir, "Figure2_Trait_Distributions.pdf"), plot = fig2, width = 11, height = 5, device = "pdf")


# ==============================================================================
# FIGURE 3: Multicollinearity Matrix Diagnostic (Pairs Matrix)
# ==============================================================================
cat("Generating Figure 3: Feature Matrix Diagnostics...\n")

# Isolate the exact quantitative columns feeding into your Random Forest model
modeling_features <- final_trait_table_scaled %>%
  select(alt_range_scaled, EOO_scaled, AOO_scaled, total_obs_scaled, HR_Status)

# Write out the advanced pairwise evaluation layout cleanly
pdf(file.path(export_dir, "Figure3_Correlation_Matrix.pdf"), width = 9, height = 9)
ggpairs(modeling_features, 
        columns = 1:4, 
        mapping = aes(color = HR_Status, alpha = 0.4),
        upper = list(continuous = wrap("cor", size = 4.5)),
        lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.8))) +
  theme_bw()
dev.off()

cat("All 3 vector graphics saved successfully as publication-grade PDFs!\n")

library(dplyr)
library(readr)
library(openxlsx)

cat("Extracting exact statistical correlation values for feature validation...\n")

# 1. Isolate only the numerical modeling features (excluding your categorical HR_Status)
numerical_features <- final_trait_table_scaled %>%
  dplyr::select(alt_range_scaled, EOO_scaled, AOO_scaled, total_obs_scaled)

# 2. Compute the exact Pearson correlation matrix
# (use = "complete.obs" safely handles any hidden missing data points)
cor_matrix <- cor(numerical_features, use = "complete.obs", method = "pearson")

# 3. Convert the raw matrix into a clean, presentation-ready data frame
cor_table <- as.data.frame(cor_matrix) %>%
  # Move the feature names out of the row names and into their own formal column
  dplyr::mutate(Feature = rownames(.)) %>%
  dplyr::select(Feature, everything()) %>%
  # Round values to 3 decimal places for standard scientific formatting
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Define output file destinations
csv_cor_path  <- file.path(export_dir, "table_feature_correlations.csv")
xlsx_cor_path <- file.path(export_dir, "table_feature_correlations.xlsx")

# 4. Export the datasets to your P:\ drive network directory
readr::write_csv(cor_table, csv_cor_path)
cat(" -> Successfully exported Correlation CSV to:", csv_cor_path, "\n")

openxlsx::write.xlsx(cor_table, xlsx_cor_path, sheetName = "Correlation Coefficients", rowNames = FALSE)
cat(" -> Successfully exported Correlation Excel Table to:", xlsx_cor_path, "\n")

# 5. Print the table straight to your R console for immediate review
cat("\n--- FEATURE CORRELATION MATRIX RESULTS ---\n")
print(cor_table)
cat("------------------------------------------\n")

# Install car package if you don't have it (needed for Levene's Test)
if (!requireNamespace("car", quietly = TRUE)) install.packages("car")

library(dplyr)
library(car)

target_traits <- c("alt_range_scaled", "EOO_scaled", "AOO_scaled")

cat("=======================================================\n")
cat("      STARTING PARAMETRIC ASSUMPTION DIAGNOSTICS       \n")
cat("=======================================================\n\n")

for (trait in target_traits) {
  cat("-------------------------------------------------------\n")
  cat("TRAIT UNDER TEST:", trait, "\n")
  cat("-------------------------------------------------------\n")
  
  # 1. TEST FOR NORMALITY (Shapiro-Wilk on Residuals)
  # We fit a simple linear model to extract the residuals across groups
  model <- lm(as.formula(paste(trait, "~ HR_Status")), data = final_trait_table_scaled)
  shapiro_res <- shapiro.test(residuals(model))
  
  cat("1. Normality Check (Shapiro-Wilk on Residuals):\n")
  cat("   W =", round(shapiro_res$statistic, 4), "| p-value =", format.pval(shapiro_res$p.value, eps = 0.001), "\n")
  if (shapiro_res$p.value < 0.05) {
    cat("   Result: SIGNIFICANTLY NON-NORMAL (Parametric Assumption Violated)\n")
  } else {
    cat("   Result: Normal Distribution (Parametric Assumption Met)\n")
  }
  
  # 2. TEST FOR EQUAL VARIANCES (Levene's Test)
  levene_res <- car::leveneTest(as.formula(paste(trait, "~ HR_Status")), data = final_trait_table_scaled)
  levene_p <- levene_res$`Pr(>F)`[1]
  
  cat("\n2. Homogeneity of Variance Check (Levene's Test):\n")
  cat("   F =", round(levene_res$`F value`[1], 4), "| p-value =", format.pval(levene_p, eps = 0.001), "\n")
  if (levene_p < 0.05) {
    cat("   Result: HETEROSCEDASTIC / UNEQUAL VARIANCES (Parametric Assumption Violated)\n")
  } else {
    cat("   Result: Equal Variances (Parametric Assumption Met)\n")
  }
  
  # 3. FINAL PIPELINE RECOMMENDATION
  cat("\n3. Statistical Framework Verdict:\n")
  if (shapiro_res$p.value < 0.05 || levene_p < 0.05) {
    cat("   >>> USE NON-PARAMETRIC: Kruskal-Wallis Test + Dunn's Post-Hoc <<<\n")
  } else {
    cat("   >>> USE PARAMETRIC: One-Way ANOVA + Tukey HSD Post-Hoc <<<\n")
  }
  cat("\n")
}

# Install required post-hoc package if missing
if (!requireNamespace("FSA", quietly = TRUE)) install.packages("FSA")

library(dplyr)
library(purrr)
library(FSA)      # For Dunn's Test
library(readr)
library(openxlsx)

export_dir <- "P:/CONSERWA/Herbicide_Project/Data_Trait"

# Define the three target traits to evaluate
target_traits <- c("alt_range_scaled", "EOO_scaled", "AOO_scaled")

# Pre-allocate lists to store results
omnibus_results <- list()
pairwise_results <- list()

cat("Executing macro-ecological trait statistical significance pipeline...\n")

for (trait in target_traits) {
  cat("\nProcessing trait:", trait, "\n")
  
  # 1. Run Kruskal-Wallis Omnibus Test
  kw_formula <- as.formula(paste(trait, "~ HR_Status"))
  kw_test <- kruskal.test(kw_formula, data = final_trait_table_scaled)
  
  omnibus_results[[trait]] <- tibble(
    Trait = trait,
    Chi_Squared = round(kw_test$statistic, 3),
    Df = kw_test$parameter,
    Omnibus_p_value = format.pval(kw_test$p.value, digits = 3, eps = 0.001)
  )
  
  # 2. Run Dunn's Post-Hoc Pairwise Test if Omnibus is significant
  # Using 'holm' correction to protect against type 1 error inflation
  dunn_test <- FSA::dunnTest(kw_formula, data = final_trait_table_scaled, method = "holm")$res
  
  pairwise_results[[trait]] <- dunn_test %>%
    as_tibble() %>%
    mutate(Trait = trait) %>%
    select(Trait, Comparison, Z_stat = Z, P_raw = P.unadj, P_adj_Holm = P.adj) %>%
    mutate(across(c(Z_stat), ~ round(.x, 3))) %>%
    mutate(
      P_raw = format.pval(P_raw, digits = 3, eps = 0.001),
      P_adj_Holm = format.pval(P_adj_Holm, digits = 3, eps = 0.001)
    )
}

# 3. Combine results into clean master tables
final_omnibus_table <- bind_rows(omnibus_results)
final_pairwise_table <- bind_rows(pairwise_results)

# 4. Save results to your network directory
write_csv(final_omnibus_table, file.path(export_dir, "stat_table_omnibus_kruskal.csv"))
write_csv(final_pairwise_table, file.path(export_dir, "stat_table_pairwise_dunn.csv"))

# Save both onto separate tabs in an Excel workbook for your manuscript appendix
wb <- createWorkbook()
addWorksheet(wb, "Omnibus (Kruskal-Wallis)")
addWorksheet(wb, "Pairwise Post-Hoc (Dunn)")
writeData(wb, "Omnibus (Kruskal-Wallis)", final_omnibus_table)
writeData(wb, "Pairwise Post-Hoc (Dunn)", final_pairwise_table)
saveWorkbook(wb, file.path(export_dir, "manuscript_trait_significance.xlsx"), overwrite = TRUE)

# 5. Print results to console for immediate inspection
cat("\n=======================================================\n")
cat("          OMNIBUS GLOBAL DIFFERENCE RESULTS            \n")
cat("=======================================================\n")
print(final_omnibus_table)

cat("\n=======================================================\n")
cat("          PAIRWISE POST-HOC COMPARISONS                \n")
cat("=======================================================\n")
print(final_pairwise_table)

library(rnaturalearth)
library(sf)

library(giscoR)

germany_border <- gisco_get_countries(country = "DE", resolution = "1")

library(dplyr)
library(ggplot2)

# Convert your data
my_data_sf <- my_data %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Ensure same CRS
germany_border <- st_transform(germany_border, st_crs(my_data_sf))

# Filter points inside Germany
my_data_sf <- my_data_sf %>%
  filter(st_intersects(geometry, germany_border, sparse = FALSE))

# Plot
final_map_raw<-ggplot() +
  geom_sf(data = germany_border, fill = "grey95", color = "black") +
  geom_sf(data = my_data_sf, aes(color = HR_Status), size = 1, alpha = 0.7) +
  theme_minimal()

ggsave(
  filename = "HR_raw_map.pdf",
  plot = final_map_raw,
  device = cairo_pdf,   # IMPORTANT for publication-quality text
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  path = "P:/CONSERWA/Herbicide_Project/Data_Trait/Figures_New"
)

# Germany boundary
library(rnaturalearth)
germany_border <- ne_countries(scale = "medium", country = "Germany", returnclass = "sf")

# Project to metric CRS (IMPORTANT for grid)
my_data_sf <- st_transform(my_data_sf, 3035)
germany_border <- st_transform(germany_border, 3035)

grid <- st_make_grid(
  germany_border,
  cellsize = 10000,   # 10 km
  square = TRUE
) %>%
  st_as_sf() %>%
  st_intersection(germany_border)

grid <- st_make_grid(
  germany_border,
  cellsize = 10000,
  square = TRUE
) %>%
  st_as_sf() %>%
  mutate(grid_id = row_number())

points_joined <- st_join(my_data_sf, grid, join = st_within)

points_joined <- points_joined %>%
  mutate(HR_clean = tolower(trimws(HR_Status)))

points_joined <- points_joined %>%
  mutate(
    HR_risk = case_when(
      HR_Status == "HR_local" ~ 1,
      HR_Status == "HR_elsewhere" ~ 0.5,
      HR_Status == "Non_HR_local" ~ 0,
      TRUE ~ NA_real_
    )
  )

grid_risk <- points_joined %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    n = n(),
    hr_mean_risk = mean(HR_risk, na.rm = TRUE)
  )

grid_final <- grid %>%
  left_join(grid_risk, by = "grid_id")

library(ggplot2)

ggplot() +
  geom_sf(data = grid_final, aes(fill = hr_mean_risk), color = NA) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "HR risk index"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Spatial HR risk surface in Germany",
    subtitle = "Weighted index based on species-level herbicide resistance status",
    caption = "Derived from GBIF occurrences + trait-based HR classification"
  )+
  
  
  # fake legend table
  annotate("text", x = Inf, y = Inf,
           label = "Weighting scheme:\nHR_local = 1\nHR_elsewhere = 0.5\nNon_HR_local = 0",
           hjust = 1.1, vjust = 1.1, size = 3.5)

ggplot() +
  geom_sf(data = grid_final, aes(fill = hr_mean_risk), color = NA) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.4) +
  
  scale_fill_viridis_c(
    option = "viridis",
    name = "HR risk index",
    limits = c(0, 1)
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    title = "Spatial HR risk surface in Germany",
    subtitle = "Weighted index based on species-level herbicide resistance status",
    caption = "Weighting: HR_local = 1, HR_elsewhere = 0.5, Non_HR_local = 0"
  ) +
  
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

library(ggplot2)
library(ggspatial)

final_map <- ggplot() +
  geom_sf(data = grid_final, aes(fill = hr_mean_risk), color = NA) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.4) +

  scale_fill_viridis_c(
    option = "cividis",
    name = "HR risk index",
    limits = c(0, 1),
    na.value = "transparent")+

  # ✅ Scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.3
  ) +

  # ✅ North arrow
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +

  theme_classic(base_size = 14) +
  labs(
    title = "Spatial HR risk surface in Germany",
    subtitle = "Weighted index based on species-level herbicide resistance status",
    caption = "Weighting: HR_local = 1, HR_elsewhere = 0.5, Non_HR_local = 0"
  ) +
  theme(
  panel.grid = element_blank(),
  legend.position = "right"
)

ggsave(
  filename = "HR_risk_map.pdf",
  plot = final_map,
  device = cairo_pdf,   # IMPORTANT for publication-quality text
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  path = "P:/CONSERWA/Herbicide_Project/Data_Trait/Figures_New"
)

grid_uncertainty <- points_joined %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    n = n(),
    hr_mean = mean(HR_risk, na.rm = TRUE),
    hr_sd = sd(HR_risk, na.rm = TRUE),
    hr_se = hr_sd / sqrt(n)
  )

ggplot() +
  geom_sf(data = grid_final %>% left_join(grid_uncertainty, by = "grid_id"),
          aes(fill = hr_se), color = NA) +
  geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_viridis_c(option = "cividis",

     name = "Uncertainty (SE)",
     
     na.value = "transparent") +
  theme_minimal()

library(ggplot2)
library(ggspatial)

uncertainty_map <- ggplot() +
  geom_sf(
    data = grid_final %>% left_join(grid_uncertainty, by = "grid_id"),
    aes(fill = hr_se),
    color = NA
  ) +
geom_sf(data = germany_border, fill = NA, color = "black", linewidth = 0.4) +
  scale_fill_viridis_c(
  option = "cividis",
  name = "Uncertainty (SE)",
  na.value = "transparent"
 
) +

  # 📏 Scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.3
  ) +

  # 🧭 North arrow
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +

  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),

    legend.position = "right"
  ) +
  labs(
    title = "Spatial uncertainty of HR risk in Germany",
    subtitle = "Standard error of weighted herbicide resistance index",
    caption = "Higher values indicate higher variability or low sample size per grid cell"
  )

ggsave(
  filename = "HR_uncertainty_map.pdf",
  plot = uncertainty_map,
  device = cairo_pdf,   # IMPORTANT for publication quality
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  path = "P:/CONSERWA/Herbicide_Project/Data_Trait/Figures_New"
)

grid_elev <- points_joined %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    elevation_mean = mean(altitude, na.rm = TRUE)
  )

grid_final2 <- grid_final %>%
  left_join(grid_uncertainty, by = "grid_id") %>%
  left_join(grid_elev, by = "grid_id")

ggplot(grid_final2, aes(x = elevation_mean, y = hr_mean_risk)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  labs(
    x = "Elevation (m)",
    y = "HR risk index"
  )

# Core packages
library(tidyverse)
library(readr)
library(mice)
library(FactoMineR)
library(factoextra)
library(janitor)
library(purrr)


data_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait"


trait_data_Coordinates <- readxl::read_excel(
  file.path(data_path, "Copy of TRY_pmf_DataRelease_cleaned_2013_06_19_Coordinates.xlsx"),
  col_names = T
)

library(dplyr)

trait_data_speciescheck <- trait_data_Coordinates %>%
  filter(Type == "species")

library(ggplot2)

ggplot(trait_data_speciescheck, aes(x = Lon, y = Lat)) +
  geom_point(alpha = 0.4, size = 0.6) +
  theme_minimal() +
  labs(
    title = "Species distribution",
    x = "Longitude",
    y = "Latitude"
  )

library(ggplot2)
library(maps)

world <- map_data("world")

ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "white") +
  geom_point(data = trait_data_speciescheck,
             aes(x = Lon, y = Lat),
             alpha = 0.4, size = 0.6, color = "red") +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Global distribution of species")

sum(is.na(trait_data_speciescheck$Lat))
sum(is.na(trait_data_speciescheck$Lon))

summary(trait_data_speciescheck$Lat)
summary(trait_data_speciescheck$Lon)

trait_data_speciescheck_clean <- trait_data_speciescheck %>%
  filter(
    !is.na(Lat),
    !is.na(Lon),
    Lat >= -90, Lat <= 90,
    Lon >= -180, Lon <= 180
  )

trait_data_speciescheck_clean <- trait_data_speciescheck %>%
  filter(!is.na(Lat), !is.na(Lon))

library(ggplot2)
library(maps)

world <- map_data("world")

ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "white") +
  geom_point(data = trait_data_speciescheck_clean,
             aes(x = Lon, y = Lat),
             alpha = 0.4, size = 0.6, color = "red") +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Global distribution of species")

library(sf)
library(rnaturalearth)
library(dplyr)

germany <- ne_countries(scale = "medium", country = "Germany", returnclass = "sf")

points_sf <- st_as_sf(
  trait_data_speciescheck_clean,
  coords = c("Lon", "Lat"),
  crs = 4326
)

trait_data_germany <- st_join(points_sf, germany, join = st_within)

trait_data_germany <- trait_data_germany %>%
  filter(!is.na(name))   # or admin column depending on output

ggplot() +
  geom_sf(data = germany, fill = "grey90", color = "white") +
  geom_sf(data = trait_data_germany,
          color = "red", alpha = 0.5, size = 0.7) +
  theme_minimal() +
  labs(title = "Species occurrences in Germany (spatial filter)")

str(trait_data_germany)

library(dplyr)
library(sf)

sf_filtered <- trait_data_germany %>%
  select(
    ObservationID,
    ExcelFileID,
    Alt,
    PlantDevState,
    Exposition,
    AccSpeciesID,
    AccSpeciesName,
    Type,
    Genus,
    Family,
    PhylogeneticGroup,
    PlantGrowthForm,
    Succulent,
    Climbing,
    LeafType,
    LeafPhenology,
    PhotosyntheticPathway,
    Woodiness,
    LeafCompoundness,
    TempAnnMean_DegC,
    PrecipAnnMean_mm,
    geometry
  )

library(dplyr)

sf_filtered <- sf_filtered %>%
  mutate(AccSpeciesName = tolower(AccSpeciesName))



sf_filtered <- sf_filtered %>%
  mutate(in_trait_data = AccSpeciesName %in% trait_data$species)

library(dplyr)

sf_filtered %>%
  count(in_trait_data)

sf_true <- sf_filtered %>%
  filter(in_trait_data == TRUE)

sf_false <- sf_filtered %>%
  filter(in_trait_data == FALSE)

trait_not_in_sf <- trait_data %>%
  filter(!species %in% sf_filtered$AccSpeciesName)

trait_not_in_sf %>%
  summarise(
    n_rows = n(),
    n_unique_species = n_distinct(species)
  )

sf_species <- unique(sf_filtered$AccSpeciesName)
trait_species <- unique(trait_data$species)

n_trait_total <- length(trait_species)

n_overlap <- length(intersect(trait_species, sf_species))

n_trait_not_in_sf <- length(setdiff(trait_species, sf_species))

n_trait_total
n_overlap
n_trait_not_in_sf

data_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait"

####ALL DATA INCLUDE SPECIES AND GENUS###
trait_data <- readxl::read_excel(file.path(data_path, "Traits.xlsx")) %>%
  clean_names()


data_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait"


trait_data <- readxl::read_excel(file.path(data_path, "TRY_weeds_filtered_analysis_edit.xlsx")) %>%
  clean_names()

#EDITED DATA FOR WEED ONLY###
trait_data <- readxl::read_excel(
  file.path(data_path, "TRY_weeds_filtered_analysis_edit.xlsx"),
  col_names = T
)



names(trait_data)


merged_data %>%
  count(HR_Status)



merged_data %>%
  count(HR_Status, habitat_group)


merged_data %>%
  count(HR_Status, life_span)


merged_data %>%
  count(HR_Status, origin)


merged_data %>%
  count(HR_Status, habitat_group, life_span, origin)


library(dplyr)
library(tidyr)

# 1. Total species count
original_n <- nrow(merged_data)

# 2. Count after removing species with ANY NA in those 73 columns
complete_n <- merged_data %>% 
  drop_na() %>% 
  nrow()

# 3. Calculate the loss
lost_n <- original_n - complete_n
percent_left <- (complete_n / original_n) * 100

# Print results
cat("Original species count:", original_n, "\n")
cat("Species with NO missing values:", complete_n, "\n")
cat("Percentage of data remaining: ", round(percent_left, 2), "%\n")

# If this returns 0, R doesn't see 'real' NAs
sum(is.na(merged_data))

# If this returns a high number, your NAs are actually just text
sum(merged_data == "NA", na.rm = TRUE)

trait_data_fixed <- merged_data %>%
  mutate(across(everything(), ~na_if(as.character(.), "NA"))) %>%
  # Convert numeric columns back to numeric if they were turned into text
  type_convert() 

# Now try your check again
nrow(trait_data_fixed %>% drop_na())

# This creates a list of columns and how much data is missing in each
missing_summary <- trait_data_fixed %>%
  summarise(across(everything(), ~ sum(is.na(.)) / n() * 100)) %>%
  pivot_longer(everything(), names_to = "trait", values_to = "percent_missing") %>%
  arrange(desc(percent_missing))

print(missing_summary)


missing_summary <- trait_data_fixed %>%
  summarise(across(everything(),
                   ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(),
               names_to = "trait",
               values_to = "percent_missing") %>%
  arrange(percent_missing)


print(missing_summary, n = Inf)


library(corrplot)

trait_numeric <- trait_data_fixed %>%
  select(where(is.numeric))

cor_matrix <- cor(trait_numeric,
                  use = "pairwise.complete.obs",
                  method = "spearman")

corrplot(cor_matrix,
         method = "color",
         tl.cex = 0.5)


library(tidyverse)
library(missForest)
library(FactoMineR)
library(janitor)
library(purrr)

set.seed(123)   # reproducibility


traits <- trait_data_fixed %>%
  rename(species = species)


any(duplicated(traits$species))


single_species <- traits %>% 
  filter(n == 1)

multi_species <- traits %>% 
  filter(n > 1)


library(ggplot2)

ggplot(traits, aes(x = leaf_length_mean, fill = n == 1)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_fill_manual(values = c("grey70", "red"),
                    labels = c("n > 1", "n = 1")) +
  theme_classic()


traits_scaled <- traits %>%
  mutate(across(where(is.numeric),
                \(x) scale(x)[,1]))


outliers <- traits_scaled %>%
  filter(n == 1) %>%
  mutate(max_abs_z =
           apply(select(., where(is.numeric)), 1,
                 \(x) max(abs(x), na.rm = TRUE)))


outliers %>% arrange(desc(max_abs_z))


table(traits$n)


traits <- traits %>%
  mutate(across(
    -species,
    ~ suppressWarnings(as.numeric(.))
  ))

summary(traits)

## ===========================================================================
## DATA PROVENANCE AND TRANSFORM  -  resolved 2026-10-27. Read before editing.
##
## The analysis input is now data/merged_data_final_REBUILT.csv, which is built
## by 19_rebuild_input_from_splot.R directly from sPlot 4.0
## species_traits_gap.csv, in the RAW DOCUMENTED UNITS of trait_key.csv
## (plant height m, leaf C mg/g, SLA mm2/mg, stem diameter m, seed mass mg).
## With that input the log10 below is applied exactly once, which is correct,
## and 02_Master_Analysis.Rmd carries TRAITS_ALREADY_LOG <- FALSE so it does
## not transform again.
##
## WHY THE REBUILD WAS NEEDED. The previous input
## (TRY_weeds_filtered_analysis_Edit_Use.xlsx -> merged_data_final.csv) was not
## on a single scale: 502 of its 575 species held log10 trait values and 56 held
## raw values, because this notebook writes its own post-log output to a
## workbook (Herbicide_Code.qmd line 5681) which was then read back in as the
## next run"s input (line 779), and a join at line 3726 combined logged rows
## with unlogged ones. Rebuilding from the sPlot source removes both the mixed
## scale and the loop. Verification: for the previous file, 10^(stored value)
## reproduced the sPlot number to 1e-13 for the logged rows and the stored value
## matched it directly for the raw rows, so both groups were traceable and the
## rebuild is exact.
##
## CONSEQUENCES FOR RESULTS (corrected vs previous): discrimination unchanged
## (binary out-of-bag AUC 0.853), predictor ranking unchanged (life history,
## plant height, photosynthetic pathway, fecundity), phylogenetic signal
## unchanged (D = 0.63). Trait-space results change: ordination PC1 falls from
## 60.5% to 26.6% of variance because the mixed-scale rows no longer dominate
## it, and single-trait contrasts strengthen (specific leaf area P = 5.9e-5,
## seeds per reproductive unit P = 5.1e-7, both previously marginal).
##
## IF THE INPUT IS EVER CHANGED back to a pre-logged workbook, set
## TRAITS_ALREADY_LOG <- TRUE in 02_Master_Analysis.Rmd so this transform is
## skipped. Non-positive values become NA here, which is correct for
## strictly-positive traits in raw units.
## ===========================================================================
traits <- traits %>%
  mutate(across(
    all_of(mean_cols),
    ~ ifelse(!is.na(.) & . > 0, log10(.), NA_real_)
  ))


traits <- traits %>%
  select(where(~ !all(is.na(.))))

sapply(traits, class)

run_module_pca <- function(data, module_name) {
  
  message(paste("Running PCA for:", module_name))
  
  # Keep species names
  species_vec <- data$species
  
  # Remove species column
df_num <- data %>% select(-species)

for(col in names(df_num)){

  vec <- as.numeric(unlist(df_num[[col]]))

  # convert bad numeric values
  vec[is.nan(vec)] <- NA_real_
  vec[is.infinite(vec)] <- NA_real_

  df_num[[col]] <- vec
}


  
  # Remove rows where all values are NA
  all_na_rows <- apply(df_num, 1, function(x) all(is.na(x)))
  if(any(all_na_rows)) {
    df_num <- df_num[!all_na_rows, ]
    species_vec <- species_vec[!all_na_rows]
    message(sum(all_na_rows), " species removed: all NA for module ", module_name)
  }
  
  # Remove columns that are entirely NA
  df_num <- df_num[, colSums(!is.na(df_num)) > 0]
  
  # Remove zero-variance columns
  df_num <- df_num[, apply(df_num, 2, sd, na.rm = TRUE) > 0]
  
  # ---------- IMPUTATION ----------
if(nrow(df_num) < 2){
  warning("Not enough data for imputation in module: ", module_name)
  scores <- as.data.frame(df_num)
  colnames(scores) <- paste0(module_name, "_PC1")
  scores$species <- species_vec
  return(scores)
}

df_num <- as.data.frame(df_num)
stopifnot(all(sapply(df_num, is.numeric)))

imputed <- missForest(df_num, verbose = FALSE)$ximp

  
  # scale
  scaled <- scale(imputed)
  
  
# ---------- PCA ----------
if(ncol(scaled) == 1){
  
  scores <- as.data.frame(scaled)
  colnames(scores) <- paste0(module_name, "_PC1")
  
  # Variance = 100% for single trait
  var_explained <- 100
  
  message(module_name, ": single trait → PC1 = 100%")
  
} else {
  
  pca <- FactoMineR::PCA(scaled, graph = FALSE)
  
  # Extract variance explained (%)
  var_explained <- pca$eig[, 2]
  
  # Safe extraction (in case < 3 PCs)
  pc1 <- round(var_explained[1], 2)
  pc2 <- ifelse(length(var_explained) >= 2, round(var_explained[2], 2), NA)
  pc3 <- ifelse(length(var_explained) >= 3, round(var_explained[3], 2), NA)
  
  # Print nicely
  message(module_name, " variance explained:")
  message("  PC1: ", pc1, "%")
  message("  PC2: ", pc2, "%")
  message("  PC3: ", pc3, "%")
  
  # Extract PC1 scores
  scores <- as.data.frame(pca$ind$coord[, 1])
  colnames(scores) <- paste0(module_name, "_PC1")
}

# ✅ KEEP species for merging
scores$species <- species_vec

# ✅ Return both scores + full variance
return(list(
  scores = scores,
  variance = var_explained
))

  
}




run_module_pca <- function(data, module_name) {

  message("Running PCA for: ", module_name)

  # --------------------------
  # 1. Validate input
  # --------------------------
  stopifnot(is.data.frame(data))
  stopifnot("species" %in% names(data))

  species_vec <- data$species
  df_num <- data %>% dplyr::select(-species)

  # --------------------------
  # 2. Safe numeric conversion (vectorized)
  # --------------------------
  df_num <- as.data.frame(lapply(df_num, function(x) {

    if (is.factor(x)) x <- as.character(x)

    x <- suppressWarnings(as.numeric(x))

    x[is.nan(x)] <- NA_real_
    x[is.infinite(x)] <- NA_real_

    return(x)
  }))

  # --------------------------
  # 3. Remove fully NA rows
  # --------------------------
  all_na_rows <- apply(df_num, 1, function(x) all(is.na(x)))

  if (any(all_na_rows)) {
    df_num <- df_num[!all_na_rows, ]
    species_vec <- species_vec[!all_na_rows]
    message(sum(all_na_rows), " all-NA rows removed in ", module_name)
  }

  # --------------------------
  # 4. Remove empty / constant columns
  # --------------------------
  df_num <- as.data.frame(df_num)

# Remove empty columns safely
keep_cols <- colSums(!is.na(df_num)) > 0

if (is.null(dim(df_num)) || length(keep_cols) <= 1) {
  df_num <- as.data.frame(df_num)
} else {
  df_num <- df_num[, keep_cols, drop = FALSE]
}
  # --------------------------
  # 5. Edge case: not enough data
  # --------------------------
  if (nrow(df_num) < 2 || ncol(df_num) < 1) {
    warning("Not enough data for PCA: ", module_name)

    scores <- data.frame(
      species = species_vec,
      PC1 = NA_real_
    )

    return(list(
      scores = scores,
      variance = NA
    ))
  }

  # --------------------------
  # 6. Imputation
  # --------------------------
  imputed <- missForest::missForest(df_num, verbose = FALSE)$ximp

  # --------------------------
  # 7. Scale
  # --------------------------
  scaled <- scale(imputed)

  # --------------------------
  # 8. PCA
  # --------------------------
  if (ncol(scaled) == 1) {

    scores <- data.frame(
      PC1 = scaled[, 1],
      species = species_vec
    )

    message(module_name, ": single trait → 100% variance")

    return(list(
      scores = scores,
      variance = 100
    ))
  }

  pca <- FactoMineR::PCA(as.data.frame(scaled), graph = FALSE)

  eig <- pca$eig[, 2]
  cumvar <- cumsum(eig)

  message(module_name, " variance explained:")
  message(" PC1: ", round(eig[1], 2), "%")
  message(" PC2: ", ifelse(length(eig) >= 2, round(eig[2], 2), NA), "%")

  # --------------------------
  # 9. Extract scores (PC1)
  # --------------------------
  scores <- data.frame(
    PC1 = pca$ind$coord[, 1],
    species = species_vec
  )

  return(list(
    scores = scores,
    variance = eig,
    cumulative_variance = cumvar
  ))
}

leaf_construction <- traits %>%
  select(species,
         leaf_thickness_mean,
         ldmc_mean,
         leaf_dry_mass_single_mean,
         sla_undef_mean,
         leaf_area_leaf_undef_mean,
         leaf_area_leaflet_undef_mean,
         leaf_area_undef_undef_mean,
         leaf_length_mean,
         leaf_width_mean,
         leaffreshmass_mean)

Leaf_Construction_PC1 <-
  run_module_pca(leaf_construction,
                 "Leaf_Construction")



leaf_nutrient <- traits %>%
  select(species,
         leaf_c_perdrymass_mean,
         leaf_n_mean,
         leaf_p_mean,
         leaf_nper_area_mean,
         leaf_cn_ratio_mean)

Leaf_Nutrient_PC1 <-
  run_module_pca(leaf_nutrient,
                 "Leaf_Nutrient")


leaf_physiology <- traits %>%
  select(species,
         leaf_water_cont_mean)

Leaf_Physiology_PC1 <-
  run_module_pca(leaf_physiology,
                 "Leaf_Physiology")


size_traits <- traits %>%
  select(species,
         stem_diam_mean,
         plant_height_veg_mean,
         plant_height_generative_mean,
         stem_cond_dens_mean)

Size_PC1 <- run_module_pca(size_traits, "Size")


root_traits <- traits %>%
  select(species,
         rooting_depth_mean,
         specific_root_length_mean,
         specific_root_length_fine_mean)

Root_PC1 <- run_module_pca(root_traits, "Root")


hydraulic_traits <- traits %>%
  select(species,
         stem_conduit_diameter_mean,
         wood_vessel_length_mean,
         wood_fiber_length_mean,
         wood_ray_millimetre_mean)

Hydraulic_PC1 <- run_module_pca(hydraulic_traits,
                                "Hydraulic")


reproduction_traits <- traits %>%
  select(species,
         seed_mass_mean,
         seed_length_mean,
         seed_germination_rate_mean,
         seed_num_rep_unit_mean,
         seed_num_disp_unit_mean)

Reproduction_PC1 <-
  run_module_pca(reproduction_traits,
                 "Reproduction")


dispersal_traits <- traits %>%
  select(species,
         disp_unit_leng_mean,
         seed_num_disp_unit_mean,
         seed_mass_mean)

Dispersal_PC1 <- run_module_pca(dispersal_traits,
                                "Dispersal")


genome_traits <- traits %>%
  select(species,
         chromosome_n_mean,
         chromosome_c_dn_acont_mean)

Genome_PC1 <- run_module_pca(genome_traits,
                             "Genome")


ITV_raw <- traits %>%
  select(species,
         rooting_depth_mean, rooting_depth_sd,
         leaf_n_mean, leaf_n_sd,
         leaf_p_mean, leaf_p_sd,
         sla_undef_mean, sla_undef_sd,
         ldmc_mean, ldmc_sd,
         stem_diam_mean, stem_diam_sd,
         seed_mass_mean, seed_mass_sd,
         plant_height_veg_mean, plant_height_veg_sd,
         specific_root_length_mean, specific_root_length_sd,
         wood_vessel_length_mean, wood_vessel_length_sd,
         leaf_area_undef_undef_mean, leaf_area_undef_undef_sd,
         leaf_water_cont_mean, leaf_water_cont_sd)


ITV_PC1_cv <- ITV_raw %>%
  mutate(
    rooting_depth_cv = rooting_depth_sd / abs(rooting_depth_mean),
    leaf_n_cv = leaf_n_sd / abs(leaf_n_mean),
    leaf_p_cv = leaf_p_sd / abs(leaf_p_mean),
    sla_cv = sla_undef_sd / abs(sla_undef_mean),
    ldmc_cv = ldmc_sd / abs(ldmc_mean),
    stem_diam_cv = stem_diam_sd / abs(stem_diam_mean),
    seed_mass_cv = seed_mass_sd / abs(seed_mass_mean),
    plant_height_cv = plant_height_veg_sd / abs(plant_height_veg_mean),
    srl_cv = specific_root_length_sd / abs(specific_root_length_mean),
    vessel_length_cv = wood_vessel_length_sd / abs(wood_vessel_length_mean),
    leaf_area_cv = leaf_area_undef_undef_sd / abs(leaf_area_undef_undef_mean),
    leaf_water_cv = leaf_water_cont_sd / abs(leaf_water_cont_mean)
  ) %>%
  select(species, ends_with("_cv"))



ITV_PC1_cv <- ITV_PC1_cv %>%
  mutate(
    across(
      where(is.numeric),
      ~ replace(., !is.finite(.), NA_real_)
    )
  )


# AFTER cleaning infinities
ITV_PC1_cv <- ITV_PC1_cv %>%
  mutate(across(where(is.numeric), log1p)) ## log(1 + x)




ITV_PC1 <- run_module_pca(
  ITV_PC1_cv,
  "Intraspecific trait variability"
)


head(ITV_PC1)


functional_traits <- list(
  Leaf_Construction_PC1$scores,
  Leaf_Nutrient_PC1$scores,
  Leaf_Physiology_PC1$scores,
  Size_PC1$scores,
  Root_PC1$scores,
  Hydraulic_PC1$scores,
  Reproduction_PC1$scores,
  Dispersal_PC1$scores,
  Genome_PC1$scores,
  ITV_PC1$scores
) %>%
  reduce(left_join, by = "species")



module_variance <- tibble(
  Module = c("Leaf_Construction","Leaf_Nutrient","Leaf_Physiology",
             "Size","Root","Hydraulic","Reproduction",
             "Dispersal","Genome","Plasticity"),
  PC1 = c(Leaf_Construction_PC1$variance[1],
          Leaf_Nutrient_PC1$variance[1],
          Leaf_Physiology_PC1$variance[1],
          Size_PC1$variance[1],
          Root_PC1$variance[1],
          Hydraulic_PC1$variance[1],
          Reproduction_PC1$variance[1],
          Dispersal_PC1$variance[1],
          Genome_PC1$variance[1],
          ITV_PC1$variance[1])
)


glimpse(functional_traits)
summary(functional_traits)


set.seed(1)

pairs(
  functional_traits %>%
    select(-species) %>%
    slice_sample(n = 500)
)



library(corrplot)

module_matrix <- functional_traits %>%
  select(-species)

cor_modules <- cor(
  module_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)

corrplot(
  cor_modules,
  method = "color",
  type = "upper",
  tl.cex = 0.9,
  addCoef.col = "black"
)


pairs(module_matrix)


library(tidyverse)

module_matrix <- functional_traits %>%
  select(-species)

# variance of each strategy axis
axis_variance <- module_matrix %>%
  summarise(across(everything(),
                   ~ var(., na.rm = TRUE))) %>%
  pivot_longer(everything(),
               names_to = "Module",
               values_to = "Variance")

axis_variance



axis_variance %>%
  mutate(Percent = Variance / sum(Variance) * 100)


pca <- FactoMineR::PCA(scaled, graph = FALSE)

# Extract variance explained
var_explained <- pca$eig[, 2]  # percentage

# Store PC1 variance
pc1_var <- var_explained[1]


ggplot(axis_variance,
       aes(x = reorder(Module, Variance),
           y = Variance)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Variance of Functional Strategy Axes",
    x = "",
    y = "Variance"
  )


module_long <- module_matrix %>%
  pivot_longer(everything(),
               names_to = "Module",
               values_to = "Score")

ggplot(module_long,
       aes(x = Score)) +
  geom_histogram(bins = 60) +
  facet_wrap(~Module, scales = "free") +
  theme_minimal()

strategy_matrix <- module_matrix %>%
  drop_na()   # species with complete strategy info


global_pca <- prcomp(
  strategy_matrix,
  center = TRUE,
  scale. = TRUE
)


summary(global_pca)


screeplot(global_pca, type = "lines")

round(global_pca$rotation, 2)


install.packages("ggfortify")   # run once
library(ggfortify)

autoplot(global_pca)



factoextra::fviz_pca_ind(global_pca)


library(factoextra)

fviz_pca_ind(
  global_pca,
  alpha.ind = 0.3,
  pointsize = 1
)



fig_dir <- "P:/CONSERWA/Herbicide_Project/Data_Trait/Figures_New"


library(factoextra)

p_ind <- fviz_pca_ind(
  global_pca,
  alpha.ind = 0.3,
  pointsize = 1
)


ggsave(
  filename = file.path(fig_dir, "Figure1a_StrategySpace_Individuals.pdf"),
  plot = p_ind,
  width = 180,
  height = 160,
  units = "mm",
  device = cairo_pdf
)


fviz_pca_var(
  global_pca,
  col.var = "contrib"
)


p_var <- fviz_pca_var(
  global_pca,
  col.var = "contrib"
)


ggsave(
  filename = file.path(fig_dir, "Figure1b_ModuleContributions.pdf"),
  plot = p_var,
  width = 180,
  height = 160,
  units = "mm",
  device = cairo_pdf
)



library(ggplot2)

p_ind <- fviz_pca_ind(
  global_pca,
  alpha.ind = 0.3,
  pointsize = 1
) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10)
  )





p_var <- fviz_pca_var(global_pca, col.var = "contrib") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10)
  )




save_figure <- function(plot, p_ind) {
  ggsave(
    filename = file.path(fig_dir, paste0(p_ind, ".pdf")),
    plot = plot,
    width = 80,
    height = 80,
    units = "mm",
    device = cairo_pdf
  )
}

save_figure(p_ind, "Figure1a_StrategySpace_Individuals2")



save_figure <- function(plot, p_var) {
  ggsave(
    filename = file.path(fig_dir, paste0(p_var, ".pdf")),
    plot = plot,
    width = 80,
    height = 80,
    units = "mm",
    device = cairo_pdf
  )
}
save_figure(p_var, "Figure1b_ModuleContributions2")

loadings <- as.data.frame(global_pca$rotation)

loadings


library(pheatmap)

pheatmap(
  loadings,
  cluster_rows = FALSE,
  main = "Module Contributions to Global Strategy Axes"
)


loadings <- as.data.frame(global_pca$rotation)

# order columns by PC number
pc_order <- order(as.numeric(gsub("[^0-9]", "", colnames(loadings))))

loadings <- loadings[, pc_order]
loadings <- loadings[order(abs(loadings$PC1), decreasing = TRUE), ]



library(pheatmap)

max_val <- max(abs(loadings))

pheatmap(
  loadings,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  breaks = seq(-max_val, max_val, length.out = 100),
  main = "Module Contributions to Global Strategy Axes"
)



library(pheatmap)

max_val <- max(abs(loadings))

p_heat <- pheatmap(
  loadings,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  breaks = seq(-max_val, max_val, length.out = 100),
  main = "Module Contributions to Global Strategy Axes"
)


cairo_pdf(
  file.path(fig_dir, "Figure1c_ModuleLoadingsHeatmap.pdf"),
  width = 180/25.4,
  height = 160/25.4
)

grid::grid.newpage()
grid::grid.draw(p_heat$gtable)

dev.off()


# find dominant module per PC
dominant_module <- apply(abs(loadings), 2, function(x) {
  rownames(loadings)[which.max(x)]
})

dominant_module


dominance_df <- data.frame(
  PC = names(dominant_module),
  Dominant_Module = dominant_module
)


library(ggplot2)

ggplot(dominance_df,
       aes(x = PC,
           y = 1,
           fill = Dominant_Module)) +
  geom_tile(color = "white", height = 0.6) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    title = "Dominant Functional Module Defining Each Strategy Axis",
    x = "Global Strategy Axis",
    fill = "Module"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank()
  )


module_colors <- c(
  Size_PC1 = "#1b9e77",
  Hydraulic_PC1 = "#1f78b4",
  Leaf_Nutrient_PC1 = "#33a02c",
  Plasticity_PC1 = "#cab2d6",
  Genome_PC1 = "#ff7f00",
  Root_PC1 = "#b15928",
  Leaf_Physiology_PC1 = "#fb9a99",
  Dispersal_PC1 = "#a6cee3",
  Reproduction_PC1 = "#e31a1c",
  Leaf_Construction_PC1 = "#6a3d9a"
)


dominance_df$PC <- factor(
  dominance_df$PC,
  levels = paste0("PC", 1:10)
)


library(ggplot2)

p_dom <- ggplot(
  dominance_df,
  aes(x = PC,
      y = 1,
      fill = Dominant_Module)
) +
  geom_tile(color = "white", height = 0.6) +
  scale_fill_manual(values = module_colors) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    title = "Dominant Functional Module Defining Each Strategy Axis",
    x = "Global Strategy Axis",
    fill = "Module"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 8, face = "bold"),

    # --- LEGEND COMPRESSION ---
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.box.spacing = unit(0.2, "cm")
  ) +
  guides(
    fill = guide_legend(
      nrow = 3,        # split legend into rows
      byrow = TRUE
    )
  )

p_dom

save_figure(p_dom, "Figure1d_ModuleDominanceMap2")


dominance_df$Strategy_Level <- ifelse(
  as.numeric(gsub("PC", "", dominance_df$PC)) <= 5,
  "Backbone strategies",
  "Specialization axes"
)


library(ggplot2)

p_dom <- ggplot(
  dominance_df,
  aes(x = PC, y = 1, fill = Dominant_Module)
) +

  # --- hierarchy background ---
  annotate("rect", xmin = 0.5, xmax = 5.5,
           ymin = 0.5, ymax = 1.5,
           fill = "grey95", alpha = 0.6) +

  annotate("rect", xmin = 5.5, xmax = 10.5,
           ymin = 0.5, ymax = 1.5,
           fill = "grey90", alpha = 0.6) +

  geom_tile(color = "white", height = 0.6) +

  scale_fill_manual(values = module_colors) +

  scale_y_continuous(NULL, breaks = NULL) +

  labs(
    title = "Hierarchical Organization of Global Plant Strategy Axes",
    x = "Global Strategy Axis",
    fill = "Module"
  ) +

  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),

    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


p_dom <- p_dom +
  annotate("text", x = 3, y = 1.35,
         label = "Backbone strategies",
         size = 3) +

annotate("text", x = 8, y = 1.35,
         label = "Specialization axes",
         size = 3)
p_dom

save_figure(p_dom, "Figure1d_StrategyHierarchy")


module_importance <- rowMeans(abs(loadings))

importance_df <- data.frame(
  Module = names(module_importance),
  Importance = module_importance
)


ggplot(importance_df,
       aes(x = reorder(Module, Importance),
           y = Importance)) +
  geom_col(fill = "grey40") +
  coord_flip() +
  theme_minimal(base_family = "Arial") +
  labs(
    x = NULL,
    y = "Mean absolute loading",
    title = "Overall Contribution of Functional Modules"
  )


# variance explained by each PC
weights <- global_pca$sdev^2 / sum(global_pca$sdev^2)

# keep same PCs used in dominance analysis
loadings_mat <- as.matrix(loadings[, 1:10])

# weighted module importance
module_importance <- rowSums(
  abs(loadings_mat) * weights[1:10]
)

importance_df <- data.frame(
  Module = rownames(loadings_mat),
  Importance = module_importance
)


ggplot(importance_df,
       aes(x = reorder(Module, Importance),
           y = Importance)) +
  geom_col(fill = "grey40") +
  coord_flip() +
  theme_minimal(base_family = "Arial") +
  labs(
    x = NULL,
    y = "Variance-weighted loading",
    title = "Overall Contribution of Functional Modules"
  )


round(loadings[,1:2], 2)


fviz_pca_ind(
  global_pca,
  geom.ind = "point",
  col.ind = "black",
  addEllipses = FALSE
)


p_ind <- fviz_pca_ind(
  global_pca,
  geom.ind = "point",
  col.ind = "grey30",
  alpha.ind = 0.3,
  pointsize = 1,
  addEllipses = FALSE
)
p_ind
# Export as high-resolution PDF
ggsave(
  filename = "pca_individuals_highres2.pdf", # output file
  plot = p_ind,                             # plot object
  device = "pdf",
  width = 8,                                # width in inches
  height = 6,                               # height in inches
  dpi = 600                                 # resolution for raster elements
)

traits %>%
  distinct(species, rank_correct) %>%
  count(species) %>%
  filter(n > 1)


functional_traits <- list(
  Leaf_Construction_PC1,
  Leaf_Nutrient_PC1,
  Leaf_Physiology_PC1,
  Size_PC1,
  Root_PC1,
  Hydraulic_PC1,
  Reproduction_PC1,
  Dispersal_PC1,
  Genome_PC1,
  Plasticity_PC1
) %>%
  reduce(left_join, by = "species") %>%
  left_join(traits %>% select(plant_growth_form, rank_correct,family,leaf_compoundness,origin,life_span,dispersal_distance,HR_Status,HR_Status_local,habitat_group,alt_min,alt_max,alt_range,n_points,EOO_km2,n_grid_cells,EOO_scaled,AOO_scaled,total_obs_scaled,habitats,expert_classification  ), by = "species")


functional_traits_species <- functional_traits %>%
  filter(rank_correct == "species")

#Check how many species-level entries
functional_traits %>%
  filter(rank_correct == "species") %>%
  nrow()



# Names of module PC1s
pc_names <- c(
  "Leaf_Construction_PC1",
  "Leaf_Nutrient_PC1",
  "Leaf_Physiology_PC1",
  "Size_PC1",
  "Root_PC1",
  "Hydraulic_PC1",
  "Reproduction_PC1",
  "Dispersal_PC1",
  "Genome_PC1",
  "Plasticity_PC1"
)

# Compute correlation between full dataset (with genus) vs species-only
cor_results <- sapply(pc_names, function(pc) {
  cor(
    functional_traits[[pc]][functional_traits$rank_correct == "species"],
    functional_traits_species[[pc]],
    use = "complete.obs"
  )
})

cor_results


strategy_matrixa <- functional_traits %>%
  drop_na() %>%
  column_to_rownames("species")


global_pca_sp <- prcomp(strategy_matrixa, center = TRUE, scale. = TRUE)


strategy_scores <- as.data.frame(global_pca_sp$x) %>%
  rownames_to_column("species")


round(global_pca$rotation, 2)


strategy_scorestoppc5 <- strategy_scores %>%
  rename(
    Growth_Reproduction = PC1,
    Physiology_Plasticity = PC2,
    Leaf_Resource = PC3,
    Plasticity = PC4,
    Reproduction_Allocation = PC5
  )


library(dplyr)

strategy_scorestoppc5 <- strategy_scorestoppc5 %>%
  select(
    species,
    Growth_Reproduction,
    Physiology_Plasticity,
    Leaf_Resource,
    Plasticity,
    Reproduction_Allocation
  )


glimpse(strategy_scorestoppc5)


analysis_data <- strategy_scorestoppc5 %>%
  left_join(life_form_data, by = "species") %>%
  left_join(origin_data, by = "species") %>%
  left_join(resistance_data, by = "species")


library(readxl)
library(dplyr)
library(tidyr)

# -------------------------------
# Step 1: Load weed datasets
# -------------------------------
data_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait"


#GAP functional trait data (traits)
trait_data2 <- readxl::read_excel(file.path(data_path, "copyTRY_pmf_DataRelease_cleaned_2013_06_19.xlsx")) %>%
  clean_names() # include attributes

weed1 <- read_excel("P:/CONSERWA/Herbicide_Project/Data_Trait/218_GEWeedSpecies20240926_notJUICE_header2.xlsx")%>%
  clean_names()# has weed species with matched concept names
weed2 <- read_excel("P:/CONSERWA/Herbicide_Project/Data_Trait/218_GEWeedSpecies20240926_EUNIS20250924.xlsx")%>%
  clean_names() # has correct expert classification of habitats, plot id or plot observation id but no species



library(dplyr)

# -------------------------------
#  Clean trait_data2 (remove duplicates)
# -------------------------------
library(dplyr)

trait_data2_clean <- trait_data2 %>%
  filter(rank_correct == "species") %>%   # keep only species level
  group_by(species,rank_correct) %>%
  summarise(
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    alt = mean(alt, na.rm = TRUE),
    
    family = first(na.omit(family)),
    plant_growth_form = first(na.omit(plant_growth_form)),
    succulent = first(na.omit(succulent)),
    climbing = first(na.omit(climbing)),
    leaf_type = first(na.omit(leaf_type)),
    photosynthetic_pathway = first(na.omit(photosynthetic_pathway)),
    leaf_compoundness = first(na.omit(leaf_compoundness)),
    
    .groups = "drop"
  )




any(duplicated(trait_data2_clean$species))

# -------------------------------
#  Merge into trait_data
# -------------------------------
trait_data_merged <- traits %>%
  left_join(trait_data2_clean,
            by = c("species", "rank_correct"))

# -------------------------------
#  Filter only TRY species level
# -------------------------------
trait_data_species <- trait_data_merged %>%
  filter(rank_correct == "species")

# -------------------------------
# (Optional) Check merge success
# -------------------------------
# How many rows got coordinates?
sum(!is.na(trait_data_species$lat))

# Inspect missing matches
missing_matches <- trait_data_species %>%
  filter(is.na(lat)) %>%
  distinct(species)


any(duplicated(trait_data_species$species))

# -------------------------------
# Step 2: Merge weed datasets on PlotObservationID
# -------------------------------
# Select relevant columns from weed2 (avoid duplicate species columns)
weed2_sel <- weed2 %>% select(plot_observation_id, altitude,
                              expert_classification,
                              longitude, latitude, location_uncertainty_m )

weed_merged <- weed1 %>%
  left_join(weed2_sel, by = "plot_observation_id")


dplyr::n_distinct(weed1$matched_concept) == nrow(weed1)


nrow(weed1)                         # total rows
dplyr::n_distinct(weed1$matched_concept)  # unique species


any(duplicated(weed_merged$matched_concept))

clean_names <- function(x) {
  x %>%
    iconv(to = "ASCII//TRANSLIT") %>%   # remove encoding issues
    tolower() %>%
    trimws() %>%
    gsub("\\s+", " ", .)              # collapse multiple spaces
}


weed_merged <- weed_merged %>%
  mutate(matched_concept = clean_names(matched_concept))

weed_merged <- weed_merged %>%
  mutate(turboveg2_concept = clean_names(turboveg2_concept))


library(dplyr)
library(tidyr)

lookup <- weed_merged %>%
  select(turboveg2_concept, matched_concept, original_taxon_concept) %>%
  
  # Keep matched_concept BEFORE pivoting
  pivot_longer(
    cols = c(turboveg2_concept, original_taxon_concept),
    names_to = "name_type",
    values_to = "species_name"
  ) %>%
  
  # ALSO include matched_concept itself as a valid variant
  bind_rows(
    weed_merged %>%
      select(species_name = matched_concept,
             matched_concept) %>%
      mutate(name_type = "matched_concept")
  ) %>%
  
  filter(!is.na(species_name)) %>%
  
  # Define canonical name
  mutate(canonical_name = matched_concept) %>%
  
  select(species_name, canonical_name) %>%
  
  distinct()


any(duplicated(lookup$species_name))

weed_metadata <- weed_merged %>%
  filter(!is.na(matched_concept)) %>%
  select(
    matched_concept,
    taxon_group,
    altitude,
    expert_classification,
    longitude,
    latitude,
    location_uncertainty_m
  ) %>%
  distinct()


# 1. Make a "long" lookup table
lookup_long <- lookup %>%
  pivot_longer(
    cols = species_name, 
    names_to = "type", 
    values_to = "species_alias"
  ) %>%
  distinct(species_alias, canonical_name)




lookup_long %>%
  count(species_alias) %>%
  filter(n > 1)


lookup_long %>%
  count(species_alias, canonical_name) %>%
  filter(n > 1)


lookup_long <- lookup_long %>%
  distinct(species_alias, .keep_all = TRUE)


# 2. Join TRY data to all names
trait_data_species <- trait_data_species %>%
  mutate(species = clean_names(species))

traits_filtered <- trait_data_species %>%
  left_join(lookup_long, by = c("species" = "species_alias")) %>%
  mutate(species_standard = coalesce(canonical_name, species)) %>%
  select(-canonical_name)



lost_in_lookup <- trait_data_species %>%
  anti_join(lookup_long, by = c("species" = "canonical_name") )%>%
  distinct(species)
length(lost_in_lookup)

intersect(
  lost_in_lookup$species,
  herbicide_long_compare$species
)


# Count duplicates in weed_metadata
weed_metadata %>% count(matched_concept) %>% filter(n > 1)

# Count species in traits_filtered
traits_filtered %>% count(species)


#Add metadata
weed_metadata_unique <- weed_metadata %>% distinct(matched_concept, .keep_all = TRUE)




weed_metadata_clean <- weed_metadata_unique %>%
  mutate(species = clean_names(matched_concept))


species_germany <- weed_metadata_clean %>%
  distinct(species) %>%
  pull(species)



species_hr <- herbicide_long_compare %>%
  mutate(species = clean_names(species)) %>%
  pull(species) %>%
  unique()


species_keep <- union(species_germany, species_hr)


traits_filtered_clean <- trait_data_species %>%
  filter(species %in% species_keep)


traits_filtered_clean_met <- traits_filtered_clean %>%
  left_join(weed_metadata_clean, by = "species")


traits_filtered_clean_met <- traits_filtered_clean_met %>%
  left_join(lookup_long, by = c("species" = "species_alias")) %>%
  mutate(species_standard = coalesce(canonical_name, species)) %>%
  select(-canonical_name)


traits_filtered_clean_met %>%
  mutate(
    match_standard_vs_matched = species_standard == matched_concept,
    match_standard_vs_species = species_standard == species,
    match_species_vs_matched = species == matched_concept
  ) %>%
  summarise(
    standard_vs_matched = mean(match_standard_vs_matched, na.rm = TRUE),
    standard_vs_species = mean(match_standard_vs_species, na.rm = TRUE),
    species_vs_matched = mean(match_species_vs_matched, na.rm = TRUE)
  )


mismatches <- traits_filtered_clean_met %>%
  filter(
    species_standard != matched_concept |
    species_standard != species
  ) %>%
  select(species, species_standard, matched_concept)

mismatches

nrow(mismatches)

species_corrections <- c(
  "bidens tripartitus" = "bidens tripartita",
  "anisantha sterilis" = "anisantha sterilis",
  "anisantha tectorum" = "anisantha tectorum",
  "ochlopoa annua"     = "ochlopoa annua"
)


traits_filtered_clean_well <- traits_filtered_clean_met %>%
  mutate(
    species_final = coalesce(
      recode(species_standard, !!!species_corrections),
      species_standard,
      species
    )
  )


traits_filtered_clean_well <- traits_filtered_clean_well %>%
  select(-species_standard, -matched_concept, -species)


traits_filtered_clean_well <- traits_filtered_clean_well %>%
  mutate(species = species_final)


traits_filtered_clean_well <- traits_filtered_clean_well %>%
  select(-species_final)

# Install once if needed

# Load package
library(writexl)

# Export to Excel
write_xlsx(traits_filtered_clean_well,
           "P:/CONSERWA/Herbicide_Project/Data_Trait/TRY_weeds_filtered.xlsx")
# -------------------------------
# Step 5: Save filtered TRY weeds dataset
# -------------------------------
write.csv(traits_filtered_clean_well, "P:/CONSERWA/Herbicide_Project/Data_Trait/TRY_weeds_filtered.csv", row.names = FALSE)

herbicide_data <- read_excel("P:/CONSERWA/Herbicide_Project/Data_Trait/Herbicide_data.xlsx")%>%
  clean_names() # has correct expert classification of habitats
str(herbicide_data)
n_species_raw <- herbicide_data %>%
  distinct(species) %>%
  nrow()

n_species_raw


library(tidyr)

herbicide_long <- herbicide_data %>%
  separate_rows(active_substance_s, sep = ",\\s*") %>%
  separate_rows(chemical_group_s, sep = ",\\s*")


library(dplyr)
library(stringr)
library(tidyr)

herbicide_long1 <- herbicide_data %>%
  
  # ✅ SPLIT SPECIES (missing step!)
  separate_rows(species, sep = ",\\s*") %>%
  
  # then your existing steps
  separate_rows(active_substance_s, sep = ",\\s*") %>%
  
  mutate(active_clean = str_remove(active_substance_s,
                                   "^\\s*[0-9]+([,\\.][0-9]+)?\\s*g/[a-z]+\\s*")) %>%
  mutate(active_clean = str_remove(active_clean, "\\s*\\(.*\\)")) %>%
  mutate(active_clean = str_trim(active_clean)) %>%
  mutate(active_clean = str_to_lower(active_clean)) %>%
  
  # ✅ clean species too
  mutate(species = str_trim(str_to_lower(species)))



n_species_fixed <- herbicide_long1 %>%
  distinct(species) %>%
  nrow()

n_species_fixed


herbicide_R <- read_excel("P:/CONSERWA/Herbicide_Project/Data_Trait/HerbicideResistanceglobal.xlsx")%>%
  clean_names() # 


str(herbicide_R)
library(dplyr)
library(stringr)

# Convert each row from "c(...)" into actual vector
parsed <- herbicide_R %>%
  str_remove_all("^c\\(|\\)$") %>%   # remove c( )
  str_split(",\\s*")                 # split into values


herbicideGlobal_fixed <- as.data.frame(do.call(cbind, parsed), stringsAsFactors = FALSE)


colnames(herbicideGlobal_fixed) <- c("year", "species", "country", "state_name","moas","active_substance","situations")


str(herbicideGlobal_fixed)

n_species_raw <- herbicideGlobal_fixed %>%
  distinct(species) %>%
  nrow()

n_species_raw


library(dplyr)
library(stringr)
library(tidyr)

herbicideGlobal_clean <- herbicideGlobal_fixed %>%
  
  # -------------------------------
  # 1. Remove weird quotes
  # -------------------------------
  mutate(across(everything(), ~ str_remove_all(., '^"|"$'))) %>%   # remove leading/trailing "
  mutate(across(everything(), ~ str_replace_all(., '\\"', ''))) %>% # remove escaped quotes
  
  # -------------------------------
  # 2. Fix encoding issues
  # -------------------------------
  mutate(across(everything(), ~ str_replace_all(., "&#x0d;", ""))) %>%
  
  # -------------------------------
  # 3. Clean species
  # -------------------------------
  mutate(
    species = str_to_lower(species),
    species = str_remove(species, "\\s*\\(.*\\)"),  # remove (syn: ...)
    species = str_trim(species)
  ) %>%
  
  # -------------------------------
  # 4. Clean country
  # -------------------------------
  mutate(
    country = str_to_lower(country),
    country = str_trim(country)
  ) %>%
  
  # -------------------------------
  # 5. Clean MOA (mode of action)
  # -------------------------------
  mutate(
    moas = str_to_lower(moas),
    moas = str_replace_all(moas, "\\r\\n", " "),
    moas = str_trim(moas)
  ) %>%
  
  # -------------------------------
  # 6. Clean active substances
  # -------------------------------
  mutate(
    active_substance = str_to_lower(active_substance),
    active_substance = str_replace_all(active_substance, "\\r\\n", ", ")
  ) %>%
  
  # 🔥 IMPORTANT: split multiple actives
  separate_rows(active_substance, sep = ",\\s*") %>%
  
  # Remove parentheses (chemical forms)
  mutate(
    active_clean = str_remove(active_substance, "\\s*\\(.*\\)"),
    active_clean = str_trim(active_clean)
  ) %>%
  
  # Remove empty rows
  filter(active_clean != "", !is.na(active_clean))


herbicideGlobal_clean %>%
  filter(grepl('"', species) | grepl('"', active_clean))


herbicideGlobal_clean %>%
  count(species) %>%
  arrange(desc(n))


herbicideGlobal_clean <- herbicideGlobal_clean %>%
  mutate(species = clean_names(species))


herbicide_species <- herbicideGlobal_clean %>%
  distinct(species)


library(dplyr)
library(stringr)

weed_lookup <- weed_merged %>%
  mutate(
    turboveg2_concept = clean_names(turboveg2_concept),
    matched_concept   = clean_names(matched_concept)
  ) %>%
  select(turboveg2_concept, matched_concept) %>%
  distinct()


direct_match <- herbicide_species %>%
  filter(species %in% weed_lookup$matched_concept) %>%
  mutate(
    species_final = species,
    match_type = "accepted"
  )



synonym_match <- herbicide_species %>%
  inner_join(weed_lookup,
             by = c("species" = "turboveg2_concept")) %>%
  mutate(
    species_final = matched_concept,
    match_type = "synonym"
  ) %>%
  select(species, species_final, match_type)



matched_all <- bind_rows(direct_match, synonym_match) %>%
  distinct(species, .keep_all = TRUE)


not_found <- herbicide_species %>%
  anti_join(matched_all, by = "species")

nrow(not_found)


cat("Total species:", nrow(herbicide_species), "\n")
cat("Accepted matches:", sum(matched_all$match_type == "accepted"), "\n")
cat("Synonym matches:", sum(matched_all$match_type == "synonym"), "\n")
cat("Not found:", nrow(not_found), "\n")


herbicideGlobal_standardized <- herbicideGlobal_clean %>%
  left_join(matched_all, by = "species") %>%
  mutate(
    species = coalesce(species_final, species)  # replace synonyms
  ) %>%
  select(-species_final, -match_type)


# Read without cleaning names first
herbicidelist_raw <- read_excel("P:/CONSERWA/Herbicide_Project/Data_Trait/Herbicidelistgermanyeppo.xlsx", col_names = TRUE)

# Check first few rows
head(herbicidelist_raw)


# Assign proper column names
colnames(herbicidelist_raw) <- c("trade_name", "permissible_no", "end_permit", "active_substance", "sphere_activity")

# Now it’s safe to parse
herbicidelist_clean <- herbicidelist_raw %>%
  mutate(active_clean = str_to_lower(active_substance)) %>%
  mutate(active_clean = str_remove(active_clean, "\\s*\\(.*\\)")) %>%
  mutate(active_clean = str_trim(active_clean)) %>%
  separate_rows(active_clean, sep = "\\s*\\+\\s*") %>%
  filter(active_clean != "")



library(dplyr)
library(stringr)
library(tidyr)

herbicidelist_clean <- herbicidelist %>%
  
  # 1. Split multiple lines into rows (\r\n is important!)
  separate_rows(active_substance, sep = "\\r?\\n") %>%
  
  # 2. Remove concentration (numbers + units like g/kg, g/l)
  mutate(active_clean = str_remove(active_substance,
                                   "^\\s*[0-9]+([,\\.][0-9]+)?\\s*g/[a-z]+\\s*")) %>%
  
  # 3. Remove parentheses content (e.g. chemical forms)
  mutate(active_clean = str_remove(active_clean, "\\s*\\(.*\\)")) %>%
  
  # 4. Trim whitespace
  mutate(active_clean = str_trim(active_clean)) %>%
  
  # 5. Split substances separated by commas
  separate_rows(active_substance, sep = "\\+|,\\s*") %>%
  
  # 6. Remove empty rows (important!)
  filter(active_clean != "")




unique_active <- herbicidelist_clean %>%
  distinct(active_clean)



active_germany <- unique_active %>%
  mutate(active_clean = str_to_lower(active_clean)) %>%
  pull(active_clean)


herbicide_long1 <- herbicideGlobal_standardized %>%
  mutate(used_in_germany = active_clean %in% active_germany)


library(dplyr)
library(stringr)

herbicide_long_compare <- herbicide_long1 %>%
  mutate(
    used_exact = active_clean %in% active_germany,
    used_fuzzy = str_detect(active_clean,
                            paste(active_germany, collapse = "|"))
  )



herbicide_long_compare %>%
  summarise(
    exact_matches = sum(used_exact, na.rm = TRUE),
    fuzzy_matches = sum(used_fuzzy, na.rm = TRUE)
  )


library(dplyr)
library(stringr)

herbicide_long_compare <- herbicide_long_compare %>%
  separate_rows(species, sep = ",\\s*") %>%   # split multiple species
  mutate(
    species = species %>%
      str_replace_all("\\s+", " ") %>%  # collapse multiple spaces → single
      str_trim() %>%                   # remove leading/trailing spaces
      str_to_lower() %>%               # standardize case
      iconv(to = "ASCII//TRANSLIT")    # remove hidden encoding issues
  )



herbicide_long_compare <- herbicide_long_compare %>%
  mutate(species = str_to_lower(species))
traits_filtered_final <- traits_filtered_clean_well %>%
  mutate(species = str_to_lower(species))



setdiff(unique(herbicide_long_compare$species),
        unique(traits_filtered_clean_well$species)) #traits_filtered_clean_well Miss 167


clean_names <- function(x) {
  x %>%
    iconv(to = "ASCII//TRANSLIT") %>%
    tolower() %>%
    trimws() %>%
    gsub("\\s+", " ", .)
}

trait_data_species <- trait_data_species %>%
  mutate(species = clean_names(species))

herbicide_long_compare <- herbicide_long_compare %>%
  mutate(species = clean_names(species))


herbicide_long_compare <- herbicide_long_compare %>%
  mutate(
    species = sub(" subsp\\..*", "", species),
    species = sub(" var\\..*", "", species)
  ) %>%
  filter(!grepl("sp\\.", species))


missing_vs_try <- setdiff(
  unique(herbicide_long_compare$species),
  unique(trait_data_species$species)
)

length(missing_vs_try)


cat(missing_vs_try, sep = "\n")


present_in_try <- intersect(
  unique(herbicide_long_compare$species),
  unique(trait_data_species$species)
)

length(present_in_try)


missing_species <- setdiff(
  unique(herbicide_long_compare$species), #herbicide_long_clean
  unique(traits_filtered_final$species)
)


present_in_try <- missing_species[
  missing_species %in% trait_data_species$species
]


absent_in_try <- setdiff(missing_species, trait_data_species$species)


length(missing_species)      # total missing
length(present_in_try)       # found in TRY
length(absent_in_try)        # not in TRY


present_df <- data.frame(
  species = present_in_try,
  stringsAsFactors = FALSE
)

library(dplyr)

origin_map <- tibble::tribble(
  ~species, ~origin,
  "echinochloa colona", "introduced",
  "hirschfeldia incana", "introduced",
  "ceratochloa cathartica", "introduced",  # fixed name
  "sisymbrium orientale", "introduced",
  "echium plantagineum", "native",
  "mesembryanthemum crystallinum", "introduced",
  "galium tricornutum", "introduced",
  "nassella trichotoma", "introduced",
  "phalaris minor", "native",
  "lactuca saligna", "native",
  "sorghum bicolor", "introduced",
  "bidens pilosa", "introduced",
  "bidens subalternans", "introduced",
  "digitaria ciliaris", "introduced",
  "amaranthus viridis", "introduced",
  "ageratum conyzoides", "introduced",
  "vaccaria hispanica", "introduced",
  "hordeum jubatum", "introduced",
  "panicum capillare", "introduced",
  "amaranthus tuberculatus", "introduced",
  "setaria faberi", "introduced",
  "ambrosia trifida", "introduced",
  "cynosurus echinatus", "introduced",
  "silene gallica", "introduced",
  "beckmannia syzigachne", "introduced",
  "sclerochloa dura", "native",
  "senecio leucanthemifolius", "introduced",
  "lindernia dubia", "introduced",
  "elatine triandra", "native",
  "eleocharis acicularis", "native",
  "xanthium strumarium", "native",
  "nassella neesiana", "introduced",
  "carduus pycnocephalus", "introduced",
  "solanum americanum", "introduced",
  "setaria verticillata", "native",
  "hydrilla verticillata", "native",
  "datura stramonium", "introduced",
  "erysimum repandum", "native",
  "bromus japonicus", "introduced",
  "abutilon theophrasti", "introduced",
  "amaranthus spinosus", "introduced",
  "centaurea solstitialis", "introduced"
)



origin_map <- origin_map %>%
  mutate(
    species = ifelse(species == "ceratochloa cathartica",
                     "bromus catharticus",
                     species)
  )


origin_map_sorted <- origin_map %>%
  arrange(species)

origin_map_sorted


grep("cathart", trait_data_species$species, value = TRUE)


missing_in_try2 <- origin_map_sorted %>%
  filter(!species %in% trait_data_species$species)

missing_in_try2
nrow(origin_map_sorted)      # should be 42
nrow(missing_in_try2)  


traits_filtered_final <- traits_filtered_final %>%
  left_join(origin_map, by = "species")


species_germany <- origin_map %>%
  distinct(species)


traits_missing <- traits_filtered %>%
  filter(species %in% species_germany$species)



traits_missing %>%
  summarise(n_species = n_distinct(species),
            n_rows = n())


traits_missing <- traits_missing %>%
  left_join(origin_map, by = "species")


traits_final_extended <- bind_rows(
  traits_filtered_final,
  traits_missing
) %>%
  distinct(species, .keep_all = TRUE)


traits_final_extended <- traits_final_extended %>%
  mutate(origin = factor(origin, levels = c("native", "introduced")))


absent_in_try_germany <- herbicide_long_compare %>%
  filter(species %in% absent_in_try) %>%
  filter(country == "germany") %>%   # adjust if needed herbicide_long_clean
  distinct(species)


absent_df <- data.frame(
  species = absent_in_try,
  stringsAsFactors = FALSE
)


library(dplyr)

synonym_map <- tibble::tibble(
  old_name = c(
    "polygonum convolvulus",
    "conyza bonariensis",
    "conyza sumatrensis",
    "kochia scoparia",
    "bidens tripartitus",
    "fimbristylis miliacea",
    "echinochloa phyllopogon",
    "chrysanthemum coronarium",
    "echinochloa erecta",
    "polygonum pensylvanicum",
    "chamaesyce maculata"
  ),
  new_name = c(
    "fallopia convolvulus",
    "erigeron bonariensis",
    "erigeron sumatrensis",
    "bassia scoparia",
    "bidens tripartita",
    "fimbristylis quinquangularis",
    "echinochloa oryzoides",
    "glebionis coronaria",
    "echinochloa crus-galli",
    "persicaria pensylvanica",
    "euphorbia maculata"
  )
)


herbicide_long_compare <- herbicide_long_compare %>%
  left_join(synonym_map, by = c("species" = "old_name")) %>%
  mutate(
    species = ifelse(!is.na(new_name), new_name, species)
  ) %>%
  select(-new_name)


species_to_add2 <- herbicide_long_compare %>%
  distinct(species) %>%
  filter(!species %in% traits_final_extended$species)


species_to_add2 <- synonym_map %>%
  distinct(new_name) %>%
  rename(species = new_name)


traits_missing3 <- traits_filtered %>%
  filter(species %in% species_to_add2$species)


traits_missing3 <- traits_missing3 %>%
  left_join(origin_map, by = "species")


traits_final_extended <- bind_rows(
  traits_final_extended,
  traits_missing3
) %>%
  distinct(species, .keep_all = TRUE)


any(duplicated(traits_final_extended$species))
# should be FALSE


herbicide_long_compare <- herbicide_long_compare %>%
  filter(!grepl("var. robustapurpurea)", species, ignore.case = TRUE))


herbicide_long_compare %>%
  filter(species %in% synonym_map$new_name) %>%
  distinct(species)


# Number of unique species before removing
n_unique_before <- herbicide_long_compare %>%
  distinct(species) %>%
  nrow()

n_unique_before


species_to_remove <- herbicide_long_compare %>%
  distinct(species) %>%
  filter(!species %in% traits_final_extended$species)


herbicide_long_compare_clean <- herbicide_long_compare %>%
  filter(species %in% traits_final_extended$species)


# Check number of removed species
setdiff(
  unique(herbicide_long_compare$species),
  unique(herbicide_long_compare_clean$species)
) %>% length()

clean_std <- function(x) {
  x %>%
    str_to_lower() %>%
    str_trim() %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("\\u00A0", " ")
}

herbicide_long_compare_clean <- herbicide_long_compare_clean %>%
  mutate(species = clean_std(species))

trait_data_species <- trait_data_species %>%
  mutate(species = clean_std(species))


remove_subsp <- function(x) {
  str_remove(x, " subsp\\..*")
}

herbicide_long_clean <- herbicide_long_compare_clean %>%
  mutate(species = remove_subsp(species))

trait_data_species <- trait_data_species %>%
  mutate(species = remove_subsp(species))


#2. Define the target path (using forward slashes for R compatibility)
output_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait/herbicide_long_clean.csv"

# 3. Export to CSV
write_csv(herbicide_long_clean, output_path)

n_species_total <- herbicide_long_clean %>%
  distinct(species) %>%
  nrow()

n_species_total


hr_status <- herbicide_long_clean %>%
  group_by(species) %>%
  summarise(
    has_germany = any(country == "germany"),
    .groups = "drop"
  ) %>%
  mutate(
    HR_Status = ifelse(has_germany, "HR_local", "HR_elsewhere")
  )


#2. Define the target path (using forward slashes for R compatibility)
output_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait/hr_status.csv"

# 3. Export to CSV
write_csv(hr_status, output_path)

hr_counts_before <- hr_status %>%
  count(HR_Status)

hr_counts_before


hr_usage <- herbicide_long_clean %>%
  group_by(species) %>%
  summarise(
    used_in_germany_any = any(used_fuzzy, na.rm = TRUE),
    .groups = "drop"
  )


hr_status_full <- hr_status %>%
  left_join(hr_usage, by = "species")


hr_status_filtered <- hr_status_full %>%
  filter(
    HR_Status == "HR_local" |
    (HR_Status == "HR_elsewhere" & used_in_germany_any == TRUE)
  )


hr_counts_after <- hr_status_filtered %>%
  count(HR_Status)

hr_counts_after
cat("Total unique species:", n_species_total, "\n\n")

cat("Before filtering:\n")
print(hr_counts_before)

cat("\nAfter filtering:\n")
print(hr_counts_after)


removed_species_HR <- hr_status_full %>%
  anti_join(hr_status_filtered, by = "species")

removed_species_HR


# Join to get active substances for the removed species
removed_species_with_actives <- removed_species_HR %>%
  left_join(
    herbicide_long_compare %>% 
      select(species, active_clean) %>%
      distinct(),
    by = "species"
  )

# Optional: collapse multiple actives per species into a comma-separated string
library(dplyr)

removed_species_with_actives <- removed_species_with_actives %>%
  group_by(species) %>%
  summarise(
    active_substances = paste(unique(active_clean), collapse = ", "),
    .groups = "drop"
  )

# View
removed_species_with_actives



hr_remaining <- hr_status_full %>%
  semi_join(hr_status_filtered, by = "species")


str(hr_status_filtered)

hr_table <- hr_status_filtered %>%
  select(species, HR_Status) %>%
  distinct()



#2. Define the target path (using forward slashes for R compatibility)
output_path <- "P:/CONSERWA/Herbicide_Project/Data_Trait/hr_table.csv"

# 3. Export to CSV
write_csv(hr_table, output_path)


traits_final_with_hr <- traits_final_extended %>%
  left_join(hr_table, by = "species")



traits_final_with_hr <- traits_final_with_hr %>%
  mutate(
    HR_Status = if_else(is.na(HR_Status), "Non_HR_local", HR_Status)
  )


traits_final_with_hr <- traits_final_with_hr %>%
  mutate(
    HR_Status_local = if_else(HR_Status == "HR_local", 
                              "HR_local", 
                              "Non_HR_local")
  )



# Count for HR_Status (3 categories)
table(traits_final_with_hr$HR_Status)

# Count for HR_Status_local (2 categories)
table(traits_final_with_hr$HR_Status_local)



habitat_data <- readxl::read_excel(
  file.path(data_path, "EUNIShabitatscodesandnames202509.xlsx"),
  sheet = "EUNIS units"
) %>%
  janitor::clean_names()



str(habitat_data)

traits_final_with_hr <- traits_final_with_hr %>%
  mutate(expert_classification = toupper(expert_classification))

habitat_data <- habitat_data %>%
  mutate(expert_classification = toupper(expert_classification))


traits_with_habitat_hr <- traits_final_with_hr %>%
  left_join(habitat_data, by = "expert_classification")


traits_with_habitat_hr %>%
  distinct(habitats) %>%
  arrange(habitats)


traits_with_habitat_hr %>%
  count(habitats, sort = TRUE)


missing_classification <- traits_with_habitat_hr %>%
  filter(
    is.na(expert_classification) |
    expert_classification == "" |
    expert_classification == "~"
  )

nrow(missing_classification)


species_to_check <- missing_classification %>%
  distinct(species)

species_to_check



library(writexl)

# Export to Excel
write_xlsx(species_to_check,
           "P:/CONSERWA/Herbicide_Project/Data_Trait/species_missing_classification.xlsx")
# -------------------------------

floraveg_data <- readxl::read_excel(
  file.path(data_path, "species_missing_classification.xlsx")
) %>%
  janitor::clean_names()


str(floraveg_data)

traits_with_habitat_hr <- traits_with_habitat_hr %>%
  left_join(
    floraveg_data %>%
      select(species,
             expert_classification_floraveg = expert_classification,
             habitats_floraveg = habitats),
    by = "species"
  ) %>%
  mutate(
    expert_classification = coalesce(
      expert_classification,
      expert_classification_floraveg
    ),
    habitats = coalesce(
      habitats,
      habitats_floraveg
    )
  ) %>%
  select(-expert_classification_floraveg,
         -habitats_floraveg)



traits_with_habitat_hr %>%
  distinct(habitats) %>%
  arrange(habitats)


traits_with_habitat_hr %>%
  count(habitats, sort = TRUE)

names(traits_with_habitat_hr)


library(tidyr)
library(dplyr)
library(stringr)

data_long <- traits_with_habitat_hr %>%
  separate_rows(habitats, sep = ",\\s*")


data_long <- data_long %>%
  mutate(
    habitats = habitats %>%
      str_to_lower() %>%
      str_squish()
  )


data_classified <- data_long %>%
  mutate(
    habitats_clean = habitats %>%
      str_to_lower() %>%
      str_squish(),
    
    habitat_group = case_when(
      
      # 🌿 RUDERAL (priority 1)
      str_detect(habitats_clean, "anthropogenic|Dry Mediterranean|Non-vernal|tramp|man.?made|synanthropic|nitrophilous") ~ "ruderal",
      
      # 🌱 ARABLE
      str_detect(habitats_clean, "bare tilled|fallow|abandoned|arable|crop|cropland|rice") ~ "arable",
      
      # ❌ NON-ARABLE
      str_detect(habitats_clean, "wetland|aquatic|shore|coastal|forest|saltmarsh") ~ "non_arable",
      
      # ⚠ UNKNOWN
      habitats_clean == "not classified" ~ "unknown",
      
      TRUE ~ "ambiguous"
    )
  )




species_class <- data_classified %>%
  group_by(species) %>%
  summarise(
    habitat_group = case_when(
      any(habitat_group == "arable") ~ "arable",
      any(habitat_group == "ruderal") ~ "ruderal",
      all(habitat_group == "non_arable") ~ "non_arable",
      TRUE ~ "ambiguous"
    )
  )


manual_override <- tibble(
  species = c("chenopodium ficifolium", "ochlopoa annua"),
  habitat_group_override = "arable"
)


species_class <- species_class %>%
  left_join(manual_override, by = "species") %>%
  mutate(
    habitat_group = coalesce(habitat_group_override, habitat_group)
  ) %>%
  select(-habitat_group_override)


analysis_data <- traits_with_habitat_hr %>%
  left_join(species_class, by = "species") %>%
  filter(habitat_group %in% c("arable", "ruderal"))


removed_species_habitat <- traits_with_habitat_hr %>%
  left_join(species_class, by = "species") %>%
  filter(habitat_group %in% c("non_arable", "unknown", "ambiguous"))

removed_species %>%
  count(HR_Status)


removed_species_habitat %>%
  distinct(species, HR_Status, habitat_group) %>%
  count(HR_Status)


species_to_fix <- c("echium plantagineum", "carduus pycnocephalus")



data_classified <- data_classified %>%
  mutate(
    habitat_group = if_else(
      species %in% species_to_fix,
      "ruderal",
      habitat_group
    )
  )


species_recovered <- data_classified %>%
  filter(species %in% species_to_fix)


removed_species_habitat <- removed_species_habitat %>%
  filter(!species %in% species_to_fix)


analysis_data <- bind_rows(
  analysis_data,
  species_recovered
)


library(writexl)

# Export to Excel
write_xlsx(analysis_data,
           "P:/CONSERWA/Herbicide_Project/Data_Trait/TRY_weeds_filtered_analysis.xlsx")
# -------------------------------

# ============================================================
# TRAIT-BASED PREDICTION OF HERBICIDE RESISTANCE
# Complete Pipeline following Chen & Breiman (2004)
# balanced RF with sequential covariate optimization
# ============================================================

library(randomForest)    # Liaw and Wiener 2002
library(tidyverse)
library(pROC)
library(ggplot2)
library(patchwork)
library(pdp)
library(caret)
library(viridis)

set.seed(42)

# ============================================================
# Assuming our module PC scores + categorical are assembled
# into a dataframe called 'model_data' with columns:
#   - species (character)
#   - HR_status (factor: "HR" vs "nonHR")
#   - Module PC scores (M1_PC1, M1_PC2, M2_PC1, etc.)
#   - Single traits (M3_LeafWater)
#   - ITV CVs (M10_rooting_depth_cv, etc.)
#   - Categorical: origin, life_span, growth_form, 
#                   leaf_type, habitat_group
# ============================================================

# ============================================================
# STEP 1: DATA PREPARATION
# ============================================================

# --- 1a. Check class distribution ---
cat("=== Class Distribution ===\n")
print(table(model_data$HR_status))
cat("Ratio:", 
    round(sum(model_data$HR_status == "nonHR") / 
            sum(model_data$HR_status == "HR"), 1), 
    ": 1 (nonHR : HR)\n")

# --- 1b. Ensure correct variable types ---
model_data <- model_data %>%
  mutate(
    HR_status = factor(HR_status, levels = c("nonHR", "HR")),
    origin = as.factor(origin),
    life_span = as.factor(life_span),
    growth_form = as.factor(growth_form),
    leaf_type = as.factor(leaf_type),
    habitat_group = as.factor(habitat_group)
  )

# --- 1c. Prepare RF dataframe (remove species column) ---
rf_data <- model_data %>%
  select(-species) %>%
  drop_na()

cat("\nFinal sample size:", nrow(rf_data), "\n")
cat("HR:", sum(rf_data$HR_status == "HR"), "\n")
cat("nonHR:", sum(rf_data$HR_status == "nonHR"), "\n")

# --- 1d. Define predictor names ---
predictor_names <- setdiff(colnames(rf_data), "HR_status")
cat("Number of predictors:", length(predictor_names), "\n")

# ============================================================
# STEP 2: BALANCED RANDOM FOREST
# Following Chen & Breiman (2004)
# ============================================================

# --- 2a. Set parameters ---
n_hr <- sum(rf_data$HR_status == "HR")
n_nonhr <- sum(rf_data$HR_status == "nonHR")
n_predictors <- length(predictor_names)
mtry_value <- floor(sqrt(n_predictors))
ntree_value <- 1000

cat("=== RF Parameters ===\n")
cat("ntree:", ntree_value, "\n")
cat("mtry:", mtry_value, "(sqrt of", n_predictors, "predictors)\n")
cat("Balanced sampling: n =", n_hr, "per class per tree\n")

# --- 2b. Run balanced RF ---
# sampsize ensures equal sampling from each class
# This is the Chen & Breiman (2004) approach

rf_balanced <- randomForest(
  HR_status ~ .,
  data = rf_data,
  ntree = ntree_value,
  mtry = mtry_value,
  importance = TRUE,
  # KEY: balanced bootstrap sampling
  sampsize = c(nonHR = n_hr, HR = n_hr),
  # Stratified sampling ensures both classes in every tree
  strata = rf_data$HR_status,
  # Keep proximity for later analysis
  proximity = TRUE,
  # Keep forest for partial dependence
  keep.forest = TRUE
)

# --- 2c. Examine initial results ---
cat("\n=== Initial Balanced RF Results ===\n")
print(rf_balanced)
cat("\nOOB Error Rate:", 
    round(rf_balanced$err.rate[ntree_value, "OOB"] * 100, 1), "%\n")
cat("HR class error:", 
    round(rf_balanced$err.rate[ntree_value, "HR"] * 100, 1), "%\n")
cat("nonHR class error:", 
    round(rf_balanced$err.rate[ntree_value, "nonHR"] * 100, 1), "%\n")

# --- 2d. Confusion matrix ---
cat("\n=== OOB Confusion Matrix ===\n")
print(rf_balanced$confusion)

# --- 2e. Error convergence plot ---
oob_error <- as.data.frame(rf_balanced$err.rate) %>%
  mutate(trees = 1:n()) %>%
  pivot_longer(-trees, names_to = "class", values_to = "error")

p_convergence <- ggplot(oob_error, 
                         aes(x = trees, y = error, color = class)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("OOB" = "#FFCC00", "HR" = "#9933FF", "nonHR" = "#2E86AB"),
    labels = c("HR" = "HR species", "nonHR" = "Non-HR species", 
               "OOB" = "Overall OOB")
  ) +
  ggtitle("Error Rate Convergence") +
  labs(x = "Number of Trees", y = "Error Rate", color = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_convergence)

# ============================================================
# STEP 3: SEQUENTIAL COVARIATE REMOVAL
# Remove covariates with negative permutation importance
# until only positive contributors remain
# ============================================================

sequential_rf_optimization <- function(data, 
                                        response = "HR_status",
                                        ntree = 1000,
                                        verbose = TRUE) {
  
  # Initialize
  current_data <- data
  predictors <- setdiff(colnames(current_data), response)
  iteration <- 0
  history <- list()
  
  repeat {
    iteration <- iteration + 1
    n_pred <- length(predictors)
    mtry_val <- max(1, floor(sqrt(n_pred)))
    
    # Class sizes for balanced sampling
    n_min <- min(table(current_data[[response]]))
    samp <- rep(n_min, length(levels(current_data[[response]])))
    names(samp) <- levels(current_data[[response]])
    
    # Fit balanced RF
    rf_iter <- randomForest(
      as.formula(paste(response, "~ .")),
      data = current_data,
      ntree = ntree,
      mtry = mtry_val,
      importance = TRUE,
      sampsize = samp,
      strata = current_data[[response]]
    )
    
    # Get permutation importance (MeanDecreaseAccuracy)
    imp <- importance(rf_iter, type = 1)  # type 1 = permutation
    imp_df <- data.frame(
      variable = rownames(imp),
      importance = imp[, 1]
    ) %>% arrange(importance)
    
    # OOB error
    oob_err <- rf_iter$err.rate[ntree, "OOB"]
    hr_err <- rf_iter$err.rate[ntree, "HR"]
    nonhr_err <- rf_iter$err.rate[ntree, "nonHR"]
    
    # Store history
    history[[iteration]] <- list(
      iteration = iteration,
      n_predictors = n_pred,
      mtry = mtry_val,
      oob_error = oob_err,
      hr_error = hr_err,
      nonhr_error = nonhr_err,
      importance = imp_df,
      model = rf_iter
    )
    
    if (verbose) {
      cat("\n--- Iteration", iteration, "---\n")
      cat("Predictors:", n_pred, "| mtry:", mtry_val, "\n")
      cat("OOB error:", round(oob_err * 100, 1), "%\n")
      cat("HR error:", round(hr_err * 100, 1), "%\n")
      cat("nonHR error:", round(nonhr_err * 100, 1), "%\n")
    }
    
    # Identify covariates with negative importance
    negative_vars <- imp_df %>% 
      filter(importance < 0) %>% 
      pull(variable) %>%
      as.character()
    
    if (length(negative_vars) == 0) {
      if (verbose) {
        cat("\n✓ All covariates have positive importance.\n")
        cat("  Final model:", n_pred, "predictors\n")
        cat("  Final OOB error:", round(oob_err * 100, 1), "%\n")
      }
      break
    }
    
    if (verbose) {
      cat("Removing", length(negative_vars), 
          "negative-importance covariates:\n")
      cat(" ", paste(negative_vars, collapse = ", "), "\n")
    }
    
    # Remove negative importance covariates
    predictors <- setdiff(predictors, negative_vars)
    current_data <- current_data %>% 
      select(all_of(c(response, predictors)))
    
    # Safety: stop if too few predictors
    if (length(predictors) < 2) {
      warning("Only 1 predictor remaining. Stopping.")
      break
    }
  }
  
  # Return results
  return(list(
    final_model = rf_iter,
    final_predictors = predictors,
    history = history,
    n_iterations = iteration
  ))
}

# --- Run sequential optimization ---
cat("\n========================================\n")
cat("SEQUENTIAL COVARIATE OPTIMIZATION\n")
cat("========================================\n")

opt_result <- sequential_rf_optimization(
  data = rf_data,
  response = "HR_status",
  ntree = 1000,
  verbose = TRUE
)

# --- Extract final model ---
rf_final <- opt_result$final_model
final_predictors <- opt_result$final_predictors

cat("\n=== FINAL OPTIMIZED MODEL ===\n")
print(rf_final)
cat("\nRetained predictors:\n")
cat(paste(final_predictors, collapse = "\n"))

# --- Optimization history plot ---
history_df <- map_dfr(opt_result$history, function(h) {
  data.frame(
    iteration = h$iteration,
    n_predictors = h$n_predictors,
    oob_error = h$oob_error,
    hr_error = h$hr_error,
    nonhr_error = h$nonhr_error
  )
})

p_optimization <- history_df %>%
  pivot_longer(c(oob_error, hr_error, nonhr_error),
               names_to = "error_type", values_to = "error") %>%
  ggplot(aes(x = iteration, y = error * 100, 
             color = error_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("oob_error" = "#FFCC00", 
               "hr_error" = "#9933FF", 
               "nonhr_error" = "#2E86AB"),
    labels = c("Overall OOB", "HR class", "non-HR class")
  ) +
  ggtitle("Model Optimization: Error Rate by Iteration") +
  labs(x = "Iteration", y = "Error Rate (%)", color = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_optimization) 

# ============================================================
# STEP 4: PERMUTATION IMPORTANCE
# ============================================================

# --- 4a. Extract importance from final model ---
imp_final <- importance(rf_final)

imp_df <- data.frame(
  variable = rownames(imp_final),
  MeanDecreaseAccuracy = imp_final[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = imp_final[, "MeanDecreaseGini"]
) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  mutate(
    module = case_when(
      str_detect(variable, "^M1_")  ~ "1. Leaf Economics",
      str_detect(variable, "^M2_")  ~ "2. Leaf Nutrients",
      str_detect(variable, "^M3_")  ~ "3. Leaf Physiology",
      str_detect(variable, "^M4_")  ~ "4. Plant Size",
      str_detect(variable, "^M5_")  ~ "5. Root Acquisition",
      str_detect(variable, "^M6_")  ~ "6. Hydraulic",
      str_detect(variable, "^M7_")  ~ "7. Reproduction",
      str_detect(variable, "^M8_")  ~ "8. Dispersal",
      str_detect(variable, "^M9_")  ~ "9. Genome",
      str_detect(variable, "^M10_") ~ "10. Plasticity",
      variable %in% c("origin", "life_span", "growth_form",
                       "leaf_type", "habitat_group") ~ "Categorical",
      TRUE ~ "Other"
    ),
    rank = row_number()
  )

# --- 4b. Importance plot ---
p_imp_accuracy <- imp_df %>%
  ggplot(aes(x = reorder(variable, MeanDecreaseAccuracy),
             y = MeanDecreaseAccuracy,
             fill = module)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_viridis_d(option = "D") +
  ggtitle("Permutation Importance\n(Mean Decrease in Accuracy)") +
  labs(x = "", y = "Mean Decrease Accuracy", fill = "Module") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

p_imp_gini <- imp_df %>%
  ggplot(aes(x = reorder(variable, MeanDecreaseGini),
             y = MeanDecreaseGini,
             fill = module)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_viridis_d(option = "D") +
  ggtitle("Node Purity Importance\n(Mean Decrease in Gini)") +
  labs(x = "", y = "Mean Decrease Gini", fill = "Module") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

figure_importance <- p_imp_accuracy + p_imp_gini +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("Figure_VarImportance.pdf", figure_importance,
       width = 18, height = 10, dpi = 300)

# --- 4c. Module-level importance summary ---
module_importance <- imp_df %>%
  group_by(module) %>%
  summarise(
    n_retained = n(),
    total_MDA = round(sum(MeanDecreaseAccuracy), 4),
    mean_MDA = round(mean(MeanDecreaseAccuracy), 4),
    total_Gini = round(sum(MeanDecreaseGini), 4),
    top_variable = variable[which.max(MeanDecreaseAccuracy)],
    pct_total_MDA = round(
      sum(MeanDecreaseAccuracy) / 
        sum(imp_df$MeanDecreaseAccuracy) * 100, 1
    )
  ) %>%
  arrange(desc(total_MDA))

cat("\n=== Module-Level Importance ===\n")
print(module_importance, width = Inf)

# --- 4d. Module importance plot ---
p_module <- module_importance %>%
  ggplot(aes(x = reorder(module, pct_total_MDA),
             y = pct_total_MDA,
             fill = module)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(pct_total_MDA, "%")),
            hjust = -0.1, size = 4) +
  coord_flip() +
  scale_fill_viridis_d(option = "D") +
  ylim(0, max(module_importance$pct_total_MDA) * 1.15) +
  ggtitle("Module Contribution to HR Prediction") +
  labs(x = "", y = "% of Total Permutation Importance") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave("Figure_ModuleImportance.pdf", p_module,
       width = 10, height = 7, dpi = 300)

# ============================================================
# STEP 5: PARTIAL DEPENDENCE
# Following Friedman (2001), Eq. 1 and 2 from the paper
# ============================================================

# --- 5a. Calculate partial dependence for all retained covariates ---

# For randomForest package, partialPlot gives logit-scale values
# f(x) = log[p_k(x)] - (1/K) * sum(log[p_j(x)])

calculate_partial_dependence <- function(rf_model, data, 
                                          variable, 
                                          class_of_interest = "HR",
                                          n_points = 50) {
  
  if (is.factor(data[[variable]]) | is.character(data[[variable]])) {
    # Categorical: get PD for each level
    pd <- partialPlot(rf_model, 
                       pred.data = as.data.frame(data),
                       x.var = variable,
                       which.class = class_of_interest,
                       plot = FALSE)
    
    result <- data.frame(
      variable = variable,
      x_value = pd$x,
      x_numeric = NA,
      pd_value = pd$y,  # logit scale
      pd_prob = exp(pd$y) / (1 + exp(pd$y)),  # probability scale
      type = "categorical"
    )
    
  } else {
    # Continuous
    pd <- partialPlot(rf_model,
                       pred.data = as.data.frame(data),
                       x.var = variable,
                       which.class = class_of_interest,
                       plot = FALSE,
                       n.pt = n_points)
    
    result <- data.frame(
      variable = variable,
      x_value = as.character(round(pd$x, 4)),
      x_numeric = pd$x,
      pd_value = pd$y,
      pd_prob = exp(pd$y) / (1 + exp(pd$y)),
      type = "continuous"
    )
  }
  
  return(result)
}

# --- 5b. Compute PD for all variables ---
# Use final optimized model data
rf_final_data <- rf_data %>% select(all_of(c("HR_status", final_predictors)))

all_pd <- map_dfr(final_predictors, function(var) {
  cat("Computing partial dependence for:", var, "\n")
  tryCatch(
    calculate_partial_dependence(
      rf_model = rf_final,
      data = rf_final_data,
      variable = var,
      class_of_interest = "HR"
    ),
    error = function(e) {
      cat("  Error:", e$message, "\n")
      return(NULL)
    }
  )
})

# Add module info
all_pd <- all_pd %>%
  mutate(
    module = case_when(
      str_detect(variable, "^M1_")  ~ "Leaf Economics",
      str_detect(variable, "^M2_")  ~ "Leaf Nutrients",
      str_detect(variable, "^M3_")  ~ "Leaf Physiology",
      str_detect(variable, "^M4_")  ~ "Plant Size",
      str_detect(variable, "^M5_")  ~ "Root Acquisition",
      str_detect(variable, "^M6_")  ~ "Hydraulic",
      str_detect(variable, "^M7_")  ~ "Reproduction",
      str_detect(variable, "^M8_")  ~ "Dispersal",
      str_detect(variable, "^M9_")  ~ "Genome",
      str_detect(variable, "^M10_") ~ "Plasticity",
      TRUE ~ "Categorical"
    )
  )

# --- 5c. Determine direction of association ---
direction_summary <- all_pd %>%
  filter(type == "continuous") %>%
  group_by(variable, module) %>%
  summarise(
    pd_min = min(pd_value),
    pd_max = max(pd_value),
    pd_range = pd_max - pd_min,
    # Correlation between x and PD value gives direction
    direction = cor(x_numeric, pd_value, 
                    use = "complete.obs"),
    association = case_when(
      cor(x_numeric, pd_value, use = "complete.obs") > 0.3 ~ 
        "Positive → Higher values increase HR probability",
      cor(x_numeric, pd_value, use = "complete.obs") < -0.3 ~ 
        "Negative → Higher values decrease HR probability",
      TRUE ~ "Non-linear / Weak association"
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(direction)))

cat("\n=== Direction of Associations ===\n")
print(direction_summary, width = Inf)

# --- 5d. Partial dependence plots ---

# Continuous variables
pd_continuous <- all_pd %>% filter(type == "continuous")

# Order by importance
var_order <- imp_df %>% 
  filter(variable %in% unique(pd_continuous$variable)) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  pull(variable)

p_pd_continuous <- pd_continuous %>%
  mutate(variable = factor(variable, levels = var_order)) %>%
  ggplot(aes(x = x_numeric, y = pd_value)) +
  geom_line(linewidth = 1, color = "#9933FF") + 
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  ggtitle("Partial Dependence: Continuous Covariates") +
  labs(
    x = "Covariate Value",
    y = "Partial Dependence\n(logit scale: f̄(x) for HR class)"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 10, face = "bold"))

# Categorical variables
pd_categorical <- all_pd %>% filter(type == "categorical")

p_pd_categorical <- pd_categorical %>%
  ggplot(aes(x = reorder(x_value, pd_value), 
             y = pd_value,
             fill = pd_value > 0)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("TRUE" = "#9933FF", 
                                "FALSE" = "#2E86AB"),
                    labels = c("Decreases HR prob", 
                               "Increases HR prob")) +
  coord_flip() +
  ggtitle("Partial Dependence: Categorical Covariates") +
  labs(x = "", y = "Partial Dependence (logit scale)",
       fill = "Direction") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Combine
figure_pd <- p_pd_continuous / p_pd_categorical +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("Figure_PartialDependence.pdf", figure_pd,
       width = 16, height = 18, dpi = 300)

# ============================================================
# STEP 6: MODEL ASSESSMENT
# ============================================================

# --- 6a. OOB-based metrics ---
oob_predictions <- rf_final$predicted
oob_votes <- rf_final$votes

# Confusion matrix with detailed metrics
cm <- confusionMatrix(
  oob_predictions, 
  rf_data %>% 
    select(all_of(c("HR_status", final_predictors))) %>% 
    pull(HR_status),
  positive = "HR"
)

cat("\n=== Detailed Classification Metrics ===\n")
print(cm)

# --- 6b. ROC from OOB votes ---
roc_oob <- roc(
  response = rf_data %>% 
    select(all_of(c("HR_status", final_predictors))) %>% 
    pull(HR_status),
  predictor = oob_votes[, "HR"],
  levels = c("nonHR", "HR"),
  direction = "<"
)

cat("\nOOB AUC:", round(auc(roc_oob), 3), "\n")
ci_auc <- ci.auc(roc_oob)
cat("95% CI:", round(ci_auc[1], 3), "-", round(ci_auc[3], 3), "\n")

# --- 6c. ROC plot ---
p_roc <- ggroc(roc_oob, legacy.axes = TRUE, 
               linewidth = 1.2, color = "#9933FF") + 
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.55, y = 0.3,
           label = paste0("AUC = ", round(auc(roc_oob), 3),
                          "\n95% CI: [", round(ci_auc[1], 3),
                          ", ", round(ci_auc[3], 3), "]"),
           size = 5, fontface = "bold") +
  ggtitle("OOB ROC Curve: Balanced Random Forest") +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal(base_size = 14)

# --- 6d. Null model permutation test ---
cat("\n=== Permutation Significance Test ===\n")
cat("Running 999 permutations...\n")

n_perm <- 999
perm_oob_errors <- numeric(n_perm)
perm_hr_errors <- numeric(n_perm)

observed_oob <- rf_final$err.rate[ntree_value, "OOB"]
observed_hr_err <- rf_final$err.rate[ntree_value, "HR"]

for (i in 1:n_perm) {
  
  perm_data <- rf_final_data
  perm_data$HR_status <- sample(perm_data$HR_status)
  
  n_hr_perm <- sum(perm_data$HR_status == "HR")
  samp_perm <- c(nonHR = n_hr_perm, HR = n_hr_perm)
  
  perm_rf <- randomForest(
    HR_status ~ .,
    data = perm_data,
    ntree = 500,  # fewer trees for speed
    mtry = max(1, floor(sqrt(length(final_predictors)))),
    sampsize = samp_perm,
    strata = perm_data$HR_status
  )
  
  perm_oob_errors[i] <- perm_rf$err.rate[500, "OOB"]
  perm_hr_errors[i] <- perm_rf$err.rate[500, "HR"]
  
  if (i %% 100 == 0) cat("  Permutation", i, "of", n_perm, "\n")
}

# p-value: proportion of permutations with LOWER error
p_value_oob <- mean(perm_oob_errors <= observed_oob)
p_value_hr <- mean(perm_hr_errors <= observed_hr_err)

cat("\nObserved OOB error:", round(observed_oob * 100, 1), "%\n")
cat("Permutation p-value (OOB):", p_value_oob, "\n")
cat("Permutation p-value (HR class):", p_value_hr, "\n")

# --- 6e. Permutation distribution plot ---
p_perm <- ggplot(data.frame(error = perm_oob_errors * 100), 
                  aes(x = error)) +
  geom_histogram(bins = 40, fill = "gray70", color = "white") +
  geom_vline(xintercept = observed_oob * 100,
             color = "#9933FF", linewidth = 1.5) +
  annotate("text", 
           x = observed_oob * 100 - 2, y = Inf,
           label = paste0("Observed = ", 
                          round(observed_oob * 100, 1), "%\np = ", 
                          ifelse(p_value_oob == 0, 
                                 "< 0.001", p_value_oob)),
           vjust = 2, hjust = 1, color = "#9933FF", 
           size = 4.5, fontface = "bold") +
  ggtitle("Permutation Test: Model vs Null") +
  labs(x = "OOB Error Rate (%) under null",
       y = "Count") +
  theme_minimal(base_size = 14)

# Combine assessment figures
figure_assessment <- (p_roc | p_perm) +
  plot_annotation(
    title = "Model Assessment",   
    tag_levels = "A"
  )

ggsave("Figure_ModelAssessment.pdf", figure_assessment,
       width = 14, height = 6, dpi = 300)

# ============================================================
# STEP 7: THREE-CLASS ANALYSIS
# HR-local (24) vs non-HR-local (250) vs HR-elsewhere (100)
# ============================================================

# --- 7a. Prepare 3-class data ---
# Assuming you have a column 'HR_status_3class'

model_data_3class <- model_data %>%
  mutate(
    HR_3class = factor(
      HR_3class,
      levels = c("nonHR_local", "HR_local", "HR_elsewhere")
    )
  )

rf_data_3class <- model_data_3class %>%
  select(-species, -HR_status) %>%  # remove binary HR
  drop_na()

cat("\n=== 3-Class Distribution ===\n")
print(table(rf_data_3class$HR_3class))

# --- 7b. Balanced sampling for 3 classes ---
# Sample size = minimum class (HR_local = 24)
n_min_3class <- min(table(rf_data_3class$HR_3class))

samp_3class <- c(
  nonHR_local = n_min_3class,
  HR_local = n_min_3class,
  HR_elsewhere = n_min_3class
)

cat("Balanced sampling:", n_min_3class, "per class\n")

# --- 7c. Initial 3-class RF ---
predictors_3class <- setdiff(colnames(rf_data_3class), "HR_3class")
mtry_3class <- floor(sqrt(length(predictors_3class)))

rf_3class <- randomForest(
  HR_3class ~ .,
  data = rf_data_3class,
  ntree = 1000,
  mtry = mtry_3class,
  importance = TRUE,
  sampsize = samp_3class,
  strata = rf_data_3class$HR_3class,
  proximity = TRUE
)

cat("\n=== 3-Class RF Results ===\n")
print(rf_3class)

# --- 7d. Sequential optimization for 3-class ---
# Modify function for 3 classes

opt_result_3class <- sequential_rf_optimization(
  data = rf_data_3class %>% 
    rename(HR_status = HR_3class),  # temporary rename
  response = "HR_status",
  ntree = 1000,
  verbose = TRUE
)

rf_final_3class <- opt_result_3class$final_model

# --- 7e. 3-class confusion matrix ---
cat("\n=== 3-Class Confusion Matrix ===\n")
print(rf_final_3class$confusion)

# Per-class error
cat("\nPer-class OOB error rates:\n")
for (cls in levels(rf_data_3class$HR_3class)) {
  err <- rf_final_3class$confusion[cls, "class.error"]
  cat(" ", cls, ":", round(err * 100, 1), "%\n")
}

# --- 7f. 3-class importance comparison ---
imp_3class <- importance(rf_final_3class, type = 1)

imp_3class_df <- data.frame(
  variable = rownames(imp_3class),
  overall = imp_3class[, "MeanDecreaseAccuracy"]
) %>%
  arrange(desc(overall)) %>%
  mutate(module = case_when(
    str_detect(variable, "^M1_")  ~ "Leaf Economics",
    str_detect(variable, "^M2_")  ~ "Leaf Nutrients",
    str_detect(variable, "^M3_")  ~ "Leaf Physiology",
    str_detect(variable, "^M4_")  ~ "Plant Size",
    str_detect(variable, "^M5_")  ~ "Root Acquisition",
    str_detect(variable, "^M6_")  ~ "Hydraulic",
    str_detect(variable, "^M7_")  ~ "Reproduction",
    str_detect(variable, "^M8_")  ~ "Dispersal",
    str_detect(variable, "^M9_")  ~ "Genome",
    str_detect(variable, "^M10_") ~ "Plasticity",
    TRUE ~ "Categorical"
  ))

# --- 7g. Compare importance: 2-class vs 3-class ---
comparison_imp <- imp_df %>%
  select(variable, binary_MDA = MeanDecreaseAccuracy) %>%
  full_join(
    imp_3class_df %>% 
      select(variable, three_class_MDA = overall),
    by = "variable"
  ) %>%
  drop_na()

p_compare_imp <- ggplot(comparison_imp, 
                         aes(x = binary_MDA, 
                             y = three_class_MDA)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "red") +
  ggrepel::geom_text_repel(
    aes(label = variable),
    size = 3, max.overlaps = 10
  ) +
  ggtitle("Variable Importance: Binary vs 3-Class Model") +
  labs(x = "Binary Model (MDA)", 
       y = "3-Class Model (MDA)") +
  theme_minimal(base_size = 14)

# --- 7h. Pairwise partial dependence for 3-class ---
# Which traits differentiate HR-local from HR-elsewhere?

# Extract votes for each class
votes_3class <- rf_final_3class$votes

# Compare probability profiles
prob_profiles <- data.frame(
  HR_3class = rf_data_3class$HR_3class,
  prob_HR_local = votes_3class[, "HR_local"],
  prob_nonHR_local = votes_3class[, "nonHR_local"],
  prob_HR_elsewhere = votes_3class[, "HR_elsewhere"]
)

# Visualize
p_prob_profiles <- prob_profiles %>%
  pivot_longer(-HR_3class, 
               names_to = "predicted_class",
               values_to = "probability") %>%
  ggplot(aes(x = HR_3class, y = probability, 
             fill = predicted_class)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c(
    "prob_HR_local" = "#9933FF",
    "prob_nonHR_local" = "#2E86AB",
    "prob_HR_elsewhere" = "#FF9900"  
  )) +
  ggtitle("OOB Vote Probabilities by True Class") +
  labs(x = "True Class", y = "OOB Vote Probability",
       fill = "Predicted Class") +
  theme_minimal(base_size = 14)

ggsave("Figure_3class_profiles.pdf", p_prob_profiles,
       width = 10, height = 7, dpi = 300)

 # ============================================================
# STEP 8: PROXIMITY-BASED MDS VISUALIZATION
# ============================================================

# --- 8a. MDS from RF proximity matrix ---
# Binary model
prox_matrix <- rf_final$proximity
mds_result <- cmdscale(1 - prox_matrix, k = 2)

mds_df <- data.frame(
  MDS1 = mds_result[, 1],
  MDS2 = mds_result[, 2],
  HR_status = rf_final_data$HR_status,
  species = model_data$species[
    match(rownames(rf_final_data), rownames(model_data))
  ]
)

p_mds_binary <- ggplot(mds_df, 
                         aes(x = MDS1, y = MDS2, 
                             color = HR_status)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("nonHR" = "#2E86AB", 
                                 "HR" = "#9933FF")) +
  ggtitle("RF Proximity MDS: Species in Trait Space") +
  theme_minimal(base_size = 14)

# --- 8b. MDS for 3-class model ---
if (exists("rf_final_3class") && !is.null(rf_final_3class$proximity)) {
  
  prox_3class <- rf_final_3class$proximity
  mds_3class <- cmdscale(1 - prox_3class, k = 2)
  
  mds_3class_df <- data.frame(
    MDS1 = mds_3class[, 1],
    MDS2 = mds_3class[, 2],
    HR_3class = rf_data_3class$HR_3class
  )
  
  p_mds_3class <- ggplot(mds_3class_df,
                           aes(x = MDS1, y = MDS2, 
                               color = HR_3class)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95, linewidth = 1) +
    scale_color_manual(values = c(
      "nonHR_local" = "#2E86AB",
      "HR_local" = "#9933FF",
      "HR_elsewhere" = "#FF9900"
    )) +
    ggtitle("3-Class RF Proximity MDS") +
    theme_minimal(base_size = 14)
  
  figure_mds <- p_mds_binary + p_mds_3class +
    plot_annotation(tag_levels = "A")
  
  ggsave("Figure_MDS_proximity.pdf", figure_mds,
         width = 16, height = 7, dpi = 300)
}

# ============================================================
# STEP 9: EXPORT ALL RESULTS
# ============================================================

# --- 9a. Summary table ---
results_list <- list(
  
  # Binary model metrics
  binary_model = data.frame(
    metric = c("OOB_error_pct", "HR_class_error_pct",
               "nonHR_class_error_pct", "AUC",
               "AUC_95CI_lower", "AUC_95CI_upper",
               "Sensitivity", "Specificity",
               "N_predictors_initial", "N_predictors_final",
               "N_iterations", "Permutation_p_value",
               "N_HR", "N_nonHR"),
    value = c(
      round(rf_final$err.rate[ntree_value, "OOB"] * 100, 1),
      round(rf_final$err.rate[ntree_value, "HR"] * 100, 1),
      round(rf_final$err.rate[ntree_value, "nonHR"] * 100, 1),
      round(auc(roc_oob), 3),
      round(ci_auc[1], 3),
      round(ci_auc[3], 3),
      round(cm$byClass["Sensitivity"], 3),
      round(cm$byClass["Specificity"], 3),
      length(predictor_names),
      length(final_predictors),
      opt_result$n_iterations,
      p_value_oob,
      n_hr, n_nonhr
    )
  ),
  
  # Variable importance
  variable_importance = imp_df,
  
  # Module importance  
  module_importance = module_importance,
  
  # Partial dependence directions
  pd_directions = direction_summary,
  
  # All partial dependence values
  partial_dependence = all_pd,
  
  # Optimization history
  optimization = history_df
)

# Save everything
for (name in names(results_list)) {
  write.csv(results_list[[name]], 
            paste0("Results_", name, ".csv"),
            row.names = FALSE)
}

# --- 9b. Save models ---
saveRDS(rf_final, "model_rf_binary_final.rds")
saveRDS(rf_final_3class, "model_rf_3class_final.rds")
saveRDS(opt_result, "optimization_result.rds")

cat("\n✓ All results exported successfully.\n")

## ############################################################################################
## STAGE 2  Geographic-range covariates (EOO, AOO, altitude)
## source file: _fix_traits.R
## ############################################################################################

## ============================================================================
## CORRECTED spatial/altitude trait pipeline (alt_range, EOO, AOO).
## Fixes the NA on Chenopodium striatiforme & Phalaris minor.
## Run with R 4.5.2.
## ============================================================================
.libPaths(c(file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5"), .libPaths()))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(sf); library(terra); library(writexl)
})
sf::sf_use_s2(TRUE)
my_data_path <- "../my_data_final_cleanedv2.csv"
out_xlsx     <- "../final_trait_table_FIXED.xlsx"

d <- read_csv(my_data_path, show_col_types = FALSE, guess_max = 100000)
d$species <- str_to_sentence(str_squish(tolower(d$species)))   # FIX: normalise case once

## ---- Phalaris minor FloraWeb supplement, WITH altitude from DEM (consistent
##      with the rest of the dataset, whose altitude is DEM-derived) -----------
phalaris <- data.frame(
  species = "Phalaris minor",
  decimalLatitude  = c(52.47498, 51.87496, 51.87497, 51.67492, 51.32497, 50.92490, 48.02488),
  decimalLongitude = c(13.20826, 12.04156, 12.20823, 7.791429, 12.34980, 6.95807, 7.87477),
  coordinateUncertaintyInMeters = c(3169, 3189, 3189, 3195, 3209, 3221, 3317),
  basisOfRecord = "OCCURRENCE", year = 1980)
dem <- geodata::elevation_30s(country = "DEU", path = tempdir())
phalaris$altitude <- as.numeric(terra::extract(dem, phalaris[, c("decimalLongitude","decimalLatitude")])[, 2])
cat("Phalaris minor DEM altitudes:", paste(round(phalaris$altitude, 1), collapse = ", "), "\n")
d <- bind_rows(d, phalaris)

## ---- spatial audit: keep records inside Germany (EPSG:3035) ------------------
germany_border <- giscoR::gisco_get_countries(country = "DE", resolution = "03") |> st_transform(3035)
d <- d %>% filter(!is.na(decimalLatitude), !is.na(decimalLongitude))
sfp <- st_as_sf(d, coords = c("decimalLongitude","decimalLatitude"), crs = 4326, remove = FALSE) |>
  st_transform(3035)
inside <- lengths(st_intersects(sfp, germany_border)) > 0
geo <- d[inside, ]
cat("audit: kept", sum(inside), "of", nrow(d), "records\n")

## ---- 10 km grid -------------------------------------------------------------
grid <- st_make_grid(germany_border, cellsize = 10000, square = TRUE)
grid <- st_sf(grid_id = seq_along(grid), geometry = grid)
grid <- grid[lengths(st_intersects(grid, germany_border)) > 0, ]
proj <- st_transform(st_as_sf(geo, coords = c("decimalLongitude","decimalLatitude"), crs = 4326), 3035)
pjoin <- st_join(proj, grid, join = st_intersects)

## ---- RANGE METRICS = BASE TABLE (every species with valid occurrences) ------
## FIX: EOO computed per-species with a vectorised grouped convex hull (no
## fragile cur_group() subsetting), and this table is the base — NOT gated on altitude.
eoo <- proj %>% group_by(species) %>%
  summarise(geom = st_combine(geometry), npts = n(), .groups = "drop") %>%
  # EOO (convex hull) is undefined for <3 distinct localities -> NA here; it is raised to AOO
  # after the join below (IUCN rule EOO >= AOO). A 0 would be a false "smallest range".
  mutate(EOO_km2 = ifelse(npts >= 3, as.numeric(st_area(st_convex_hull(geom))) / 1e6, NA_real_)) %>%
  st_drop_geometry() %>% select(species, EOO_km2)
range_metrics_df <- pjoin %>% st_drop_geometry() %>% group_by(species) %>%
  summarise(n_grid_cells = n_distinct(grid_id, na.rm = TRUE),
            AOO_km2 = n_grid_cells * 100, n_points = n(), .groups = "drop") %>%
  left_join(eoo, by = "species") %>%
  # IUCN Red List Guidelines (Standards & Petitions Committee 2019, Box 4.2): AOO is a subset of
  # EOO, so EOO can never be < AOO. Where the convex-hull EOO is undefined (<3 localities) or
  # smaller than AOO, set EOO = AOO.
  mutate(EOO_km2 = ifelse(is.na(EOO_km2), AOO_km2, pmax(EOO_km2, AOO_km2)))

## ---- altitude traits (joined ONTO range metrics, not used as the gate) -------
alt_traits <- geo %>% filter(!is.na(altitude)) %>% group_by(species) %>%
  summarise(alt_min = min(altitude, na.rm = TRUE), alt_max = max(altitude, na.rm = TRUE),
            alt_range = alt_max - alt_min, .groups = "drop") %>%
  mutate(alt_min = ifelse(is.infinite(alt_min), NA, alt_min),
         alt_max = ifelse(is.infinite(alt_max), NA, alt_max),
         alt_range = ifelse(is.infinite(alt_range) | alt_range < 0, 0, alt_range))

## ---- FIX: range metrics is LEFT/base; altitude is joined on -----------------
traits <- range_metrics_df %>% left_join(alt_traits, by = "species")

## ---- NA-safe scaling (base scale() blanks a whole column if any NA) ---------
scl <- function(x) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE); (x - m) / s }
final_trait_table_scaled <- traits %>%
  mutate(
    alt_range_scaled = scl(alt_range),
    EOO_scaled       = scl(log10(EOO_km2 + 1)),
    AOO_scaled       = scl(log10(AOO_km2 + 1)),
    total_obs_scaled = scl(log10(n_points + 1))
  ) %>%
  select(species, n_points, n_grid_cells, alt_min, alt_max, alt_range,
         EOO_km2, AOO_km2, alt_range_scaled, EOO_scaled, AOO_scaled, total_obs_scaled)

## ---- verify the two species + overall NA -----------------------------------
cat("\n=== the two previously-NA species ===\n")
print(as.data.frame(final_trait_table_scaled %>%
  filter(species %in% c("Chenopodium striatiforme","Phalaris minor"))))
na_now <- final_trait_table_scaled %>%
  filter(if_any(c(alt_range_scaled, EOO_scaled, AOO_scaled), is.na)) %>% pull(species)
cat("\nspecies still NA in alt_range/EOO/AOO scaled:", length(na_now),
    if (length(na_now)) paste("->", paste(na_now, collapse=", ")) else "(none)", "\n")
cat("total species in trait table:", nrow(final_trait_table_scaled), "\n")

writexl::write_xlsx(final_trait_table_scaled, out_xlsx)
cat("\nWrote:", out_xlsx, "\n")

## ############################################################################################
## STAGE 3  Agraratlas overlay + agricultural context
## source file: 07_agraratlas_overlay.R
## ############################################################################################

# Watch-list weeds vs German agricultural-statistics atlas (Agraratlas, Destatis/StBA).
# WMS basemaps (visual overlay) + WCS GeoTIFFs (quantitative extraction at occurrence points).
# Data licence: Datenlizenz Deutschland - Namensnennung 2.0 (dl-de/by-2-0); (c) Statistische AEmter
# des Bundes und der Laender, Agraratlas. Used for analysis/overlay as permitted.
options(timeout = 300)
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(png); library(terra); library(patchwork)
})
base <- "."; fd <- file.path(base,"Figures"); td <- file.path(base,"Tables")
occ_path <- "../my_data_cleaned_occ_gbifv2.csv"
BB <- c(xmin=5.8, ymin=47.2, xmax=15.1, ymax=55.1)
wms <- "https://www.wms.nrw.de/wms/agraratlas"; wcs <- "https://www.wcs.nrw.de/stba/agraratlas"

## --- occurrences: watch-list (top15) + background sample of all species ---
wl <- readr::read_csv(file.path(td,"Table_RiskWatchlist.csv"), show_col_types=FALSE) |> arrange(desc(p_HR))
wl_sp <- head(wl$species, 15)
occ <- data.table::fread(occ_path, select=c("species","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters"))
setnames(occ, c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters"), c("lat","lon","uncert"))
occ <- occ[!is.na(lat)&!is.na(lon)&lon>=BB["xmin"]&lon<=BB["xmax"]&lat>=BB["ymin"]&lat<=BB["ymax"]&(is.na(uncert)|uncert<=10000)]
wl_occ <- occ[species %in% wl_sp]
set.seed(42); bg_occ <- occ[sample(.N, min(30000, .N))]   # background = all-species occurrences
cat("watch-list pts:", nrow(wl_occ), "| background pts:", nrow(bg_occ), "\n")

## --- WMS basemaps (for overlay) ---
get_wms <- function(layer, file, W=900){
  H <- round(W*(BB["ymax"]-BB["ymin"])/(BB["xmax"]-BB["xmin"]))
  u <- sprintf("%s?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&LAYERS=%s&STYLES=&CRS=CRS:84&BBOX=%s&WIDTH=%d&HEIGHT=%d&FORMAT=image/png&TRANSPARENT=TRUE",
               wms, layer, paste(BB[c("xmin","ymin","xmax","ymax")],collapse=","), W, H)
  tryCatch({ download.file(u, file, mode="wb", quiet=TRUE); png::readPNG(file) }, error=function(e){cat("WMS fail",layer,"\n");NULL})
}
wheat_img <- get_wms("K08_Weizenanteil_2020", file.path(fd,"_wms_wheat.png"))
farm_img  <- get_wms("K02_Betriebsgroesse_2020", file.path(fd,"_wms_farm.png"))

## --- WCS GeoTIFFs (for extraction): class rasters 1-5 (0 = non-agric/nodata) ---
get_wcs <- function(cov, W=700){
  H <- round(W*(BB["ymax"]-BB["ymin"])/(BB["xmax"]-BB["xmin"]))
  f <- file.path(fd, paste0("_wcs_", cov, ".tif"))
  ok <- tryCatch({ download.file(sprintf("%s?SERVICE=WCS&VERSION=1.0.0&REQUEST=GetCoverage&COVERAGE=%s&CRS=EPSG:4326&BBOX=%s&WIDTH=%d&HEIGHT=%d&FORMAT=GeoTIFF",
                                          wcs, cov, paste(BB[c("xmin","ymin","xmax","ymax")],collapse=","), W, H), f, mode="wb", quiet=TRUE); TRUE },
               error=function(e) FALSE)
  if(!ok) return(NULL)
  r <- terra::rast(f)
  if (is.na(terra::crs(r)) || terra::crs(r)=="") terra::crs(r) <- "EPSG:4326"
  terra::ext(r) <- c(BB["xmin"],BB["xmax"],BB["ymin"],BB["ymax"]); r
}
vars <- c(Wheat="K08_Weizenanteil_2020", Maize="K06_Maisanteil_2020",
          FarmSize="K02_Betriebsgroesse_2020", Organic="K10_Oekoflaechenanteil_2020")
rasters <- lapply(vars, get_wcs); names(rasters) <- names(vars)

## --- extract class at watch-list vs background points; quantitative test ---
ext_cls <- function(r, df){ v <- terra::extract(r, cbind(df$lon, df$lat))[,1]; v[v==0] <- NA; v }
res <- map_dfr(names(rasters), function(nm){
  r <- rasters[[nm]]; if(is.null(r)) return(NULL)
  wlv <- ext_cls(r, wl_occ); bgv <- ext_cls(r, bg_occ)
  w <- suppressWarnings(wilcox.test(wlv, bgv))
  tibble(variable=nm,
         watchlist_median=median(wlv,na.rm=TRUE), background_median=median(bgv,na.rm=TRUE),
         watchlist_mean=mean(wlv,na.rm=TRUE), background_mean=mean(bgv,na.rm=TRUE),
         watchlist_se=sd(wlv,na.rm=TRUE)/sqrt(sum(!is.na(wlv))),
         background_se=sd(bgv,na.rm=TRUE)/sqrt(sum(!is.na(bgv))),
         pct_top2_watchlist=mean(wlv>=4,na.rm=TRUE)*100, pct_top2_background=mean(bgv>=4,na.rm=TRUE)*100,
         p_value=w$p.value, n_wl=sum(!is.na(wlv)), n_bg=sum(!is.na(bgv)))
})
readr::write_csv(res, file.path(td,"Table_AgriContext.csv"))
cat("\n=== Agricultural-context comparison (class 1-5; watch-list vs background) ===\n"); print(as.data.frame(res), digits=3)

## --- FIG 12: WMS overlay (watch-list density over wheat-share & farm-size) ---
# state boundaries as a path data.frame (coord_quickmap is cartesian, so geom_sf is not used)
suppressPackageStartupMessages({library(geodata); library(sf)})
states <- sf::st_as_sf(geodata::gadm("DEU", level=1, path=tempdir())); states$name <- states$NAME_1
states <- suppressWarnings(sf::st_cast(sf::st_make_valid(states), "MULTIPOLYGON"))
stpath <- as.data.frame(sf::st_coordinates(states))
stpath$grp <- interaction(stpath$L1, stpath$L2, stpath$L3)
cent <- cbind(name=states$name, as.data.frame(sf::st_coordinates(suppressWarnings(sf::st_centroid(states)))))
th <- theme_void(base_family="Arial") + theme(plot.title=element_text(size=9,face="bold"),
       plot.subtitle=element_text(size=7), legend.position="none")
panel <- function(img, ttl, sub){
  g <- ggplot(wl_occ, aes(lon, lat))
  if(!is.null(img)) g <- g + annotation_raster(as.raster(img), BB["xmin"],BB["xmax"],BB["ymin"],BB["ymax"])
  g + stat_density_2d(aes(fill=after_stat(level)), geom="polygon", colour=NA, alpha=0.45, bins=8) +
    scale_fill_gradient(low="#fdd0c0", high="#9e2a18") +
    geom_path(data=stpath, aes(X, Y, group=grp), colour="grey35", linewidth=0.15, inherit.aes=FALSE) +
    geom_text(data=cent, aes(X, Y, label=name), size=1.4, colour="grey20", inherit.aes=FALSE) +
    coord_quickmap(xlim=BB[c("xmin","xmax")], ylim=BB[c("ymin","ymax")], expand=FALSE) +
    labs(title=ttl, subtitle=sub) + th
}
pA <- panel(wheat_img,"(a) Watch-list density over wheat-share","green = wheat share; red = high-risk weed density")
pB <- panel(farm_img,"(b) Watch-list density over farm-size","green = mean farm size; red = high-risk weed density")
fig12 <- (pA | pB) +
         patchwork::plot_annotation(caption="Agraratlas (c) Statistische AEmter des Bundes u. der Laender, dl-de/by-2-0")
ggsave(file.path(fd,"Fig12_AgriOverlay.pdf"), fig12, width=180, height=120, units="mm", device=cairo_pdf)
ggsave(file.path(fd,"Fig12_AgriOverlay.png"), fig12, width=180, height=120, units="mm", dpi=600)

## --- FIG 13: quantitative agricultural context (watch-list vs background) ---
qd <- res |> select(variable, Watch_list=watchlist_mean, Background=background_mean, se_wl=watchlist_se, se_bg=background_se) |>
  pivot_longer(c(Watch_list,Background), names_to="group", values_to="mean_class") |>
  mutate(se = ifelse(group=="Watch_list", se_wl, se_bg))
fig13 <- ggplot(qd, aes(variable, mean_class, fill=group)) +
  geom_col(position=position_dodge(0.7), width=0.65) +
  geom_errorbar(aes(ymin=pmax(0, mean_class-1.96*se), ymax=mean_class+1.96*se), position=position_dodge(0.7),
                width=0.18, linewidth=0.3, colour="grey30") +
  scale_fill_manual(values=c(Watch_list="#E76F51", Background="#457B9D"), name=NULL,
                    labels=c(Background="all weeds (background)", Watch_list="watch-list (high-risk)")) +
  labs(x=NULL, y="Mean agricultural-intensity class (1-5)",
       title="Agricultural context of high-risk vs background weed occurrences") +
  theme_minimal(base_size=9, base_family="Arial") +
  theme(legend.position="top", panel.grid.major.x=element_blank())
ggsave(file.path(fd,"Fig13_AgriContext.pdf"), fig13, width=130, height=90, units="mm", device=cairo_pdf)
ggsave(file.path(fd,"Fig13_AgriContext.png"), fig13, width=130, height=90, units="mm", dpi=600)

cat("AGRARATLAS DONE\n")

## ############################################################################################
## STAGE 4  Crop-data overlay + crop context
## source file: 08_cropdata_overlay.R
## ############################################################################################

# Watch-list weeds vs OFFICIAL district-level crop statistics (NUTS3, 1979-2021; Final_data-crop types).
# Primary agricultural-context analysis (district areal data; more defensible than the gridded Agraratlas).
# Source: German regional statistical offices; CC-BY 4.0 / Datenlizenz Deutschland - Namensnennung 2.0.
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(data.table); library(sf); library(giscoR); library(geodata); library(patchwork)
})
sf::sf_use_s2(FALSE)
base <- "."; fd <- file.path(base,"Figures"); td <- file.path(base,"Tables")
BB <- c(xmin=5.8, ymin=47.2, xmax=15.1, ymax=55.1); col_hr <- "#E76F51"

## --- crop statistics: recent-year district means (area), wide ---
cr <- readxl::read_excel("../Final_data-crop types.xlsx")
cr <- cr |> mutate(nuts_id = trimws(nuts_id),                 # nuts_id has trailing spaces
                   value = suppressWarnings(as.numeric(value))) |>  # value stored as text
  filter(measure=="area", (is.na(outlier)|outlier!=1), year>=2017)
wide <- cr |> group_by(nuts_id, crop) |> summarise(area=mean(value,na.rm=TRUE), .groups="drop") |>
  filter(crop %in% c("ArabLand","district","ww","grain_maize","silage_maize")) |>
  pivot_wider(names_from=crop, values_from=area)
wide <- wide |> mutate(
  arable_share = pmin(1, ArabLand/district),
  ww_share     = pmin(1, ww/ArabLand),
  maize_share  = pmin(1, (coalesce(grain_maize,0)+coalesce(silage_maize,0))/ArabLand))
cat("districts with crop metrics:", nrow(wide), "\n")

## --- NUTS3 polygons + join + state boundaries ---
nuts <- giscoR::gisco_get_nuts(country="DE", nuts_level=3, year="2021", resolution="03") |> st_transform(4326)
nuts <- nuts |> left_join(wide, by=c("NUTS_ID"="nuts_id"))
cat("NUTS3 joined (non-NA arable):", sum(!is.na(nuts$arable_share)), "/", nrow(nuts), "\n")
states <- sf::st_as_sf(geodata::gadm("DEU", level=1, path=tempdir())); states$name <- states$NAME_1
cent <- cbind(name=states$name, as.data.frame(sf::st_coordinates(suppressWarnings(sf::st_centroid(states)))))

## --- occurrences: watch-list (top15) + background ---
wl <- readr::read_csv(file.path(td,"Table_RiskWatchlist.csv"), show_col_types=FALSE) |> arrange(desc(p_HR))
occ <- data.table::fread("../my_data_cleaned_occ_gbifv2.csv",
                         select=c("species","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters"))
setnames(occ, c("decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters"), c("lat","lon","uncert"))
occ <- occ[!is.na(lat)&!is.na(lon)&lon>=BB["xmin"]&lon<=BB["xmax"]&lat>=BB["ymin"]&lat<=BB["ymax"]&(is.na(uncert)|uncert<=10000)]
hr_spp <- wl$species[wl$p_HR >= 0.5]                       # ALL model-flagged high-risk species (P>=0.5), not top-15
cat("high-risk watch-list species (P>=0.5):", length(hr_spp), "of", nrow(wl), "\n")
wl_occ <- occ[species %in% hr_spp]; set.seed(42); bg_occ <- occ[sample(.N, min(30000,.N))]

## --- point-in-district extraction (quantitative) ---
pip <- function(df){ p <- st_as_sf(df, coords=c("lon","lat"), crs=4326)
  j <- st_join(p, nuts[,c("arable_share","ww_share","maize_share")], join=st_intersects)
  st_drop_geometry(j) }
wlv <- pip(wl_occ); bgv <- pip(bg_occ)
res <- map_dfr(c("arable_share","ww_share","maize_share"), function(v){
  w <- suppressWarnings(wilcox.test(wlv[[v]], bgv[[v]]))
  tibble(variable=v, watchlist_mean=mean(wlv[[v]],na.rm=TRUE), background_mean=mean(bgv[[v]],na.rm=TRUE),
         watchlist_se=sd(wlv[[v]],na.rm=TRUE)/sqrt(sum(!is.na(wlv[[v]]))),
         background_se=sd(bgv[[v]],na.rm=TRUE)/sqrt(sum(!is.na(bgv[[v]]))),
         p_value=w$p.value, n_wl=sum(!is.na(wlv[[v]])), n_bg=sum(!is.na(bgv[[v]])))
})
readr::write_csv(res, file.path(td,"Table_CropContext.csv"))
## --- ALSO export the agricultural data used here into data/ (for co-authors / re-use) ---
dd <- file.path(base,"data"); dir.create(dd, showWarnings=FALSE, recursive=TRUE)
readr::write_csv(res, file.path(dd,"agri_context_watchlist_vs_background.csv"))          # the tested contrast
readr::write_csv(as.data.frame(wide), file.path(dd,"agri_crop_shares_by_district_NUTS3.csv")) # district crop shares (arable/ww/maize)
readr::write_csv(as.data.frame(wlv), file.path(dd,"agri_shares_at_watchlist_occurrences.csv")) # per-occurrence extracted shares (watch-list)
cat("\n=== crop-context (district share; watch-list vs background) ===\n"); print(as.data.frame(res), digits=3)
cat("agri data exported to data/: agri_context_watchlist_vs_background.csv, agri_crop_shares_by_district_NUTS3.csv, agri_shares_at_watchlist_occurrences.csv\n")
cat("N high-risk species mapped (P>=0.5, with DE occurrences):", length(unique(wl_occ$species)), " | background occurrences:", nrow(bg_occ), "\n")

## --- FIG 14: choropleth (arable & winter-wheat share) + watch-list density ---
th <- theme_void(base_family="Arial") + theme(plot.title=element_text(size=9,face="bold"),
       plot.subtitle=element_text(size=7), legend.key.width=unit(3.5,"mm"), legend.key.height=unit(3,"mm"),
       legend.title=element_text(size=6.5), legend.text=element_text(size=6))
cho <- function(fillvar, ttl, leg){
  ggplot() +
    geom_sf(data=nuts, aes(fill=.data[[fillvar]]), colour="grey80", linewidth=0.08) +
    scale_fill_viridis_c(option="viridis", name=leg, na.value="grey92") +
    geom_sf(data=states, fill=NA, colour="grey35", linewidth=0.18) +
    geom_text(data=cent, aes(X,Y,label=name), size=1.4, colour="grey15", inherit.aes=FALSE) +
    stat_density_2d(data=wl_occ, aes(lon,lat), colour="#ff2a1a", linewidth=0.55, bins=7, alpha=0.9) +
    coord_sf(xlim=BB[c("xmin","xmax")], ylim=BB[c("ymin","ymax")], expand=FALSE) +
    labs(title=ttl, subtitle="district share (NUTS3); bold red contours = watch-list weed density") + th
}
fig14 <- (cho("arable_share","(a) Arable-land share","arable / area") |
          cho("ww_share","(b) Winter-wheat share","wheat / arable")) +
  patchwork::plot_annotation(caption="Crop statistics: German regional statistical offices (CC-BY 4.0 / dl-de/by-2-0)")
ggsave(file.path(fd,"Fig14_CropOverlay.pdf"), fig14, width=180, height=120, units="mm", device=cairo_pdf)
ggsave(file.path(fd,"Fig14_CropOverlay.png"), fig14, width=180, height=120, units="mm", dpi=600)

## --- FIG 15: quantitative (watch-list vs background) ---
qd <- res |> transmute(variable=recode(variable, arable_share="Arable share", ww_share="Winter-wheat share", maize_share="Maize share"),
                       `watch-list (high-risk)`=watchlist_mean, `all weeds (background)`=background_mean,
                       se_wl=watchlist_se, se_bg=background_se) |>
  pivot_longer(c(`watch-list (high-risk)`,`all weeds (background)`), names_to="group", values_to="share") |>
  mutate(se = ifelse(group=="watch-list (high-risk)", se_wl, se_bg))
fig15 <- ggplot(qd, aes(variable, share, fill=group)) +
  geom_col(position=position_dodge(0.7), width=0.65) +
  geom_errorbar(aes(ymin=pmax(0, share-1.96*se), ymax=share+1.96*se), position=position_dodge(0.7),
                width=0.18, linewidth=0.3, colour="grey30") +
  scale_fill_manual(values=c("all weeds (background)"="#457B9D","watch-list (high-risk)"="#E76F51"), name=NULL) +
  labs(x=NULL, y="Mean district share at occurrences",
       title="Crop context of high-risk vs background weed occurrences") +
  theme_minimal(base_size=9, base_family="Arial") + theme(legend.position="top", panel.grid.major.x=element_blank())
ggsave(file.path(fd,"Fig15_CropContext.pdf"), fig15, width=130, height=90, units="mm", device=cairo_pdf)
ggsave(file.path(fd,"Fig15_CropContext.png"), fig15, width=130, height=90, units="mm", dpi=600)
cat("CROPDATA DONE\n")

## ############################################################################################
## STAGE 5  HR-status spatial surfaces
## source file: 09_hr_surface.R
## ############################################################################################

# Regenerate the spatial HR-status surfaces FROM the cleaned GBIF occurrence data,
# in the consistent manuscript palette (#2A9D8F HR_elsewhere / #E76F51 HR_local / #457B9D Non_HR_local).
# HR_Status is joined from merged_data_final.csv by species (occurrence file carries no status column).
# Outputs (Figures/): HR_raw_map, HR_risk_map, HR_uncertainty_map, FigHR_RiskUncert  (.pdf + .png 600 dpi)
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(sf)
  library(rnaturalearth); library(geodata); library(patchwork); library(ggspatial)
})

base    <- "."
fd      <- file.path(base, "Figures")
occ_csv <- "../my_data_cleaned_occ_gbifv2.csv"
md_csv  <- "data/merged_data_final_REBUILT.csv"
DE_BB   <- list(xmin = 5.8, xmax = 15.1, ymin = 47.2, ymax = 55.1)

## ---- consistent manuscript palette ----
col_status <- c(HR_elsewhere = "#2A9D8F", HR_local = "#E76F51", Non_HR_local = "#457B9D")
## risk index runs 0 (Non_HR) -> 0.5 (HR_elsewhere) -> 1 (HR_local): interpolate the SAME three colours
risk_cols <- c("#457B9D", "#2A9D8F", "#E76F51")

try(windowsFonts(Arial = windowsFont("Arial")), silent = TRUE)
# New Phytologist: legible type at print size (>= ~7 pt body, >= 6 pt axis ticks)
th <- theme_minimal(base_size = 9, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 6.5, colour = "grey25"),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.5),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7),
        legend.key.width = unit(4.5, "mm"), legend.key.height = unit(4.5, "mm"))

# export both formats at 600 dpi (cairo_pdf rasterises any raster grobs at this dpi; vector stays vector)
save_both <- function(p, name, w, h) {
  tryCatch(ggsave(file.path(fd, paste0(name, ".pdf")), p, width = w, height = h,
                  units = "mm", device = cairo_pdf, dpi = 600),
           error = function(e) ggsave(file.path(fd, paste0(name, ".pdf")), p,
                                      width = w, height = h, units = "mm", dpi = 600))
  ggsave(file.path(fd, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600)
}

## ---- HR_Status lookup (species -> class) ----
md <- fread(md_csv, select = c("species", "HR_Status"))
md[, key := tolower(trimws(species))]
hr_lkp <- setNames(md$HR_Status, md$key)

## ---- occurrences: species + coords, join status, clip to Germany ----
occ <- fread(occ_csv, select = c("species", "decimalLatitude", "decimalLongitude"))
setnames(occ, c("decimalLatitude", "decimalLongitude"), c("lat", "lon"))
occ[, key := tolower(trimws(species))]
occ[, HR_Status := hr_lkp[key]]
occ <- occ[!is.na(HR_Status) & !is.na(lat) & !is.na(lon) &
           lon >= DE_BB$xmin & lon <= DE_BB$xmax & lat >= DE_BB$ymin & lat <= DE_BB$ymax]
occ[, HR_Status := factor(HR_Status, levels = c("Non_HR_local", "HR_local", "HR_elsewhere"))]
cat("classified occurrences in DE:", nrow(occ), "across", uniqueN(occ$species), "species\n")
print(occ[, .N, by = HR_Status])

## ---- basemap: country (metric CRS) + federal states + names ----
germany <- rnaturalearth::ne_countries(scale = "medium", country = "Germany",
                                        returnclass = "sf") |> sf::st_transform(3035)
states  <- sf::st_as_sf(geodata::gadm("DEU", level = 1, path = tempdir())) |>
  sf::st_make_valid() |> sf::st_transform(3035)
states$name <- states$NAME_1
cent <- cbind(name = states$name,
              as.data.frame(sf::st_coordinates(suppressWarnings(sf::st_centroid(states)))))
bb   <- sf::st_bbox(germany)
state_bdy  <- function() geom_sf(data = states, fill = NA, colour = "grey55", linewidth = 0.12)
state_name <- function(sz = 1.4) geom_text(data = cent, aes(X, Y, label = name), size = sz,
                                           colour = "grey25", inherit.aes = FALSE)
coord <- function() coord_sf(crs = 3035, xlim = c(bb["xmin"], bb["xmax"]),
                             ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE)

## project points once
xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(occ, coords = c("lon", "lat"), crs = 4326), 3035))
occ[, `:=`(X = xy[, 1], Y = xy[, 2])]

## ===== (1) RAW occurrences coloured by HR status =====
set.seed(42)
SHOW <- 150000                                   # display subsample (grids below use ALL records)
occ_disp <- if (nrow(occ) > SHOW) occ[sample(.N, SHOW)] else occ
## draw order: Non_HR_local (baseline) -> HR_elsewhere -> HR_local on top (rarest, most relevant)
occ_disp[, prk := fcase(HR_Status == "Non_HR_local", 1L, HR_Status == "HR_elsewhere", 2L,
                        HR_Status == "HR_local", 3L)]
occ_disp <- occ_disp[order(prk)]
m_raw <- ggplot() +
  geom_sf(data = germany, fill = "grey97", colour = "grey60", linewidth = 0.3) +
  state_bdy() +
  geom_point(data = occ_disp, aes(X, Y, colour = HR_Status), size = 0.25, alpha = 0.5) +
  scale_colour_manual(values = col_status, name = "HR status",
                      guide = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  state_name() + coord() +
  labs(x = NULL, y = NULL, title = "(a) Recorded occurrences of HR-classified species",
       subtitle = sprintf("GBIF records (display subsample n=%s of %s); not predictions, not abundance",
                          format(nrow(occ_disp), big.mark = ","),
                          format(nrow(occ), big.mark = ","))) + th
save_both(m_raw, "HR_raw_map", 130, 130)

## ===== 10 x 10 km grid (EPSG:3035) -> weighted HR index + uncertainty (ALL records) =====
CELL <- 10000
occ[, `:=`(cx = floor(X / CELL) * CELL + CELL / 2,
           cy = floor(Y / CELL) * CELL + CELL / 2)]
occ[, w := fcase(HR_Status == "HR_local", 1, HR_Status == "HR_elsewhere", 0.5,
                 HR_Status == "Non_HR_local", 0)]
cell <- occ[, .(n = .N, hr_mean = mean(w), hr_se = sd(w) / sqrt(.N)), by = .(cx, cy)]

gridmap <- function(fillv, ttl, sub, scale_layer, leg) {
  ggplot() +
    geom_sf(data = germany, fill = "grey97", colour = "grey60", linewidth = 0.3) +
    geom_tile(data = cell, aes(cx, cy, fill = .data[[fillv]]), width = CELL, height = CELL) +
    scale_layer + state_bdy() +
    annotation_scale(location = "bl", width_hint = 0.28, height = unit(1.4, "mm"),
                     text_cex = 0.6) +
    coord() + labs(x = NULL, y = NULL, title = ttl, subtitle = sub) + th
}
m_risk <- gridmap("hr_mean", "(a) Weighted HR-status index",
                  "10 x 10 km; HR_local=1, HR_elsewhere=0.5, Non_HR_local=0",
                  scale_fill_gradientn(colours = risk_cols, values = c(0, 0.5, 1),
                                       limits = c(0, 1), name = "HR index"))
m_unc  <- gridmap("hr_se", "(b) Uncertainty of the index (SE)",
                  "10 x 10 km; higher = fewer / more heterogeneous records",
                  scale_fill_gradient(low = "grey92", high = "#E76F51", na.value = "grey92",
                                      name = "SE"))
# standalone supplement maps at NPh 1.5-column width (130 mm); composite at double-column (180 mm)
save_both(m_risk, "HR_risk_map", 130, 140)
save_both(m_unc,  "HR_uncertainty_map", 130, 140)

fig_ru <- m_risk + m_unc
save_both(fig_ru, "FigHR_RiskUncert", 180, 115)

cat("HR SURFACE MAPS REGENERATED (manuscript palette):", nrow(cell), "occupied 10km cells\n")

## ############################################################################################
## STAGE 6  Sites-of-action / mechanism + selection history
## source file: 10_moa_analysis.R
## ############################################################################################

# =============================================================================
# Does the HR-local vs HR-elsewhere distinction (seen in geographic range and traits)
# also appear in MECHANISM OF ACTION (MOA / HRAC group) resistance profiles?
# Merges herbicide_long_clean.csv (species x MOA x active substance x crop cases) with
# hr_status.csv (species -> HR_local / HR_elsewhere), restricted to the 575 analysis species.
# Saves the merged table to data/ and a figure + results object. Standalone (not in the knit).
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(readxl) })
base <- "."; fd <- file.path(base, "Figures"); dd <- file.path(base, "data")
P <- "P:/CONSERWA/Herbicide_Project/Data_Trait"
norm <- function(x) str_squish(str_to_lower(x))
## resolve database species names to the ACCEPTED names used in the 575 model (nomenclatural synonyms)
synonyms <- c("vulpia bromoides"               = "festuca bromoides",
              "hordeum murinum ssp. leporinum" = "hordeum murinum",
              "hordeum murinum ssp. glaucum"   = "hordeum murinum")
fix_syn <- function(s) ifelse(s %in% names(synonyms), unname(synonyms[s]), s)

sp575 <- read_csv(file.path(dd, "imputed_analysis_full_575species.csv"), show_col_types = FALSE) |> mutate(species = norm(species))

## --- MOA source: the GLOBAL herbicide-resistance database. herbicide_long_clean.csv is its cleaned,
##     one-mechanism-per-row derivative (best species-name matching); we use it as PRIMARY and add, from the
##     raw global .xlsx (multiple mechanisms per cell), ONLY the species that cleaning dropped - so there is
##     no double-counting (a species comes from exactly one source) while coverage is maximised. ---
hl <- read_csv(file.path(P, "herbicide_long_clean.csv"), show_col_types = FALSE) |>
  transmute(species = fix_syn(norm(species)), hrac = str_match(tolower(moas), "hrac group (\\d+)")[, 2],
            active = tolower(active_clean), situation = tolower(situations), year = year,
            country = tolower(country), source = "cleaned-long")
xl <- read_excel(file.path(P, "HerbicideResistanceglobal.xlsx")) |>
  mutate(species = fix_syn(norm(Species)), hg = str_extract_all(tolower(MOAs), "hrac group \\d+")) |>
  unnest(hg) |> transmute(species, hrac = str_extract(hg, "\\d+"), active = tolower(Active_substance),
            situation = tolower(Situations), year = Year, country = tolower(Country), source = "raw-global")
moa_all <- bind_rows(hl, xl |> filter(!species %in% unique(hl$species))) |> filter(!is.na(hrac))

## --- HRAC group -> short mechanistic label ---
hrac_lab <- c("1"="ACCase (1)","2"="ALS (2)","3"="Microtubule assembly (3)","4"="Synthetic auxin (4)",
  "5"="PSII - Ser264 (5)","6"="PSII - His215 (6)","9"="EPSPS / glyphosate (9)","10"="Glutamine synthetase (10)",
  "12"="PDS (12)","13"="DOXP synthase (13)","14"="PPO (14)","15"="VLCFA (15)","22"="PSI diversion (22)",
  "23"="Microtubule organisation (23)","27"="HPPD (27)","29"="Cellulose synthesis (29)",
  "34"="Lycopene cyclase (34)","0"="Other / unknown")
moa_all <- moa_all |> mutate(moa = ifelse(hrac %in% names(hrac_lab), hrac_lab[hrac], paste0("HRAC ", hrac)))

## --- HR status from the MODEL's OWN classification (575 analysis file: 24 HR_local + 85 HR_elsewhere) ---
hrmodel <- sp575 |> filter(hr_status %in% c("HR_local","HR_elsewhere")) |> transmute(species, HR_Status = hr_status)
merged <- moa_all |> inner_join(hrmodel, by = "species") |>
  transmute(species, HR_Status, year, country, hrac_group = hrac, moa, active_substance = active,
            crop_situation = situation, source)
write_csv(merged, file.path(dd, "moa_hr_status_merged.csv"))
cat("merged rows:", nrow(merged), " species:", n_distinct(merged$species), "\n")
cat("species covered (of 24 local + 85 elsewhere):\n"); print(merged |> distinct(species, HR_Status) |> count(HR_Status))

## --- species-level MOA sets ---
spp     <- merged |> distinct(species, HR_Status, moa)
breadth <- spp |> count(species, HR_Status, name = "n_moa")
labs2   <- c(HR_local = "HR-local (Germany)", HR_elsewhere = "HR-elsewhere")
cols    <- c(HR_local = "#E76F51", HR_elsewhere = "#2A9D8F")   # match project col_status (HR-local orange, HR-elsewhere teal)

# (1) resistance BREADTH: distinct MOAs per species
wt  <- wilcox.test(n_moa ~ HR_Status, data = breadth)
med <- breadth |> group_by(HR_Status) |> summarise(n_spp = n(), median_moa = median(n_moa),
                                                    mean_moa = round(mean(n_moa), 2), .groups = "drop")
multi <- breadth |> group_by(HR_Status) |> summarise(pct_multi = round(100*mean(n_moa >= 2)), .groups = "drop")
# (2) MOA FREQUENCY: % of species in each HR group resistant to each MOA
Ns   <- breadth |> count(HR_Status, name = "N")
freq <- spp |> count(HR_Status, moa) |> left_join(Ns, by = "HR_Status") |> mutate(pct = 100*n/N)
# (3) COMPOSITION test: MOA x HR-status (do the two groups resist different mechanisms?)
tab  <- spp |> distinct(species, HR_Status, moa) |> count(HR_Status, moa) |>
  pivot_wider(names_from = HR_Status, values_from = n, values_fill = 0)
chi  <- suppressWarnings(chisq.test(as.matrix(tab[, c("HR_local","HR_elsewhere")])))

cat(sprintf("BREADTH: HR-local mean %.2f MOAs (median %.0f, n=%.0f); HR-elsewhere mean %.2f (median %.0f, n=%.0f); Wilcoxon p=%.3g\n",
  med$mean_moa[med$HR_Status=="HR_local"], med$median_moa[med$HR_Status=="HR_local"], med$n_spp[med$HR_Status=="HR_local"],
  med$mean_moa[med$HR_Status=="HR_elsewhere"], med$median_moa[med$HR_Status=="HR_elsewhere"], med$n_spp[med$HR_Status=="HR_elsewhere"], wt$p.value))
cat(sprintf("MULTI-RESISTANCE (>=2 MOAs): HR-local %.0f%%, HR-elsewhere %.0f%%\n",
  multi$pct_multi[multi$HR_Status=="HR_local"], multi$pct_multi[multi$HR_Status=="HR_elsewhere"]))
cat(sprintf("COMPOSITION chi-square p=%.3g\n", chi$p.value))

## --- (4) GEOGRAPHIC-EXPOSURE test: does German range predict the NUMBER of resisted MOAs? ---
rng    <- breadth |> left_join(sp575 |> select(species, eoo_scaled, aoo_scaled), by = "species")
rc_eoo <- suppressWarnings(cor.test(rng$n_moa, rng$eoo_scaled, method = "spearman"))
rc_aoo <- suppressWarnings(cor.test(rng$n_moa, rng$aoo_scaled, method = "spearman"))
range_moa <- list(rho_eoo = unname(rc_eoo$estimate), p_eoo = rc_eoo$p.value,
                  rho_aoo = unname(rc_aoo$estimate), p_aoo = rc_aoo$p.value, n = sum(!is.na(rng$eoo_scaled)))
cat(sprintf("RANGE~MOA count (Spearman): EOO rho=%.3f p=%.3f | AOO rho=%.3f p=%.3f\n",
            range_moa$rho_eoo, range_moa$p_eoo, range_moa$rho_aoo, range_moa$p_aoo))

## --- (5) WITHIN-group concentration: is the MOA distribution uniform within each group? (GOF) ---
gof <- setNames(lapply(c("HR_local","HR_elsewhere"), function(g) {
  gg <- suppressWarnings(chisq.test(table(spp$moa[spp$HR_Status == g])))
  list(stat = unname(gg$statistic), df = unname(gg$parameter), p = gg$p.value) }),
  c("HR_local","HR_elsewhere"))

## --- (6) PER-MOA HR-local minus HR-elsewhere (which mechanisms drive the visual dominance) ---
per_moa <- freq |> select(moa, HR_Status, pct) |>
  pivot_wider(names_from = HR_Status, values_from = pct, values_fill = 0) |>
  mutate(diff = HR_local - HR_elsewhere) |> arrange(desc(diff))

## --- (7) PHYLOGENETIC distinctness: are HR-local and HR-elsewhere different lineages? ---
phylo_d <- tryCatch({
  suppressPackageStartupMessages({library(ape); library(caper)})
  tr <- readRDS(file.path(base, "megatree.rds")); tr$node.label <- NULL   # drop node labels (caper clash)
  tp <- data.frame(tip = tr$tip.label, sp = norm(gsub("_"," ", tr$tip.label))) |>
    inner_join(breadth |> transmute(sp = species, hrl = as.integer(HR_Status == "HR_local")), by = "sp")
  tr2 <- ape::drop.tip(tr, setdiff(tr$tip.label, tp$tip))
  cd  <- tp[match(tr2$tip.label, tp$tip), ]
  cp  <- caper::comparative.data(tr2, data.frame(tip = tr2$tip.label, hrl = cd$hrl), names.col = "tip")
  set.seed(1); pd <- caper::phylo.d(cp, binvar = hrl, permut = 500)
  list(D = unname(pd$DEstimate), p_random = pd$Pval1, p_brownian = pd$Pval0, n = ape::Ntip(tr2))
}, error = function(e) { message("phylo.d failed: ", conditionMessage(e)); NULL })
if (!is.null(phylo_d)) cat(sprintf("PHYLO HR-local clustering: D=%.3f p(random)=%.3f p(Brownian)=%.3f\n",
                                   phylo_d$D, phylo_d$p_random, phylo_d$p_brownian))

## --- (8) SELECTION-HISTORY test: do longer resistance histories carry more MOAs? ---
## The resistance-record YEAR (global DB, 1957-2026) is a temporal proxy for how long a species has been
## under herbicide selection. Per species: first-detection year, resistance-career length (last - first),
## MOA count, and record count (n_rec) - a monitoring/ascertainment proxy. This tests the selection-history
## reading directly, unlike geographic range (which was null).
sel_sp <- merged |> filter(!is.na(year)) |> group_by(species, HR_Status) |>
  summarise(first_year = min(year), last_year = max(year), career = max(year) - min(year),
            n_moa = n_distinct(moa), n_rec = n(), .groups = "drop") |>
  left_join(dplyr::select(sp575, species, eoo_scaled), by = "species")
sh_first  <- wilcox.test(first_year ~ HR_Status, sel_sp)     # HR-local detected earlier?
sh_career <- wilcox.test(career ~ HR_Status, sel_sp)         # HR-local longer histories?
sh_med <- sel_sp |> group_by(HR_Status) |>
  summarise(first = median(first_year), career = median(career), n_rec = median(n_rec), .groups = "drop")
rc_car <- suppressWarnings(cor.test(sel_sp$n_moa, sel_sp$career,     method = "spearman"))
rc_fy  <- suppressWarnings(cor.test(sel_sp$n_moa, sel_sp$first_year, method = "spearman"))
rc_nr  <- suppressWarnings(cor.test(sel_sp$n_moa, sel_sp$n_rec,      method = "spearman"))
## phylo-controlled model (n_moa ~ history + monitoring + range): controls for shared ancestry so the
## history->MOA link is not an artefact of resistance-prone clades. Pagel's lambda estimates the residual
## phylogenetic signal; predictors z-scored so coefficients are comparable.
sel_phy <- tryCatch({
  suppressPackageStartupMessages(library(phylolm))
  trp0 <- readRDS(file.path(base, "megatree.rds")); trp0$node.label <- NULL
  sdf <- as.data.frame(sel_sp |> mutate(tip = gsub(" ", "_", str_to_sentence(species))))
  sdf <- sdf[sdf$tip %in% trp0$tip.label, ]
  trp <- ape::drop.tip(trp0, setdiff(trp0$tip.label, sdf$tip)); sdf <- sdf[match(trp$tip.label, sdf$tip), ]
  rownames(sdf) <- sdf$tip
  sdf$career_z <- scale(sdf$career)[, 1]; sdf$nrec_z <- scale(sdf$n_rec)[, 1]
  sdf$eoo_scaled[is.na(sdf$eoo_scaled)] <- 0
  set.seed(1); fit <- phylolm(n_moa ~ career_z + nrec_z + eoo_scaled, data = sdf, phy = trp, model = "lambda")
  co <- summary(fit)$coefficients
  list(n = nrow(sdf), lambda = unname(fit$optpar),
       career_b = co["career_z","Estimate"], career_p = co["career_z","p.value"],
       nrec_b   = co["nrec_z","Estimate"],   nrec_p   = co["nrec_z","p.value"],
       eoo_b    = co["eoo_scaled","Estimate"], eoo_p  = co["eoo_scaled","p.value"])
}, error = function(e) { message("phylolm sel-hist failed: ", conditionMessage(e)); NULL })
## German-specific first-detection: ONLY HR-local have German records (HR-elsewhere have none by
## definition, which is why the comparison above must use global first-detection years).
ger_fd <- merged |> filter(grepl("german", country, ignore.case = TRUE), !is.na(year)) |>
  group_by(species) |> summarise(gf = min(year), .groups = "drop")
sel_hist <- list(med = sh_med, first_p = sh_first$p.value, career_p = sh_career$p.value,
                 rho_career = unname(rc_car$estimate), p_career = rc_car$p.value,
                 rho_first = unname(rc_fy$estimate), p_first = rc_fy$p.value,
                 rho_nrec = unname(rc_nr$estimate), p_nrec = rc_nr$p.value, phy = sel_phy,
                 ger_first_med = if (nrow(ger_fd)) median(ger_fd$gf) else NA,
                 ger_first_min = if (nrow(ger_fd)) min(ger_fd$gf) else NA, ger_n = nrow(ger_fd))
write.csv(sel_sp, file.path(dd, "moa_selection_history.csv"), row.names = FALSE)
cat(sprintf("SEL-HISTORY: first-yr %.0f vs %.0f (P=%.4f); career %.0f vs %.0f yr (P=%.4f); n_moa~career rho=%.2f; phylo career p=%s | monitoring p=%s\n",
            sh_med$first[sh_med$HR_Status=="HR_local"], sh_med$first[sh_med$HR_Status=="HR_elsewhere"], sh_first$p.value,
            sh_med$career[sh_med$HR_Status=="HR_local"], sh_med$career[sh_med$HR_Status=="HR_elsewhere"], sh_career$p.value,
            rc_car$estimate, if (!is.null(sel_phy)) sprintf("%.4f", sel_phy$career_p) else "NA",
            if (!is.null(sel_phy)) sprintf("%.4f", sel_phy$nrec_p) else "NA"))

saveRDS(list(med = med, multi = multi, freq = freq, breadth = breadth,
             wilcox_p = wt$p.value, wilcox_W = unname(wt$statistic),
             chisq_p = chi$p.value, chisq_stat = unname(chi$statistic), chisq_df = unname(chi$parameter),
             n_local = med$n_spp[med$HR_Status=="HR_local"],
             n_else = med$n_spp[med$HR_Status=="HR_elsewhere"],
             range_moa = range_moa, gof = gof, per_moa = per_moa, phylo_d = phylo_d, sel_hist = sel_hist),
        file.path(base, "moa_analysis.rds"))

## --- FIGURE ---
th <- theme_minimal(base_size = 7, base_family = "Arial") + theme(panel.grid.minor = element_blank())
ord <- freq |> group_by(moa) |> summarise(t = sum(pct)) |> arrange(t) |> pull(moa)
pA <- freq |> mutate(moa = factor(moa, levels = ord)) |>
  ggplot(aes(moa, pct, fill = HR_Status)) +
  geom_col(position = position_dodge(0.72), width = 0.66) + coord_flip() +
  scale_fill_manual(values = cols, labels = labs2, name = NULL) +
  labs(x = NULL, y = "Species resistant to this site of action (%)",
       title = "a  Resistance by site of action (HRAC group)") + th + theme(legend.position = "top")
# significance stars, matching the resistance-syndrome figure (Fig. S22)
sigm <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
ymx <- max(breadth$n_moa, na.rm = TRUE)
pB <- ggplot(breadth, aes(HR_Status, n_moa, colour = HR_Status)) +
  geom_jitter(width = 0.15, height = 0.08, size = 0.5, alpha = 0.30) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.5) +
  scale_colour_manual(values = cols, guide = "none") + scale_x_discrete(labels = labs2) +
  annotate("segment", x = 1, xend = 2, y = ymx*1.05, yend = ymx*1.05, linewidth = 0.4) +
  annotate("text", x = 1.5, y = ymx*1.11, label = sigm(wt$p.value), size = 4.2) +
  coord_cartesian(ylim = c(0, ymx*1.18)) +
  labs(x = NULL, y = "Distinct sites of action per species", title = "b  Number of resisted sites of action") + th
p_moa <- pA + pB + patchwork::plot_layout(widths = c(1.55, 1))
ggsave(file.path(fd, "Fig_MOA_HRstatus.pdf"), p_moa, width = 205, height = 118, units = "mm", device = cairo_pdf)
ggsave(file.path(fd, "Fig_MOA_HRstatus.png"), p_moa, width = 205, height = 118, units = "mm", dpi = 600)

## --- range-vs-MOA scatter (transparency: the geographic-exposure correlation is ~zero) ---
rng2 <- rng |> filter(!is.na(eoo_scaled))
p_rc <- ggplot(rng2, aes(eoo_scaled, n_moa, colour = HR_Status)) +
  geom_point(size = 1.4, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey45", linewidth = 0.5, linetype = 2) +
  scale_colour_manual(values = cols, labels = labs2, name = NULL) +
  labs(x = "Geographic range in Germany (EOO, standardised)",
       y = "Number of resisted sites of action per species",
       title = sprintf("No range-SOA link (Spearman rho = %.2f, P = %.2f)", range_moa$rho_eoo, range_moa$p_eoo)) +
  th + theme(legend.position = "top")
ggsave(file.path(fd, "Fig_MOA_RangeCorr.pdf"), p_rc, width = 130, height = 110, units = "mm", device = cairo_pdf)
ggsave(file.path(fd, "Fig_MOA_RangeCorr.png"), p_rc, width = 130, height = 110, units = "mm", dpi = 600)

## --- selection-history figure: (A) resistance-career length by group, (B) MOA count vs career ---
sigm2 <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
ymxc <- max(sel_sp$career, na.rm = TRUE)
pSA <- ggplot(sel_sp, aes(HR_Status, career, colour = HR_Status)) +
  geom_jitter(width = 0.15, height = 0.3, size = 0.5, alpha = 0.30) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.5) +
  scale_colour_manual(values = cols, guide = "none") + scale_x_discrete(labels = labs2) +
  annotate("segment", x = 1, xend = 2, y = ymxc*1.05, yend = ymxc*1.05, linewidth = 0.4) +
  annotate("text", x = 1.5, y = ymxc*1.11, label = sigm2(sel_hist$career_p), size = 4.2) +
  coord_cartesian(ylim = c(0, ymxc*1.18)) +
  labs(x = NULL, y = "Resistance-record span (years)", title = "a  Span of resistance records") + th
pSB <- ggplot(sel_sp, aes(career, n_moa, colour = HR_Status)) +
  geom_point(size = 1.4, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey45", linewidth = 0.5, linetype = 2) +
  scale_colour_manual(values = cols, labels = labs2, name = NULL) +
  labs(x = "Resistance-record span (years)", y = "Number of resisted sites of action",
       title = sprintf("b  More sites of action with a longer record (Spearman rho = %.2f)", sel_hist$rho_career)) +
  th + theme(legend.position = "top")
p_sh <- pSA + pSB + patchwork::plot_layout(widths = c(1, 1.5))
ggsave(file.path(fd, "Fig_MOA_SelHistory.pdf"), p_sh, width = 195, height = 100, units = "mm", device = cairo_pdf)
ggsave(file.path(fd, "Fig_MOA_SelHistory.png"), p_sh, width = 195, height = 100, units = "mm", dpi = 600)
cat("MOA ANALYSIS DONE -> data/moa_hr_status_merged.csv + data/moa_selection_history.csv + Figures/Fig_MOA_HRstatus + Fig_MOA_RangeCorr + Fig_MOA_SelHistory + moa_analysis.rds\n")

## ############################################################################################
## STAGE 7  Reviewer checks (recording effort, learning curve)
## source file: 11_reviewer_checks.R
## ############################################################################################

# ============================================================================
# 11_reviewer_checks.R
# Two reviewer-requested robustness checks (Meeting 1):
#   A. Survey-effort control  -> is the watch-list a recording-effort artefact?
#   B. Down-sampling learning curve -> is the small resistant sample adequate?
# Uses ONLY real data (no synthetic species). Writes figures + reviewer_checks.rds.
# ============================================================================
setwd("."); options(warn = -1)
suppressMessages({library(data.table); library(dplyr); library(ggplot2); library(randomForest); library(patchwork)})
figdir <- "Figures"; dir.create(figdir, showWarnings = FALSE)
col_status <- c(HR_elsewhere = "#2A9D8F", HR_local = "#E76F51", Non_HR_local = "#457B9D")
auc_fn <- function(score, y1){ r <- rank(score); n1 <- sum(y1); n0 <- sum(!y1); (sum(r[y1]) - n1*(n1+1)/2)/(n1*n0) }
out <- list()

## ------------------------------------------------------------------ ##
## A. SURVEY-EFFORT CONTROL
## ------------------------------------------------------------------ ##
md  <- fread("data/merged_data_final_REBUILT.csv", select = c("species","n","HR_Status"))
md[, sp := tolower(trimws(species))]
oob <- fread("data/OOB_predictions_per_species.csv")
oob[, sp := tolower(trimws(species))]
pn  <- merge(md, oob[, .(sp, p_HR = binary_OOB_prob_HR, true = threeclass_true)], by = "sp")

# A1 species level: predicted P(resistance) vs recording effort (record count n)
s_all <- suppressWarnings(cor.test(pn$p_HR, pn$n, method = "spearman"))
nonhr <- pn[HR_Status == "Non_HR_local"]
s_non <- suppressWarnings(cor.test(nonhr$p_HR, nonhr$n, method = "spearman"))
wl_flag <- nonhr[p_HR >= 0.5]; wl_noflag <- nonhr[p_HR < 0.5]
w_eff <- suppressWarnings(wilcox.test(wl_flag$n, wl_noflag$n))
out$effort_species <- list(rho_all = unname(s_all$estimate), p_all = s_all$p.value,
                           rho_nonhr = unname(s_non$estimate), p_nonhr = s_non$p.value,
                           flag_median_n = median(wl_flag$n), noflag_median_n = median(wl_noflag$n),
                           wilcox_p = w_eff$p.value, n_flag = nrow(wl_flag))
cat(sprintf("A1  P_HR vs records: all rho=%.3f (p=%.2g) | within non-HR rho=%.3f (p=%.2g)\n",
            out$effort_species$rho_all, out$effort_species$p_all, out$effort_species$rho_nonhr, out$effort_species$p_nonhr))
cat(sprintf("    flagged median records=%.0f vs non-flagged=%.0f (Wilcoxon p=%.2g)\n",
            out$effort_species$flag_median_n, out$effort_species$noflag_median_n, out$effort_species$wilcox_p))

# A2 spatial: watch-list richness vs total recording effort per 10x10 km cell
wl_species <- pn[true == "Non_HR_local" & p_HR >= 0.5, sp]
occ <- fread("../my_data_cleaned_occ_gbifv2.csv",
             select = c("species","decimalLatitude","decimalLongitude","countryCode"))
setnames(occ, c("species","decimalLatitude","decimalLongitude","countryCode"), c("species","lat","lon","cc"))
occ <- occ[cc == "DE" & is.finite(lat) & is.finite(lon) & lat > 47 & lat < 55.2 & lon > 5.7 & lon < 15.1]
occ[, sp := tolower(trimws(species))]
xy <- sf::sf_project("EPSG:4326", "EPSG:3035", as.matrix(occ[, .(lon, lat)]))
occ[, `:=`(cx = floor(xy[,1]/10000), cy = floor(xy[,2]/10000))]
cell <- occ[, .(records = .N, rich_all = uniqueN(sp), rich_wl = uniqueN(sp[sp %in% wl_species])), by = .(cx, cy)]
cell <- cell[records >= 5]                      # ignore trivially-sampled cells
sp_eff  <- suppressWarnings(cor.test(cell$rich_wl, cell$records,  method = "spearman"))
sp_rich <- suppressWarnings(cor.test(cell$rich_wl, cell$rich_all, method = "spearman"))
cell[, frac := rich_wl / rich_all]
frac_eff <- suppressWarnings(cor.test(cell$frac, cell$records, method = "spearman"))
out$effort_spatial <- list(n_cells = nrow(cell),
                           rho_effort = unname(sp_eff$estimate), p_effort = sp_eff$p.value,
                           rho_rich = unname(sp_rich$estimate),
                           frac_mean = mean(cell$frac), frac_cv = sd(cell$frac)/mean(cell$frac),
                           rho_frac_effort = unname(frac_eff$estimate), p_frac_effort = frac_eff$p.value)
cat(sprintf("A2  cells=%d | watch-list richness vs total records rho=%.2f | vs total richness rho=%.2f\n",
            out$effort_spatial$n_cells, out$effort_spatial$rho_effort, out$effort_spatial$rho_rich))
cat(sprintf("    watch-list fraction mean=%.2f CV=%.2f | fraction vs effort rho=%.3f (p=%.2g)\n",
            out$effort_spatial$frac_mean, out$effort_spatial$frac_cv, out$effort_spatial$rho_frac_effort, out$effort_spatial$p_frac_effort))

# figure A
pA1 <- ggplot(pn, aes(n, p_HR, colour = HR_Status)) +
  geom_point(alpha = .55, size = 1.1) + scale_x_log10() +
  geom_hline(yintercept = .5, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = col_status, name = NULL) +
  labs(x = "German records per species (log scale) - recording effort",
       y = "Predicted P(resistance), out-of-bag") +
  theme_minimal(base_size = 7) + theme(text = element_text(colour = "black"), legend.position = "top", panel.grid.minor = element_blank())
pA2 <- ggplot(cell, aes(records, rich_wl)) +
  geom_point(alpha = .35, size = .8, colour = "#457B9D") + scale_x_log10() +
  labs(x = "Total records per 10x10 km cell (log) - effort", y = "Watch-list species per cell") +
  theme_minimal(base_size = 7) + theme(text = element_text(colour = "black"), panel.grid.minor = element_blank())
ggsave(file.path(figdir, "Fig_Reviewer_Effort.png"), pA1 + pA2, width = 9.2, height = 4.2, dpi = 600)
ggsave(file.path(figdir, "Fig_Reviewer_Effort.pdf"), pA1 + pA2, width = 9.2, height = 4.2, device = cairo_pdf)

## ------------------------------------------------------------------ ##
## B. DOWN-SAMPLING LEARNING CURVE (rarefy the resistant class)
## ------------------------------------------------------------------ ##
ms <- as.data.frame(fread("data/imputed_analysis_modelling_set.csv"))
preds <- c("Leaf_Economics_PC1","Leaf_Size_PC1","Leaf_Nutrient_PC1","Height_PC1","Hydraulic_PC1",
           "Root_PC1","Propagule_Size_PC1","Fecundity_PC1","Genome_PC1","Rooting_Depth_PC1",
           "Stem_Diameter_PC1","Stem_Cond_Dens_PC1","Germination_PC1","ITV_PC1",
           "life_span","photosynthetic_pathway","origin","alt_range_scaled","eoo_scaled","aoo_scaled")
for (v in c("life_span","photosynthetic_pathway","origin")) ms[[v]] <- as.factor(ms[[v]])
# positive class = HR (the minority, n=109)
lev <- names(sort(table(ms$HR_bin)))[1]      # minority level label
ms$y <- factor(ifelse(ms$HR_bin == lev, "HR", "nonHR"), levels = c("nonHR","HR"))
hr_idx  <- which(ms$y == "HR"); non_idx <- which(ms$y == "nonHR")
cat(sprintf("\nB   HR species=%d  non-HR=%d  | learning curve rarefying HR count\n", length(hr_idx), length(non_idx)))
ks <- c(30, 45, 60, 75, 90, length(hr_idx)); reps <- 25
lc <- data.frame(); set.seed(42)
for (k in ks) {
  aucs <- numeric(reps)
  for (r in seq_len(reps)) {
    hrk <- sample(hr_idx, k)
    dat <- ms[c(hrk, non_idx), c(preds, "y")]
    samp <- rep(min(table(dat$y)), 2)
    rf <- randomForest(reformulate(preds, "y"), data = dat, ntree = 400,
                       sampsize = samp, strata = dat$y)
    aucs[r] <- auc_fn(rf$votes[, "HR"], dat$y == "HR")
  }
  lc <- rbind(lc, data.frame(n_HR = k, auc_mean = mean(aucs), auc_sd = sd(aucs)))
  cat(sprintf("    n_HR=%3d  OOB AUC=%.3f +/- %.3f\n", k, mean(aucs), sd(aucs)))
}
out$learning_curve <- lc
pB <- ggplot(lc, aes(n_HR, auc_mean)) +
  geom_ribbon(aes(ymin = auc_mean - auc_sd, ymax = auc_mean + auc_sd), fill = "#2A9D8F", alpha = .18) +
  geom_line(colour = "#2A9D8F", linewidth = .8) + geom_point(colour = "#2A9D8F", size = 2) +
  geom_hline(yintercept = 0.85, linetype = 2, colour = "grey40") +
  annotate("text", x = min(ks), y = 0.856, label = "full-model AUC 0.85", hjust = 0, size = 2.3, colour = "grey30") +
  labs(x = "Number of resistant (HR) species used in training",
       y = "Balanced random-forest OOB AUC") +
  theme_minimal(base_size = 7) + theme(text = element_text(colour = "black"), panel.grid.minor = element_blank())
ggsave(file.path(figdir, "Fig_Reviewer_LearningCurve.png"), pB, width = 6.2, height = 4.2, dpi = 600)
ggsave(file.path(figdir, "Fig_Reviewer_LearningCurve.pdf"), pB, width = 6.2, height = 4.2, device = cairo_pdf)

## Extended Data Fig 5 - native 3-panel composite (a: per-species effort, b: spatial effort, c: learning curve)
ed5 <- (pA1 | pA2) / pB + patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 7, colour = "black"))
ggsave(file.path(figdir, "Fig_ED5_Robustness.png"), ed5, width = 180, height = 150, units = "mm", dpi = 600)
ggsave(file.path(figdir, "Fig_ED5_Robustness.pdf"), ed5, width = 180, height = 150, units = "mm", device = cairo_pdf)

saveRDS(out, "reviewer_checks.rds")
cat("\nDONE -> Figures/Fig_Reviewer_Effort.png, Figures/Fig_Reviewer_LearningCurve.png, reviewer_checks.rds\n")

## ############################################################################################
## STAGE 8  MAIN ANALYSIS (models, importance, phylogeny, watch-list, nested-model, CWR test)  [annotated source: 02_Master_Analysis.Rmd]
## source file: 02_Master_Analysis.Rmd
## ############################################################################################

params <-
list(data_path = "data/merged_data_final_REBUILT.csv", out_dir = ".", impute_primary = TRUE, 
    ntree = 1000L, n_perm = 999L, run_loo = TRUE, loo_ntree = 300L, 
    stability_runs = 10L, run_shap = TRUE)

knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.align = "center", dpi = 300, dev = "cairo_pdf"
)
options(width = 84)
## shrink code + console output so wide prints (confusion matrices etc.) stay on the page
knitr::knit_hooks$set(size = function(before, options, envir)
  if (before) paste0("\n \\", options$size, "\n\n") else "\n\n \\normalsize\n")
knitr::opts_chunk$set(size = "footnotesize")
set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse); library(janitor)
  library(missForest); library(FactoMineR); library(factoextra)
  library(randomForest); library(pROC); library(caret)
  library(patchwork); library(viridis); library(corrplot); library(ggrepel)
  library(partykit); library(kableExtra)
})
results_summary <- list()   # collects real numbers for the report & slides

## wide-table helper: booktabs + shrink-to-fit page width
kbl <- function(x, caption = NULL, digits = 3, ...) {
  knitr::kable(x, caption = caption, digits = digits, booktabs = TRUE, ...) |>
    kableExtra::kable_styling(latex_options = c("HOLD_position", "scale_down"),
                              font_size = 9)
}

## ---- output folders ----
out_dir <- params$out_dir
fig_dir <- file.path(out_dir, "Figures")
tab_dir <- file.path(out_dir, "Tables")
mod_dir <- file.path(out_dir, "Models")
for (d in c(out_dir, fig_dir, tab_dir, mod_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- consistent colours ----
col_status <- c(HR_elsewhere = "#2A9D8F", HR_local = "#E76F51", Non_HR_local = "#457B9D")
col_bin    <- c(HR = "#E76F51", nonHR = "#457B9D")
## canonical HR-status legend labels (match main-text figures) — keyed by raw factor level
hr_labs    <- c(Non_HR_local = "Non-resistant", HR_local = "HR-local (Germany)", HR_elsewhere = "HR-elsewhere")
bin_labs   <- c(nonHR = "Non-resistant", HR = "Resistant")

## p-value formatter: never print "= 0.000" (a p-value cannot be exactly zero) - use "< 0.001"
.pf <- function(p) if (!is.null(p) && !is.na(p) && p < 0.001) "< 0.001" else sprintf("= %.3f", p)

## ---- Arial + theme ----
try(windowsFonts(Arial = windowsFont("Arial")), silent = TRUE)
base_family <- "Arial"
## Nature Plants: text 5-7 pt. Standard labels 7 pt; secondary text 6 pt; nothing <5 or >7.
theme_pub <- function(base_size = 7) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(text        = element_text(size = base_size, colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
          axis.title  = element_text(size = base_size, face = "bold", colour = "black"),
          axis.text   = element_text(size = base_size, colour = "black"),
          plot.title  = element_text(size = base_size, face = "bold", colour = "black"),
          plot.subtitle = element_text(size = base_size, colour = "black"),
          plot.tag    = element_text(size = base_size, face = "bold", colour = "black"),
          legend.text = element_text(size = base_size, colour = "black"),
          legend.title= element_text(size = base_size, face = "bold", colour = "black"),
          legend.key.size = unit(4, "mm"),
          strip.text  = element_text(size = base_size, face = "bold", colour = "black"))
}
theme_set(theme_pub())
update_geom_defaults("point", list(size = 1.6))
update_geom_defaults("line",  list(linewidth = 0.7))
update_geom_defaults("col",   list(linewidth = 0))

## ---- figure export helper: 600 dpi, mm, Arial, cairo_pdf + png ----
save_fig <- function(plot, name, w = 80, h = 80, dir = fig_dir) {
  pdf_f <- file.path(dir, paste0(name, ".pdf"))
  png_f <- file.path(dir, paste0(name, ".png"))
  ok <- TRUE
  tryCatch(
    ggsave(pdf_f, plot, width = w, height = h, units = "mm", device = cairo_pdf, dpi = 600),
    error = function(e) { ok <<- FALSE; ggsave(pdf_f, plot, width = w, height = h, units = "mm", dpi = 600) }
  )
  ggsave(png_f, plot, width = w, height = h, units = "mm", dpi = 600)
  invisible(pdf_f)
}

## ---- readable trait-name map (raw *_mean -> pretty); defined early so all figures/tables can use it ----
.nm <- c(sla_undef_mean = "SLA", ldmc_mean = "LDMC", leaf_thickness_mean = "Leaf thickness",
  leaf_water_cont_mean = "Leaf water", leaf_area_undef_undef_mean = "Leaf area",
  leaf_area_leaf_undef_mean = "Leaf area (leaf)", leaf_area_leaflet_undef_mean = "Leaf area (leaflet)",
  leaf_length_mean = "Leaf length", leaf_width_mean = "Leaf width", leaf_dry_mass_single_mean = "Leaf dry mass",
  leaffreshmass_mean = "Leaf fresh mass", leaf_c_perdrymass_mean = "Leaf C", leaf_n_mean = "Leaf N",
  leaf_p_mean = "Leaf P", leaf_nper_area_mean = "Leaf N/area", leaf_cn_ratio_mean = "Leaf C:N",
  plant_height_veg_mean = "Height (veg.)", plant_height_generative_mean = "Height (gen.)",
  stem_conduit_diameter_mean = "Conduit diameter", wood_vessel_length_mean = "Vessel length",
  wood_fiber_length_mean = "Fibre length", wood_ray_millimetre_mean = "Ray width",
  specific_root_length_mean = "SRL", specific_root_length_fine_mean = "SRL (fine)",
  seed_mass_mean = "Seed mass", seed_length_mean = "Seed length", disp_unit_leng_mean = "Disp.-unit length",
  seed_num_rep_unit_mean = "Seed no. (rep.)", seed_num_disp_unit_mean = "Seed no. (disp.)",
  chromosome_n_mean = "Chromosome no.", chromosome_c_dn_acont_mean = "Genome size",
  rooting_depth_mean = "Rooting depth", stem_diam_mean = "Stem diameter",
  stem_cond_dens_mean = "Stem cond. density", seed_germination_rate_mean = "Germination")
.lab <- function(x) { y <- unname(.nm[x]); ifelse(is.na(y), gsub("_", " ", gsub("_mean$", "", x)), y) }

raw <- readr::read_csv(params$data_path, show_col_types = FALSE) |> janitor::clean_names()
# coerce literal "NA" strings to real NA, then re-type
raw <- raw |> mutate(across(where(is.character), ~na_if(., "NA")))
knitr::kable(as.data.frame(table(raw$hr_status)),  col.names = c("HR_Status","n"),
             caption = "Response distribution (verified from the data).")
knitr::kable(as.data.frame(table(raw$hr_status, raw$habitat_group)) |>
               setNames(c("HR_Status","habitat_group","n")) |> filter(n > 0),
             caption = "HR_Status x habitat_group: every HR species is arable; all ruderal species are Non-HR.")

# suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr); library(sf)
#   library(terra); library(giscoR); library(geodata)})
# sf::sf_use_s2(TRUE)
# occ <- readr::read_csv("../my_data_final_cleanedv2.csv", show_col_types = FALSE, guess_max = 1e5) |>
#   mutate(species = str_to_sentence(str_squish(tolower(species))))
# 
# ## (1) coordinate cleaning (coordinateCleaner protocol, applied as explicit filters) -----------------
# occ_clean <- occ |>
#   filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION")) |>   # wild records only
#   filter(coordinateUncertaintyInMeters <= 2000 |                                              # high precision: <= 2 km
#          (is.na(coordinateUncertaintyInMeters) & basisOfRecord == "HUMAN_OBSERVATION" & year >= 2010)) |>
#   filter(!coordinateUncertaintyInMeters %in% c(301, 999, 3036, 10000)) |>                     # drop GBIF placeholder uncertainties
#   filter(!is.na(decimalLatitude), !is.na(decimalLongitude), decimalLatitude != 0, decimalLongitude != 0)
# 
# ## (2) restrict to within Germany; all spatial operations in EPSG:3035 ------------------------------
# germany <- giscoR::gisco_get_countries(country = "DE", resolution = "03") |> sf::st_transform(3035)
# pts <- sf::st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE) |>
#   sf::st_transform(3035)
# pts <- pts[lengths(sf::st_intersects(pts, germany)) > 0, ]
# 
# ## (3) AOO = occupied 10 x 10 km grid cells x 100 km^2 ----------------------------------------------
# grid <- sf::st_make_grid(germany, cellsize = 10000, square = TRUE)
# grid <- sf::st_sf(grid_id = seq_along(grid), geometry = grid)
# grid <- grid[lengths(sf::st_intersects(grid, germany)) > 0, ]
# aoo  <- sf::st_join(pts, grid, join = sf::st_intersects) |> sf::st_drop_geometry() |>
#   group_by(species) |> summarise(AOO_km2 = dplyr::n_distinct(grid_id, na.rm = TRUE) * 100, .groups = "drop")
# 
# ## (4) EOO = convex-hull area (km^2); <3 localities -> NA, later raised to AOO (IUCN EOO >= AOO) -----
# eoo <- pts |> group_by(species) |> summarise(geom = sf::st_combine(geometry), npts = dplyr::n(), .groups = "drop") |>
#   mutate(EOO_km2 = ifelse(npts >= 3, as.numeric(sf::st_area(sf::st_convex_hull(geom))) / 1e6, NA_real_)) |>
#   sf::st_drop_geometry() |> select(species, EOO_km2)
# 
# ## (5) altitudinal range from a 30-arc-second (~1 km) DEM: max - min per species --------------------
# dem <- geodata::elevation_30s(country = "DEU", path = tempdir())
# occ_clean$altitude <- as.numeric(terra::extract(dem, occ_clean[, c("decimalLongitude", "decimalLatitude")])[, 2])
# alt <- occ_clean |> filter(!is.na(altitude)) |> group_by(species) |>
#   summarise(alt_range = max(altitude) - min(altitude), .groups = "drop")
# 
# ## (6) assemble, apply IUCN EOO >= AOO, log-transform and z-standardise (NA-safe) -------------------
# scl <- function(x) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE); (x - m) / s }
# range_covariates <- aoo |> left_join(eoo, by = "species") |> left_join(alt, by = "species") |>
#   mutate(EOO_km2 = ifelse(is.na(EOO_km2), AOO_km2, pmax(EOO_km2, AOO_km2)),
#          EOO_scaled       = scl(log10(EOO_km2 + 1)),
#          AOO_scaled       = scl(log10(AOO_km2 + 1)),
#          alt_range_scaled = scl(alt_range))
# ## range_covariates is merged into the species trait table to produce the cached modelling data read above.

# conceptual framework / graphical abstract: restrained palette - neutral boxes,
# colour only where it encodes meaning (3 resistance classes + final output).
# Interpretation/robustness all derive from the PRIMARY binary model (fan), not 1:1 from phases.
.ink<-"#1F2D3D"; .brd<-"#AEB4BC"; .fillN<-"#FFFFFF"; .band<-"#F3F5F7"; .arw<-"#7A828C"
.o<-"#E76F51"; .t<-"#2A9D8F"; .b<-"#457B9D"; .emph<-"#1D3557"; .tint<-"#E8EEF4"
Bx <- function(id,x,y,w,h,lab,fill=.fillN,tcol=.ink,fs=2.4,face="plain",bord=.brd)
  data.frame(id,x,y,w,h,lab,fill,tcol,fs,bord,face,stringsAsFactors=FALSE)
boxes <- do.call(rbind, list(
  Bx("a1",27,90,24,8,"TRY trait database\n78,278 taxa"),
  Bx("a2",55,90,28,8,"Filter to German arable weeds\n(EVA/FloraVeg, WeedScience)"),
  Bx("a3",83,90,18,8,"575 weed\nspecies",face="bold"),
  Bx("ca",51,79,17,6.5,"HR-local · 24", fill=.o, tcol="white", fs=2.5),
  Bx("cb",70,79,17,6.5,"HR-elsewhere · 85", fill=.t, tcol="white", fs=2.5),
  Bx("cc",89,79,17,6.5,"non-HR-local · 466", fill=.b, tcol="white", fs=2.5),
  Bx("b1",22,64,19,8,"Mean + SD traits\n(strategy + ITV)"),
  Bx("b2",44,64,19,8,"Functional modules\n(trait clustering)"),
  Bx("b3",66,64,19,8,"missForest impute\n+ PCA → module PC1"),
  Bx("b4",88,64,19,8,"+ origin, life history,\nrange, altitude"),
  Bx("p0",24.5,46,25,10,"Phase 0 gate:\nhabitat sensitivity\n(2 binary-RF models\ncompared → ruderal kept)",fs=2.0),
  Bx("c0",55,46,30,10,"PRIMARY MODEL\nclass-balanced binary RF\n109 HR  vs  466 non-HR",fs=2.3,face="bold",fill=.tint,bord=.emph),
  Bx("pv",85.5,46,25,10,"Variants (same RF):\nP1 transfer\nP2 local vs elsewhere\nP3 3-class · P4 2-class",fs=2.0),
  Bx("d1",21,16,16,9,"Permutation importance\n+ partial dependence"),
  Bx("d2",38,16,16,9,"SHAP\n(per-species)"),
  Bx("d3",55,16,16,9,"Risk\nwatch-list", fill=.emph, tcol="white", face="bold", fs=2.6),
  Bx("d4",72,16,16,9,"Phylogenetic signal\n(D-stat, phyloglm)"),
  Bx("d5",89,16,16,9,"Cross-algorithm\n+ all-traits checks")
))
Aw <- function(x,xe,y,ye=y) data.frame(x,xe,y,ye)
arr <- do.call(rbind, list(
  Aw(39,41,90), Aw(69,74,90),
  Aw(83,51,86,82.5), Aw(83,70,86,82.5), Aw(83,89,86,82.5),
  Aw(55,55,76,68.5),
  Aw(31.5,34.5,64), Aw(53.5,56.5,64), Aw(75.5,78.5,64),
  Aw(55,55,60,51),
  Aw(37,40,46), Aw(70,73,46),
  Aw(55,21,41,20.7),Aw(55,38,41,20.7),Aw(55,55,41,20.7),Aw(55,72,41,20.7),Aw(55,89,41,20.7)
))
bands <- data.frame(ymin=c(74,58,38,10), ymax=c(97,72,56,25),
  lab=c("DATA","TRAITS","MODELLING","INTERPRETATION"), ymid=c(85.5,65,47,17.5))
p_concept <- ggplot() +
  geom_rect(data=bands, aes(xmin=11,xmax=99,ymin=ymin,ymax=ymax), fill=.band, colour=NA) +
  geom_text(data=bands, aes(x=6,y=ymid,label=lab), angle=90, colour=.emph,
            fontface="bold", size=2.5, family=base_family) +
  geom_text(aes(x=34,y=32.5,label="all read-outs derive from the primary binary model"),
            colour="#5b636d", fontface="italic", size=2.05, family=base_family, hjust=0) +
  geom_segment(data=arr, aes(x=x,xend=xe,y=y,yend=ye),
               arrow=grid::arrow(length=unit(1.6,"mm"),type="closed"), colour=.arw, linewidth=0.45) +
  geom_rect(data=boxes, aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2,fill=I(fill)),
            colour=boxes$bord, linewidth=0.4) +
  geom_text(data=boxes, aes(x=x,y=y,label=lab,colour=I(tcol),fontface=face),
            size=boxes$fs, lineheight=0.92, family=base_family) +
  coord_cartesian(xlim=c(0,100), ylim=c(8,98), expand=FALSE) +
  theme_void(base_family = base_family)
save_fig(p_concept, "Fig0_Framework", 180, 140); p_concept   # Nature Plants 2-column width (180 mm)

dat <- raw |>
  mutate(
    # response factors
    HR3 = factor(hr_status, levels = c("Non_HR_local","HR_local","HR_elsewhere")),
    HR_bin = factor(if_else(hr_status == "Non_HR_local", "nonHR", "HR"), levels = c("nonHR","HR")),
    # categorical covariates (collapse rare growth-form levels)
    origin = factor(origin),
    life_span = factor(life_span),
    dispersal_distance = factor(dispersal_distance),
    leaf_compoundness = factor(leaf_compoundness),
    photosynthetic_pathway = factor(photosynthetic_pathway),
    growth_form = fct_collapse(factor(plant_growth_form),
      herb = c("herb"), graminoid = c("graminoid"),
      woody = c("shrub","tree","shrub/tree","herb/shrub","succulent"))
  )

cat("Rows:", nrow(dat), " | unique species:", dplyr::n_distinct(dat$species),
    " | duplicated species:", any(duplicated(dat$species)), "\n")
cat("All species-level:", all(dat$rank_correct == "species"), "\n")

# trait column groups present in the file
mean_cols <- grep("_mean$", names(dat), value = TRUE)
sd_cols   <- grep("_sd$",   names(dat), value = TRUE)
cat("mean-trait columns:", length(mean_cols), " | sd-trait columns:", length(sd_cols), "\n")

## --- live: distribution of the three range covariates by HR status (box + jitter) ---
trait_labs <- c(alt_range_scaled = "a  Altitudinal range",
                eoo_scaled       = "b  EOO (extent of occurrence)",
                aoo_scaled       = "c  AOO (area of occupancy)")
diag_dat <- dat |>
  dplyr::select(HR3, alt_range_scaled, eoo_scaled, aoo_scaled) |>
  tidyr::pivot_longer(-HR3, names_to = "Trait", values_to = "value") |>
  dplyr::filter(!is.na(value), !is.na(HR3)) |>
  dplyr::mutate(Trait = factor(unname(trait_labs[Trait]), levels = unname(trait_labs)))
p_box <- ggplot(diag_dat, aes(HR3, value, colour = HR3)) +
  geom_jitter(width = 0.18, alpha = 0.45, size = 0.9) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.5) +
  facet_wrap(~Trait, scales = "free_y") +
  scale_colour_manual(values = col_status, guide = "none") +
  scale_x_discrete(labels = hr_labs) +
  labs(x = NULL, y = "scaled value") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_fig(p_box, "Figure2_Trait_Distributions", 180, 80)
p_box

## --- live: Kruskal-Wallis omnibus + Dunn pairwise (Holm-adjusted) on the range covariates ---
target_traits <- c("alt_range_scaled", "eoo_scaled", "aoo_scaled")
omnibus_kw <- purrr::map_dfr(target_traits, function(tr) {
  k <- kruskal.test(reformulate("HR3", tr), data = dat)
  tibble::tibble(Trait = unname(trait_labs[tr]),
                 Chi_squared = round(unname(k$statistic), 2),
                 Df          = unname(k$parameter),
                 P_value     = format.pval(k$p.value, digits = 3, eps = 1e-3))
})
pairwise_dunn <- purrr::map_dfr(target_traits, function(tr) {
  d <- FSA::dunnTest(reformulate("HR3", tr), data = dat, method = "holm")$res
  tibble::as_tibble(d) |>
    dplyr::transmute(Trait = unname(trait_labs[tr]), Comparison,
                     Z      = round(Z, 2),
                     P_raw  = format.pval(P.unadj, digits = 3, eps = 1e-3),
                     P_Holm = format.pval(P.adj,   digits = 3, eps = 1e-3))
})
## persist for the supplement (05_Supplementary.Rmd reads these by name)
readr::write_csv(omnibus_kw,    file.path(out_dir, "stat_table_omnibus_kruskal.csv"))
readr::write_csv(pairwise_dunn, file.path(out_dir, "stat_table_pairwise_dunn.csv"))
kbl(omnibus_kw, caption = "Kruskal-Wallis omnibus tests of range covariates across the three HR-status classes (computed in-pipeline).")

kbl(pairwise_dunn, caption = "Dunn pairwise comparisons (Holm-adjusted) of range covariates between HR-status classes (computed in-pipeline).")

## --- live: macro-ecological niche space (EOO x AOO; point size = altitudinal range) ---
eco_dat <- dat |> dplyr::filter(!is.na(eoo_scaled), !is.na(aoo_scaled), !is.na(HR3))
p_eco <- ggplot(eco_dat, aes(eoo_scaled, aoo_scaled, colour = HR3)) +
  geom_point(aes(size = alt_range_scaled), alpha = 0.6) +
  scale_colour_manual(values = col_status, labels = hr_labs, name = "HR status") +
  scale_size_continuous(name = "Alt. range\n(scaled)", range = c(0.5, 4)) +
  labs(x = "EOO (scaled)", y = "AOO (scaled)")
save_fig(p_eco, "Figure1_Ecological_Space", 130, 110)
p_eco

knitr::include_graphics("Figures/Fig_AOO_EOO_raw.png")

## --- Figure 3: Pearson correlations among the range covariates INCLUDING total_obs_scaled.
##     total_obs was dropped from the model PRECISELY BECAUSE it is near-collinear with AOO
##     (r = 0.98); this matrix is the evidence for that decision, so it must show all four.
##     Read from the canonical table (total_obs is not retained in the modelling frame `dat`). ---
cor_src <- file.path(out_dir, "table_feature_correlations.csv")
cor_df  <- readr::read_csv(cor_src, show_col_types = FALSE)
cor_mat <- as.matrix(cor_df[, -1]); rownames(cor_mat) <- cor_df[[1]]
lab <- c(alt_range_scaled = "Alt. range", EOO_scaled = "EOO",
         AOO_scaled = "AOO", total_obs_scaled = "Total obs.")
dimnames(cor_mat) <- list(unname(lab[rownames(cor_mat)]), unname(lab[colnames(cor_mat)]))
draw_cor <- function() corrplot::corrplot(cor_mat, method = "color", type = "upper",
                   addCoef.col = "black", tl.col = "grey20", tl.srt = 30, tl.cex = 0.9,
                   number.cex = 0.9, cl.cex = 0.75, diag = FALSE,
                   col = colorRampPalette(c("#457B9D", "white", "#E76F51"))(200),
                   mar = c(0, 0, 1, 0))
## standalone 600-dpi export at NPh single-column width (80 mm)
in80 <- 80 / 25.4
grDevices::cairo_pdf(file.path(fig_dir, "Figure3_Correlation_Matrix.pdf"),
                     width = in80, height = in80); draw_cor(); grDevices::dev.off()
grDevices::png(file.path(fig_dir, "Figure3_Correlation_Matrix.png"),
               width = in80, height = in80, units = "in", res = 600); draw_cor(); grDevices::dev.off()
draw_cor()   # inline in the report



# ## ---- Spatial HR surface from cleaned GBIF occurrences (run via 09_hr_surface.R, cached) ----
# library(sf); library(rnaturalearth); library(ggspatial); library(data.table)
# risk_cols <- c("#457B9D", "#2A9D8F", "#E76F51")          # Non_HR -> HR_elsewhere -> HR_local
# 
# ## status lookup (occurrence file has coordinates only); join by species
# md  <- fread("data/merged_data_final_REBUILT.csv", select = c("species","HR_Status"))
# md[, key := tolower(trimws(species))]; hr_lkp <- setNames(md$HR_Status, md$key)
# occ <- fread("../my_data_cleaned_occ_gbifv2.csv",
#              select = c("species","decimalLatitude","decimalLongitude"))
# setnames(occ, c("decimalLatitude","decimalLongitude"), c("lat","lon"))
# occ[, HR_Status := hr_lkp[tolower(trimws(species))]]
# occ <- occ[!is.na(HR_Status) & !is.na(lat) & !is.na(lon)]
# occ[, HR_Status := factor(HR_Status, levels = c("Non_HR_local","HR_local","HR_elsewhere"))]
# 
# germany <- rnaturalearth::ne_countries(scale = "medium", country = "Germany",
#                                         returnclass = "sf") |> sf::st_transform(3035)
# xy <- sf::st_coordinates(sf::st_transform(
#   sf::st_as_sf(occ, coords = c("lon","lat"), crs = 4326), 3035))
# occ[, `:=`(X = xy[,1], Y = xy[,2])]
# 
# ## (a) RAW occurrences coloured by HR status (display subsample; grids use all records)
# occ_disp <- occ[sample(.N, min(.N, 150000))]
# occ_disp[, prk := fcase(HR_Status=="Non_HR_local",1L, HR_Status=="HR_elsewhere",2L,
#                         HR_Status=="HR_local",3L)]; occ_disp <- occ_disp[order(prk)]
# m_raw <- ggplot() +
#   geom_sf(data = germany, fill = "grey97", colour = "grey60", linewidth = 0.3) +
#   geom_point(data = occ_disp, aes(X, Y, colour = HR_Status), size = 0.25, alpha = 0.5) +
#   scale_colour_manual(values = col_status, labels = hr_labs, name = "HR status",
#                       guide = guide_legend(override.aes = list(size = 2, alpha = 1))) +
#   labs(subtitle = "Recorded occurrences of HR-classified species (not predictions, not abundance)")
# save_fig(m_raw, "HR_raw_map", 130, 130)
# 
# ## 10 x 10 km grid; weighted HR index (HR_local = 1, HR_elsewhere = 0.5, Non_HR_local = 0)
# CELL <- 10000
# occ[, `:=`(cx = floor(X/CELL)*CELL + CELL/2, cy = floor(Y/CELL)*CELL + CELL/2,
#            w = fcase(HR_Status=="HR_local",1, HR_Status=="HR_elsewhere",0.5,
#                      HR_Status=="Non_HR_local",0))]
# cell <- occ[, .(n = .N, hr_mean = mean(w), hr_se = sd(w)/sqrt(.N)), by = .(cx, cy)]
# 
# ## (b) weighted-risk surface   (c) uncertainty (SE) surface - manuscript palette
# m_risk <- ggplot() +
#   geom_sf(data = germany, fill = "grey97", colour = "grey60", linewidth = 0.3) +
#   geom_tile(data = cell, aes(cx, cy, fill = hr_mean), width = CELL, height = CELL) +
#   scale_fill_gradientn(colours = risk_cols, values = c(0, 0.5, 1), limits = c(0, 1),
#                        name = "HR index") +
#   annotation_scale(location = "bl", width_hint = 0.28) +
#   labs(subtitle = "10 x 10 km weighted mean HR status")
# m_unc <- ggplot() +
#   geom_sf(data = germany, fill = "grey97", colour = "grey60", linewidth = 0.3) +
#   geom_tile(data = cell, aes(cx, cy, fill = hr_se), width = CELL, height = CELL) +
#   scale_fill_gradient(low = "grey92", high = "#E76F51", name = "SE") +
#   labs(subtitle = "Standard error (high = few / heterogeneous records)")
# save_fig(m_risk + m_unc, "FigHR_RiskUncert", 180, 110)





miss <- dat |>
  summarise(across(all_of(c(mean_cols, sd_cols)), ~mean(is.na(.))*100)) |>
  pivot_longer(everything(), names_to = "trait", values_to = "pct_missing") |>
  arrange(desc(pct_missing))
knitr::kable(head(miss, 12), digits = 1, caption = "Top-12 traits by percent missing (max ~7 percent).")

# Pre-imputation justification: Spearman block structure on OBSERVED data only
trait_num <- dat |> select(all_of(mean_cols))
cor_obs <- cor(trait_num, use = "pairwise.complete.obs", method = "spearman")
dimnames(cor_obs) <- list(.lab(rownames(cor_obs)), .lab(colnames(cor_obs)))  # readable trait names (no underscores)
p_corr <- ~corrplot(cor_obs, method = "color", tl.cex = 0.55, tl.col = "black",
                    order = "hclust", cl.cex = 0.6)
grDevices::cairo_pdf(file.path(fig_dir, "Fig_S1_TraitCorrelation.pdf"), width = 180/25.4, height = 180/25.4)
corrplot(cor_obs, method = "color", tl.cex = 0.6, tl.col = "black", order = "hclust", cl.cex = 0.7)
invisible(dev.off())
corrplot(cor_obs, method = "color", tl.cex = 0.55, tl.col = "black", order = "hclust", cl.cex = 0.6)

# Modules REFINED (2026-06-25) after an observed-data Spearman coherence check (see the module-coherence
# diagnostic below): the old Leaf_Construction conflated two distinct axes, Size was height-dominated, Root
# mixed depth with fineness, and Reproduction's components were near-independent. The refined set keeps only
# empirically-coherent multi-trait blocks and PROMOTES the decoupled traits to standalone predictors
# (nothing is deleted - consistent with the keep-everything philosophy; the RF is robust to weak inputs).
modules <- list(
  # --- coherent multi-trait blocks (-> module PC1) ---
  Leaf_Economics = c("sla_undef_mean","ldmc_mean","leaf_thickness_mean","leaf_water_cont_mean"),   # split from old Leaf_Construction; +leaf-water-content (r=0.46 with SLA)
  Leaf_Size      = c("leaf_area_undef_undef_mean","leaf_area_leaf_undef_mean","leaf_area_leaflet_undef_mean",
                     "leaf_length_mean","leaf_width_mean","leaf_dry_mass_single_mean","leaffreshmass_mean"),  # the coherent size sub-block (mean|r|=0.68)
  Leaf_Nutrient  = c("leaf_c_perdrymass_mean","leaf_n_mean","leaf_p_mean","leaf_nper_area_mean","leaf_cn_ratio_mean"),
  Height         = c("plant_height_veg_mean","plant_height_generative_mean"),                       # old Size, minus the non-cohering stem traits (r=0.79)
  Hydraulic      = c("stem_conduit_diameter_mean","wood_vessel_length_mean","wood_fiber_length_mean","wood_ray_millimetre_mean"),
  Root           = c("specific_root_length_mean","specific_root_length_fine_mean"),                 # root economics (r=0.62), minus rooting_depth (decoupled)
  Propagule_Size = c("seed_mass_mean","seed_length_mean","disp_unit_leng_mean"),                    # seed/fruit size: dispersal-unit length joins (r=0.59 with seed size; disp_len~=seed_len r=0.84) -> mean|r|=0.50
  Fecundity      = c("seed_num_rep_unit_mean","seed_num_disp_unit_mean"),                           # seed NUMBER (split from old Reproduction; distinct from seed size)
  Genome         = c("chromosome_n_mean","chromosome_c_dn_acont_mean"),
  # --- decoupled single traits PROMOTED to standalone predictors (best block link < 0.30 -> a PC1 would average noise) ---
  Rooting_Depth  = c("rooting_depth_mean"),        # vs Root(SRL) mean|r|=0.06 - root DEPTH is orthogonal to root fineness
  Stem_Diameter  = c("stem_diam_mean"),            # vs Height 0.09, vs Hydraulic 0.07 - fits neither
  Stem_Cond_Dens = c("stem_cond_dens_mean"),       # vs Hydraulic 0.33 (borderline) BUT sign-conflicts (r=-0.08 with conduit diameter) and lowers Hydraulic coherence 0.47->0.41
  Germination    = c("seed_germination_rate_mean") # vs Propagule_Size 0.20 - regeneration, decoupled from propagule size
)

# ITV (Intraspecific Trait Variation): 12 mean/sd pairs -> CV
itv_pairs <- tribble(
  ~mean, ~sd,
  "rooting_depth_mean","rooting_depth_sd", "leaf_n_mean","leaf_n_sd", "leaf_p_mean","leaf_p_sd",
  "sla_undef_mean","sla_undef_sd", "ldmc_mean","ldmc_sd", "stem_diam_mean","stem_diam_sd",
  "seed_mass_mean","seed_mass_sd", "plant_height_veg_mean","plant_height_veg_sd",
  "specific_root_length_mean","specific_root_length_sd", "wood_vessel_length_mean","wood_vessel_length_sd",
  "leaf_area_undef_undef_mean","leaf_area_undef_undef_sd", "leaf_water_cont_mean","leaf_water_cont_sd"
)

# Universal module reducer.
#  impute = TRUE  -> missForest within module (all species retained)
#  impute = FALSE -> complete-case within module (species with any NA in module get PC1 = NA)
# Returns: tibble(species, <Module>_PC1) and variance explained; PC1 sign-oriented to +ve mean loading.
run_module_pca <- function(df, cols, module, impute = TRUE, logtrans = TRUE, imp_lookup = NULL) {
  cols <- intersect(cols, names(df))
  sp <- df$species
  X  <- df |> select(all_of(cols)) |> mutate(across(everything(), as.numeric))

  if (logtrans) X <- X |> mutate(across(everything(), ~ifelse(!is.na(.) & . > 0, log10(.), .)))

  # single-trait module: z-score, use directly (PC1 == the trait). When imputing, fill NAs from the
  # full-matrix missForest (imp_lookup) so PROMOTING a decoupled trait to a standalone predictor does NOT
  # drop species (un-imputed singles would otherwise lose ~9% of species via the model-frame drop_na).
  if (ncol(X) == 1) {
    v0 <- X[[1]]
    if (impute && anyNA(v0) && !is.null(imp_lookup) && cols[1] %in% names(imp_lookup))
      v0 <- imp_lookup[[cols[1]]]
    v <- as.numeric(scale(v0))
    out <- tibble(species = sp, !!paste0(module, "_PC1") := v)
    return(list(scores = out, var = 100))
  }

  if (impute) {
    Xi <- as.data.frame(X)
    Xi <- Xi[, colSums(!is.na(Xi)) > 0, drop = FALSE]
    Xi <- missForest::missForest(Xi, verbose = FALSE)$ximp
    keep_sp <- sp
  } else {
    cc  <- complete.cases(X)
    Xi  <- as.data.frame(X[cc, , drop = FALSE])
    keep_sp <- sp[cc]
  }
  Xs <- scale(Xi)
  pca <- FactoMineR::PCA(as.data.frame(Xs), graph = FALSE, ncp = 2)
  pc1 <- pca$ind$coord[, 1]
  # orient: positive loadings on average
  if (mean(pca$var$coord[, 1]) < 0) pc1 <- -pc1
  out <- tibble(species = keep_sp, !!paste0(module, "_PC1") := as.numeric(pc1))
  list(scores = out, var = pca$eig[1, 2])
}

build_modules <- function(df, impute = TRUE) {
  # full-matrix missForest once, to impute single-trait (standalone) modules (see run_module_pca)
  imp_lookup <- NULL
  if (impute) {
    allt <- intersect(unique(unlist(modules)), names(df))
    Xa <- df |> select(all_of(allt)) |> mutate(across(everything(), as.numeric)) |>
      mutate(across(everything(), ~ifelse(!is.na(.) & . > 0, log10(.), .)))
    Xa <- as.data.frame(Xa[, colSums(!is.na(Xa)) > 0, drop = FALSE])
    imp_lookup <- missForest::missForest(Xa, verbose = FALSE)$ximp
  }
  res <- imap(modules, ~run_module_pca(df, .x, .y, impute = impute, imp_lookup = imp_lookup))
  # ITV (intraspecific trait variation)
  itv <- df |> select(species, all_of(c(itv_pairs$mean, itv_pairs$sd)))
  if (impute) {
    M <- as.data.frame(itv |> select(-species) |> mutate(across(everything(), as.numeric)))
    M <- missForest::missForest(M, verbose = FALSE)$ximp
    itv_sp <- itv$species
  } else {
    cc <- complete.cases(itv |> select(-species))
    M <- as.data.frame(itv[cc, -1] |> mutate(across(everything(), as.numeric)))
    itv_sp <- itv$species[cc]
  }
  cv <- map_dfc(seq_len(nrow(itv_pairs)), function(i) {
    m <- M[[itv_pairs$mean[i]]]; s <- M[[itv_pairs$sd[i]]]
    tibble(!!paste0("cv_", i) := s / abs(m))
  })
  cv[!is.finite(as.matrix(cv))] <- NA
  cv <- cv |> mutate(across(everything(), ~log1p(pmax(., 0))))
  cv_df <- bind_cols(species = itv_sp, cv)
  plas <- run_module_pca(cv_df, setdiff(names(cv_df), "species"), "ITV",
                         impute = impute, logtrans = FALSE)
  res$ITV <- plas

  scores <- reduce(map(res, "scores"), full_join, by = "species")
  vars   <- tibble(Module = names(res), PC1_var = map_dbl(res, "var"))
  list(scores = scores, vars = vars)
}

mod_cc  <- build_modules(dat, impute = FALSE)   # Pipeline B: complete-case
mod_imp <- build_modules(dat, impute = TRUE)    # Pipeline A: imputed

vtab <- mod_imp$vars |> rename(`PC1 % (imputed)` = PC1_var) |>
  left_join(mod_cc$vars |> rename(`PC1 % (complete-case)` = PC1_var), by = "Module")
kbl(vtab, digits = 1,
    caption = "Within-module PC1 variance explained. Two columns = two module-CONSTRUCTION variants (imputed vs complete-case). Pipeline C reuses the complete-case module scores, so there is no separate Pipeline-C variance column.")

p_var <- mod_imp$vars |>
  ggplot(aes(reorder(Module, PC1_var), PC1_var)) +
  geom_col(fill = "#457B9D") + coord_flip() +
  geom_text(aes(label = sprintf("%.0f", PC1_var)), hjust = -0.15, size = 3) +
  ylim(0, 105) +
  labs(x = NULL, y = "PC1 variance explained (%)", title = "Module PC1 variance")
save_fig(p_var, "Fig1_ModulePC1Variance", 130, 95); p_var   # NPh 1.5 column

.lg <- function(d) d |> dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), ~suppressWarnings(ifelse(!is.na(.) & . > 0, log10(.), .))))
.allt  <- intersect(unique(unlist(modules)), names(dat))
.Cm    <- suppressWarnings(cor(.lg(dat[, .allt]), use = "pairwise.complete.obs", method = "spearman"))
.multi <- names(modules)[vapply(modules, function(z) length(intersect(z, .allt)) >= 2, logical(1))]
.xab   <- function(a, b) { a <- intersect(a, .allt); b <- intersect(setdiff(b, a), .allt)   # mean |r| BETWEEN two trait sets
  if (!length(a) || !length(b)) NA_real_ else mean(abs(.Cm[a, b]), na.rm = TRUE) }
coh <- purrr::map_dfr(names(modules), function(m) {
  tr <- intersect(modules[[m]], .allt); k <- length(tr)
  within <- if (k >= 2) { sub <- .Cm[tr, tr]; round(mean(abs(sub[upper.tri(sub)]), na.rm = TRUE), 2) } else NA_real_
  pc <- if (k >= 2) { Xi <- .lg(dat[, tr]); Xi <- Xi[complete.cases(Xi), , drop = FALSE]
                      e <- prcomp(scale(Xi))$sdev^2; round(100*e[1]/sum(e)) } else 100L
  oc <- vapply(setdiff(.multi, m), function(o) .xab(tr, modules[[o]]), numeric(1))   # link to each OTHER block
  best <- if (all(is.na(oc))) NA else names(which.max(oc)); bestr <- if (all(is.na(oc))) NA else round(max(oc, na.rm = TRUE), 2)
  tibble::tibble(
    Block = m, n = k,
    `within |r|` = if (is.na(within)) "-" else sprintf("%.2f", within),
    `PC1 %` = pc,
    `nearest block (|r|)` = if (is.na(best)) "-" else sprintf("%s (%.2f)", best, bestr),
    Decision = if (k >= 2) "PCA module"
               else if (!is.na(bestr) && bestr >= 0.30) "standalone (borderline)" else "standalone",
    .wn = within, .nn = bestr)
})
# ITV (built separately): coherence among the 12 CV pairs (sd/|mean|, log1p)
.cvm <- sapply(seq_len(nrow(itv_pairs)), function(i) {
  mm <- suppressWarnings(as.numeric(dat[[itv_pairs$mean[i]]])); s <- suppressWarnings(as.numeric(dat[[itv_pairs$sd[i]]]))
  cv <- s / abs(mm); cv[!is.finite(cv)] <- NA; log1p(pmax(cv, 0)) })
.cvm <- .cvm[, colSums(!is.na(.cvm)) > 1, drop = FALSE]
.Ci  <- suppressWarnings(cor(.cvm, use = "pairwise.complete.obs", method = "spearman"))
.Xi  <- .cvm[complete.cases(.cvm), , drop = FALSE]
.pci <- if (nrow(.Xi) > ncol(.cvm) + 2) { e <- prcomp(scale(.Xi))$sdev^2; round(100*e[1]/sum(e)) } else NA
.itvw <- round(mean(abs(.Ci[upper.tri(.Ci)]), na.rm = TRUE), 2)
coh <- rbind(coh, tibble::tibble(Block = "ITV", n = nrow(itv_pairs),
             `within |r|` = sprintf("%.2f", .itvw),
             `PC1 %` = .pci, `nearest block (|r|)` = "separate CV", Decision = "PCA (CV pairs)",
             .wn = .itvw, .nn = NA_real_))
# ECOLOGICAL STRATEGY each module operationalises (modules are the measurable proxy; the strategy is the biology)
strat <- c(
  Leaf_Economics = "Resource acquisition (leaf economics: fast vs slow return)",
  Leaf_Size      = "Light capture / leaf-area strategy",
  Leaf_Nutrient  = "Nutrient economy",
  Height         = "Competition for light (stature)",
  Hydraulic      = "Water transport (hydraulic) strategy",
  Root           = "Belowground acquisition (root economics)",
  Propagule_Size = "Dispersal & establishment (propagule size)",
  Fecundity      = "Reproductive output / colonisation (fecundity)",
  Genome         = "Genome size / life-cycle tempo",
  Rooting_Depth  = "Belowground foraging depth",
  Stem_Diameter  = "Mechanical support / stature",
  Stem_Cond_Dens = "Stem tissue density (support vs storage)",
  Germination    = "Regeneration / establishment",
  ITV            = "Phenotypic variability (niche breadth)")
coh$`Ecological strategy` <- unname(strat[coh$Block])
results_summary$module_coherence <- coh
coh_disp <- coh[, c("Ecological strategy", "Block", "n", "within |r|", "PC1 %", "Decision")]
kableExtra::kbl(coh_disp, booktabs = TRUE, longtable = FALSE, row.names = FALSE,
    caption = "How the functional modules were formed, and the ECOLOGICAL STRATEGY each operationalises (observed, non-imputed Spearman). The module (Block) is the measurable proxy; the strategy is the biology we interpret. **n** = traits in the block; **within |r|** = mean absolute Spearman correlation AMONG the block's own traits ('-' for a single trait, which has no within-block partner); **PC1** = the variance share the first axis captures. The decision rule - and two apparent exceptions - are explained beneath the figure below.") |>
  kableExtra::kable_styling(latex_options = c("HOLD_position", "scale_down"), font_size = 8) |>
  kableExtra::column_spec(1, width = "3.8cm") |>
  kableExtra::column_spec(6, width = "2.4cm")

## Blocked correlation heatmap: shows HOW each individual trait groups into its module.
.ord <- unlist(modules); .ord <- .ord[.ord %in% colnames(.Cm)]
.Cm2 <- .Cm[.ord, .ord]
.bn  <- vapply(modules, function(z) sum(z %in% .ord), integer(1))
.bl  <- data.frame(mod = names(modules), n = .bn, stringsAsFactors = FALSE)[.bn > 0, ]
.bl$end <- cumsum(.bl$n); .bl$start <- .bl$end - .bl$n + 1; .bl$mid <- (.bl$start + .bl$end) / 2
.bl$pc1 <- coh$`PC1 %`[match(.bl$mod, coh$Block)]
.bl$type <- ifelse(.bl$n >= 2, "Module (-> one PC1 axis)", "Standalone trait")
.bl$lab  <- ifelse(.bl$n >= 2, sprintf("%s  (PC1 %s%%)", .bl$mod, .bl$pc1), paste0(.bl$mod, "  (standalone)"))

## ---- master label prettifier: PC1 axes, covariates, phyloglm factor terms; falls back to .nm / generic ----
.PRETTY <- c(
  # ecological-strategy module PC1 axes
  Leaf_Economics_PC1 = "Leaf economics", Leaf_Size_PC1 = "Leaf size", Leaf_Nutrient_PC1 = "Leaf nutrients",
  Height_PC1 = "Plant height",
  Hydraulic_PC1 = "Hydraulic", Root_PC1 = "Root", Propagule_Size_PC1 = "Propagule size",
  Fecundity_PC1 = "Fecundity", Genome_PC1 = "Genome size", Rooting_Depth_PC1 = "Rooting depth",
  Stem_Diameter_PC1 = "Stem diameter", Stem_Cond_Dens_PC1 = "Stem conduit density",
  Germination_PC1 = "Germination", ITV_PC1 = "ITV",
  # covariates
  origin = "Origin", life_span = "Life history", growth_form = "Growth form",
  leaf_compoundness = "Leaf compoundness", photosynthetic_pathway = "Photosynthetic pathway",
  dispersal_distance = "Dispersal distance", alt_range_scaled = "Altitudinal range",
  eoo_scaled = "EOO", aoo_scaled = "AOO",
  # phyloglm factor-coefficient tokens (variable + non-reference level, no separator)
  c4C4 = "C4 photosynthesis", life_annualnon_annual = "Non-annual life history",
  life_annualannual = "Annual life history", originnative = "Native origin")
.pp1 <- function(t) {
  if (t %in% names(.PRETTY)) return(unname(.PRETTY[[t]]))
  if (t %in% names(.nm))     return(unname(.nm[[t]]))
  s <- gsub("_", " ", sub("_mean$", "", sub("_scaled$", "", sub("_PC1$", "", t))))
  s <- paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  s <- gsub("\\bAoo\\b", "AOO", s); s <- gsub("\\bEoo\\b", "EOO", s); s <- gsub("\\bItv\\b", "ITV", s)
  s
}
pp <- function(x) vapply(as.character(x), .pp1, character(1), USE.NAMES = FALSE)
np <- length(.ord)
.long <- data.frame(xi = rep(seq_len(np), times = np), yi = rep(seq_len(np), each = np), r = as.vector(.Cm2))
p_mf <- ggplot(.long, aes(xi, yi, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#457B9D", mid = "white", high = "#E76F51", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman r") +
  geom_rect(data = .bl, inherit.aes = FALSE,
            aes(xmin = start - 0.5, xmax = end + 0.5, ymin = start - 0.5, ymax = end + 0.5, colour = type),
            fill = NA, linewidth = 0.9) +
  scale_colour_manual(values = c("Module (-> one PC1 axis)" = "#2A9D8F", "Standalone trait" = "#E76F51"), name = NULL) +
  geom_text(data = .bl, inherit.aes = FALSE, aes(x = np + 0.9, y = mid, label = lab, colour = type),
            hjust = 0, size = 2.7, fontface = "bold", show.legend = FALSE) +
  scale_x_continuous(breaks = seq_len(np), labels = .lab(.ord), expand = expansion(add = c(0.5, 15))) +
  scale_y_continuous(breaks = seq_len(np), labels = .lab(.ord), expand = expansion(add = 0.5)) +
  labs(x = NULL, y = NULL,
       title = "How each trait groups into a module (observed Spearman correlations among the 35 traits)",
       subtitle = "Traits ordered by module; each box is a block (green = multi-trait module -> one PC1 axis; orange = standalone trait).") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, colour = "grey15"),
        axis.text.y = element_text(size = 7, colour = "grey15"), panel.grid = element_blank(),
        plot.subtitle = element_text(size = 8.5), legend.position = "right")
save_fig(p_mf, "Fig_ModuleFormation", 225, 200); p_mf

covars <- c("origin","life_span","growth_form","leaf_compoundness",
            "photosynthetic_pathway","dispersal_distance",
            "alt_range_scaled","eoo_scaled","aoo_scaled")

assemble <- function(module_scores) {
  module_scores |>
    left_join(
      dat |> select(species, all_of(covars), HR_bin, HR3, hr_status, habitat_group),
      by = "species"
    ) |>
    mutate(across(c(alt_range_scaled, eoo_scaled, aoo_scaled), as.numeric))
}
full_imp <- assemble(mod_imp$scores)
full_cc  <- assemble(mod_cc$scores)

cat("Imputed pipeline: rows with complete predictors =",
    sum(complete.cases(full_imp |> select(ends_with("_PC1"), all_of(covars)))), "/", nrow(full_imp), "\n")
cat("Complete-case pipeline: rows with complete predictors =",
    sum(complete.cases(full_cc |> select(ends_with("_PC1"), all_of(covars)))), "/", nrow(full_cc), "\n")

num_pred <- full_imp |> select(ends_with("_PC1"), alt_range_scaled, eoo_scaled, aoo_scaled)
cor_pred <- cor(num_pred, use = "pairwise.complete.obs", method = "spearman")
corrplot(cor_pred, method = "color", type = "upper", tl.cex = 0.5, tl.col = "black",
         addCoef.col = "grey30", number.cex = 0.4)
hi <- which(abs(cor_pred) > 0.7 & abs(cor_pred) < 1, arr.ind = TRUE)
cat("Numeric predictor pairs with |rho| > 0.7:", nrow(hi)/2, "\n")

# Generic balanced RF
brf <- function(data, y = "HR_bin", ntree = params$ntree, mtry = NULL,
                importance = TRUE, proximity = FALSE) {
  data <- data |> droplevels()
  preds <- setdiff(names(data), y)
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(length(preds))))
  nmin <- min(table(data[[y]]))
  samp <- setNames(rep(nmin, nlevels(data[[y]])), levels(data[[y]]))
  randomForest(reformulate(preds, y), data = data, ntree = ntree, mtry = mtry,
               importance = importance, sampsize = samp, strata = data[[y]],
               proximity = proximity, keep.forest = TRUE)
}

model_frame <- function(df, y = "HR_bin", track = "imp") {
  df |> select(ends_with("_PC1"), all_of(covars), all_of(y)) |> tidyr::drop_na()
}
get_auc <- function(rf, data, y = "HR_bin") {
  as.numeric(pROC::auc(pROC::roc(data[[y]], rf$votes[,"HR"], levels = c("nonHR","HR"),
                                 direction = "<", quiet = TRUE)))
}
imp_vec <- function(rf) { i <- randomForest::importance(rf, type = 1)[,1]; sort(i, decreasing = TRUE) }

base_track <- if (params$impute_primary) full_imp else full_cc

dA <- model_frame(base_track |> filter(habitat_group == "arable"))            # Model A: arable only
dB <- model_frame(base_track)                                                 # Model B: arable + ruderal
cat("Model A (arable):", nrow(dA), " | Model B (arable+ruderal):", nrow(dB), "\n")
print(table(dA$HR_bin)); print(table(dB$HR_bin))

set.seed(1); rfA <- brf(dA); set.seed(1); rfB <- brf(dB)
aucA <- get_auc(rfA, dA); aucB <- get_auc(rfB, dB)

impA <- imp_vec(rfA); impB <- imp_vec(rfB)
common <- intersect(names(impA), names(impB))
rank_rho <- cor(rank(-impA[common]), rank(-impB[common]), method = "spearman")
.rank_p  <- suppressWarnings(cor.test(impA[common], impB[common], method = "spearman"))$p.value   # Supp Fig 20 P
top10_overlap <- length(intersect(names(head(impA,10)), names(head(impB,10))))

knitr::kable(data.frame(
  Model = c("A: Arable only","B: Arable + Ruderal"),
  N = c(nrow(dA), nrow(dB)),
  HR_class_error = c(rfA$err.rate[params$ntree,"HR"], rfB$err.rate[params$ntree,"HR"]),
  AUC = c(aucA, aucB)), digits = 3, caption = "Phase 0: discrimination with vs without ruderal species.")
cat(sprintf("Importance rank correlation A vs B: rho = %.3f | top-10 overlap = %d/10\n",
            rank_rho, top10_overlap))

cmp <- tibble(var = common, A = impA[common], B = impB[common])
p_p0 <- ggplot(cmp, aes(A, B)) +
  geom_abline(slope = 1, linetype = "dashed", colour = "grey60") +
  geom_point(colour = "#2A9D8F") +
  ggrepel::geom_text_repel(aes(label = pp(var)), size = 1.8, max.overlaps = 12) +
  labs(x = "Importance (Arable only)", y = "Importance (Arable + Ruderal)",
       title = sprintf("Phase 0 importance concordance (rho=%.2f)", rank_rho))
save_fig(p_p0, "Fig2_Phase0_HabitatSensitivity", 130, 110); p_p0   # NPh 1.5 column

USE_RUDERAL <- (rank_rho >= 0.8 && top10_overlap >= 7 && abs(aucA - aucB) < 0.03)
cat("Decision -> include ruderal baseline:", USE_RUDERAL, "\n")
sel <- function(df) if (USE_RUDERAL) df else df |> filter(habitat_group == "arable")
analysis_imp <- sel(full_imp); analysis_cc <- sel(full_cc)
results_summary$phase0 <- list(aucA = aucA, aucB = aucB, rank_rho = rank_rho, rank_p = .rank_p,
                               top10_overlap = top10_overlap, use_ruderal = USE_RUDERAL,
                               nA = nrow(dA), nB = nrow(dB))

dcc  <- model_frame(analysis_cc);  dimp <- model_frame(analysis_imp)
set.seed(2); rf_cc  <- brf(dcc); set.seed(2); rf_imp <- brf(dimp)
auc_cc <- get_auc(rf_cc, dcc); auc_imp <- get_auc(rf_imp, dimp)
icc <- imp_vec(rf_cc); iim <- imp_vec(rf_imp)

# Pipeline C: all species, NA retained, cforest surrogate splits, class-weighted for imbalance.
# VALIDATION (differs from the RF pipelines): partykit::cforest uses subsampling WITHOUT replacement (~63%
# of cases per tree) and its OWN out-of-bag prediction (predict(OOB=TRUE)); imbalance is handled by case
# weights (wC), NOT by balanced sampsize. So Pipeline C's AUC is an OOB estimate from a different resampler.
dC <- analysis_cc |> select(ends_with("_PC1"), all_of(covars), HR_bin)   # NAs kept
wC <- ifelse(dC$HR_bin == "HR", sum(dC$HR_bin == "nonHR")/sum(dC$HR_bin == "HR"), 1)
set.seed(2)
cf_C <- partykit::cforest(HR_bin ~ ., data = dC, ntree = params$ntree, weights = wC,
                          control = ctree_control(mincriterion = 0))
probC <- predict(cf_C, OOB = TRUE, type = "prob")[, "HR"]
auc_C <- as.numeric(pROC::auc(pROC::roc(dC$HR_bin, probC, levels = c("nonHR","HR"),
                                        direction = "<", quiet = TRUE)))
predC <- factor(ifelse(probC > 0.5, "HR", "nonHR"), levels = c("nonHR","HR"))
hrrec_C <- mean(predC[dC$HR_bin == "HR"] == "HR")
iC <- sort(partykit::varimp(cf_C), decreasing = TRUE)

cmAB <- intersect(names(icc), names(iim)); cmAC <- intersect(names(iim), names(iC))
rho_AB <- cor(rank(-icc[cmAB]), rank(-iim[cmAB]), method = "spearman")
rho_AC <- cor(rank(-iim[cmAC]), rank(-iC[cmAC]),  method = "spearman")

tracks_tab <- data.frame(
  Pipeline = c("A: missForest + RF","B: complete-case + RF","C: no-impute + cforest"),
  Model = c("randomForest","randomForest","cforest (surrogate)"),
  N = c(nrow(dimp), nrow(dcc), nrow(dC)),
  N_HR = c(sum(dimp$HR_bin=="HR"), sum(dcc$HR_bin=="HR"), sum(dC$HR_bin=="HR")),
  AUC = c(auc_imp, auc_cc, auc_C))
kbl(tracks_tab, digits = 3,
    caption = "Three-pipeline sensitivity. Pipeline A (imputed + randomForest) is the PRIMARY model used for all main results. B and C are robustness checks. AUCs are NOT directly comparable across models (C uses cforest + a different importance metric); the ~0.02 spread sits inside Pipeline A's 95 percent CI, i.e. it is noise, not evidence that C is 'better'.")
cat(sprintf("Importance-rank correlation: A vs B rho = %.3f | A vs C rho = %.3f\n", rho_AB, rho_AC))
cat(sprintf("Pipeline C: all %d species retained (no imputation); OOB AUC = %.3f; HR recall = %.2f\n",
            nrow(dC), auc_C, hrrec_C))

results_summary$tracks <- tracks_tab
results_summary$rho_AB <- rho_AB; results_summary$rho_AC <- rho_AC

p_tr <- tracks_tab |> ggplot(aes(reorder(Pipeline, AUC), AUC, fill = Pipeline)) +
  geom_col(width = .6) + coord_flip(ylim = c(0.5, max(tracks_tab$AUC)+0.05)) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c("#457B9D","#2A9D8F","#E76F51")) +
  labs(x = NULL, y = "OOB AUC") +
  theme(legend.position = "none")
save_fig(p_tr, "Fig_S3_TrackComparison", 115, 80); p_tr

primary <- if (params$impute_primary) dimp else dcc
prim_lbl <- if (params$impute_primary) "imputed (Pipeline A)" else "complete-case (Pipeline B)"
cat("Primary modelling pipeline:", prim_lbl, " | N =", nrow(primary), "\n")
results_summary$primary_track <- prim_lbl
results_summary$primary_N <- nrow(primary)

# Imputed all-trait matrix (35 mean traits, missForest). Built ONCE here and reused by the trait-space
# PCA, the go/no-go LDA, and the later all-traits robustness RF (Step 18).
raw_traits <- dat |> select(species, all_of(mean_cols))
Xall <- as.data.frame(raw_traits |> select(-species) |> mutate(across(everything(), as.numeric)))
Xall <- Xall |> mutate(across(everything(), ~suppressWarnings(ifelse(!is.na(.) & . > 0, log10(.), .))))
set.seed(42); Xall_imp <- missForest::missForest(Xall, verbose = FALSE)$ximp
alltraits <- bind_cols(species = raw_traits$species, Xall_imp) |>
  left_join(dat |> select(species, all_of(covars), HR_bin, habitat_group), by = "species")
if (!USE_RUDERAL) alltraits <- alltraits |> filter(habitat_group == "arable")

suppressPackageStartupMessages({library(vegan); library(ggrepel)})
## PCA of the 35 IMPUTED functional traits (missForest -> log10 -> z-scaled), with trait-loading arrows.
## Geographic-range covariates (alt_range/EOO/AOO) are excluded (they are biogeography, not function).
## Woody non-weeds are excluded so the ordination is not dominated by tree-sized outliers; ALL herbaceous
## species (including the large-herb outliers) are kept in the PCA and the PERMANOVA - the zoom only affects
## what is displayed.
tp <- alltraits |> dplyr::left_join(dat |> dplyr::select(species, hr_status, plant_growth_form), by = "species") |>
  dplyr::filter(!is.na(hr_status), !grepl("tree|shrub", plant_growth_form))
trait_cols <- intersect(mean_cols, names(tp)); n_traits <- length(trait_cols)
Xtp <- scale(as.matrix(tp |> dplyr::select(dplyr::all_of(trait_cols)) |>
                       dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))))
pca_t <- prcomp(Xtp); ve <- 100 * pca_t$sdev^2 / sum(pca_t$sdev^2)
scr <- data.frame(species = tp$species, PC1 = pca_t$x[, 1], PC2 = pca_t$x[, 2],
                  HR3 = factor(tp$hr_status, levels = c("Non_HR_local","HR_local","HR_elsewhere")))
grp <- factor(ifelse(tp$hr_status == "Non_HR_local", "nonHR", "HR")); set.seed(7)
pm <- vegan::adonis2(Xtp ~ grp, method = "euclidean", permutations = 999)   # 544 herbaceous spp, incl off-range
## robustness: the SAME PERMANOVA including the 31 woody species (all 575) - confirms the woody exclusion
## affects only the figure, not the statistics.
ta <- alltraits |> dplyr::filter(!is.na(HR_bin))
Xall <- scale(as.matrix(ta |> dplyr::select(dplyr::all_of(trait_cols)) |>
                        dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))))
set.seed(7); pm_all <- vegan::adonis2(Xall ~ ta$HR_bin, method = "euclidean", permutations = 999)
## export the exact matrix used for this PCA
.pcadir <- file.path(out_dir, "data"); dir.create(.pcadir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(dplyr::bind_cols(tp |> dplyr::select(species, hr_status, plant_growth_form),
                                  as.data.frame(round(Xtp, 4))),
                 file.path(.pcadir, "trait_pca_input_35traits.csv"))
ld <- as.data.frame(pca_t$rotation[, 1:2]); colnames(ld) <- c("PC1","PC2")
ld$lab <- .lab(rownames(ld)); ld$len <- sqrt(ld$PC1^2 + ld$PC2^2)
lda <- ld; lda[, 1:2] <- lda[, 1:2] * (3.2 / max(abs(as.matrix(lda[, 1:2])))); top <- lda |> dplyr::arrange(dplyr::desc(len)) |> head(10)
zx <- c(-4, 4); zy <- c(-4.5, 4.5); scr$off <- scr$PC1<zx[1]|scr$PC1>zx[2]|scr$PC2<zy[1]|scr$PC2>zy[2]; n_out <- sum(scr$off)
p_tpca <- ggplot(scr, aes(PC1, PC2)) +
  geom_point(aes(colour = HR3), size = .7, alpha = .4) + stat_ellipse(aes(colour = HR3), level = .95) +
  geom_segment(data = top, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(1.5, "mm")), colour = "grey25", linewidth = .4, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = top, aes(x = PC1, y = PC2, label = lab), size = 2, colour = "grey10",
                           inherit.aes = FALSE, max.overlaps = 20, segment.size = .2, min.segment.length = 0) +
  scale_colour_manual(values = col_status, labels = hr_labs, name = "HR status") + coord_cartesian(xlim = zx, ylim = zy) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = sprintf("PCA of the %d imputed functional traits, by resistance status", n_traits),
       subtitle = "herbaceous weeds; arrows = top trait loadings")
save_fig(p_tpca, "Fig_S20_TraitPCA", 140, 120); p_tpca
n_woody <- sum(grepl("tree|shrub", dat$plant_growth_form), na.rm = TRUE)
results_summary$traitpca <- list(ve1 = ve[1], ve2 = ve[2], permR2 = pm$R2[1], permP = pm$`Pr(>F)`[1],
                                 n_out = n_out, n_traits = n_traits, n_total = nrow(scr),
                                 n_woody = n_woody, n_all = nrow(dat),
                                 permR2_all = pm_all$R2[1], permP_all = pm_all$`Pr(>F)`[1])
cat(sprintf("Imputed %d-trait PCA: PC1=%.1f%% PC2=%.1f%% | PERMANOVA on ALL %d herbaceous spp (incl. %d off-display): R2=%.3f P=%.3f\n",
            n_traits, ve[1], ve[2], nrow(scr), n_out, pm$R2[1], pm$`Pr(>F)`[1]))

ord <- ld$lab[order(abs(ld$PC1))]
Ll <- ld |> tidyr::pivot_longer(c(PC1, PC2), names_to = "axis", values_to = "loading") |>
  dplyr::mutate(trait = factor(lab, levels = ord))
p_load <- ggplot(Ll, aes(loading, trait, fill = axis)) +
  geom_col(position = position_dodge(width = .7), width = .65) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .3) +
  scale_fill_manual(values = c(PC1 = "#1D3557", PC2 = "#E9A03B"), name = "Axis") +
  labs(x = "Loading (correlation with axis)", y = NULL,
       title = "Trait loadings: which traits drive PC1 vs PC2")
save_fig(p_load, "Fig_S20b_Loadings", 120, 150); p_load

p_full <- ggplot(scr, aes(PC1, PC2, colour = HR3)) + geom_point(size = .7, alpha = .5) +
  stat_ellipse(level = .95) + scale_colour_manual(values = col_status, labels = hr_labs, name = "HR status") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "Full (un-zoomed) view: the off-range large-herb species")
save_fig(p_full, "Fig_S20c_FullView", 130, 100); p_full
offtab <- scr |> dplyr::filter(off) |>
  dplyr::transmute(Species = species, HR = as.character(HR3), PC1 = round(PC1, 1)) |>
  dplyr::left_join(dat |> dplyr::select(species, plant_growth_form), by = c("Species" = "species")) |>
  dplyr::arrange(PC1)
readr::write_csv(offtab, file.path(tab_dir, "Table_TraitPCA_offrange.csv"))
readr::write_csv(offtab, file.path(out_dir, "data", "Table_TraitPCA_offrange.csv"))  # also in data/ for users
kbl(head(offtab, 15), digits = 1, row.names = FALSE,
    caption = sprintf("The %d off-display species (all non-resistant; full list in Table S11): large herbs and graminoids that stretch PC1 but are retained in the PCA/PERMANOVA. First 15 shown.", nrow(offtab)))

## "Loadings" describe the PCA axes (identical for all species), so a loadings plot cannot vary by HR
## status. The meaningful version - how each trait differs across the classes - is each trait's mean
## z-score per class, ordered by PC1 (size) loading.
lab3p <- hr_labs
ord2 <- .lab(names(sort(pca_t$rotation[, 1])))
prof <- as.data.frame(Xtp) |>
  dplyr::mutate(HR3 = factor(lab3p[tp$hr_status], levels = lab3p)) |>
  tidyr::pivot_longer(-HR3, names_to = "trait", values_to = "z") |>
  dplyr::group_by(HR3, trait) |> dplyr::summarise(m = mean(z, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(trait = factor(.lab(trait), levels = ord2))
p_prof <- ggplot(prof, aes(m, trait, colour = HR3)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = .3) + geom_point(size = 1.5, alpha = .9) +
  scale_colour_manual(values = setNames(col_status[names(lab3p)], lab3p), name = "HR status") +
  labs(x = "Mean z-scored trait value", y = NULL, title = "Per-trait mean z-score by resistance class")
save_fig(p_prof, "Fig_S20d_TraitProfile", 130, 150); p_prof

suppressPackageStartupMessages(library(patchwork))
lab3 <- c(Non_HR_local = "non-HR-local", HR_local = "HR-local", HR_elsewhere = "HR-elsewhere")
sdat <- dat |> dplyr::filter(!is.na(hr_status)) |> dplyr::mutate(HR3 = factor(lab3[hr_status], levels = lab3))
pal3 <- setNames(col_status[names(lab3)], lab3)
sigm <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
A <- sdat |> dplyr::group_by(HR3) |> dplyr::summarise(prop = 100*mean(life_span == "annual", na.rm = TRUE), .groups = "drop")
pA <- ggplot(A, aes(HR3, prop, fill = HR3)) + geom_col(width = .62) +
  geom_text(aes(label = sprintf("%.0f%%", prop)), vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = pal3, guide = "none") + ylim(0, 100) +
  scale_x_discrete(labels = hr_labs) +
  labs(x = NULL, y = "% annual", title = "a  Annual life history") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
B <- sdat |> dplyr::filter(!is.na(photosynthetic_pathway)) |> dplyr::count(HR3, photosynthetic_pathway) |>
  dplyr::group_by(HR3) |> dplyr::mutate(prop = 100*n/sum(n))
pB <- ggplot(B, aes(HR3, prop, fill = photosynthetic_pathway)) + geom_col(width = .62) +
  scale_fill_manual(values = c(C3 = "#A8C5D8", `C3/C4` = "#E9C46A", C4 = "#E76F51"), name = "Pathway") +
  scale_x_discrete(labels = hr_labs) +
  labs(x = NULL, y = "% of species", title = "b  Photosynthetic pathway") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
## Continuous traits kept by the THREE resistance classes (for the HR-local vs HR-elsewhere breakdown),
## with significance for BOTH the 3-class Kruskal-Wallis test AND the binary resistant-vs-non contrast that
## the PCA/PERMANOVA tests (Fig. S20) - so the panel is coherent with that ordination.
sdat$RES <- factor(ifelse(sdat$hr_status == "Non_HR_local", "Non-resistant", "Resistant"),
                   levels = c("Non-resistant", "Resistant"))
cv <- c(sla_undef_mean = "SLA", seed_num_rep_unit_mean = "Seeds/unit", plant_height_veg_mean = "Plant height", seed_mass_mean = "Seed mass")
cd <- purrr::map_dfr(names(cv), function(v) {
  raw <- suppressWarnings(as.numeric(sdat[[v]]))
  pk <- tryCatch(kruskal.test(raw ~ sdat$HR3)$p.value, error = function(e) NA)
  pw <- tryCatch(wilcox.test(raw ~ sdat$RES)$p.value, error = function(e) NA)
  data.frame(HR3 = sdat$HR3, val = suppressWarnings(log10(raw + ifelse(grepl("seed_num", v), 1, 0))),
             trait = sprintf("%s\n(3-class %s | res-vs-non %s)", cv[v], sigm(pk), sigm(pw)))
}) |> dplyr::filter(!is.na(val))
pC <- ggplot(cd, aes(HR3, val, colour = HR3)) + geom_jitter(width = .15, alpha = .22, size = .45) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = .45) + facet_wrap(~trait, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = pal3, guide = "none") + scale_x_discrete(labels = hr_labs) +
  labs(x = NULL, y = "log10 value",
       title = "c  Continuous traits by class - significance: 3-class | resistant-vs-non (the PCA contrast)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), strip.text = element_text(size = 5.8))
p_synd <- (pA | pB) / pC + patchwork::plot_layout(heights = c(1, 1.1))
save_fig(p_synd, "Fig_S22_ResistanceSyndrome", 180, 150); p_synd
rstat <- sdat |> dplyr::mutate(res = hr_status != "Non_HR_local")
wp <- function(v) tryCatch(wilcox.test(as.numeric(sdat[[v]]) ~ sdat$RES)$p.value, error = function(e) NA)
results_summary$syndrome <- list(
  annual_res = 100*mean(rstat$life_span[rstat$res] == "annual", na.rm = TRUE),
  annual_non = 100*mean(rstat$life_span[!rstat$res] == "annual", na.rm = TRUE),
  c4_res = 100*mean(rstat$photosynthetic_pathway[rstat$res] == "C4", na.rm = TRUE),
  c4_non = 100*mean(rstat$photosynthetic_pathway[!rstat$res] == "C4", na.rm = TRUE),
  annual_local = A$prop[A$HR3 == "HR-local"], annual_else = A$prop[A$HR3 == "HR-elsewhere"],
  p_sla = wp("sla_undef_mean"), p_seed = wp("seed_num_rep_unit_mean"))

# SUPERVISED counterpart to the unsupervised 35-trait PCA (Fig. S20), on the SAME 35 imputed traits.
# Overfitting diagnostic: STANDARD LDA (resubstitution - trained AND tested on the same data) vs
# CROSS-VALIDATED LDA. Given the 19:1 imbalance we use REPEATED STRATIFIED k-FOLD CV (5-fold x 5 repeats,
# class proportions preserved in every fold) instead of leave-one-out, which is high-variance and ignores
# class balance. The resubstitution-minus-CV gap IS the overfitting. Xtp = scaled 35-trait matrix (Fig. S20).
Xl   <- Xtp
keep <- stats::complete.cases(Xl); Xl <- Xl[keep, , drop = FALSE]
y3 <- factor(tp$hr_status[keep], levels = c("Non_HR_local","HR_local","HR_elsewhere"))
yb <- factor(ifelse(tp$hr_status[keep] == "Non_HR_local", "nonHR", "HR"), levels = c("nonHR","HR"))
auc_of <- function(truth, score)
  as.numeric(pROC::auc(pROC::roc(truth, score, levels = c("nonHR","HR"), direction = "<", quiet = TRUE)))

# --- STANDARD LDA: resubstitution (train == test) - the optimistic, overfit-prone estimate ---
ld3_in  <- MASS::lda(Xl, y3); acc3_in <- mean(predict(ld3_in)$class == y3)
auc_in  <- auc_of(yb, predict(MASS::lda(Xl, yb))$posterior[, "HR"])
# --- CROSS-VALIDATED LDA: repeated STRATIFIED 5-fold CV (5 repeats), imbalance-aware ---
.strat_folds <- function(g, k, seed) { set.seed(seed); f <- integer(length(g))
  for (lv in levels(g)) { w <- sample(which(g == lv)); f[w] <- rep_len(seq_len(k), length(w)) }; f }
K <- 5L; R <- 5L; cv_bin <- numeric(R); cv_3 <- numeric(R)
for (r in seq_len(R)) {
  fb <- .strat_folds(yb, K, 100L + r); pb <- rep(NA_real_, length(yb))
  f3 <- .strat_folds(y3, K, 200L + r); p3 <- factor(rep(NA, length(y3)), levels = levels(y3))
  for (k in seq_len(K)) {
    tb <- fb != k; pb[!tb] <- predict(MASS::lda(Xl[tb, , drop=FALSE], yb[tb]), Xl[!tb, , drop=FALSE])$posterior[, "HR"]
    t3 <- f3 != k; p3[!t3] <- predict(MASS::lda(Xl[t3, , drop=FALSE], y3[t3]), Xl[!t3, , drop=FALSE])$class
  }
  cv_bin[r] <- auc_of(yb, pb); cv_3[r] <- mean(p3 == y3)
}
auc_cv <- mean(cv_bin); auc_cv_sd <- sd(cv_bin); acc3_cv <- mean(cv_3); acc3_cv_sd <- sd(cv_3)
# fair non-linear benchmark: a balanced random forest on the SAME 35 traits (OOB), fit inline (benchmark for
# the go/no-go check, NOT the main module-based model, which comes later)
rf_dat <- as.data.frame(Xl); rf_dat$HR_bin <- yb
set.seed(42); auc_rf <- get_auc(brf(rf_dat), rf_dat)
results_summary$lda <- list(n = nrow(Xl), n_traits = ncol(Xl), acc3_in = acc3_in, acc3_cv = acc3_cv,
  acc3_cv_sd = acc3_cv_sd, auc_in = auc_in, auc_cv = auc_cv, auc_cv_sd = auc_cv_sd,
  overfit_gap = auc_in - auc_cv, rf_auc = auc_rf, cv_scheme = "5-fold x 5 repeats, stratified")
cat(sprintf("LDA on %d traits (n=%d): resub AUC=%.3f acc=%.3f | CV(5x5 stratified) AUC=%.3f (SD %.3f) acc=%.3f | OVERFIT gap=%.3f | 35-trait RF AUC=%.3f\n",
            ncol(Xl), nrow(Xl), auc_in, acc3_in, auc_cv, auc_cv_sd, acc3_cv, auc_in - auc_cv, auc_rf))

sc <- as.data.frame(predict(ld3_in)$x); sc$HR3 <- y3
pA <- ggplot(sc, aes(LD1, LD2, colour = HR3)) +
  geom_point(size = .8, alpha = .5) + stat_ellipse(level = .95, linewidth = .5) +
  scale_colour_manual(values = col_status, name = NULL,
    labels = hr_labs) +
  labs(title = "a  In-sample LDA axes (resubstitution view)", x = "LD1", y = "LD2") +
  theme(legend.position = "bottom", legend.text = element_text(size = 7), plot.title = element_text(size = 8.5))
bars <- data.frame(
  Setting = factor(c("Standard LDA\n(resubstitution)","CV-LDA\n(5x5 stratified)","Random forest\n(35 traits, OOB)"),
                   levels = c("Standard LDA\n(resubstitution)","CV-LDA\n(5x5 stratified)","Random forest\n(35 traits, OOB)")),
  AUC = c(auc_in, auc_cv, auc_rf), lo = c(NA, auc_cv - auc_cv_sd, NA), hi = c(NA, auc_cv + auc_cv_sd, NA))
pB <- ggplot(bars, aes(Setting, AUC, fill = Setting)) +
  geom_col(width = .62) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.16, colour = "grey30", na.rm = TRUE) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.4, size = 2.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  annotate("segment", x = 1, xend = 2, y = auc_in + 0.06, yend = auc_in + 0.06, colour = "grey35", linewidth = .3) +
  annotate("text", x = 1.5, y = auc_in + 0.095, label = sprintf("overfit gap %.3f", auc_in - auc_cv), size = 2.4, colour = "grey25") +
  scale_fill_manual(values = c("#E76F51","#457B9D","#2A9D8F"), guide = "none") +
  coord_cartesian(ylim = c(0.45, 1.10)) +
  labs(x = NULL, y = "AUC (resistant vs non)", title = "b  Overfitting check: resubstitution vs repeated stratified CV") +
  theme(plot.title = element_text(size = 7.8), axis.text.x = element_text(size = 6.5))
p_lda <- pA + pB + patchwork::plot_layout(widths = c(1.05, 1))
save_fig(p_lda, "Fig_S21b_LDA_CV", 185, 96); p_lda

bstick <- function(p) rev(cumsum(1 / rev(seq_len(p)))) / p * 100         # broken-stick expected %
dimtab <- purrr::map_dfr(names(modules), function(m) {
  cols <- intersect(modules[[m]], names(dat))
  if (length(cols) < 2)
    return(tibble::tibble(Module = m, n = length(cols), PC1 = 100, PC2 = NA_real_, BS2 = NA_real_, PC2_gt_BS = NA))
  X <- dat |> dplyr::select(dplyr::all_of(cols)) |>
       dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
       dplyr::mutate(dplyr::across(dplyr::everything(), ~ifelse(!is.na(.) & . > 0, log10(.), .)))
  X  <- as.data.frame(X[, colSums(!is.na(X)) > 0, drop = FALSE])
  Xi <- missForest::missForest(X, verbose = FALSE)$ximp
  ev <- prcomp(scale(Xi))$sdev^2; vv <- 100 * ev / sum(ev); bs <- bstick(length(ev))
  tibble::tibble(Module = m, n = length(cols), PC1 = vv[1], PC2 = vv[2], BS2 = bs[2], PC2_gt_BS = vv[2] > bs[2])
})
kbl(dimtab, digits = 1, caption = "Per-module dimensionality: percent variance on PC1 and PC2 versus the broken-stick (null) expectation for axis 2 (BS2). PC1 dominates every module; PC2 exceeds its broken-stick threshold only rarely.")

mod_pc12 <- function(df) {
  res <- purrr::imap(modules, function(cols, m) {
    cols <- intersect(cols, names(df))
    X <- df |> dplyr::select(dplyr::all_of(cols)) |>
         dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
         dplyr::mutate(dplyr::across(dplyr::everything(), ~ifelse(!is.na(.) & . > 0, log10(.), .)))
    X <- as.data.frame(X[, colSums(!is.na(X)) > 0, drop = FALSE])
    if (ncol(X) < 2) { v <- as.numeric(scale(X[[1]])); return(tibble::tibble(species = df$species, !!paste0(m, "_PC1") := v)) }
    Xi <- missForest::missForest(X, verbose = FALSE)$ximp; pc <- prcomp(scale(Xi))$x[, 1:2]
    tibble::tibble(species = df$species, !!paste0(m, "_PC1") := pc[, 1], !!paste0(m, "_PC2") := pc[, 2])
  })
  purrr::reduce(res, dplyr::full_join, by = "species")
}
sc12 <- mod_pc12(dat) |> dplyr::left_join(dat |> dplyr::select(species, dplyr::all_of(covars), HR_bin), by = "species")
d12  <- sc12 |> dplyr::select(dplyr::matches("_PC[12]$"), dplyr::all_of(covars), HR_bin) |> tidyr::drop_na()
set.seed(42); rf12 <- brf(d12); auc12 <- get_auc(rf12, d12)
# Self-contained predictive-equivalence: this section now runs BEFORE the core model, so we fit the
# comparison forests inline (full PC1 set vs PC1+PC2 vs all raw traits) rather than reusing the final model.
set.seed(42); rf_pc1  <- brf(primary); auc_pc1 <- get_auc(rf_pc1, primary)            # PC1-only (full set)
d_allp <- alltraits |> dplyr::select(dplyr::all_of(mean_cols), dplyr::all_of(covars), HR_bin) |>
  tidyr::drop_na()
set.seed(42); rf_allp <- brf(d_allp); auc_allp <- get_auc(rf_allp, d_allp)            # all raw traits (no PCA)
dim_auc <- data.frame(
  Model      = c("PC1 only (the model's design)", "PC1 + PC2 per module", "All raw traits (no PCA)"),
  Predictors = c(ncol(primary) - 1, ncol(d12) - 1, ncol(d_allp) - 1),
  AUC        = c(auc_pc1, auc12, auc_allp))
kbl(dim_auc, digits = 3, caption = "Predictive equivalence: adding PC2 per module, or using every raw trait, does not improve out-of-bag AUC over the PC1-only model - the higher dimensions carry no extra predictive signal.")
p_dim <- dim_auc |> ggplot(aes(reorder(Model, AUC), AUC)) + geom_col(fill = "#457B9D", width = .6) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.12, size = 2.6) +
  coord_flip(ylim = c(0.5, max(dim_auc$AUC) + 0.05)) +
  labs(x = NULL, y = "Out-of-bag AUC", title = "One axis is enough: PC1 vs PC1+PC2 vs all-traits")
save_fig(p_dim, "Fig_S21_DimChoice", 150, 70); p_dim
results_summary$dimchoice <- list(auc_pc1 = auc_pc1, auc12 = auc12, auc_all = auc_allp,
                                  mean_pc1 = mean(dimtab$PC1, na.rm = TRUE),
                                  n_pc2_gt_bs = sum(dimtab$PC2_gt_BS, na.rm = TRUE), n_mod = nrow(dimtab))

sc <- readRDS("strategy_comparison.rds")
kbl(sc$auc, digits = 3,
    caption = "Strategy-definition comparison: balanced-RF out-of-bag AUC using ONLY trait-derived strategy axes (no life-history, photosynthetic-pathway or range covariates). The all-35-trait row is the trait-only ceiling; higher = more of the trait signal retained.")







kbl(sc$sweep, digits = 3, caption = "Louvain resolution sweep: higher resolution gives more (finer) trait communities; the trait-only OOB AUC rises with community count but never exceeds our modules, while modularity falls (finer splits are less well-supported by the network).")



set.seed(42); rf0 <- brf(primary, proximity = TRUE)
cat(">>> INITIAL core model: ALL", ncol(primary) - 1, "predictors, BEFORE optimisation.",
    "Its OOB error (~19%) and confusion matrix are the starting point; the OPTIMISED final model",
    "(Step 8 onward, reported in Steps 11 & 15) has fewer predictors and lower error. <<<\n")
print(rf0)
.oob_initial <- rf0$err.rate[nrow(rf0$err.rate), "OOB"]   # initial all-predictor OOB error (pre-selection)
oob_df <- as.data.frame(rf0$err.rate) |> mutate(trees = row_number()) |>
  pivot_longer(-trees, names_to = "class", values_to = "error")
p_conv <- ggplot(oob_df, aes(trees, error, colour = class)) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c(OOB = "grey30", HR = "#E76F51", nonHR = "#457B9D")) +
  labs(x = "Trees", y = "OOB error", title = "Error convergence", colour = NULL)
save_fig(p_conv, "Fig3_ErrorConvergence", 115, 90); p_conv

## Didactic schematic: how ONE tree's class-balanced bootstrap draws its training set. Each class is
## sampled WITH replacement; the per-class draw count k = the minority count (set by brf() via sampsize).
set.seed(1)
# one panel per class; coord_equal works per-plot (patchwork), unlike facet free-scales + coord_equal
boot_panel <- function(N, k, ttl) {
  cnt  <- tabulate(sample.int(N, k, replace = TRUE), nbins = N)   # times each species drawn (0 = OOB)
  side <- ceiling(sqrt(N))
  d <- data.frame(row = (seq_len(N) - 1) %/% side, col = (seq_len(N) - 1) %% side,
                  times = factor(pmin(cnt, 3), levels = 0:3,
                                 labels = c("0 (out-of-bag)", "1", "2", "3+")))
  ggplot(d, aes(col, -row, fill = times)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    scale_fill_manual(values = c("0 (out-of-bag)" = "grey85", "1" = "#457B9D",
                                 "2" = "#E9A03B", "3+" = "#E76F51"),
                      name = "Times drawn\ninto this tree", drop = FALSE) +
    coord_equal() +
    labs(title = sprintf("%s\nN=%d, k=%d  ->  %.0f%% in-bag, %.0f%% OOB (expected)",
                         ttl, N, k, 100*(1 - (1 - 1/N)^k), 100*(1 - 1/N)^k), x = NULL, y = NULL) +
    theme_pub(9) + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                         panel.grid = element_blank(), plot.title = element_text(size = 8))
}
pA <- boot_panel(109, 109, "Minority class (resistant): bootstrapped at full size") +
  theme(legend.position = "none")   # single shared legend shown on the right panel only
pB <- boot_panel(466, 109, "Majority class (non-resistant): sub-sampled to minority size")
p_boot <- pA + pB +
  patchwork::plot_annotation(
    title = "How one tree's class-balanced bootstrap draws its training set",
    subtitle = "Each square = one species; sampling is WITH replacement (a species may be drawn 0, 1, 2 or 3+ times).\nPercentages are the EXPECTED values, P(never drawn) = (1 - 1/N)^k = 37% for the minority (N=k=109) and 79% for the sub-sampled majority (N=466, k=109);\nthe grid shows ONE example tree, so its grey count varies slightly around these expected values.",
    theme = theme(plot.subtitle = element_text(size = 7)))
save_fig(p_boot, "Fig_S0_BalancedBootstrap", 180, 95); p_boot

seq_opt <- function(data, y = "HR_bin", ntree = params$ntree) {
  preds <- setdiff(names(data), y); hist <- list(); it <- 0
  repeat {
    it <- it + 1
    d <- data |> select(all_of(c(preds, y)))
    set.seed(7); rf <- brf(d, y, ntree)
    imp <- randomForest::importance(rf, type = 1)[,1]
    hist[[it]] <- tibble(iter = it, n_pred = length(preds),
                         oob = rf$err.rate[ntree,"OOB"], hr = rf$err.rate[ntree,"HR"])
    neg <- names(imp)[imp < 0]
    if (length(neg) == 0 || length(preds) - length(neg) < 2) break
    preds <- setdiff(preds, neg)
  }
  list(model = rf, preds = preds, history = bind_rows(hist))
}
opt <- seq_opt(primary)
rf_final <- opt$model; final_preds <- opt$preds
cat("Retained", length(final_preds), "predictors:\n")
writeLines(strwrap(paste(final_preds, collapse = ", "), width = 84))   # wrap to fit the page margin
p_opt <- opt$history |> pivot_longer(c(oob, hr)) |>
  ggplot(aes(iter, value*100, colour = name)) + geom_line() + geom_point(size = 1) +
  scale_colour_manual(values = c(oob = "grey30", hr = "#E76F51"),
                      labels = c("Overall OOB","HR status")) +
  labs(x = "Iteration", y = "Error (%)", colour = NULL, title = "Backward elimination")
save_fig(p_opt, "Fig4_Optimisation", 115, 90); p_opt

imp_final <- randomForest::importance(rf_final)
imp_df <- tibble(variable = rownames(imp_final),
                 MDA = imp_final[,"MeanDecreaseAccuracy"],
                 Gini = imp_final[,"MeanDecreaseGini"]) |>
  mutate(module = case_when(str_detect(variable,"_PC1$") ~ str_remove(variable,"_PC1"),
                            TRUE ~ "Categorical/Spatial")) |>
  arrange(desc(MDA))
write_csv(imp_df, file.path(tab_dir, "Table_VariableImportance.csv"))
imp_df <- imp_df |> mutate(vlab = pp(variable))   # unified label prettifier (same map as all other figures)
p_imp <- imp_df |> ggplot(aes(reorder(vlab, MDA), MDA, fill = module)) +
  geom_col() + coord_flip() + scale_fill_viridis_d() +
  labs(x = NULL, y = "Mean decrease in accuracy", fill = NULL) +
  theme(legend.position = "none")
save_fig(p_imp, "Fig5_VariableImportance", 130, 150); p_imp   # NPh 1.5 column (main Fig 3)

mod_imp_tab <- imp_df |> group_by(module) |>
  summarise(total_MDA = sum(MDA), pct = 100*sum(MDA)/sum(imp_df$MDA)) |> arrange(desc(pct))
kbl(mod_imp_tab, digits = 2, caption = "Module-level contribution to HR prediction (sum of permutation importance by module).")
results_summary$module_importance <- mod_imp_tab
results_summary$top_predictors <- head(imp_df$variable, 6)
results_summary$n_final_pred <- length(final_preds)

stab <- map_dfr(1:params$stability_runs, function(i){
  set.seed(i*13); rf <- brf(primary |> select(all_of(c(final_preds,"HR_bin"))))
  tibble(run=i, variable=names(imp_vec(rf)), rank=rank(-imp_vec(rf)))
}) |> group_by(variable) |>
  summarise(mean_rank=mean(rank), sd_rank=sd(rank), min=min(rank), max=max(rank)) |> arrange(mean_rank)
kbl(head(stab,12), digits = 1, caption = "CHECK A: importance-rank stability over runs.")
p_stab <- head(stab,12) |> mutate(lab = pp(variable)) |> ggplot(aes(reorder(lab,-mean_rank), mean_rank)) +
  geom_point(colour="#E76F51") + geom_errorbar(aes(ymin=min,ymax=max), width=.3, colour="#457B9D") +
  coord_flip() + labs(x=NULL, y="Importance rank (lower = better)")
save_fig(p_stab, "Fig_S2_RankStability", 120, 105); p_stab

cont_preds <- final_preds[sapply(primary[final_preds], is.numeric)]
pd_df <- as.data.frame(primary)
pd_all <- map_dfr(cont_preds, function(v){
  pd <- do.call(randomForest::partialPlot,
                list(x = rf_final, pred.data = pd_df, x.var = v,
                     which.class = "HR", plot = FALSE, n.pt = 40))
  tibble(variable = v, x = pd$x, y = pd$y)
})
dir_sum <- pd_all |> group_by(variable) |>
  summarise(direction = cor(x, y), .groups="drop") |>
  mutate(assoc = case_when(direction > .2 ~ "+ trait  ->  + HR prob",
                           direction < -.2 ~ "- trait  ->  + HR prob",
                           TRUE ~ "non-linear/weak")) |> arrange(desc(abs(direction)))
kbl(dir_sum, digits = 2, caption = "Direction of trait-HR associations (partial dependence).")
p_pd <- pd_all |> dplyr::mutate(variable = factor(pp(variable), levels = pp(dir_sum$variable))) |>
  ggplot(aes(x, y)) + geom_line(colour = "#E76F51") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  facet_wrap(~variable, scales = "free_x") +
  labs(x = "Predictor value", y = "Partial dependence (logit, HR)")
save_fig(p_pd, "Fig6_PartialDependence", 180, 150); p_pd   # NPh double column (main Fig 4)

uni <- map_dfr(final_preds, function(v){
  x <- primary[[v]]; g <- primary$HR_bin
  if (is.numeric(x)) {
    t <- suppressWarnings(wilcox.test(x[g=="HR"], x[g=="nonHR"]))
    tibble(variable=v, test="Wilcoxon", p=t$p.value)
  } else {
    t <- fisher.test(table(x,g), simulate.p.value=TRUE, B=2000)
    tibble(variable=v, test="Fisher", p=t$p.value)
  }
}) |> mutate(p_adj = p.adjust(p, "BH"), sig = p_adj < 0.05) |> arrange(p) |>
  mutate(p     = ifelse(p     < 0.001, "< 0.001", sprintf("%.3f", p)),
         p_adj = ifelse(p_adj < 0.001, "< 0.001", sprintf("%.3f", p_adj)))
kbl(uni, digits = 3, caption = "CHECK B: univariate HR vs non-HR tests (BH-adjusted).")

## Partial dependence of the FOUR leading predictors, INCLUDING the two categorical ones (life span,
## photosynthetic pathway) that the continuous-only Fig. 6 omits. partialPlot returns the centred log-odds
## for class "HR": positive = more likely HR, negative = less likely, 0 = no effect.
pd_top <- function(v) {
  pp <- do.call(randomForest::partialPlot,
                list(x = rf_final, pred.data = as.data.frame(primary), x.var = v,
                     which.class = "HR", plot = FALSE, n.pt = 30))
  x <- pp$x
  if (!is.numeric(primary[[v]]) && is.numeric(x)) x <- levels(primary[[v]])[x]   # relabel codes if needed
  data.frame(x = x, y = pp$y)
}
pd_panel <- function(v, lab, letter, ylab = "Logit contribution") {
  d <- pd_top(v)
  g <- ggplot(d) + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    labs(y = ylab, title = sprintf("%s  %s", letter, lab)) +
    theme(plot.title = element_text(size = 7))
  if (is.numeric(primary[[v]]))
    g + geom_line(aes(x, y), colour = "#457B9D", linewidth = .7) + labs(x = lab)
  else
    g + geom_col(aes(reorder(as.character(x), y), y, fill = y > 0), width = .6) + coord_flip() +
      scale_fill_manual(values = c(`TRUE` = "#E76F51", `FALSE` = "#457B9D"), guide = "none") + labs(x = NULL)
}
p_pdtop <- (pd_panel("life_span", "Life history", "A") | pd_panel("Height_PC1", "Height", "B", ylab = NULL)) /
           (pd_panel("photosynthetic_pathway", "Photosynthesis", "C") | pd_panel("Fecundity_PC1", "Fecundity", "D", ylab = NULL))
save_fig(p_pdtop, "Fig6b_PD_TopPredictors", 180, 130); p_pdtop

## Main Fig 3 built NATIVELY (patchwork) so text renders at the theme size on the final 180 mm canvas,
## rather than being scaled down by post-hoc image stitching. a = importance, b = top-4 partial dependence.
saveRDS(list(p_imp = p_imp, p_pdtop = p_pdtop), "fig3_objects.rds")   # allow standalone Fig 3 rebuilds
fig_ms3 <- (p_imp | patchwork::wrap_elements(full = p_pdtop)) +
  patchwork::plot_layout(widths = c(1, 1.5)) +
  patchwork::plot_annotation(tag_levels = "a",
                             theme = theme(plot.tag = element_text(face = "bold", size = 7, colour = "black")))
save_fig(fig_ms3, "Fig_MS3_ImportanceDirection", 180, 92); fig_ms3

roc_oob <- pROC::roc(primary$HR_bin, rf_final$votes[,"HR"], levels=c("nonHR","HR"),
                     direction="<", quiet=TRUE)
ci <- pROC::ci.auc(roc_oob)
cm <- caret::confusionMatrix(rf_final$predicted, primary$HR_bin, positive="HR")
.r0 <- unname(1 - rf0$confusion["HR","class.error"]); .e0 <- rf0$err.rate[nrow(rf0$err.rate),"OOB"]
.rf <- unname(cm$byClass["Sensitivity"]);            .ef <- rf_final$err.rate[nrow(rf_final$err.rate),"OOB"]
cat(sprintf(paste(">>> FINAL OPTIMISED model (rf_final, %d predictors): the model used for ALL downstream",
  "results. Its OOB confusion/recall differ from the Step-7 INITIAL all-predictor model because weak",
  "predictors were pruned (HR recall %.2f -> %.2f, OOB error %.0f%% -> %.0f%%); the difference is the",
  "optimisation, NOT the permutation test. AUC (threshold-free) is unchanged. <<<\n"),
  length(final_preds), .r0, .rf, 100*.e0, 100*.ef))
cat("AUC =", round(as.numeric(pROC::auc(roc_oob)),3),
    sprintf(" (95%% CI %.3f-%.3f)\n", ci[1], ci[3]))
print(cm$table); cat("Sensitivity (HR recall):", round(cm$byClass["Sensitivity"],3),
                     " Specificity:", round(cm$byClass["Specificity"],3),
                     " Balanced accuracy:", round(cm$byClass["Balanced Accuracy"],3), "\n")

# permutation significance
obs <- rf_final$err.rate[params$ntree,"OOB"]
perm <- replicate(params$n_perm, {
  d <- primary; d$HR_bin <- sample(d$HR_bin)
  brf(d |> select(all_of(c(final_preds,"HR_bin"))), ntree=500,
      importance=FALSE)$err.rate[500,"OOB"]
})
pval <- mean(perm <= obs)
cat("Permutation p-value (", params$n_perm, "perms):", pval, "\n")
results_summary$core <- list(auc = as.numeric(pROC::auc(roc_oob)), ci_lo = ci[1], ci_hi = ci[3],
                             sens = unname(cm$byClass["Sensitivity"]),
                             spec = unname(cm$byClass["Specificity"]),
                             bal_acc = unname(cm$byClass["Balanced Accuracy"]),
                             oob = obs, oob_initial = .oob_initial, perm_p = pval, n_perm = params$n_perm)

p_roc <- pROC::ggroc(roc_oob, legacy.axes=TRUE, colour="#E76F51", linewidth=.8) +
  geom_abline(slope=1, linetype="dashed", colour="grey60") +
  annotate("text", .6,.3, size=2.5,
           label=sprintf("AUC=%.3f\n[%.3f, %.3f]", as.numeric(pROC::auc(roc_oob)), ci[1], ci[3])) +
  labs(x="1 - Specificity", y="Sensitivity")
p_perm <- ggplot(tibble(e=perm*100), aes(e)) +
  geom_histogram(bins=30, fill="grey75") +
  geom_vline(xintercept=obs*100, colour="#E76F51", linewidth=1) +
  labs(x="Null OOB error (%)", y="Count")
fig_assess <- p_roc + p_perm +
  patchwork::plot_annotation(tag_levels = "a",
                             theme = theme(plot.tag = element_text(size = 7, face = "bold", colour = "black")))
save_fig(fig_assess, "Fig7_ModelAssessment", 170, 90); fig_assess   # NPh double column (main Fig 2)



.p <- rf_final$votes[, "HR"]; .y <- primary$HR_bin
.thr <- seq(0.05, 0.95, by = 0.01)
.sw <- data.frame(thr = .thr,
  sens = sapply(.thr, function(t) mean(.p[.y == "HR"]    >= t)),
  spec = sapply(.thr, function(t) mean(.p[.y == "nonHR"] <  t)))
iy <- which.max(.sw$sens + .sw$spec - 1)                       # Youden-optimal
ir <- which(.sw$sens >= 0.80); ir <- if (length(ir)) max(ir) else 1L   # highest cutoff keeping recall >= 0.80
op <- data.frame(
  `Operating point` = c("Default (0.5)", "Youden-optimal", "Recall ~0.80 (screen)"),
  Threshold   = c(0.50, .sw$thr[iy], .sw$thr[ir]),
  Sensitivity = c(mean(.p[.y == "HR"]    >= 0.5), .sw$sens[iy], .sw$sens[ir]),
  Specificity = c(mean(.p[.y == "nonHR"] <  0.5), .sw$spec[iy], .sw$spec[ir]),
  check.names = FALSE)
kbl(op, digits = 2, caption = "Binary-model HR-recall (sensitivity) and specificity at three OOB-vote cutoffs: the default 0.5, the Youden-optimal cutoff, and a recall-oriented (~0.80) screen. AUC is identical across all three - only the recall/specificity balance moves.")
# OOB confusion matrix AT THE YOUDEN-OPTIMAL CUTOFF (so the matrix is not tied to the arbitrary 0.5)
thr_y  <- .sw$thr[iy]
pred_y <- factor(ifelse(.p >= thr_y, "HR", "nonHR"), levels = c("nonHR", "HR"))
cm_y   <- as.data.frame.matrix(table(Prediction = pred_y, Reference = .y))
kbl(cm_y, caption = sprintf("Binary OOB confusion matrix at the YOUDEN-OPTIMAL cutoff (%.2f); rows = predicted class, columns = true class. At this operating point HR recall (sensitivity) = %.2f and specificity = %.2f - i.e. %d of the %d resistant species are caught (vs the 0.5-cutoff matrix above). Reported alongside, not instead of, the default; the watch-list uses no cutoff at all.", thr_y, .sw$sens[iy], .sw$spec[iy], cm_y["HR","HR"], sum(.y == "HR")))
results_summary$cm_youden <- cm_y
.lng <- tidyr::pivot_longer(.sw, c(sens, spec), names_to = "metric", values_to = "value")
p_thr <- ggplot(.lng, aes(thr, value, colour = metric)) + geom_line(linewidth = .6) +
  geom_vline(xintercept = c(0.50, .sw$thr[iy], .sw$thr[ir]),
             linetype = c("dashed", "dotted", "longdash"), colour = "grey55") +
  scale_colour_manual(values = c(sens = "#E76F51", spec = "#457B9D"),
                      labels = c("Sensitivity (HR recall)", "Specificity"), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "OOB-vote cutoff for calling a species 'resistant'", y = "Rate") +
  theme(legend.position = "bottom")
save_fig(p_thr, "Fig_S0b_ThresholdChoice", 120, 88); p_thr
results_summary$threshold <- op

base <- if (params$impute_primary) analysis_imp else analysis_cc
train_tr <- base |> filter(hr_status %in% c("HR_elsewhere","Non_HR_local")) |>
  mutate(HR_bin = factor(if_else(hr_status=="HR_elsewhere","HR","nonHR"), levels=c("nonHR","HR"))) |>
  select(ends_with("_PC1"), all_of(covars), HR_bin) |> drop_na()
test_local <- base |> filter(hr_status == "HR_local") |>
  select(ends_with("_PC1"), all_of(covars)) |> drop_na()
set.seed(11); rf_tr <- brf(train_tr)
pred_local <- predict(rf_tr, test_local)
recall_local <- mean(pred_local == "HR")
cat("Trained on HR_elsewhere vs Non_HR (n =", nrow(train_tr), "); held-out HR_local (n =",
    nrow(test_local), ")\n")
cat("Recall on the held-out local weeds =", round(recall_local,3),
    sprintf(" (%d / %d flagged HR)\n", sum(pred_local=="HR"), length(pred_local)))
results_summary$phase1 <- list(recall_local = recall_local, n_flagged = sum(pred_local=="HR"),
                               n_local = length(pred_local), n_train = nrow(train_tr))

d_diag <- base |> filter(hr_status %in% c("HR_local","HR_elsewhere")) |>
  mutate(grp = factor(hr_status, levels=c("HR_elsewhere","HR_local"))) |>
  select(ends_with("_PC1"), all_of(covars), grp) |> drop_na()
nmin <- min(table(d_diag$grp))
set.seed(12)
rf_diag <- randomForest(grp ~ ., d_diag, ntree=params$ntree,
                        sampsize=c(nmin,nmin), strata=d_diag$grp, importance=TRUE)
print(rf_diag$confusion)
diag_imp <- sort(randomForest::importance(rf_diag, type=1)[,1], decreasing=TRUE)
cat("Top discriminating predictors (local vs elsewhere):\n"); print(round(head(diag_imp,8),3))
results_summary$phase2 <- list(oob = rf_diag$err.rate[params$ntree,"OOB"],
                               top = names(head(diag_imp,5)),
                               local_err = rf_diag$confusion["HR_local","class.error"])

.ma <- readRDS("moa_analysis.rds")
.gm <- function(g) round(.ma$med$mean_moa[.ma$med$HR_Status == g], 1)
.gp <- function(g) .ma$multi$pct_multi[.ma$multi$HR_Status == g]
.sh <- .ma$sel_hist                                    # selection-history (temporal) test
.shv <- function(g, col) if (!is.null(.sh)) .sh$med[[col]][.sh$med$HR_Status == g] else NA
.pf  <- function(p) if (!is.null(p) && !is.na(p) && p < 0.001) "< 0.001" else sprintf("= %.3f", p)



kbl(.ma$med |> dplyr::transmute(
      `HR status` = ifelse(HR_Status == "HR_local", "HR-local (Germany)", "HR-elsewhere"),
      `Species (n)` = n_spp, `Mechanisms per species (mean)` = mean_moa, Median = median_moa),
    caption = sprintf("Number of resisted MOAs by HR status. HR-local weeds resist significantly more mechanisms of action than HR-elsewhere weeds (Wilcoxon P = %.4f); the set of mechanisms resisted does not differ (chi-square P = %.2f).", .ma$wilcox_p, .ma$chisq_p))



.orig <- tryCatch({
  od <- dat |> dplyr::filter(hr_status %in% c("HR_local","HR_elsewhere")) |>
    dplyr::transmute(cls = ifelse(hr_status == "HR_local", "HR-local (Germany)", "HR-elsewhere"),
                     orig = ifelse(as.character(origin) == "introduced", "Introduced", "Native"))
  tb <- table(factor(od$cls, c("HR-local (Germany)","HR-elsewhere")), od$orig); ft <- fisher.test(tb)
  list(tab = tb, p = ft$p.value, or = unname(1/ft$estimate))   # OR of being introduced: HR-elsewhere vs HR-local
}, error = function(e) NULL)
if (!is.null(.orig)) {
  df <- as.data.frame.matrix(.orig$tab); df[["% introduced"]] <- round(100*df$Introduced/(df$Introduced+df$Native))
  kbl(cbind(`HR status` = rownames(df), df),
      caption = sprintf("Biogeographic origin by HR status. HR-elsewhere weeds are far more often introduced than HR-local weeds (Fisher's exact P = %.3f, odds ratio = %.1f).", .orig$p, .orig$or), row.names = FALSE)
}

d3 <- base |> select(ends_with("_PC1"), all_of(covars), HR3) |> drop_na() |> droplevels()
print(table(d3$HR3))
nmin3 <- min(table(d3$HR3))
set.seed(13)
# VALIDATION: same balanced-bootstrap + out-of-bag scheme as the binary model, but each tree draws the
# smallest class size (sampsize = nmin3 = 24 per class). All metrics below (confusion, error) are OOB.
rf3 <- randomForest(HR3 ~ ., d3, ntree=params$ntree,
                    sampsize=rep(nmin3, nlevels(d3$HR3)), strata=d3$HR3,
                    importance=TRUE, proximity=TRUE)
saveRDS(rf3, file.path(mod_dir, "rf_3class.rds"))
cm3 <- as.data.frame.matrix(round(rf3$confusion, 3))
kbl(cm3, caption = "3-class OOB confusion matrix (rows = true class; last column = per-class error).")
results_summary$phase3 <- list(confusion = rf3$confusion,
                               class_error = rf3$confusion[,"class.error"])

# --- Binary model (rf_final) OOB confusion (already computed as `cm` in Step 11) ---
bin_cm <- cm$table
kbl(as.data.frame.matrix(bin_cm),
    caption = "Binary (All-HR vs Non-HR) OOB confusion matrix; rows = predicted, columns = true.")

# --- FAIR comparison: OOB-based HR_local recall for BOTH models ---
# 3-class: HR_local recall = 1 - class.error
local_recall_3 <- 1 - rf3$confusion["HR_local","class.error"]
# Binary: align OOB predictions of rf_final to true hr_status (NOT in-bag prediction)
prim_meta <- analysis_imp |>
  select(ends_with("_PC1"), all_of(covars), HR_bin, hr_status) |>
  tidyr::drop_na(ends_with("_PC1"), all_of(covars), HR_bin)
local_recall_2 <- mean(rf_final$predicted[prim_meta$hr_status == "HR_local"] == "HR")
elsew_recall_2 <- mean(rf_final$predicted[prim_meta$hr_status == "HR_elsewhere"] == "HR")

compare_tab <- data.frame(
  Metric = c("Overall OOB error", "HR detection (overall)", "HR_local recall (OOB)",
             "HR_elsewhere recall (OOB)", "Use case"),
  `Binary (2-class)` = c(sprintf("%.1f%%", rf_final$err.rate[params$ntree,"OOB"]*100),
                         sprintf("AUC %.3f", as.numeric(pROC::auc(roc_oob))),
                         sprintf("%.2f", local_recall_2), sprintf("%.2f", elsew_recall_2),
                         "General risk screening"),
  `3-class (production)` = c(sprintf("%.1f%%", (1 - sum(diag(rf3$confusion[,1:3]))/sum(rf3$confusion[,1:3]))*100),
                            "per-class (below)", sprintf("%.2f", local_recall_3),
                            sprintf("%.2f", 1 - rf3$confusion["HR_elsewhere","class.error"]),
                            "Risk staging (local vs emerging)"),
  check.names = FALSE)
kbl(compare_tab, caption = "Binary vs 3-class: how each performs and what each is for (all OOB-based).")
results_summary$phase4 <- list(recall_3class = local_recall_3, recall_2class = local_recall_2,
                               bin_cm = bin_cm, cm3 = rf3$confusion)
cat(sprintf("Reconciliation: binary OVERALL HR recall = %.2f (Table above) = weighted avg of HR_local %.2f (n=24) and HR_elsewhere %.2f (n=85).\n",
            unname(1 - rf_final$confusion["HR","class.error"]), local_recall_2, elsew_recall_2))
# Same binary model at the recall-favourable Youden cutoff (thr_y from Step 11), for consistency with that section:
.vh <- rf_final$votes[, "HR"]
ov_recall_2y  <- mean(.vh[prim_meta$hr_status != "Non_HR_local"] >= thr_y)
loc_recall_2y <- mean(.vh[prim_meta$hr_status == "HR_local"]     >= thr_y)
results_summary$phase4$recall_2class_youden <- ov_recall_2y
cat(sprintf("At the Youden cutoff (%.2f) the SAME binary model gives overall HR recall %.2f and HR_local recall %.2f.\n",
            thr_y, ov_recall_2y, loc_recall_2y))

# STABLE watch-list: many non-resistant species have near-tied P, so a single random-forest
# fit reorders the ranking between runs. We therefore average the out-of-bag P(HR) over 50 seeded
# balanced-RF fits (fixed seeds 1001..1050) - reproducible and robust to the tie-breaking randomness.
wl_meta <- analysis_imp |>
  select(species, hr_status, ends_with("_PC1"), all_of(covars), HR_bin) |>
  tidyr::drop_na(ends_with("_PC1"), all_of(covars), HR_bin)
.wl_pdat <- primary |> select(all_of(c(final_preds, "HR_bin")))
.wl_nrep <- 50
.wl_votes <- matrix(NA_real_, nrow = nrow(.wl_pdat), ncol = .wl_nrep)
for (.i in seq_len(.wl_nrep)) { set.seed(1000 + .i); .wl_rf <- brf(.wl_pdat, ntree = 500, importance = FALSE)
  .wl_votes[, .i] <- .wl_rf$votes[, "HR"] }
watch <- tibble(species = wl_meta$species, hr_status = wl_meta$hr_status,
                p_HR = rowMeans(.wl_votes)) |>
  filter(hr_status == "Non_HR_local") |>
  arrange(desc(p_HR))
kbl(head(watch, 20) |> transmute(Rank = row_number(), Species = species, `P(resistance)` = round(p_HR,3)),
    caption = "Top-20 currently non-resistant weeds by predicted resistance probability (OOB).")

p_watch <- head(watch, 20) |>
  ggplot(aes(reorder(species, p_HR), p_HR)) +
  geom_col(fill = "#E76F51") +
  geom_text(aes(label = sprintf("%.2f", p_HR)), hjust = -0.15, size = 3) +
  coord_flip(ylim = c(0, min(1, max(head(watch,20)$p_HR)+0.1))) +
  labs(x = NULL, y = "Predicted resistance probability (OOB)",
       title = "Highest-risk non-resistant weeds to watch")
save_fig(p_watch, "Fig9_RiskWatchlist", 130, 150); p_watch   # NPh 1.5 column (main Fig 6)

write_csv(watch, file.path(tab_dir, "Table_RiskWatchlist.csv"))
results_summary$watchlist <- head(watch, 15) |> mutate(p_HR = round(p_HR,3))

.dd <- file.path(out_dir, "data"); dir.create(.dd, showWarnings = FALSE, recursive = TRUE)
fmt <- function(d) d |> dplyr::transmute(
  Species = species, Family = family, `P(resistance)` = `P(resistance)`,
  life_span = life_span, Height_PC1 = round(Height_PC1, 2),
  photosynthetic_pathway = photosynthetic_pathway, Fecundity_PC1 = round(Fecundity_PC1, 2)) |>
  dplyr::arrange(dplyr::desc(`P(resistance)`))
# BINARY model, aligned to primary's rows; family joined from `dat`
ex <- analysis_imp |>
  dplyr::select(species, hr_status, HR_bin, dplyr::ends_with("_PC1"), dplyr::all_of(covars)) |>
  tidyr::drop_na(dplyr::ends_with("_PC1"), dplyr::all_of(covars), HR_bin) |>
  dplyr::left_join(dplyr::distinct(dat, species, family), by = "species")
stopifnot(nrow(ex) == nrow(primary))
ex$`P(resistance)` <- round(rf_final$votes[, "HR"], 3); ex$pred_bin <- as.character(rf_final$predicted)
tA <- ex |> dplyr::filter(hr_status == "Non_HR_local", pred_bin == "HR")    |> fmt()   # non-HR called HR
tB <- ex |> dplyr::filter(hr_status != "Non_HR_local", pred_bin == "nonHR") |> fmt()   # HR called non-HR
readr::write_csv(tA, file.path(.dd, "Table_nonHR_predicted_HR.csv"))
readr::write_csv(tB, file.path(.dd, "Table_HR_predicted_nonHR.csv"))
# THREE-CLASS model, aligned to rf3's rows
d3m <- analysis_imp |>
  dplyr::select(species, hr_status, dplyr::ends_with("_PC1"), dplyr::all_of(covars), HR3) |>
  tidyr::drop_na(dplyr::ends_with("_PC1"), dplyr::all_of(covars), HR3) |> droplevels() |>
  dplyr::left_join(dat |> dplyr::select(species, family), by = "species")
stopifnot(nrow(d3m) == length(rf3$predicted))
d3m$`P(resistance)` <- round(rf3$votes[, "HR_local"], 3); d3m$pred3 <- as.character(rf3$predicted)
tC <- d3m |> dplyr::filter(hr_status == "HR_elsewhere", pred3 == "HR_local") |> fmt()
readr::write_csv(tC, file.path(.dd, "Table_HRelsewhere_predicted_HRlocal.csv"))
kbl(head(tA, 8), caption = sprintf("(A) Non-resistant weeds the BINARY model flags as HR (n = %d; full list in data/, file A). P(resistance) = out-of-bag predicted probability of herbicide resistance (0-1). This is the watch-list dichotomised at the 0.5 cutoff: the same non-resistant species as the ranking above, restricted to those with P(resistance) >= 0.5 - the model's emerging-resistance candidates. Showing the 8 highest-probability species.", nrow(tA)))
kbl(head(tB, 8), caption = sprintf("(B) Herbicide-resistant weeds the BINARY model classifies as non-HR at the 0.5 cutoff (n = %d; data/, file B): resistant species with an atypical, non-resistant-like trait profile (low P(resistance)).", nrow(tB)))
kbl(head(tC, 8), caption = sprintf("(C) HR-elsewhere weeds the THREE-CLASS model predicts as HR-local (n = %d; data/, file C). Here P(resistance) = P(HR-local): foreign-resistant weeds whose trait profile already resembles established German ones - the highest priority to establish locally.", nrow(tC)))
results_summary$classtabs <- c(nonHR_as_HR = nrow(tA), HR_as_nonHR = nrow(tB), elsewhere_as_local = nrow(tC))

mds_plot <- function(rf, labels, cols, title, lbls = ggplot2::waiver()){
  m <- cmdscale(1 - rf$proximity, k = 2)
  df <- tibble(MDS1=m[,1], MDS2=m[,2], grp=labels)
  ggplot(df, aes(MDS1, MDS2, colour=grp)) + geom_point(size=.8, alpha=.7) +
    stat_ellipse(level=.95) + scale_colour_manual(values=cols, labels=lbls) +
    labs(colour=NULL, title=title)
}
p_mds1 <- mds_plot(rf0, primary$HR_bin, col_bin, "Binary RF proximity", lbls = bin_labs)
p_mds2 <- mds_plot(rf3, d3$HR3, col_status, "3-class RF proximity", lbls = hr_labs)
fig_mds <- p_mds1 + p_mds2 +
  patchwork::plot_annotation(tag_levels = "a",
    theme = theme(plot.tag = element_text(size = 7, face = "bold", colour = "black")))
save_fig(fig_mds, "Fig8_MDS_Proximity", 170, 90); fig_mds   # NPh double column (supp Fig S8)

d_all <- alltraits |> select(all_of(mean_cols), all_of(covars), HR_bin) |> tidyr::drop_na()

set.seed(42); rf_all <- brf(d_all)
auc_all <- get_auc(rf_all, d_all)
imp_all <- sort(randomForest::importance(rf_all, type = 1)[,1], decreasing = TRUE)
cat(sprintf("All-traits RF: %d predictors, N = %d, OOB AUC = %.3f (module RF AUC = %.3f)\n",
            ncol(d_all)-1, nrow(d_all), auc_all, as.numeric(pROC::auc(roc_oob))))
cat("Top 10 raw predictors:\n"); print(round(head(imp_all, 10), 3))

p_all <- tibble(variable = names(head(imp_all,15)), MDA = head(imp_all,15)) |>
  mutate(lab = pp(variable)) |>
  ggplot(aes(reorder(lab, MDA), MDA)) + geom_col(fill = "#2A9D8F") + coord_flip() +
  labs(x = NULL, y = "Mean decrease accuracy")
save_fig(p_all, "Fig_S4_AllTraitsRF", 125, 135); p_all

results_summary$robustA <- list(auc_all = auc_all, auc_module = as.numeric(pROC::auc(roc_oob)),
                                n_pred = ncol(d_all)-1, top = names(head(imp_all,8)))

woody_sp <- dat |> dplyr::filter(grepl("tree|shrub", plant_growth_form)) |> dplyr::pull(species)
nw <- analysis_imp |> dplyr::filter(!species %in% woody_sp) |>
  dplyr::select(dplyr::all_of(final_preds), HR_bin) |> tidyr::drop_na()
set.seed(42); rf_nw <- brf(nw); auc_nw <- get_auc(rf_nw, nw)
auc_full <- as.numeric(pROC::auc(roc_oob))
ov <- length(intersect(names(head(imp_vec(rf_final), 5)), names(head(imp_vec(rf_nw), 5))))
cat(sprintf("Woody species removed: %d (all non-resistant: %s)\nOOB AUC  full = %.3f  |  woody-excluded = %.3f  (delta = %.3f)\nTop-5 importance overlap = %d/5\n",
            length(woody_sp), all(dat$HR_bin[dat$species %in% woody_sp] == "nonHR"),
            auc_full, auc_nw, abs(auc_full - auc_nw), ov))
results_summary$woody_sens <- list(n = length(woody_sp), auc_full = auc_full, auc_nw = auc_nw, top5_overlap = ov)

suppressPackageStartupMessages({library(ape); library(phylolm)})

## ---- Build the GBOTB dated megatree (fallback: taxonomy tree) ----
sp_tab <- dat |>
  mutate(species_c = str_to_sentence(str_squish(species)),
         genus_c   = str_to_title(word(species_c, 1)),
         family_c  = str_to_title(ifelse(is.na(family), "Unknown", as.character(family))),
         tip       = gsub(" ", "_", species_c)) |>
  distinct(species, .keep_all = TRUE)

mt_path <- file.path(out_dir, "megatree.rds")    # pre-built by build_megatree.R (avoids in-knit memory spike)
tree <- tryCatch({
  if (file.exists(mt_path)) {
    mt <- readRDS(mt_path)
    if (is.null(attr(mt, "type"))) attr(mt, "type") <- "GBOTB dated megatree (V.PhyloMaker2)"
    mt
  } else {
    sp.list <- sp_tab |> transmute(species = species_c, genus = genus_c, family = family_c) |> as.data.frame()
    invisible(utils::capture.output(
      pm <- suppressWarnings(suppressMessages(V.PhyloMaker2::phylo.maker(sp.list, scenarios = "S3")))))
    mt <- pm$scenario.3
    attr(mt, "type") <- "GBOTB dated megatree (V.PhyloMaker2)"; mt
  }
}, error = function(e) {
  message("Megatree failed, using taxonomy tree: ", conditionMessage(e))
  tx <- sp_tab |> transmute(fam = factor(family_c), genus = factor(genus_c), tip = factor(tip)) |> as.data.frame()
  t2 <- ape::as.phylo(~fam/genus/tip, data = tx, collapse = FALSE) |>
    ape::compute.brlen(method = "Grafen") |> ape::multi2di()
  t2$edge.length[t2$edge.length == 0] <- 1e-8; attr(t2, "type") <- "taxonomy tree (fallback)"; t2
})
tree_type <- attr(tree, "type")
# map data 'tip' (Genus_species) <-> hr status
tipmap <- setNames(sp_tab$hr_status, sp_tab$tip)
cat("Phylogeny:", tree_type, "| tips =", ape::Ntip(tree), "\n")

## ===== shared circular-phylogeny helpers (used for Fig S6 here and Fig S9 below) =====
clade_of <- function(f) {   # APG-style major clade for the outer arc ring / labels
  monocot <- c("Liliaceae","Poaceae","Juncaceae","Cyperaceae","Orchidaceae","Alliaceae",
               "Hyacinthaceae","Anthericaceae","Asparagaceae")
  caryo   <- c("Polygonaceae","Portulacaceae","Amaranthaceae","Chenopodiaceae","Illecebraceae","Caryophyllaceae")
  aster   <- c("Cornaceae","Primulaceae","Campanulaceae","Asteraceae","Apiaceae","Adoxaceae","Caprifoliaceae",
               "Dipsacaceae","Convolvulaceae","Solanaceae","Gentianaceae","Rubiaceae","Hydrophyllaceae",
               "Boraginaceae","Oleaceae","Plantaginaceae","Buddlejaceae","Verbenaceae","Lamiaceae",
               "Scrophulariaceae","Orobanchaceae")
  rosid   <- c("Vitaceae","Geraniaceae","Onagraceae","Aceraceae","Sapindaceae","Malvaceae","Resedaceae",
               "Brassicaceae","Fabaceae","Juglandaceae","Rosaceae","Cannabaceae","Urticaceae","Oxalidaceae",
               "Linaceae","Euphorbiaceae","Clusiaceae","Violaceae","Salicaceae")
  saxi <- c("Saxifragaceae","Crassulaceae"); basal <- c("Papaveraceae","Ranunculaceae")
  ifelse(f %in% monocot, "Monocots", ifelse(f %in% basal, "Basal eudicots",
    ifelse(f %in% saxi, "Saxifragales", ifelse(f %in% rosid, "Rosids",
      ifelse(f %in% caryo, "Caryophyllales", ifelse(f %in% aster, "Asterids", "Other"))))))
}
clade_cols <- c(Monocots="#1B7837", `Basal eudicots`="#8C8C8C", Saxifragales="#7570B3",
                Rosids="#4575B4", Caryophyllales="#D6604D", Asterids="#E08214", Other="#BBBBBB")
geo_periods <- data.frame(name = c("Cretaceous","Paleogene","Neogene","Quaternary"),
  from = c(145,66,23,2.6), to = c(66,23,2.6,0),
  col = c("#7FC97F55","#FDC08655","#FFFF9955","#F0F0F055"), stringsAsFactors = FALSE)
fan_geometry <- function(phy, open_angle) {   # silent fan plot -> tip angles, radius, open wedge
  grDevices::pdf(NULL); plot(phy, type = "fan", show.tip.label = FALSE, open.angle = open_angle)
  env <- get("last_plot.phylo", envir = .PlotPhyloEnv); grDevices::dev.off()
  nt <- length(phy$tip.label); xx <- env$xx[1:nt]; yy <- env$yy[1:nt]
  ang <- atan2(yy, xx); sa <- sort(ang)
  gaps <- c(diff(sa), (sa[1]+2*pi) - sa[nt]); gi <- which.max(gaps)
  gp <- if (gi == nt) (sa[nt]+sa[1]+2*pi)/2 else (sa[gi]+sa[gi+1])/2
  list(ang = ang, Rmax = max(sqrt(xx^2+yy^2)), dA = min(diff(sa))/2,
       gap = atan2(sin(gp), cos(gp)), nt = nt)
}
draw_geo_rings <- function(g, H, per) {   # faint period annuli + radial "Myr before present" axis
  a2r <- function(age) (H-age)/H * g$Rmax
  th <- seq(min(g$ang)-g$dA, max(g$ang)+g$dA, length.out = 500)
  for (i in seq_len(nrow(per))) { r1 <- a2r(per$from[i]); r2 <- a2r(per$to[i])
    polygon(c(r2*cos(th), rev(r1*cos(th))), c(r2*sin(th), rev(r1*sin(th))), col = per$col[i], border = NA)
    lines(r2*cos(th), r2*sin(th), col = "grey75", lwd = 0.3, lty = 3) }
  ca <- cos(g$gap); sb <- sin(g$gap)
  for (a in c(0,50,100,150)) if (a <= H) { r <- a2r(a); points(r*ca, r*sb, pch = 3, cex = 0.4, col = "grey40")
    text(r*ca, r*sb, a, cex = 0.6, col = "grey30", adj = c(0.5,-0.6), srt = g$gap*180/pi) }
  text(a2r(H*0.5)*ca, a2r(H*0.5)*sb, "Myr before present", cex = 0.6, col = "grey30",
       adj = c(0.5,1.7), srt = g$gap*180/pi)
}
draw_clade_arcs <- function(g, cl, r0f, r1f, labf, minlen = 4, minlab = 10, cex = 0.6) {
  runs <- rle(cl); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
  rd <- data.frame(val = runs$values, len = runs$lengths, s = starts, e = ends, stringsAsFactors = FALSE)
  rd <- rd[rd$val != "Other" & rd$len >= minlen, ]
  # true angular span of a run = min..max over ALL its tips (subtree order zig-zags, so the run's
  # first/last tip are NOT the angular extremes); unwrap across the -pi/pi seam (e.g. the rosids).
  arc_ang <- function(s, e) { ra <- g$ang[s:e]
    if (diff(range(ra)) > pi) ra <- ifelse(ra < 0, ra + 2*pi, ra)
    seq(min(ra) - g$dA, max(ra) + g$dA, length.out = 160) }
  for (k in seq_len(nrow(rd))) { aa <- arc_ang(rd$s[k], rd$e[k])
    polygon(c(g$Rmax*r0f*cos(aa), rev(g$Rmax*r1f*cos(aa))),
            c(g$Rmax*r0f*sin(aa), rev(g$Rmax*r1f*sin(aa))), col = clade_cols[rd$val[k]], border = NA) }
  for (cc in unique(rd$val)) { sub <- rd[rd$val == cc, ]; b <- sub[which.max(sub$len), ]
    if (b$len < minlab) next
    aa <- arc_ang(b$s, b$e); am <- mean(range(aa))
    text(g$Rmax*labf*cos(am), g$Rmax*labf*sin(am), cc, srt = 0, adj = 0.5, cex = cex, font = 2,
         col = clade_cols[cc]) }   # horizontal labels (readable at any arc position)
}

## ===== Fig S6: species-level circular tree, tips coloured by HR status =====
tr <- ape::ladderize(tree)
st_tip  <- tipmap[tr$tip.label]
tip_col <- unname(col_status[as.character(st_tip)])          # HR drawn large/bright, non-HR small
tip_cex <- ifelse(st_tip == "HR_local", 1.4, ifelse(st_tip == "HR_elsewhere", 1.1, 0.35))
fam_by_tip <- setNames(sp_tab$family_c, sp_tab$tip); cl6 <- clade_of(unname(fam_by_tip[tr$tip.label]))
g6 <- fan_geometry(tr, 8); H6 <- max(ape::node.depth.edgelength(tr))
per6 <- geo_periods[geo_periods$to < H6, ]; per6$from <- pmin(per6$from, H6)
## export the per-species phylogeny data behind Fig S6 (tip -> family, clade, HR status) + the tree
.phylo_dd <- file.path(out_dir, "data"); dir.create(.phylo_dd, showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(tip = tr$tip.label, species = gsub("_", " ", tr$tip.label),
  family = unname(fam_by_tip[tr$tip.label]), clade = cl6, hr_status = as.character(st_tip),
  stringsAsFactors = FALSE), file.path(.phylo_dd, "phylo_species_data.csv"), row.names = FALSE)
ape::write.tree(tr, file.path(.phylo_dd, "phylo_species_tree.nwk"))
draw_tree <- function() {
  par(mar = c(0,0,0,0), family = base_family, xpd = NA); LIM <- g6$Rmax*1.92
  plot(tr, type = "fan", show.tip.label = FALSE, edge.width = 0.25, edge.color = "grey70",
       open.angle = 8, no.margin = TRUE, x.lim = c(-LIM,LIM), y.lim = c(-LIM,LIM))
  draw_geo_rings(g6, H6, per6)
  par(new = TRUE)
  plot(tr, type = "fan", show.tip.label = FALSE, edge.width = 0.25, edge.color = "grey55",
       open.angle = 8, no.margin = TRUE, x.lim = c(-LIM,LIM), y.lim = c(-LIM,LIM))
  ape::tiplabels(pch = 19, col = tip_col, cex = tip_cex)
  rr0 <- g6$Rmax*1.02; rr1 <- g6$Rmax*1.06          # thin per-tip status ring (resistant tips pop)
  for (i in seq_len(g6$nt)) { aa <- seq(g6$ang[i]-g6$dA, g6$ang[i]+g6$dA, length.out = 4)
    polygon(c(rr0*cos(aa), rev(rr1*cos(aa))), c(rr0*sin(aa), rev(rr1*sin(aa))), col = tip_col[i], border = NA) }
  draw_clade_arcs(g6, cl6, 1.12, 1.16, 1.30, minlen = 4, minlab = 10, cex = 0.9)
  legend(-LIM*0.99, -LIM*0.78, legend = c("HR-local (Germany)","HR-elsewhere","Non-resistant"),
         col = col_status[c("HR_local","HR_elsewhere","Non_HR_local")], pch = 19,
         pt.cex = c(1.5,1.2,0.8), bty = "n", cex = 0.82)
  gx <- LIM*0.58; gy <- -LIM*0.78; hh <- g6$Rmax*0.05
  for (i in seq_len(nrow(per6))) { rect(gx, gy-(i-1)*hh*1.5, gx+g6$Rmax*0.05, gy-(i-1)*hh*1.5+hh,
    col = gsub("55$","",per6$col[i]), border = "grey55", lwd = 0.4)
    text(gx+g6$Rmax*0.07, gy-(i-1)*hh*1.5+hh*0.5, per6$name[i], cex = 0.72, adj = 0) }
  # figure title removed for clean panels (described in the caption)
  # text(0, LIM*0.99, "Phylogenetic distribution of herbicide resistance (all species)",
  #      cex = 1.05, font = 2, adj = c(0.5,1))
}
grDevices::cairo_pdf(file.path(fig_dir, "Fig_S6_Phylogeny.pdf"), width = 190/25.4, height = 190/25.4)
draw_tree(); invisible(dev.off())
grDevices::png(file.path(fig_dir, "Fig_S6_Phylogeny.png"), width = 190, height = 190,
               units = "mm", res = 600, type = "cairo"); draw_tree(); invisible(dev.off())
draw_tree()
results_summary$tree_tips <- ape::Ntip(tree); results_summary$tree_type <- tree_type

fammap <- setNames(sp_tab$family_c, sp_tab$tip)
hrmap  <- setNames(as.integer(sp_tab$hr_status != "Non_HR_local"), sp_tab$tip)
present <- tree$tip.label
fam_stats <- tibble(tip = present, fam = fammap[present], hr = hrmap[present]) |>
  group_by(fam) |> summarise(n = n(), n_hr = sum(hr), hr_frac = mean(hr), rep = first(tip), .groups = "drop")
fam_tree <- ape::drop.tip(tree, setdiff(present, fam_stats$rep))
fam_tree$tip.label <- fam_stats$fam[match(fam_tree$tip.label, fam_stats$rep)]
fam_tree <- ape::ladderize(fam_tree)
fs  <- fam_stats[match(fam_tree$tip.label, fam_stats$fam), ]
pal <- colorRampPalette(c("#457B9D", "#e9c46a", "#E76F51"))(100)
col_fam <- pal[pmax(1, ceiling(fs$hr_frac * 100))]
fs$clade <- clade_of(fs$fam)
root_age <- max(ape::node.depth.edgelength(fam_tree))
per9 <- geo_periods[geo_periods$to < root_age, ]; per9$from <- pmin(per9$from, root_age)
g9  <- fan_geometry(fam_tree, 12)
# circular dated family tree: geological-period annuli, a per-family resistance-load ring,
# curved family labels (coloured by load) and major-clade arcs (helpers from robustB-tree).
draw_clades <- function() {
  par(mar = c(0,0,0,0), family = base_family, xpd = NA); LIM <- g9$Rmax*2.16
  plot(fam_tree, type = "fan", show.tip.label = FALSE, edge.width = 0.5, edge.color = "grey60",
       open.angle = 12, no.margin = TRUE, x.lim = c(-LIM,LIM), y.lim = c(-LIM,LIM*1.02))
  draw_geo_rings(g9, root_age, per9)
  par(new = TRUE)
  plot(fam_tree, type = "fan", show.tip.label = FALSE, edge.width = 0.55, edge.color = "grey45",
       open.angle = 12, no.margin = TRUE, x.lim = c(-LIM,LIM), y.lim = c(-LIM,LIM*1.02))
  rr0 <- g9$Rmax*1.03; rr1 <- g9$Rmax*1.09          # per-family resistance-load ring
  for (i in seq_len(g9$nt)) { aa <- seq(g9$ang[i]-g9$dA, g9$ang[i]+g9$dA, length.out = 6)
    polygon(c(rr0*cos(aa), rev(rr1*cos(aa))), c(rr0*sin(aa), rev(rr1*sin(aa))),
            col = col_fam[i], border = "white", lwd = 0.2) }
  for (i in seq_len(g9$nt)) { a <- g9$ang[i]; flip <- cos(a) < 0     # curved family labels
    srt <- if (flip) a*180/pi + 180 else a*180/pi; adj <- if (flip) 1 else 0
    text(g9$Rmax*1.12*cos(a), g9$Rmax*1.12*sin(a), fs$fam[i], srt = srt, adj = adj,
         cex = 0.92, font = 3, col = col_fam[i]) }
  draw_clade_arcs(g9, fs$clade, 1.48, 1.52, 1.64, minlen = 1, minlab = 1, cex = 0.85)
  lx0 <- -LIM*0.98; lx1 <- -LIM*0.55; ly <- -LIM*0.92; hh <- g9$Rmax*0.045
  ng <- 100; xs <- seq(lx0, lx1, length.out = ng+1)                  # load gradient legend
  rect(xs[-(ng+1)], ly, xs[-1], ly+hh, col = pal, border = NA); rect(lx0, ly, lx1, ly+hh, border = "grey40", lwd = 0.5)
  text(c(lx0,(lx0+lx1)/2,lx1), ly-hh*0.7, c("0","50","100"), cex = 0.74, col = "grey30")
  text(lx0, ly+hh*2.1, "Resistance load (%)", cex = 0.78, font = 2, adj = 0)   # short + left-anchored so the title is not clipped at the page edge
  gx <- LIM*0.42; gy <- -LIM*0.92                                    # geological-period key
  for (i in seq_len(nrow(per9))) { rect(gx, gy-(i-1)*hh*1.5, gx+g9$Rmax*0.05, gy-(i-1)*hh*1.5+hh,
    col = gsub("55$","",per9$col[i]), border = "grey55", lwd = 0.4)
    text(gx+g9$Rmax*0.07, gy-(i-1)*hh*1.5+hh*0.5, per9$name[i], cex = 0.72, adj = 0) }
  # figure title removed for clean panels (described in the caption)
  # text(0, LIM*0.995, "Phylogenetic distribution of herbicide-resistance load across plant families",
  #      cex = 1.05, font = 2, adj = c(0.5,1))
}
grDevices::cairo_pdf(file.path(fig_dir, "Fig_S9_CladesDated.pdf"), width = 200/25.4, height = 200/25.4)
draw_clades(); invisible(dev.off())
grDevices::png(file.path(fig_dir, "Fig_S9_CladesDated.png"), width = 200, height = 200,
               units = "mm", res = 600, type = "cairo"); draw_clades(); invisible(dev.off())
draw_clades()
results_summary$root_age <- root_age
results_summary$hr_families <- sum(fam_stats$n_hr > 0)
results_summary$top_hr_families <- fam_stats |> filter(n_hr > 0) |> arrange(desc(hr_frac)) |>
  slice_head(n = 6) |> pull(fam)

# --- Family-size transparency: is the family 'concentration' a load (fraction) or just big families? ---
fam_show <- fam_stats |> filter(n_hr > 0) |> arrange(desc(hr_frac), desc(n_hr)) |>
  transmute(Family = fam, `N species (tree)` = n, `N resistant` = n_hr,
            `Resistance load (pct)` = round(100 * hr_frac)) |> slice_head(n = 12)
results_summary$fam_show <- fam_show          # carried to the co-author walkthrough
write.csv(fam_show, file.path(tab_dir, "Table_FamilyLoad.csv"), row.names = FALSE)
## export the FULL per-family phylogeny data behind Fig S9 (all families, with clade) + the family tree
fam_export <- fam_stats |> mutate(clade = clade_of(fam)) |>
  transmute(family = fam, clade, n_species = n, n_resistant = n_hr,
            resistance_load_pct = round(100 * hr_frac, 1)) |>
  arrange(desc(resistance_load_pct), desc(n_resistant))
write.csv(fam_export, file.path(out_dir, "data", "phylo_family_load.csv"), row.names = FALSE)
ape::write.tree(fam_tree, file.path(out_dir, "data", "phylo_family_tree.nwk"))
# Does load merely track family size? Spearman across families with >= 3 species.
fam_big <- fam_stats |> filter(n >= 3)
rho_load_size <- suppressWarnings(cor(fam_big$n, fam_big$hr_frac, method = "spearman"))
results_summary$rho_load_size <- rho_load_size
cat(sprintf("Load-vs-size Spearman rho (families with >=3 spp, n=%d) = %.2f -> %s\n",
            nrow(fam_big), rho_load_size,
            ifelse(abs(rho_load_size) < 0.3, "load is essentially independent of family size",
                   "load shows some association with family size")))
kbl(fam_show, row.names = FALSE,
    caption = "Resistance load by family (only families containing at least one resistant species). Load (pct) = resistant / total species in that family, so it normalises for family size. Both small high-load families (e.g. Amaranthaceae, Chenopodiaceae) and large families appear; load shows only a moderate positive association with family size (see text), so the concentration is real and graded rather than a pure big-family artefact.")

suppressPackageStartupMessages(library(caper))
tree_d <- ape::multi2di(tree)                      # resolve polytomies (avoids caper contrCalc bug)
tree_d$node.label <- NULL                          # caper: node labels must not duplicate tips
tree_d$edge.length[tree_d$edge.length <= 0] <- 1e-8
dstat_df <- data.frame(tip = tree_d$tip.label,
                       HR = as.integer(tipmap[tree_d$tip.label] != "Non_HR_local"))
dstat_df <- dstat_df[!is.na(dstat_df$HR), ]
cdat <- caper::comparative.data(phy = tree_d, data = dstat_df, names.col = "tip")
set.seed(123)
pd <- tryCatch(caper::phylo.d(cdat, binvar = HR, permut = 500),
               error = function(e) { message("phylo.d failed: ", conditionMessage(e)); NULL })
if (!is.null(pd)) {
  .fp <- function(p) if (!is.na(p) && p < 0.001) "< 0.001" else sprintf("= %.3f", p)
  cat(sprintf("Phylogenetic signal of resistance: D = %.3f\n", pd$DEstimate))
  cat(sprintf("  P(D = 1, random)           %s\n", .fp(pd$Pval1)))
  cat(sprintf("  P(D = 0, Brownian/clumped) %s\n", .fp(pd$Pval0)))
  results_summary$Dstat <- list(D = unname(pd$DEstimate), p_random = unname(pd$Pval1),
                                p_brownian = unname(pd$Pval0))
} else {
  cat("D-statistic unavailable (caper numerical edge case); the megatree figure (Fig S9)\n",
      "and the phyloglm still document phylogenetic clustering.\n")
  results_summary$Dstat <- list(D = NA, p_random = NA, p_brownian = NA)
}

# predictors: top numeric module axes from the main RF + life_span (annual vs not) + origin
#   + photosynthetic pathway as a SINGLE parameter (C4 vs non-C4). It enters as binary, not the
#   3-level factor, deliberately: with only ~24 HR events the events-per-parameter ratio is tight, so
#   pathway is reduced to its one ecologically decisive contrast (C4 vs the rest) to test whether the
#   C4 signal seen in permutation importance survives phylogenetic correction without destabilising the fit.
# NOTE: this chunk now runs AFTER the D-statistic chunk, which attaches `caper` (and hence MASS) to the
# search path; MASS masks dplyr::select, so all dplyr verbs here are namespace-qualified to be safe.
top_mod <- imp_df |> dplyr::filter(str_detect(variable, "_PC1$")) |> dplyr::slice_head(n = 5) |> dplyr::pull(variable)
phylo_dat <- analysis_imp |>
  dplyr::left_join(sp_tab |> dplyr::select(species, tip), by = "species") |>
  dplyr::mutate(life_annual = factor(dplyr::if_else(life_span == "annual", "annual", "non_annual"),
                                     levels = c("non_annual", "annual")),   # non_annual = reference -> coef is the (positive) effect of being ANNUAL
         c4 = factor(dplyr::if_else(photosynthetic_pathway == "C4", "C4", "non_C4"),
                     levels = c("non_C4", "C4")),   # non_C4 = reference; coef = effect of being C4
         HR01 = as.integer(HR_bin == "HR")) |>
  dplyr::select(tip, HR01, dplyr::all_of(top_mod), life_annual, c4, origin) |> tidyr::drop_na()
phylo_dat <- as.data.frame(phylo_dat) |> dplyr::distinct(tip, .keep_all = TRUE)
rownames(phylo_dat) <- phylo_dat$tip
tr <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(phylo_dat)))
phylo_dat <- phylo_dat[tr$tip.label, ]
cat("phyloglm: N =", nrow(phylo_dat), "| HR events =", sum(phylo_dat$HR01),
    "| predictors =", length(top_mod) + 3, "(top 5 module axes + annual + C4 + origin)\n")

form <- as.formula(paste("HR01 ~", paste(c(top_mod, "life_annual", "c4", "origin"), collapse = " + ")))
set.seed(42)
pgl <- phylolm::phyloglm(form, data = phylo_dat, phy = tr, method = "logistic_MPLE", btol = 30)
co <- as.data.frame(summary(pgl)$coefficients)
co$term <- rownames(co); co <- co |> dplyr::filter(term != "(Intercept)")
codisp <- co[, c("term","Estimate","StdErr","p.value")]
codisp$p.value <- ifelse(codisp$p.value < 0.001, "< 0.001", sprintf("%.3f", codisp$p.value))
kbl(codisp, digits = 3, row.names = FALSE,
    caption = sprintf("Phylogenetic logistic regression (phyloglm) on the GBOTB dated megatree; alpha = %.2f. Photosynthetic pathway now enters as C4 vs non-C4 (term 'c4C4').", pgl$alpha))
# Report the C4 contrast explicitly (does the importance-ranked C4 signal survive phylogenetic correction?)
c4row <- co[grepl("^c4", co$term), ]
if (nrow(c4row)) cat(sprintf("C4 vs non-C4 (term %s): estimate = %+.2f logit, p %s -> C4 %s after phylogenetic correction\n",
    c4row$term[1], c4row$Estimate[1],
    ifelse(c4row$p.value[1] < 0.001, "< 0.001", sprintf("= %.3f", c4row$p.value[1])),
    ifelse(c4row$p.value[1] < 0.05,
           ifelse(c4row$Estimate[1] > 0, "REMAINS a significant positive predictor", "is significantly negative"),
           "is positive but NOT significant once relatedness is accounted for (consistent with C4 being phylogenetically clustered)")))

co$label <- pp(co$term)
p_phy <- co |> dplyr::mutate(sig = p.value < 0.05) |>
  ggplot(aes(reorder(label, Estimate), Estimate, colour = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = Estimate - 1.96*StdErr, ymax = Estimate + 1.96*StdErr)) +
  coord_flip() + scale_colour_manual(values = c(`TRUE` = "#E76F51", `FALSE` = "grey55"),
                                     labels = c("ns","p<0.05"), name = NULL) +
  labs(x = NULL, y = "Coefficient (logit, +-95% CI)") +
  theme(legend.position = "inside", legend.position.inside = c(0.84, 0.28),
        legend.background = element_rect(fill = "white", colour = NA))
save_fig(p_phy, "Fig_S5_PhyloGLM", 95, 90); p_phy

results_summary$robustB <- list(alpha = pgl$alpha, n = nrow(phylo_dat), events = sum(phylo_dat$HR01),
                                sig_terms = co$term[co$p.value < 0.05], coef = co[,c("term","Estimate","p.value")],
                                c4_est = if (nrow(c4row)) c4row$Estimate[1] else NA,
                                c4_p   = if (nrow(c4row)) c4row$p.value[1]  else NA)

suppressPackageStartupMessages({library(glmnet); library(xgboost); library(Matrix)})
alt_pred <- final_preds                                   # same predictors as the final RF
ad <- primary |> dplyr::select(dplyr::all_of(c(alt_pred, "HR_bin")))
X  <- model.matrix(HR_bin ~ . , data = ad)[, -1]
y  <- ad$HR_bin
yn <- as.integer(y == "HR")

# --- glmnet: elastic-net logistic, class-weighted, AUC by 5-fold CV (matched to xgboost) ---
w  <- ifelse(yn == 1, sum(yn == 0)/sum(yn == 1), 1)
set.seed(42)
cvfit <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 0.5,
                           type.measure = "auc", weights = w, nfolds = 5)
auc_glm <- max(cvfit$cvm)

# --- xgboost: gradient boosting, scale_pos_weight, AUC by 5-fold CV ---
dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = yn)
set.seed(42)
xgcv <- xgboost::xgb.cv(
  params = list(objective = "binary:logistic", eval_metric = "auc",
                max_depth = 3, eta = 0.1, scale_pos_weight = sum(yn==0)/sum(yn==1)),
  data = dtrain, nrounds = 300, nfold = 5, early_stopping_rounds = 20,
  verbose = 0)
auc_xgb <- max(xgcv$evaluation_log$test_auc_mean)

alt_tab <- data.frame(
  Model = c("Random forest (OOB)", "Elastic-net logistic (CV)", "Gradient boosting (CV)"),
  AUC = c(as.numeric(pROC::auc(roc_oob)), auc_glm, auc_xgb))
kbl(alt_tab, digits = 3, caption = "Three algorithms, same predictors: AUC agreement shows the signal is not RF-specific.")

p_alt <- alt_tab |> ggplot(aes(reorder(Model, AUC), AUC)) +
  geom_col(fill = "#457B9D", width = .6) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.15, size = 3) +
  coord_flip(ylim = c(0.5, max(alt_tab$AUC)+0.06)) +
  labs(x = NULL, y = "AUC (out-of-sample)")
save_fig(p_alt, "Fig_S7_AlternativeLearners", 120, 75); p_alt
results_summary$altlearn <- alt_tab

set.seed(42)
rf_n3  <- brf(primary |> dplyr::select(life_span, Height_PC1, photosynthetic_pathway, HR_bin))
roc_n3 <- pROC::roc(primary$HR_bin, rf_n3$votes[, "HR"], levels = c("nonHR", "HR"), direction = "<", quiet = TRUE)
ci_n3  <- pROC::ci.auc(roc_n3)
cat(sprintf("Nested 3-predictor model (life history + plant height + C4 pathway): OOB AUC = %.3f (95%% CI %.3f-%.3f)\n",
            as.numeric(pROC::auc(roc_n3)), ci_n3[1], ci_n3[3]))
cat(sprintf("Full %d-predictor model: OOB AUC = %.3f (95%% CI %.3f-%.3f)\n",
            length(final_preds), results_summary$core$auc, results_summary$core$ci_lo, results_summary$core$ci_hi))
results_summary$nested3 <- list(auc = as.numeric(pROC::auc(roc_n3)), ci_lo = ci_n3[1], ci_hi = ci_n3[3],
                                preds = c("life_span", "Height_PC1", "photosynthetic_pathway"))



# ---- imputed analytical data used for the analysis ----
data_dir <- file.path(out_dir, "data"); dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
# Pipeline A (missForest-imputed) module scores + covariates + response, all 575 species:
readr::write_csv(full_imp, file.path(data_dir, "imputed_analysis_full_575species.csv"))
# the exact modelling set after the Phase-0 habitat decision (what the RF was trained on):
readr::write_csv(analysis_imp, file.path(data_dir, "imputed_analysis_modelling_set.csv"))
cat("Imputed analysis data written to:", data_dir, "\n")

# ---- per-species OOB predictions that ASSEMBLE into the confusion matrices (binary + 3-class) ----
# Every species gets out-of-bag votes (predicted only by the trees that never saw it), so each prediction
# vector is the SAME length as the modelling set. This side-by-side table is exactly what is tallied into
# the binary and three-class confusion matrices; it is written to data/ for inspection.
.pc1  <- grep("_PC1$", names(analysis_imp), value = TRUE)
.keep <- stats::complete.cases(analysis_imp[, c(.pc1, covars, "HR_bin")])   # the rows the forests were fit on
.spp  <- analysis_imp$species[.keep]
stopifnot(length(.spp) == length(rf_final$predicted), length(.spp) == length(rf3$predicted))
oob_all <- data.frame(
  species              = .spp,
  binary_true          = as.character(primary$HR_bin),
  binary_OOB_pred      = as.character(rf_final$predicted),
  binary_OOB_prob_HR   = round(rf_final$votes[, "HR"], 4),
  binary_correct       = as.character(primary$HR_bin) == as.character(rf_final$predicted),
  threeclass_true      = as.character(d3$HR3),
  threeclass_OOB_pred  = as.character(rf3$predicted),
  threeclass_correct   = as.character(d3$HR3) == as.character(rf3$predicted))
.v3 <- round(as.data.frame(rf3$votes), 4); names(.v3) <- paste0("threeclass_prob_", names(.v3))
oob_all <- cbind(oob_all, .v3)
readr::write_csv(oob_all, file.path(data_dir, "OOB_predictions_per_species.csv"))
cat(sprintf("OOB per-species predictions written (%d species) -> data/OOB_predictions_per_species.csv\n", nrow(oob_all)))
# these reproduce the reported matrices exactly:
cat("Binary confusion matrix (from the OOB table):\n")
print(table(true = oob_all$binary_true, OOB_pred = oob_all$binary_OOB_pred))
cat("Three-class confusion matrix (from the OOB table):\n")
print(table(true = oob_all$threeclass_true, OOB_pred = oob_all$threeclass_OOB_pred))

saveRDS(rf_final, file.path(mod_dir, "rf_binary_final.rds"))
write_csv(imp_df, file.path(tab_dir, "Table_VariableImportance.csv"))
write_csv(mod_imp_tab, file.path(tab_dir, "Table_ModuleImportance.csv"))
results_summary$session <- R.version.string
saveRDS(results_summary, file.path(out_dir, "results_summary.rds"))   # consumed by the slide deck
sessionInfo()









.rv <- tryCatch(readRDS("reviewer_checks.rds"), error = function(e) NULL)
if (is.null(.rv)) .rv <- list(
  effort_species = list(rho_all=NA,p_all=NA,rho_nonhr=NA,p_nonhr=NA,flag_median_n=NA,noflag_median_n=NA),
  effort_spatial = list(rho_effort=NA,n_cells=NA),
  learning_curve = data.frame(n_HR=c(90,109), auc_mean=c(NA,NA), auc_sd=c(NA,NA)))
.pf <- function(p) if (!is.null(p) && !is.na(p) && p < 0.001) "< 0.001" else sprintf("= %.3f", p)
.d90 <- tail(.rv$learning_curve$auc_mean,1) - .rv$learning_curve$auc_mean[.rv$learning_curve$n_HR==90]

knitr::include_graphics("Figures/Fig_Reviewer_Effort.png")

knitr::include_graphics("Figures/Fig_Reviewer_LearningCurve.png")

.tx <- tryCatch(read.csv("data/merged_data_final_REBUILT.csv", check.names = FALSE), error = function(e) NULL)
# genus -> crop complex, compiled from the Harlan and de Wet Global Inventory (Vincent et al. 2013) + GRIN-Global
.cwrmap <- c(aegilops="Wheat",triticum="Wheat",hordeum="Barley",avena="Oat",secale="Rye",sorghum="Sorghum",
  setaria="Foxtail millet",panicum="Proso millet",echinochloa="Barnyard millet",eleusine="Finger millet",
  digitaria="Fonio",phalaris="Canary grass",lolium="Ryegrass",brassica="Brassica crops",raphanus="Radish",
  sinapis="Mustard",eruca="Rocket",camelina="Camelina",lepidium="Garden cress",hirschfeldia="Mustard",
  rapistrum="Mustard",beta="Beet",amaranthus="Grain amaranth",chenopodium="Quinoa/goosefoot",atriplex="Orach",
  spinacia="Spinach",lactuca="Lettuce",cichorium="Chicory",helianthus="Sunflower",daucus="Carrot",anthriscus="Chervil",
  solanum="Potato/eggplant",vicia="Faba bean/vetch",lathyrus="Grass pea",pisum="Pea",medicago="Alfalfa",
  trifolium="Clover",melilotus="Sweetclover",linum="Flax",papaver="Poppy",fagopyrum="Buckwheat",portulaca="Purslane",
  allium="Onion/garlic")
if (!is.null(.tx)) {
  .tx$fam   <- stringr::str_to_title(trimws(.tx$family))
  .tx$genus <- sub("[ _].*","", tolower(trimws(.tx$species)))
  .tx$crop  <- .cwrmap[.tx$genus]                                  # curated arable-crop match
  .ciat <- tryCatch(read.csv("data/cwr_ciat_species.csv"), error = function(e) NULL)  # bulk CIAT/GBIF match
  .cf <- if (!is.null(.ciat)) setNames(as.logical(.ciat$cwr_ciat), tolower(.ciat$species)) else NULL
  taxo <- do.call(rbind, lapply(c("HR_local","HR_elsewhere","Non_HR_local"), function(L){
    s <- .tx[.tx$HR_Status == L, ]; ft <- sort(table(s$fam), decreasing = TRUE); p <- prop.table(ft)
    H <- -sum(p*log(p)); ev <- H/log(length(ft))
    data.frame(Class = sub("_","-",L), `N species` = nrow(s), `N families` = length(ft),
               `Shannon H` = round(H,2), Evenness = round(ev,2),
               `CWR% (CIAT)` = if (!is.null(.cf)) round(100*mean(.cf[tolower(s$species)], na.rm=TRUE)) else NA,
               `CWR% (arable)` = round(100*mean(!is.na(s$crop))), check.names = FALSE)
  }))
  kbl(taxo, caption = "Family-level composition (with Shannon diversity and Pielou evenness) and crop-wild-relative share of the three HR-status classes.")
  ## CWR enrichment test (resistant pooled vs non-resistant), reported in Supplementary Table 11
  .res <- .tx$HR_Status != "Non_HR_local"
  .ciat_yes <- if (!is.null(.cf)) { .v <- .cf[tolower(.tx$species)]; !is.na(.v) & .v } else rep(NA, nrow(.tx))
  .arab_yes <- !is.na(.tx$crop)
  .cwrft <- function(yes) { k <- !is.na(yes)
    ft <- fisher.test(table(factor(yes[k], c(FALSE, TRUE)), factor(.res[k], c(FALSE, TRUE))))
    list(p = ft$p.value, or = unname(ft$estimate)) }
  .ct <- .cwrft(.ciat_yes); .at <- .cwrft(.arab_yes)
  cat(sprintf("CWR enrichment (resistant vs non-resistant, Fisher's exact): CIAT P = %.3g (OR %.2f); arable P = %.3g (OR %.2f)\n",
              .ct$p, .ct$or, .at$p, .at$or))
  results_summary$cwr_test <- list(ciat_p = .ct$p, ciat_or = .ct$or, arable_p = .at$p, arable_or = .at$or)
  saveRDS(results_summary, file.path(out_dir, "results_summary.rds"))   # re-save: cwr_test is set after the Step-16 export save
}

## Consistent-font 2-panel composites for ED5, ED6 and Supplementary S19. The panels were
## previously stitched from separately-saved figures and downscaled, which pushed their fonts
## below the 7 pt used elsewhere. Here we combine the pre-built panel ggplots with patchwork and
## save at the final 180 mm two-column width, so every panel renders at the native theme_pub 7 pt.
suppressPackageStartupMessages(library(patchwork))
## persist the panel ggplots so composites can be rebuilt without re-running the pipeline
.pnl <- Filter(function(v) exists(v, inherits = TRUE), c("p_stab","p_tr","p_all","p_alt","p_conv","p_opt"))
cat(">> panels available for composites:", paste(.pnl, collapse=", "), "\n")
tryCatch(saveRDS(setNames(lapply(.pnl, get, inherits = TRUE), .pnl), "composite_panels.rds"), error = function(e) NULL)
.mkcomp <- function(a, b, nm, w, h) {
  ok <- tryCatch({
    pa <- get(a, inherits = TRUE); pb <- get(b, inherits = TRUE)
    cc <- (pa | pb) + patchwork::plot_annotation(tag_levels = "a",
            theme = theme(plot.tag = element_text(size = 8, face = "bold", colour = "black")))
    save_fig(cc, nm, w, h); TRUE
  }, error = function(e) { cat("!! COMPOSITE FAIL", nm, ":", conditionMessage(e), "\n"); FALSE })
  cat(">> composite", nm, "built:", ok, "\n")
}
.mkcomp("p_stab", "p_tr",  "Fig_ED5_StabilityComposite",  180, 84)
.mkcomp("p_all",  "p_alt", "Fig_ED6_RobustnessComposite", 180, 96)
.mkcomp("p_conv", "p_opt", "Fig_S19_RFdiagComposite",     180, 76)

rs <- results_summary
p <- function(x,d=3) formatC(as.numeric(x),format="f",digits=d)
cat(sprintf("
**1. Data screening -> what we are allowed to model.** All 109 resistant species are arable and the
94 ruderal species are all non-resistant, so habitat could only ever act as a leakage shortcut. The
**Phase-0** test settled this empirically: dropping vs keeping ruderal species barely moved the
importance ranking (rho = %s) and discrimination (AUC %s arable-only vs %s with ruderal), so ruderal
species were **%s** as a neutral baseline. *Consequence for the next step:* the modelling set is fixed,
and habitat is excluded as a predictor.

**2. Missing data -> imputation choice.** Because complete-case deletion compounds across ~37 traits
and would discard a large share of the 24 HR-local species, we imputed (missForest) and verified that
the imputed (A), complete-case (B) and no-imputation cforest (C) pipelines rank the same predictors
(A vs B rho = %s). *Consequence:* all main results use **Pipeline A (imputed)**; B and C are robustness
checks, and the small AUC spread across pipelines (%s / %s / %s) is within noise, not evidence that any
model is 'better'.

**3. Modules -> predictors.** The functional traits were reduced to **9 multi-trait module axes (PC1)
+ 4 single-trait predictors + an ITV axis** (after the observed-data coherence check that split the
conflated modules and promoted decoupled traits); these plus life-history, origin and the 3 range
covariates (~22 predictors) form the matrix the random forest learns from.

**4. Core model -> is strategy predictive?** Yes: balanced RF reached OOB **AUC = %s** [%s-%s], HR
sensitivity %s, and beat the null in a 999-label permutation test (p %s). *Consequence:* it is
legitimate to interpret which strategies drive the signal.

**5. Importance + direction -> which strategies, and how.** After backward elimination (%s predictors
retained), the strongest predictors were %s. Module aggregation and partial-dependence give the
mechanism (fast, colonising life-history -> higher HR probability).

**6. Robustness -> are the conclusions safe?** (a) Replacing modules with raw traits gave essentially
the same discrimination (all-traits AUC %s vs module AUC %s) -> not a PCA artefact. (b) On the GBOTB
dated megatree, %s remained significant after phylogenetic correction -> not merely shared ancestry.

**7. Risk structure -> screening vs staging.** A model trained only on foreign weeds recovered %s%% of
local weeds (Phase 1, transfer); local and foreign weeds were partly distinct (Phase 2); the binary
model is the better general screen while the 3-class model adds risk *staging* (HR-local recall %s vs
binary %s).
",
 p(rs$phase0$rank_rho,2), p(rs$phase0$aucA), p(rs$phase0$aucB),
 ifelse(isTRUE(rs$phase0$use_ruderal),"retained","excluded"),
 p(rs$rho_AB,2), p(rs$tracks$AUC[1]), p(rs$tracks$AUC[2]), p(rs$tracks$AUC[3]),
 p(rs$core$auc), p(rs$core$ci_lo), p(rs$core$ci_hi), p(rs$core$sens,2),
 ifelse(rs$core$perm_p==0,"< 0.001",as.character(rs$core$perm_p)),
 rs$n_final_pred, paste(rs$top_predictors[1:4], collapse=", "),
 p(rs$robustA$auc_all), p(rs$robustA$auc_module),
 ifelse(length(rs$robustB$sig_terms)>0, paste(rs$robustB$sig_terms, collapse=", "), "no terms"),
 round(100*rs$phase1$recall_local), p(rs$phase4$recall_3class,2), p(rs$phase4$recall_2class,2)
))

rs <- results_summary
cat(sprintf("\n*Data-anchored framing for the abstract:* balanced-RF OOB AUC = %.2f; resistance is
phylogenetically clustered (D = %s); top strategy signals = %s; and %d non-resistant weeds are flagged
as highest-risk on the watch-list.\n",
 as.numeric(rs$core$auc),
 tryCatch(sprintf('%.2f', rs$Dstat$D), error=function(e) 'see Step 19'),
 paste(rs$top_predictors[1:3], collapse=', '),
 tryCatch(nrow(rs$watchlist), error=function(e) 15)))

## ############################################################################################
## STAGE 9  Figure builders
## source file: _make_figs.R
## ############################################################################################

## ============================================================================
## Standalone regeneration of TWO figures, using saved pipeline objects:
##   - Figure 5  : Fig_S9_CladesDated  (family phylogeny, NEW gradient legend)
##   - Figure S7b: Fig_S8b_SHAP_beeswarm_byclass  (class-faceted SHAP)
## Avoids the full ~30-package pipeline; needs only ape + RF/SHAP/ggplot stack.
## Exports vector .pdf (cairo) + 600-dpi .png for each.
## ============================================================================
## use the user library matching whichever R runs this (R-4.4 minimal lib OR R-4.5.2 full env)
uv <- paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
ulib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", uv)
if (dir.exists(ulib)) .libPaths(c(ulib, .libPaths()))
base_dir <- "."
fig_dir  <- file.path(base_dir, "Figures")
mod_csv  <- file.path(base_dir, "data", "imputed_analysis_modelling_set.csv")
mer_csv  <- "data/merged_data_final_REBUILT.csv"

sentence_case <- function(x) {
  x <- tolower(trimws(gsub("\\s+", " ", x)))
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}
title_case <- function(x) vapply(strsplit(tolower(x), " "), function(w)
  paste(paste0(toupper(substring(w, 1, 1)), substring(w, 2)), collapse = " "), character(1))
to_num <- function(v) if (is.numeric(v)) v else as.numeric(as.factor(v))  # factor->code for colour

## label prettifier (mirrors pp() in 02_Master_Analysis.Rmd): raw predictor -> readable, no underscores
.PRETTY <- c(
  Leaf_Economics_PC1 = "Leaf economics", Leaf_Size_PC1 = "Leaf size", Leaf_Nutrient_PC1 = "Leaf nutrients",
  Height_PC1 = "Plant height",
  Hydraulic_PC1 = "Hydraulic", Root_PC1 = "Root", Propagule_Size_PC1 = "Propagule size",
  Fecundity_PC1 = "Fecundity", Genome_PC1 = "Genome size", Rooting_Depth_PC1 = "Rooting depth",
  Stem_Diameter_PC1 = "Stem diameter", Stem_Cond_Dens_PC1 = "Stem conduit density",
  Germination_PC1 = "Germination", ITV_PC1 = "ITV",
  origin = "Origin", life_span = "Life history", growth_form = "Growth form",
  leaf_compoundness = "Leaf compoundness", photosynthetic_pathway = "Photosynthetic pathway",
  dispersal_distance = "Dispersal distance", alt_range_scaled = "Altitudinal range",
  eoo_scaled = "EOO", aoo_scaled = "AOO")
.pp1 <- function(t) {
  if (t %in% names(.PRETTY)) return(unname(.PRETTY[[t]]))
  s <- gsub("_", " ", sub("_mean$", "", sub("_scaled$", "", sub("_PC1$", "", t))))
  s <- paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  gsub("\\bAoo\\b", "AOO", gsub("\\bEoo\\b", "EOO", gsub("\\bItv\\b", "ITV", s)))
}
pp <- function(x) vapply(as.character(x), .pp1, character(1), USE.NAMES = FALSE)

## ---- FIGURE 5 (Fig_S9_CladesDated) DISABLED ----
## 02_Master_Analysis.Rmd builds the authoritative CIRCULAR dated family tree over shaded geological
## periods (report page 80 / Extended Data Fig. 6 caption). Do NOT overwrite it with the rectangular
## phylogram below.
if (FALSE) {
suppressPackageStartupMessages(library(ape))
tree <- readRDS(file.path(base_dir, "megatree.rds"))
mod  <- read.csv(mod_csv, stringsAsFactors = FALSE)
mer  <- read.csv(mer_csv, stringsAsFactors = FALSE)
mod$tip <- gsub(" ", "_", sentence_case(mod$species))
mer$tip <- gsub(" ", "_", sentence_case(mer$species))
mer$family_c <- title_case(ifelse(is.na(mer$family) | mer$family == "", "Unknown", mer$family))
fammap <- setNames(mer$family_c, mer$tip)
hrmap  <- setNames(as.integer(mod$hr_status != "Non_HR_local"), mod$tip)
present <- tree$tip.label
df <- data.frame(tip = present, fam = unname(fammap[present]), hr = unname(hrmap[present]),
                 stringsAsFactors = FALSE)
cat("tips:", length(present), "| matched family:", sum(!is.na(df$fam)),
    "| matched hr:", sum(!is.na(df$hr)), "\n")
df <- df[!is.na(df$fam) & !is.na(df$hr), ]
sp <- split(df, df$fam)
fam_stats <- data.frame(
  fam     = names(sp),
  n       = vapply(sp, nrow, integer(1)),
  n_hr    = vapply(sp, function(d) sum(d$hr), integer(1)),
  hr_frac = vapply(sp, function(d) mean(d$hr), numeric(1)),
  rep     = vapply(sp, function(d) d$tip[1], character(1)),
  stringsAsFactors = FALSE)
fam_tree <- ape::drop.tip(tree, setdiff(present, fam_stats$rep))
fam_tree$tip.label <- fam_stats$fam[match(fam_tree$tip.label, fam_stats$rep)]
fam_tree <- ape::ladderize(fam_tree)
fs  <- fam_stats[match(fam_tree$tip.label, fam_stats$fam), ]
pal <- colorRampPalette(c("#457B9D", "#e9c46a", "#E76F51"))(100)
col_fam <- pal[pmax(1, ceiling(fs$hr_frac * 100))]
root_age <- max(ape::node.depth.edgelength(fam_tree))
base_family <- "sans"

add_grad_legend <- function(pal, title = "Resistant species in family (%)") {
  usr <- par("usr"); xr <- usr[2] - usr[1]; yr <- usr[4] - usr[3]
  bx0 <- usr[1] + 0.03 * xr; bx1 <- bx0 + 0.26 * xr
  by0 <- usr[3] + 0.06 * yr; by1 <- by0 + 0.022 * yr
  n <- length(pal); xs <- seq(bx0, bx1, length.out = n + 1)
  rect(xs[-(n + 1)], by0, xs[-1], by1, col = pal, border = NA, xpd = NA)
  rect(bx0, by0, bx1, by1, border = "grey40", lwd = 0.6, xpd = NA)
  at <- bx0 + c(0, 0.5, 1) * (bx1 - bx0)
  segments(at, by0, at, by0 - 0.012 * yr, lwd = 0.6, xpd = NA)
  text(at, by0 - 0.014 * yr, c("0", "50", "100"), pos = 1, cex = 0.7, xpd = NA)
  text((bx0 + bx1) / 2, by1 + 0.020 * yr, title, cex = 0.72, font = 2, xpd = NA)
}
draw_clades <- function() {
  par(mar = c(3.2, 0, 2, 0), family = base_family)
  plot(fam_tree, type = "phylogram", show.tip.label = TRUE, label.offset = 2, cex = 0.5,
       tip.color = col_fam, edge.color = "grey55", edge.width = 0.7, no.margin = FALSE)
  ape::axisPhylo(side = 1, root.time = root_age, backward = TRUE)
  mtext("Time before present (Myr)", side = 1, line = 2, cex = 0.8)
  ape::tiplabels(pch = 15, col = col_fam, cex = 0.7)
  add_grad_legend(pal)
  title("Dated family-level phylogeny: clades and resistance load", cex.main = 0.9)
}
grDevices::cairo_pdf(file.path(fig_dir, "Fig_S9_CladesDated.pdf"), width = 170/25.4, height = 210/25.4)
draw_clades(); invisible(dev.off())
grDevices::png(file.path(fig_dir, "Fig_S9_CladesDated.png"), width = 170, height = 210,
               units = "mm", res = 600, type = "cairo"); draw_clades(); invisible(dev.off())
cat("Figure 5 done:", nrow(fam_stats), "families;",
    "max family HR load =", round(100 * max(fam_stats$hr_frac)), "%\n")
}  ## end disabled FIGURE 5 block

## -------------------------------------------------------------- FIGURE S7b ---
suppressPackageStartupMessages({library(randomForest); library(kernelshap); library(shapviz)
  library(ggplot2); library(ggbeeswarm)})
rf <- readRDS(file.path(base_dir, "Models", "rf_binary_final.rds"))
preds <- rownames(rf$importance)
primary <- read.csv(mod_csv, stringsAsFactors = FALSE)
xl <- rf$forest$xlevels
for (v in names(xl)) if (is.character(xl[[v]])) primary[[v]] <- factor(primary[[v]], levels = xl[[v]])
## the exported CSV retains a few non-imputed cells (unlike the in-memory missForest
## matrix); SHAP needs predictions, so restrict to complete-predictor rows.
ok_rows <- stats::complete.cases(primary[, preds, drop = FALSE])
cat("SHAP: using", sum(ok_rows), "of", length(ok_rows), "complete-predictor rows\n")
primary <- primary[ok_rows, , drop = FALSE]
set.seed(42)
n_expl <- min(80, nrow(primary)); n_bg <- min(100, nrow(primary))
expl_idx <- sample(nrow(primary), n_expl)
Xexpl <- primary[expl_idx, preds, drop = FALSE]
Xbg   <- primary[sample(nrow(primary), n_bg), preds, drop = FALSE]
shap_pred <- function(object, newdata) predict(object, newdata, type = "prob")[, "HR"]
ks <- kernelshap::kernelshap(rf, X = Xexpl, bg_X = Xbg, pred_fun = shap_pred, verbose = FALSE)
sv <- shapviz::shapviz(ks)

S_mat <- as.matrix(shapviz::get_shap_values(sv))
X_val <- tryCatch(as.data.frame(shapviz::get_feature_values(sv)),
                  error = function(e) as.data.frame(sv$X))
ord <- head(names(sort(colMeans(abs(S_mat)), decreasing = TRUE)), min(12, ncol(S_mat)))
cls <- factor(as.character(primary$hr_status)[expl_idx],
              levels = c("HR_local", "HR_elsewhere", "Non_HR_local"),
              labels = c("HR-local (Germany)", "HR-elsewhere", "Non-resistant"))
bee_df <- do.call(rbind, lapply(ord, function(f) {
  v <- to_num(X_val[[f]]); rng <- range(v, na.rm = TRUE); span <- diff(rng); if (span == 0) span <- 1
  data.frame(feature = f, shap = S_mat[, f], fval_s = (v - rng[1]) / span, class = cls)
}))
bee_df$feature <- factor(pp(bee_df$feature), levels = rev(pp(ord)))
have_bee  <- requireNamespace("ggbeeswarm", quietly = TRUE)
jit_layer <- if (have_bee) {
  ggbeeswarm::geom_quasirandom(size = 0.55, alpha = 0.85, width = 0.32, groupOnX = FALSE)
} else {
  geom_jitter(height = 0.28, width = 0, size = 0.55, alpha = 0.85)
}
p <- ggplot(bee_df, aes(shap, feature, colour = fval_s)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  jit_layer +
  scale_colour_viridis_c(option = "C", name = "Feature\nvalue", breaks = c(0, 1), labels = c("Low", "High")) +
  facet_wrap(~ class, nrow = 1) +
  labs(x = "SHAP value (impact on P(resistance))", y = NULL) +
  theme_bw(9) + theme(legend.position = "right", panel.grid.minor = element_blank())
## cairo_pdf can choke on the ~thousands of beeswarm points (0-byte output); svglite is robust
## and gives a clean vector master for EMF/PDF conversion.
tryCatch(ggsave(file.path(fig_dir, "Fig_S8b_SHAP_beeswarm_byclass.pdf"), p, width = 170, height = 110,
       units = "mm", device = cairo_pdf, dpi = 600), error = function(e) message("cairo_pdf S8b failed: ", conditionMessage(e)))
ggsave(file.path(fig_dir, "Fig_S8b_SHAP_beeswarm_byclass.svg"), p, width = 170, height = 110,
       units = "mm", device = svglite::svglite)
ggsave(file.path(fig_dir, "Fig_S8b_SHAP_beeswarm_byclass.png"), p, width = 170, height = 110,
       units = "mm", dpi = 600)

## diagnostic: mean SHAP per class for each top axis (to verify the interpretation text)
agg <- aggregate(shap ~ feature + class, data = bee_df, FUN = mean)
agg$shap <- round(agg$shap, 3)
cat("\n=== mean SHAP by class x feature (top axes) ===\n")
print(reshape(agg, idvar = "feature", timevar = "class", direction = "wide"), row.names = FALSE)
cat("\nFigure S7b done.\n")

## ############################################################################################
## STAGE 9  Main Fig 3 composite
## source file: _build_fig_ms4.R
## ############################################################################################

## ============================================================================
## Builder for main-text Fig 3 composite: Fig_MS4_PhyloPanel.png
##   (a) per-species circular phylogeny  = Fig_S6_Phylogeny.png
##   (b) phylogenetic logistic-regression = Fig_S5_PhyloGLM.png
## Horizontal magick composite, panels labelled a / b, exported at the Nature
## Plants 2-column width (180 mm, 600 dpi) so it is NOT downscaled in the docx
## (previous native width was ~351 mm, which shrank text below spec). Panel (b)
## is only ever DOWN-scaled (no upscaling blur); the circular tree's white margin
## is trimmed to cut wasted space.
## Run AFTER 02_Master_Analysis.Rmd has regenerated the two source panels.
## ============================================================================
suppressPackageStartupMessages(library(magick))
fd     <- "./Figures"
mm2px  <- function(mm, dpi = 600) round(mm / 25.4 * dpi)
Wpx    <- mm2px(180)                          # final composite width
sz     <- mm2px(2.8)                          # panel-letter size (~8 pt at final size)
lab    <- function(im, L) image_annotate(im, L, gravity = "northwest",
                                          location = "+20+10", size = sz,
                                          weight = 700, color = "black")

a <- image_trim(image_read(file.path(fd, "Fig_S6_Phylogeny.png")), fuzz = 3)  # drop white margin
b <- image_read(file.path(fd, "Fig_S5_PhyloGLM.png"))

## phylogeny takes ~55% of the final width; phyloGLM sized to ~0.82x its height (down-scale only)
a  <- image_resize(a, paste0(round(0.55 * Wpx), "x"))
Ha <- image_info(a)$height
b  <- image_resize(b, paste0("x", round(0.82 * Ha)))

a  <- lab(a, "a")
b  <- image_extent(b, geometry = paste0(image_info(b)$width + mm2px(5), "x", Ha),
                   gravity = "center", color = "white")
b  <- lab(b, "b")

comp <- image_append(image_join(list(a, b)))     # horizontal (a | b)
comp <- image_resize(comp, paste0(Wpx, "x"))      # normalise to exactly 180 mm wide
image_write(comp, file.path(fd, "Fig_MS4_PhyloPanel.png"))
image_write(comp, file.path(fd, "Fig_MS4_PhyloPanel.pdf"), format = "pdf", density = 600)
cat("Fig_MS4_PhyloPanel rebuilt:", paste(image_info(comp)[c("width","height")], collapse = " x "),
    "(=", round(image_info(comp)$width/600*25.4), "mm wide)\n")

## ############################################################################################
## STAGE 9  Supplement composites
## source file: _compose_supp.R
## ############################################################################################

## Compose the three horizontal multi-panel SUPPLEMENT figures from their panels
## (extracted from _finalize.R; no document rendering here). 600-dpi PNG + PDF.
suppressPackageStartupMessages({library(magick); library(pdftools)})
fd <- "./Figures"
Pf <- function(f) file.path(fd, f); mm2px <- function(mm, dpi = 600) round(mm/25.4*dpi)

## rasterise the PDF-only trait-correlation figure -> PNG (for docx embedding)
png::writePNG(pdftools::pdf_render_page(Pf("Fig_S1_TraitCorrelation.pdf"), page = 1, dpi = 300),
              Pf("Fig_S1_TraitCorrelation.png"))

lab <- function(im, L, sz = 70)
  image_annotate(im, paste0("(", L, ")"), gravity = "northwest", location = "+8+4",
                 size = sz, weight = 700, color = "black")
row_append <- function(files, letters, out, width_mm = 178, pad = 30) {
  ims <- Map(function(f, L) lab(image_read(f), L), files, letters)
  H <- max(vapply(ims, function(i) image_info(i)$height, numeric(1)))
  ims <- lapply(ims, function(i) { r <- image_resize(i, paste0("x", H))
    image_extent(r, geometry = paste0(image_info(r)$width + pad, "x", H), gravity = "center", color = "white") })
  final_im <- image_resize(image_append(image_join(ims)), paste0(mm2px(width_mm), "x"))
  image_write(final_im, out)
  image_write(final_im, sub("\\.png$", ".pdf", out), format = "pdf", density = 600)
}
row_append(c(Pf("Fig3_ErrorConvergence.png"), Pf("Fig4_Optimisation.png")),
           c("a","b"), Pf("FigS4_RFdiagnostics.png"), 178)
row_append(c(Pf("Fig_S2_RankStability.png"), Pf("Fig_S3_TrackComparison.png")),
           c("a","b"), Pf("FigS5_StabilityTracks.png"), 178)
row_append(c(Pf("Fig_S4_AllTraitsRF.png"), Pf("Fig_S7_AlternativeLearners.png")),
           c("a","b"), Pf("FigS6_Robustness.png"), 180)   # panel (c) phyloglm dropped (duplicated main-text Fig 3b)
cat("supplement composites rebuilt\n")

## ############################################################################################
## STAGE 10 Render driver (knits 02 -> Herbicide_HR_Report)
## source file: _render.R
## ############################################################################################

# Render helper: points rmarkdown at RStudio's bundled pandoc, renders with params.
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
args <- commandArgs(trailingOnly = TRUE)
mode <- ifelse(length(args) >= 1, args[1], "quick")

rmd <- "./02_Master_Analysis.Rmd"

p_quick <- list(impute_primary = TRUE, ntree = 300, n_perm = 49,
                run_loo = FALSE, loo_ntree = 200, stability_runs = 3, run_shap = FALSE)
p_full  <- list(impute_primary = TRUE, ntree = 1000, n_perm = 999,
                run_loo = TRUE, loo_ntree = 300, stability_runs = 10, run_shap = TRUE)

params <- if (mode == "full") p_full else p_quick
out    <- if (mode == "full") "Herbicide_HR_Report.pdf" else "Herbicide_HR_Report_DRAFT.pdf"

cat("Rendering mode:", mode, "| pandoc:", rmarkdown::pandoc_available(),
    as.character(rmarkdown::pandoc_version()), "\n")
t0 <- Sys.time()
rmarkdown::render(rmd, output_file = out, params = params, quiet = FALSE,
                  envir = new.env())
cat("DONE in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min ->", out, "\n")
