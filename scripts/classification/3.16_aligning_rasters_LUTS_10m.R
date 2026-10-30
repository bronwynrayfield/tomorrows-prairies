
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS TO REFERENCE                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 3: align all grassland layers to classification study area file from step 1
#         to prep for analysis

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# AAFC LUTS directory
luts_tile_root <- file.path(raw_grass_data_dir, "Landcover/AAFC Land Use 2020")
luts_merged_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/merged")
luts_reproj_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/reproj")
luts_cropped_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/cropped")
luts_aligned_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/aligned")
luts_pilot_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/pilot")
luts_masked_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_10m_v1/masked")

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# LUTS rasters
luts_w_file <- file.path(luts_tile_root, "AAFC_10m_Land_Use_2020_West_version1_30Apr2025/LU2020_10m_u07to15__v1_30Apr2025.tif")
luts_e_file <- file.path(luts_tile_root, "AAFC_10m_Land_Use_2020_East_provisional_2May2025/LU2020_10m_u16to22__provisional_2May2025.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
res(sa)
crs(sa)
sa

####--- ALGIN RASTERS ---#### 

#----- 2020 LUTS 10m  ----#

# Load rasters
r_w <- rast(luts_w_file)
r_e <- rast(luts_e_file)

plot(r_w) # warning: [plot] unknown categories in raster values 
plot(r_e)

origin(r_w) # 0 0
origin(r_e) # 0 0

cats(r_w)
cats(r_e)

r_w
#class       : SpatRaster 
#dimensions  : 201397, 313221, 1  (nrow, ncol, nlyr)
#resolution  : 10, 10  (x, y)
#extent      : -2685510, 446700, 891220, 2905190  (xmin, xmax, ymin, ymax)
#coord. ref. : Albers_Conic_Equal_Area 
#source      : LU2020_10m_u07to15__v1_30Apr2025.tif 
#categories  : LU2020, Definition, IPCC_Class, EUDR 
#name        :     LU2020 
#min value   : Settlement 
#max value   :         NA

r_e
#class       : SpatRaster 
#dimensions  : 244024, 293165, 1  (nrow, ncol, nlyr)
#resolution  : 10, 10  (x, y)
#extent      : 337570, 3269220, 264140, 2704380  (xmin, xmax, ymin, ymax)
#coord. ref. : Albers_Conic_Equal_Area 
#source      : LU2020_10m_u16to22__provisional_2May2025.tif 
#color table : 1 
#categories  : LU2020, Definition, IPCC_Class 
#name        :       LU2020 
#min value   :   Settlement 
#max value   : Snow and Ice 

# print frequency table
freq_table_w <- freq(r_w) 
print(freq_table_w) 

#   layer                                                     value       count
#1      1                                                Settlement    57360993
#2      1                               High Reflectance Settlement    47449753
#3      1                                         Settlement Forest    52200980
#4      1                                      Vegetated Settlement   101069297
#5      1                          Very High Reflectance Settlement    30707669
#6      1                                                     Water  4568437305
#7      1                                                    Forest 12295927102
#8      1                                            Forest Wetland  3621175785
#9      1               Forest Regenerating after Harvest <20 years   797426784
#10     1       Forest Wetland Regenerating after Harvest <20 years    16348322
#11     1             Forest Regenerating after Harvest 20-29 years   347235606
#12     1     Forest Wetland Regenerating after Harvest 20-29 years     5362631
#13     1                  Forest Regenerating after Fire <20 years  1640725146
#14     1                                                  Cropland  1190327676
#15     1                                           Annual Cropland  2897346403
#16     1                                Land Converted to Cropland    18439582
#17     1                         Land Converted to Annual Cropland     8985993
#18     1                                         Prairie Grassland   722442363
#19     1                                       Temperate Grassland    84295660
#20     1                            Temperate Rainforest Grassland    13539964
#21     1                                          Boreal Grassland   190226888
#22     1                                         Montane Grassland   213431867
#23     1                                   Alpine Tundra Grassland   297352163
#24     1                                   Arctic Tundra Grassland     1218187
#25     1                                           Prairie Wetland   244768534
#26     1                                         Temperate Wetland    44168022
#27     1                              Temperate Rainforest Wetland    12291126
#28     1                                            Boreal Wetland  2467682879
#29     1                                     Arctic Tundra Wetland    12033960
#30     1                       Newly-Detected Settlement <10 years     6419321
#31     1      Newly-Detected High Reflectance Settlement <10 years    10040098
#32     1                Newly-Detected Settlement Forest <10 years      793684
#33     1             Newly-Detected Vegetated Settlement <10 years    13475599
#34     1 Newly-Detected Very High Reflectance Settlement <10 years     3110866
#35     1                                                Other Land   896835222
#36     1                                              Snow and Ice   276783518
#37     1                                                      <NA> 29874332789

freq_table_e <- freq(r_e) 
print(freq_table_e)

#   layer                                                     value       count
#1      1                                                Settlement    56116437
#2      1                               High Reflectance Settlement    23333417
#3      1                                         Settlement Forest   105274706
#4      1                                      Vegetated Settlement   138216592
#5      1                          Very High Reflectance Settlement    10492299
#6      1                                                     Water 14187349212
#7      1                                                    Forest 12217405752
#8      1                                            Forest Wetland  2790931819
#9      1               Forest Regenerating after Harvest <20 years   728514836
#10     1       Forest Wetland Regenerating after Harvest <20 years    23945618
#11     1             Forest Regenerating after Harvest 20-29 years   506497196
#12     1     Forest Wetland Regenerating after Harvest 20-29 years     7401602
#13     1                  Forest Regenerating after Fire <20 years   760161590
#14     1                                                  Cropland   350699911
#15     1                                           Annual Cropland   401538605
#16     1                                Land Converted to Cropland     2853738
#17     1                         Land Converted to Annual Cropland     1359762
#18     1                                       Temperate Grassland     8823457
#19     1                                          Boreal Grassland  1637985847
#20     1                                   Alpine Tundra Grassland    24422484
#21     1                                   Arctic Tundra Grassland   633391300
#22     1                                         Temperate Wetland    95280265
#23     1                                            Boreal Wetland  2920681167
#24     1                                     Alpine Tundra Wetland    35756380
#25     1                                     Arctic Tundra Wetland   114222070
#26     1                       Newly-Detected Settlement <10 years     6853463
#27     1      Newly-Detected High Reflectance Settlement <10 years     3817349
#28     1                Newly-Detected Settlement Forest <10 years     5155317
#29     1             Newly-Detected Vegetated Settlement <10 years    14444266
#30     1 Newly-Detected Very High Reflectance Settlement <10 years     1993687
#31     1                                                Other Land   212166362
#32     1                                              Snow and Ice     1601217

# Set pixels with value 0 to NA - NO! need to preserve zeroes
#r_ab[r_ab == 0] <- NA
#r_mb[r_mb == 0] <- NA
#r_sk[r_sk == 0] <- NA

# Basic checks
crs(r_w) == crs(r_e)                     # TRUE; only checks string equality
crs(r_w, proj=TRUE)                      # print full WKT to eyeball differences (SAME)
crs(r_e, proj=TRUE)

res(r_w); res(r_e)                       # should both be 10 10
origin(r_w); origin(r_e)                 # MUST be identical for seamless merge - both 0,0
nlyr(r_w); nlyr(r_e)

# Check if different
compareGeom(r_w, r_e, stopOnError = FALSE) # FALSE

# Check nodata (if set inconsistently, can also get visible seams)
NAflag(r_w); NAflag(r_e)

############ 

### align each raster separately, then merge

luts10_out_path <- file.path(luts_reproj_dir, "luts10_e_reproj.tif")

# Conservative and robust GDAL options
wopt <- list(
  datatype = "INT1U",            # or "INT1U"/"INT2U" 
  gdal     = c(
    "COMPRESS=LZW",
    "PREDICTOR=2",
    "TILED=YES",
    "BIGTIFF=YES",
    "BLOCKXSIZE=512",
    "BLOCKYSIZE=512"
  )
)

# Reproject to 'sa' with nearest neighbor - try without alignment to avoid memory issues
luts_e_reproj <- project(r_e, sa, method = "near", #align = TRUE,
                         filename = luts10_out_path, overwrite = TRUE, wopt = wopt)
plot(luts_e_reproj)

# Ensure both are now on the same CRS (sa's)
stopifnot(same.crs(luts_e_reproj, sa))

# Force exact grid alignment to sa
luts_e_aligned <- resample(
  luts_e_reproj, sa, method = "near",
  filename = file.path(luts_aligned_dir, "luts10_e_aligned.tif"),
  overwrite = TRUE, wopt = wopt
)

plot(luts_e_aligned)

# same for western raster

# Reproject to 'sa' with nearest neighbor - try without alignment to avoid memory issues
luts_w_reproj <- project(r_w, sa, method = "near", #align = TRUE,
                         filename = file.path(luts_reproj_dir, "luts10_w_reproj.tif"), 
                         overwrite = TRUE, wopt = wopt)
plot(luts_w_reproj)

# Ensure both are now on the same CRS (sa's)
stopifnot(same.crs(luts_w_reproj, sa))

# Force exact grid alignment to sa
luts_w_aligned <- resample(
  luts_w_reproj, sa, method = "near",
  filename = file.path(luts_aligned_dir, "luts10_w_aligned.tif"),
  overwrite = TRUE, wopt = wopt
)

plot(luts_w_aligned)

luts_e_aligned
#class       : SpatRaster 
#dimensions  : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : luts10_e_aligned.tif 
#color table : 1 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : LU2020, Definition, IPCC_Class 
#name        :       LU2020 
#min value   :   Settlement 
#max value   : Snow and Ice 

luts_w_aligned
#class       : SpatRaster 
#dimensions  : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : luts10_w_aligned.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : LU2020, Definition, IPCC_Class, EUDR 
#name        : LU2020 
#min value   :     NA 
#max value   :     NA 

freq_table_e <- freq(luts_e_aligned) 
print(freq_table_e)

#   layer                                                     value     count
#1      1                                                Settlement   3457101  21
#2      1                               High Reflectance Settlement   1533668  22
#3      1                                         Settlement Forest   5355388  24
#4      1                                      Vegetated Settlement   7593490  28
#5      1                          Very High Reflectance Settlement    656847  29
#6      1                                                     Water 833815545  31 
#7      1                                                    Forest 619265365  41
#8      1                                            Forest Wetland 244114677  42
#9      1               Forest Regenerating after Harvest <20 years  35127417  43
#10     1       Forest Wetland Regenerating after Harvest <20 years   1779158  44
#11     1             Forest Regenerating after Harvest 20-29 years  31865497  47
#12     1     Forest Wetland Regenerating after Harvest 20-29 years    571760  48
#13     1                  Forest Regenerating after Fire <20 years  56583335  49
#14     1                                                  Cropland  22846422  51
#15     1                                           Annual Cropland  31599704  52
#16     1                                Land Converted to Cropland    163640  55
#17     1                         Land Converted to Annual Cropland     71312  56
#18     1                                       Temperate Grassland    560880  64
#19     1                                          Boreal Grassland  97565140  66
#20     1                                   Arctic Tundra Grassland  61925199  69
#21     1                                         Temperate Wetland   5316897  74
#22     1                                            Boreal Wetland 243803964  76
#23     1                                     Arctic Tundra Wetland  11154165  79
#24     1                       Newly-Detected Settlement <10 years    439394  81
#25     1      Newly-Detected High Reflectance Settlement <10 years    286426  82
#26     1                Newly-Detected Settlement Forest <10 years    318315  84
#27     1             Newly-Detected Vegetated Settlement <10 years   1092295  88
#28     1 Newly-Detected Very High Reflectance Settlement <10 years    167707  89
#29     1                                                Other Land   6169797  91
#30     1                                              Snow and Ice        61  92

freq_table_w <- freq(luts_w_aligned) 
print(freq_table_w)

#   layer                                                     value      count
#1      1                                                      <NA> 5254292037
#2      1                                                Settlement    6375289  21
#3      1                               High Reflectance Settlement    5273028  22
#4      1                                         Settlement Forest    5799402  24
#5      1                                      Vegetated Settlement   11230136  28
#6      1                          Very High Reflectance Settlement    3412015  29
#7      1                                                     Water  507361918  31
#8      1                                                    Forest 1366215351  41
#9      1                                            Forest Wetland  402353283  42
#10     1               Forest Regenerating after Harvest <20 years   88602313  43
#11     1       Forest Wetland Regenerating after Harvest <20 years    1816966  44
#12     1             Forest Regenerating after Harvest 20-29 years   38579815  47
#13     1     Forest Wetland Regenerating after Harvest 20-29 years     595685  48
#14     1                  Forest Regenerating after Fire <20 years  182299150  49
#15     1                                                  Cropland  132258946  51
#16     1                                           Annual Cropland  321923912  52
#17     1                                Land Converted to Cropland    2047724  55
#18     1                         Land Converted to Annual Cropland     998316  56
#19     1                                         Prairie Grassland   80274323  63
#20     1                                       Temperate Grassland    9365631  64
#21     1                            Temperate Rainforest Grassland    1505217  65
#22     1                                          Boreal Grassland   21134759  66
#23     1                                         Montane Grassland   23714987  67
#24     1                                   Alpine Tundra Grassland   33037024  68
#25     1                                   Arctic Tundra Grassland     135304  69
#26     1                                           Prairie Wetland   27196780  73
#27     1                                         Temperate Wetland    4907608  74
#28     1                              Temperate Rainforest Wetland    1365985  75
#29     1                                            Boreal Wetland  274191084  76
#30     1                                     Arctic Tundra Wetland    1337406  78
#31     1                       Newly-Detected Settlement <10 years     713336  79
#32     1      Newly-Detected High Reflectance Settlement <10 years    1115506  81
#33     1                Newly-Detected Settlement Forest <10 years      88232  82
#34     1             Newly-Detected Vegetated Settlement <10 years    1496887  88
#35     1 Newly-Detected Very High Reflectance Settlement <10 years     345304  89
#36     1                                                Other Land   99650901  91
#37     1                                              Snow and Ice   30752521  92
#38     1                                                      <NA> 2325604116

# ---- Merge two aligned categorical rasters (luts_w, luts_e) into a single coded map ----

# load binary raster
luts_e_file <- file.path(luts_aligned_dir, "luts10_e_aligned.tif")
luts_w_file <- file.path(luts_aligned_dir, "luts10_w_aligned.tif")
luts_e <- rast(luts_e_file)
luts_w <- rast(luts_w_file)

plot(luts_e)
plot(luts_w)

# 0) Robust write options (can tweak block sizes if needed)
wopt <- list(
  datatype = "INT1U",  # fits 0..255
  gdal = c(
    "COMPRESS=LZW",
    "PREDICTOR=2",
    "TILED=YES",
    "BIGTIFF=YES",
    "BLOCKXSIZE=512",
    "BLOCKYSIZE=512"
  )
)

# Normalization function for labels
normalize_label <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("&lt;", "<") %>%
    str_replace_all("&gt;", ">") %>%
    str_replace_all("[()]", "") %>%     # remove parentheses
    str_replace_all("[\u2010-\u2015]", "-") %>%  # normalize various hyphens to ASCII hyphen
    str_squish()                         # collapse inner whitespace
  # For case-insensitive matching, could add: tolower()
}

# 1) Target class code domain (code -> label)
lut <- tribble(
  ~code, ~label,
  21, "Settlement",
  22, "High Reflectance Settlement",
  24, "Settlement Forest",
  #25, "Roads",
  28, "Vegetated Settlement",
  29, "Very High Reflectance Settlement",
  31, "Water",
  41, "Forest",
  42, "Forest Wetland",
  43, "Forest Regenerating after Harvest (<20 years)",
  44, "Forest Wetland Regenerating after Harvest (<20 years)",
  47, "Forest Regenerating after Harvest (20-29 years)",
  48, "Forest Wetland Regenerating after Harvest (20-29 years)",
  49, "Forest Regenerating after Fire (<20 years)",
  51, "Cropland",
  52, "Annual Cropland",
  55, "Land Converted to Cropland",
  56, "Land Converted to Annual Cropland",
  63, "Prairie Grassland",
  64, "Temperate Grassland",
  65, "Temperate Rainforest Grassland",
  66, "Boreal Grassland",
  67, "Montane Grassland",
  68, "Alpine Tundra Grassland",
  69, "Arctic Tundra Grassland",
  73, "Prairie Wetland",
  74, "Temperate Wetland",
  75, "Temperate Rainforest Wetland",
  76, "Boreal Wetland",
  78, "Alpine Tundra Wetland",
  79, "Arctic Tundra Wetland",
  81, "Newly-Detected Settlement (<10 years)",
  82, "Newly-Detected High Reflectance Settlement (<10 years)",
  84, "Newly-Detected Settlement Forest (<10 years)",
  #85, "Newly-detected roads (<10 years)",
  88, "Newly-Detected Vegetated Settlement (<10 years)",
  89, "Newly-Detected Very High Reflectance Settlement (<10 years)",
  91, "Other Land",
  92, "Snow and Ice"
)

# 2) Ensure grids match exactly
stopifnot(compareGeom(luts_w, luts_e, crs=TRUE, res=TRUE, rowcol=TRUE, ext=TRUE, stopOnError=FALSE))

# --- Diagnostics: see what the RAT looks like on each raster ---
cat("RAT columns (west): ", paste(names(cats(luts_w)[[1]]), collapse=", "), "\n")
cat("RAT columns (east): ", paste(names(cats(luts_e)[[1]]), collapse=", "), "\n")

# Robust helper: detect ID column (handles 'Value'), normalize, and build reclass
mk_rcl <- function(r, lut) {
  cx <- cats(r)[[1]]
  
  if (is.null(cx) || !is.data.frame(cx)) {
    stop("Raster has no categories/RAT. You may need to supply a value->label mapping another way.")
  }
  
  # Find ID column
  id_candidates <- c("Value", "ID", "value", "id")
  id_col <- id_candidates[id_candidates %in% names(cx)][1]
  if (is.na(id_col)) {
    stop("Could not find an ID column among: ", paste(id_candidates, collapse=", "),
         ". RAT columns are: ", paste(names(cx), collapse=", "))
  }
  
  # Ensure LU2020 exists
  if (!"LU2020" %in% names(cx)) {
    stop("RAT does not have a 'LU2020' column. RAT columns: ", paste(names(cx), collapse=", "))
  }
  
  cx2 <- cx %>%
    dplyr::select(all_of(id_col), LU2020) %>%
    mutate(
      LU2020 = as.character(LU2020),
      norm_label = normalize_label(LU2020)
    )
  
  lut2 <- lut %>%
    mutate(
      label = as.character(label),
      norm_label = normalize_label(label)
    )
  
  # Join on normalized label
  jn <- cx2 %>%
    left_join(lut2, by = "norm_label")
  
  # Warn for labels that still didn't map
  miss <- jn %>%
    filter(is.na(code)) %>%
    distinct(LU2020) %>%
    arrange(LU2020)
  if (nrow(miss) > 0) {
    warning(
      "Unmapped labels in this raster (will become NA unless handled):\n  - ",
      paste(miss$LU2020, collapse = "\n  - ")
    )
  }
  
  # Reclass matrix: old ID -> new code
  rcl <- as.matrix(jn[, c(id_col, "code")])
  
  # Preserve 0 explicitly
  rcl <- rbind(rcl, c(0, 0))
  rcl <- unique(rcl)
  
  # Drop rows with NA in old ID
  rcl[!is.na(rcl[,1]), , drop = FALSE]
}


# 4) Reclassify to the shared numeric code domain
w_rcl <- mk_rcl(luts_w, lut)
e_rcl <- mk_rcl(luts_e, lut)

luts_w_code <- classify(luts_w, rcl = w_rcl, others = NA, right = TRUE,
                        filename = "", overwrite = TRUE, wopt = wopt)
luts_e_code <- classify(luts_e, rcl = e_rcl, others = NA, right = TRUE,
                        filename = "", overwrite = TRUE, wopt = wopt)

plot(luts_w_code)
plot(luts_e_code)

# 5) Merge (prefer west, fill with east)
out_file <- file.path(luts_merged_dir, "luts10_merged_aligned_take5.tif")
#merged <- cover(luts_w_code, luts_e_code, filename = out_file, overwrite = TRUE, wopt = wopt)
merged <- mosaic(luts_w_code, luts_e_code, fun = "max", filename = out_file, overwrite = TRUE, wopt = wopt)

plot(merged)

# 6) Attach clean categories table
cats(merged) <- data.frame(ID = lut$code, LU2020 = lut$label)

# 7) QA
print(merged)
qa <- freq(merged)
cx <- cats(merged)[[1]]
qa %>%
  as.data.frame() %>%
  left_join(cx, by = c("value" = "ID")) %>%
  arrange(value) %>%
  print(n = Inf)

# load merged raster
#luts_merged_file <- file.path(luts_merged_dir, "luts10_merged_aligned.tif")
#luts_merged <- rast(luts_merged_file)
#plot(luts_merged)

# reclassify to binary (1 = grass, 0 = not grass)
luts10_algn_binary <- classify(merged, 
                             rcl = matrix(c(0, NA,
                                            21, 0, # Settlement
                                            22, 0, # High Reflectance Settlement
                                            24, 0, # Settlement Forest
                                            25, 0, # Roads
                                            28, 0, # Vegetated Settlement
                                            29, 0, # Very High Reflectance Settlement
                                            31, 0, # Water
                                            41, 0, # Forest
                                            42, 0, # Forest Wetland
                                            43, 0, # Forest Regenerating after Harvest (<20 years)
                                            44, 0, # Forest Wetland Regenerating after Harvest (<20 years)
                                            47, 0, # Forest Regenerating after Harvest (20-29 years)
                                            48, 0, # Forest Wetland Regenerating after Harvest (20-29 years)
                                            49, 0, # Forest Regenerating after Fire (<20 years)
                                            51, 0, # Cropland
                                            52, 0, # Annual Cropland
                                            55, 0, # Land Converted to Cropland
                                            56, 0, # Land Converted to Annual Cropland
                                            63, 1, # Prairie Grassland
                                            64, 1, # Temperate Grassland
                                            65, 1, # Temperate Rainforest Grassland
                                            66, 1, # Boreal Grassland
                                            67, 1, # Montane Grassland
                                            68, 1, # Alpine Tundra Grassland
                                            69, 1, # Arctic Tundra Grassland
                                            73, 0, # Prairie Wetland
                                            74, 0, # Temperate Wetland
                                            75, 0, # Temperate Rainforest Wetland
                                            76, 0, # Boreal Wetland
                                            78, 0, # Alpine Tundra Wetland
                                            79, 0, # Arctic Tundra Wetland
                                            81, 0, # Newly-Detected Settlement (<10 years)
                                            82, 0, # Newly-Detected High Reflectance Settlement (<10 years)
                                            84, 0, # Newly-Detected Settlement Forest (<10 years)
                                            85, 0, # Newly-detected roads (<10 years)
                                            88, 0, # Newly-Detected Vegetated Settlement (<10 years)
                                            89, 0, # Newly-Detected Very High Reflectance Settlement (<10 years)
                                            91, 0, # Other Land
                                            92, 0), # Snow and Ice
                                          ncol = 2, byrow = TRUE))

# print raster details and plot
luts10_algn_binary
plot(luts10_algn_binary)

# Save the raster as 8-bit unsigned integer
luts_binary_out_path <- file.path(luts_aligned_dir, "aafc_luts_10m_2020_merged_reproj_algn_binary.tif")
writeRaster(luts10_algn_binary, luts_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load luts raster
#luts <- file.path(luts_aligned_dir, "cfga_luts_merged_reproj_algn.tif")
#luts_algn <- rast(luts)
#plot(luts_algn)

# reclassify to isolate each class
luts10_algn_crop <- classify(merged, 
                               rcl = matrix(c(0, NA,
                                              21, 0, # Settlement
                                              22, 0, # High Reflectance Settlement
                                              24, 0, # Settlement Forest
                                              25, 0, # Roads
                                              28, 0, # Vegetated Settlement
                                              29, 0, # Very High Reflectance Settlement
                                              31, 0, # Water
                                              41, 0, # Forest
                                              42, 0, # Forest Wetland
                                              43, 0, # Forest Regenerating after Harvest (<20 years)
                                              44, 0, # Forest Wetland Regenerating after Harvest (<20 years)
                                              47, 0, # Forest Regenerating after Harvest (20-29 years)
                                              48, 0, # Forest Wetland Regenerating after Harvest (20-29 years)
                                              49, 0, # Forest Regenerating after Fire (<20 years)
                                              51, 1, # Cropland
                                              52, 0, # Annual Cropland
                                              55, 1, # Land Converted to Cropland
                                              56, 0, # Land Converted to Annual Cropland
                                              63, 0, # Prairie Grassland
                                              64, 0, # Temperate Grassland
                                              65, 0, # Temperate Rainforest Grassland
                                              66, 0, # Boreal Grassland
                                              67, 0, # Montane Grassland
                                              68, 0, # Alpine Tundra Grassland
                                              69, 0, # Arctic Tundra Grassland
                                              73, 0, # Prairie Wetland
                                              74, 0, # Temperate Wetland
                                              75, 0, # Temperate Rainforest Wetland
                                              76, 0, # Boreal Wetland
                                              78, 0, # Alpine Tundra Wetland
                                              79, 0, # Arctic Tundra Wetland
                                              81, 0, # Newly-Detected Settlement (<10 years)
                                              82, 0, # Newly-Detected High Reflectance Settlement (<10 years)
                                              84, 0, # Newly-Detected Settlement Forest (<10 years)
                                              85, 0, # Newly-detected roads (<10 years)
                                              88, 0, # Newly-Detected Vegetated Settlement (<10 years)
                                              89, 0, # Newly-Detected Very High Reflectance Settlement (<10 years)
                                              91, 0, # Other Land
                                              92, 0), # Snow and Ice
                                            ncol = 2, byrow = TRUE))

# print raster details and plot
luts10_algn_crop
plot(luts10_algn_crop)

# Save the raster as 8-bit unsigned integer
luts_crop_out_path <- file.path(luts_aligned_dir, "aafc_luts_10m_2020_merged_reproj_algn_cropland.tif")
writeRaster(luts10_algn_crop, luts_crop_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa)
pa

# crop to pilot area
crop_ext <- ext(pa)
luts_pa <- crop(merged, crop_ext)
luts_pa
plot (luts_pa)

# mask to pilot area
luts_masked <- mask(luts_pa, pa)
luts_masked
plot(luts_masked)

# print raster details 
origin(luts_masked) # 6.696121 -2.550334
res(luts_masked)
crs(luts_masked)

# Save the raster as 8-bit unsigned integer
luts_pa_out_path <- file.path(luts_pilot_dir, "luts_10m_merged_reproj_algn_crop_pa.tif")
writeRaster(luts_masked, luts_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
luts_bin_pa <- crop(luts10_algn_binary, crop_ext)
luts_bin_pa
plot (luts_bin_pa)

# mask to pilot area
luts_bin_masked <- mask(luts_bin_pa, pa)
luts_bin_masked
plot(luts_bin_masked)

# Save the raster as 8-bit unsigned integer
luts_bin_pa_out_path <- file.path(luts_pilot_dir, "luts_10m_merged_reproj_algn_crop_pa_binary.tif")
writeRaster(luts_bin_masked, luts_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop cropland raster to pilot area
luts_cr_pa <- crop(luts10_algn_crop, crop_ext)
luts_cr_pa
plot (luts_cr_pa)

# mask to pilot area
luts_cr_masked <- mask(luts_cr_pa, pa)
luts_cr_masked
plot(luts_cr_masked)

# Save the raster as 8-bit unsigned integer
luts_cr_pa_out_path <- file.path(luts_pilot_dir, "luts_10m_merged_reproj_algn_crop_pa_cropland.tif")
writeRaster(luts_cr_masked, luts_cr_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load binary raster
#luts_binary_file <- file.path(luts_merged_dir, "luts_merged_aligned_bad.tif")
#luts_binary <- rast(luts_binary_file)

crop_ext <- ext(sa)

# crop raster to study area
luts_bin_sa <- crop(luts10_algn_binary, crop_ext)
luts_bin_sa
plot (luts_bin_sa)

# mask to study area
luts_bin_masked <- mask(luts_bin_sa, sa)
luts_bin_masked
plot(luts_bin_masked)

# Save the raster as 8-bit unsigned integer
luts_bin_sa_out_path <- file.path(luts_mask_dir, "luts10_reproj_algn_binary_masked.tif")
writeRaster(luts_bin_masked, luts_bin_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop raster to study area
luts_crop_sa <- crop(luts10_algn_crop, crop_ext)
luts_crop_sa
plot(luts_crop_sa)

# mask to study area
luts_crop_masked <- mask(luts_crop_sa, sa)
luts_crop_masked
plot(luts_crop_masked)

# Save the raster as 8-bit unsigned integer
luts_crop_sa_out_path <- file.path(luts_mask_dir, "luts10_reproj_algn_cropland_masked.tif")
writeRaster(luts_crop_masked, luts_crop_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load raster
luts_file <- file.path(luts_merged_dir, "luts10_merged_aligned_take5.tif")
luts <- rast(luts_file)

crop_ext <- ext(sa)

# crop raster to study area
luts_sa <- crop(luts, crop_ext)
luts_sa
plot (luts_sa)

# mask to study area
luts_masked <- mask(luts_sa, sa)
luts_masked
plot(luts_masked)

# Save the raster as 8-bit unsigned integer
luts_sa_out_path <- file.path(luts_masked_dir, "luts10_reproj_algn_merged_masked.tif")
writeRaster(luts_masked, luts_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
