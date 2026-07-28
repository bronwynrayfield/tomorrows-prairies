
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
library(readxl)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# Prairies Ecozone
pe_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecozones"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# CNVC directory
cnvc_root <- file.path(raw_grass_data_dir, "Grassland/NatureServe_CNVC/updated_version_23Oct2025/NS_CNVC_TIFF_COMPRESSED")
cnvc_aligned_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/aligned")
cnvc_cropped_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/cropped")
cnvc_masked_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
cnvc_pilot_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/pilot")
cnvc_cat_dir <- file.path(raw_grass_data_dir, "./Grassland/NatureServe_CNVC/updated_version_23Oct2025")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# Prairies Ecozone
pe_file <- file.path(pe_dir, "ecozone_clean_prairies.shp")

# CNVC raster
cnvc_file <- file.path(cnvc_root, "NatureServe_CNVC_COMPRESSED_20251023.tif")

# CNVC broad categories
cnvc_cat_file <- file.path(cnvc_cat_dir, "CNVC_CAgrasslands_BroadCat_CANonly_DFL.xlsx")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- NatureServe CNVC ----#

# Load rasters
r_cnvc <- rast(cnvc_file)

origin(r_cnvc) # 6.696121 -2.550334
r_cnvc

#class       : SpatRaster 
#dimensions  : 273158, 243836, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2676983, 4638097, -3677433, 4517307  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : NatureServe_CNVC_COMPRESSED_20251023.tif 
#color table : 1 
#categories  : ELCODE, GNAME, GNAME2, Grank, RndGrank, ClassifCD, NATURAL_, WETLAND, NSX_LINK, NATIONS, SUBNATIONS, BiomeKey, BiomeName, SubBiKey, SubBiNam, EcoBiKey, EcoBiName, DivKey, DivName, MGKey, MGName, GroupKey, Group_, ClassName, HexCode, Red, Green, Blue 
#name        : ELCODE 
#min value   :   G002 
#max value   :  G9901

# print frequency table
#freq_table <- freq(r_cnvc) 
#print(freq_table) 

# Mark as categorical (factor)
#r_cnvc <- as.factor(r_cnvc)

# crop before reprojecting to make manageable size 

# Define the bbox in lon/lat degrees (WGS84)
canada_ext_ll  <- ext(-140, -70, 30, 60)
canada_bbox_ll <- as.polygons(canada_ext_ll, crs = "EPSG:4326")

# Get CRS from the raster 
r_crs <- crs(r_cnvc)

# Project bbox polygon to the raster's CRS
canada_bbox_rcrs <- project(canada_bbox_ll, r_crs)

# Meter/CRS-based extent derived from projected polygon
canada_ext_rcrs <- ext(canada_bbox_rcrs)

# polygon for Canada
#canada_sf <- ne_countries(country = "Canada", returnclass = "sf")
#canada_v  <- vect(canada_sf)

# plot
#plot(canada_v, col = "grey90", border = "grey40", main = "Canada + bounding box")
#plot(canada_bbox, add = TRUE, border = "red", lwd = 2, lty = 2)

# crop and mask to bounding box
r_cnvc_crop <- crop(r_cnvc, canada_ext_rcrs)
r_cnvc_mask <- mask(r_cnvc_crop, canada_bbox_rcrs)

r_cnvc_mask
plot(r_cnvc_mask)
cats(r_cnvc_mask)

#class       : SpatRaster 
#dimensions  : 124332, 168576, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2676983, 2380297, -840242.6, 2889717  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_27d8390d5e4d_10200.tif 
#color table : 1 
#varname     : NatureServe_CNVC_COMPRESSED_20251023 
#categories  : ELCODE, GNAME, GNAME2, Grank, RndGrank, ClassifCD, NATURAL_, WETLAND, NSX_LINK, NATIONS, SUBNATIONS, BiomeKey, BiomeName, SubBiKey, SubBiNam, EcoBiKey, EcoBiName, DivKey, DivName, MGKey, MGName, GroupKey, Group_, ClassName, HexCode, Red, Green, Blue 
#name        : ELCODE 
#min value   :   G007 
#max value   :  G9901 

# Fix missing ELCODE
stopifnot(length(levels(r_cnvc_mask)) > 0)
names(cats)
cnvc_levels <- levels(r_cnvc_mask)[[1]]
cats$ELCODE[cats$Value == 359 & (is.na(cats$ELCODE) | cats$ELCODE == "")]  <- "G359"
cats$ELCODE[cats$Value == 639 & (is.na(cats$ELCODE) | cats$ELCODE == "")]  <- "G639"
cats$ELCODE[cats$Value == 843 & (is.na(cats$ELCODE) | cats$ELCODE == "")]  <- "G843"
cats$ELCODE[cats$Value == 1215 & (is.na(cats$ELCODE) | cats$ELCODE == "")] <- "G1215"
cats$ELCODE[cats$Value == 1216 & (is.na(cats$ELCODE) | cats$ELCODE == "")] <- "G1216"
levels(r_cnvc_mask) <- cnvc_levels

sum(is.na(cats$ELCODE) | trimws(cats$ELCODE) == "")

# Save the raster as 8-bit unsigned integer
cnvc_crop_out_path <- file.path(cnvc_cropped_dir, "cnvc_crop_bbox.tif")
writeRaster(r_cnvc_mask, cnvc_crop_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# load cropped raster
#cnvc_crop_file <- file.path(cnvc_cropped_dir, "cnvc_crop_bbox.tif")
#cnvc_crop <- rast(cnvc_crop_file)
#plot(cnvc_crop)

# reproject before resampling
#cnvc_reproj <- project(r_cnvc, crs(sa)) 

#freq_table2 <- freq(cnvc_reproj) 
#print(freq_table2) 

#origin(cnvc_reproj) # 1.8392815 0.4291859
#cnvc_reproj

# align to reference raster 
cnvc_algn <- resample(r_cnvc_mask, sa, method = "near")
plot(cnvc_algn)
cnvc_algn

#class       : SpatRaster
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008)
#source      : spat_ae288953647_44584_DvggOkjkTOQeAAI.tif
#color table : 1
#varname     : class_study_area_aligned_cnvc_Dec2025
#categories  : ELCODE
#name        : ELCODE
#min value   :   G016
#max value   :  G9901

# Get values and counts
freq_table3 <- freq(cnvc_algn) 
print(freq_table3) 

#    layer value     count
#1       1  G016     22021
#2       1  G020      4310
#3       1  G021   7772304
#4       1  G030   3793862
#5       1  G032  12740535
#6       1  G045    685189
#7       1  G046  44322990
#8       1  G048  43416248
#9       1  G059  25938287
#10      1  G061     95427
#11      1  G075   2973883
#12      1  G089    107937
#13      1  G101    566426
#14      1  G103   1360058
#15      1  G105    643847
#16      1  G107      4572
#17      1  G113    819217
#18      1  G125   4427779
#19      1  G133  14867219
#20      1  G136        96
#21      1  G141  36156015
#22      1  G144    453115
#23      1  G145   4399116
#24      1  G146  10318706
#25      1  G147  24886575
#26      1  G160   1427837
#27      1  G161     29347
#28      1  G183    212311
#29      1  G184         1
#30      1  G195   3845190
#31      1  G198     12890
#32      1  G205      1439
#33      1  G206   2401001
#34      1  G208  21823307
#35      1  G209   2936622
#36      1  G210 125801316
#37      1  G211  17914993
#38      1  G212   7158703
#39      1  G213  32974636
#40      1  G215  31535634
#41      1  G216  23138785
#42      1  G217  31942018
#43      1  G218  81844257
#44      1  G219  53160465
#45      1  G220  41357027
#46      1  G222   5971662
#47      1  G223  51797423
#48      1  G225   3561711
#49      1  G226    144024
#50      1  G235   2229219
#51      1  G237  14247228
#52      1  G240  95220363
#53      1  G241  32703860
#54      1  G243   1820593
#55      1  G248  16843255
#56      1  G249   3228113
#57      1  G257   4839078
#58      1  G258      7901
#59      1  G261     15391
#60      1  G267   6435838
#61      1  G268   1044883
#62      1  G271  18442642
#63      1  G272  20047679
#64      1  G273  32621888
#65      1  G276    802763
#66      1  G277    180528
#67      1  G282   3168807
#68      1  G285     17764
#69      1  G300   5822471
#70      1  G301   8818938
#71      1  G302 152094545
#72      1  G303 120589819
#73      1  G304  88889586
#74      1  G305  16815640
#75      1  G307   5557618
#76      1  G308  27317761
#77      1  G310   5680069
#78      1  G311  11497182
#79      1  G314   1331588
#80      1  G316   2510588
#81      1  G317   8885012
#82      1  G318   3541308
#83      1  G319   6494469
#84      1  G320     32972
#85      1  G322     95544
#86      1  G328  19750064
#87      1  G329   1389576
#88      1  G331 278264400
#89      1  G332  10302944
#90      1  G333   1865354
#91      1  G334        29
#92      1  G340       141
#93      1  G342       621
#94      1  G344  36404808
#95      1  G345    122067
#96      1  G359     34485
#97      1  G488    865028
#98      1  G496    120460
#99      1  G497   3222431
#100     1  G498    235482
#101     1  G499    986508
#102     1  G505     38749
#103     1  G506   9094000
#104     1  G507   1651487
#105     1  G510    838468
#106     1  G517   1214790
#107     1  G520     68210
#108     1  G521   4271399
#109     1  G524   1957471
#110     1  G526   3782950
#111     1  G527   1080444
#112     1  G531   3649599
#113     1  G537   6748148
#114     1  G538   8488750
#115     1  G552    333610
#116     1  G556   2723639
#117     1  G559    875140
#118     1  G563    100185
#119     1  G565  12303815
#120     1  G566   7087237
#121     1  G567    968124
#122     1  G570   8348150
#123     1  G571  12852855
#124     1  G573    672814
#125     1  G600  39171673
#126     1  G624  13305322
#127     1  G638     22081
#128     1  G639     72101
#129     1  G647     13006
#130     1  G648   7788016
#131     1  G650  10615564
#132     1  G652   4494441
#133     1  G653   3945240
#134     1  G655     32350
#135     1  G662     70990
#136     1  G663     74395
#137     1  G673   1037361
#138     1  G678    159894
#139     1  G679  26039442
#140     1  G681      1334
#141     1  G742  37365562
#142     1  G744   1026742
#143     1  G749   6567328
#144     1  G751  11554959
#145     1  G767     15541
#146     1  G770   5831822
#147     1  G775    468331
#148     1  G788     15527
#149     1  G793      1863
#150     1  G800   1218025
#151     1  G804   1704028
#152     1  G805     30254
#153     1  G806  29528098
#154     1  G807     11558
#155     1  G809    125278
#156     1  G831  15599777
#157     1  G838     71280
#158     1  G839   1815317
#159     1  G841      3002
#160     1  G843      1905
#161     1  G844      1256
#162     1  G849  17649480
#163     1  G851   8831581
#164     1  G853    410058
#165     1  G888  19867421
#166     1  G889  14004894
#167     1  G891    160408
#168     1  G899   1125759
#169     1  G900     45904
#170     1  G902    850547
#171     1  G903    495965
#172     1  G904  11921923
#173     1  G906    619580
#174     1  G907  26783612
#175     1  G908    129539
#176     1  G918   2279657
#177     1  G919  17079791
#178     1  G920   1316682
#179     1  G921 125011825
#180     1  G922   7353214
#181     1  G962        11
#182     1  G984   1129272
#183     1  G990  22572425
#184     1  G991   8248041
#185     1  G999   1001288
#186     1 G1172   4211604
#187     1 G1176    129299
#188     1 G1177   5324896
#189     1 G1213   9288524
#190     1 G1214    709365
#191     1 G1215   7823152
#192     1 G1216    138128
#193     1 G8801 132156345
#194     1 G8802 987330411
#195     1 G8803   8834375
#196     1 G8804 380388505
#197     1 G8805  41561339
#198     1 G8806   2335353
#199     1 G8807 233339511
#200     1 G8809   5118800
#201     1 G9901    240608

# Save the raster as 8-bit unsigned integer
cnvc_out_path <- file.path(cnvc_aligned_dir, "cnvc_algn_good.tif")
writeRaster(cnvc_algn, cnvc_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# need broad categories to reclassify to binary

# load cropped raster
#cnvc_algn_file <- file.path(cnvc_aligned_dir, "cnvc_algn.tif")
#cnvc_algn <- rast(cnvc_algn_file)
#plot(cnvc_algn)

# print frequency table
#freq_table <- freq(cnvc_algn) 
#print(freq_table) 

# load spreadsheet with broad categories
cnvc_cat <- read_excel(cnvc_cat_file)
head(cnvc_cat)

# Inspect current levels (RAT)
lev <- levels(cnvc_algn)[[1]]   # data.frame or NULL if not categorical
if (is.null(lev)) {
  stop("cnvc_algn does not have a levels table (RAT)")
}

# Peek at column names to find the correct join key in the RAT
print(names(lev))
# The join key in the RAT might be called "Value", "value", or similar.

# Add missing ELCODE values
cnvc_cat <- cnvc_cat %>%
  mutate(
    ELCODE = if_else(
      is.na(ELCODE) | ELCODE == "",
        case_when(
          Value == "359"  ~ "G359",
          Value == "639"  ~ "G639",
          Value == "843"  ~ "G843",
          Value == "1215" ~ "G1215",
          Value == "1216" ~ "G1216",
        TRUE ~ ELCODE),
      ELCODE ))

# Value 359 also missing a broad category - need to confirm with Don
cnvc_cat <- cnvc_cat %>%
  mutate(
    Broad_Cat = if_else(
      ELCODE == "G359",
      "?",
      Broad_Cat))

# Harmonize a join key named "Value" in the RAT
if (!"Value" %in% names(lev)) {
  if ("value" %in% names(lev)) {
    names(lev)[names(lev) == "value"] <- "Value"
  } else if ("ID" %in% names(lev)) {
    # In some rasters, the RAT only has ID and no explicit Value column.
    # If cnvc_cat maps ID -> Broad_Cat, can join on ID instead.
    warning("No 'Value' column found in raster levels; found only 'ID'. If cnvc_cat maps original pixel values, use Option B to ratify.")
  }
}

# Join: keep all raster categories; attach Broad_Cat by matching Value
lev_joined <- merge(
  lev,
  cnvc_cat[, c("Value", "Broad_Cat")],
  by = "Value",
  all.x = TRUE
)

# Write the joined table back to the raster
levels(cnvc_algn)[[1]] <- lev_joined

# check that new field is present
cats(cnvc_algn)

# Define the mapping from Broad_Cat -> binary (0/1/NA)
bcat_to_bin <- function(bc) {
  if (is.na(bc)) return(NA_integer_)  # keep NA as NA
  
  switch(bc,
         "Other Native Upland Vegetation"           = 0L,
         "Wetland & Water"                          = 0L,
         "Ruderal Grassland"                        = 1L,
         "Alvar Native Grassland"                   = 1L,
         "Tallgrass Prairie Native Grassland"       = 1L,
         "Mixedgrass Prairie Native Grassland"      = 1L,
         "Rocky Mountain Montane Native Grassland"  = 1L,
         "Rocky Mountain Lowland Native Grassland"  = 1L,
         "Intermountain Steppe Native Grassland"    = 1L,
         "?"                                        = 0L,
         "Pacific Coastal Prairie Native Grassland" = 1L,
         "Atlantic Coastal Native Grassland"        = 1L,
         "Eastern Acidic Rocky Native Grassland"    = 1L,
         "Developed & Other Ag"                     = 0L,
         "Ag Pasture/Hay"                           = 1L,
         "Ruderal Vegetation"                       = 0L,
         # default if a Broad_Cat shows up that isn't listed:
         { warning(sprintf("Unmapped Broad_Cat encountered: '%s' -> set to NA", bc)); NA_integer_ }
  )
}

# Build the lookup table: Value -> binary
stopifnot(all(c("Value", "Broad_Cat") %in% names(cnvc_cat)))

# Normalize spaces just in case (optional but helpful)
cnvc_cat$Broad_Cat <- trimws(cnvc_cat$Broad_Cat)

# Vectorized mapping
bin_vec <- vapply(cnvc_cat$Broad_Cat, bcat_to_bin, integer(1))

lut <- data.frame(
  from = cnvc_cat$Value,  # original raster codes
  to   = bin_vec          # 0/1/NA per Broad_Cat
)

# Use classify with a two-column data.frame for exact value mapping
# 'others=NA' means any raster value not present in 'from' becomes NA 
cnvc_algn_binary <- classify(cnvc_algn, rcl = lut, others = NA)

# Verify results and plot
print(cnvc_algn_binary)
plot(cnvc_algn_binary, main = "Binary reclass by Broad_Cat (1 = grass)")

# Check unique output values (should be 0, 1, and NA)
unique_vals <- unique(values(cnvc_algn_binary, mat = FALSE))
print(unique_vals)

# Save the raster as 8-bit unsigned integer
cnvc_binary_out_path <- file.path(cnvc_aligned_dir, "cnvc_algn_binary.tif")
writeRaster(cnvc_algn_binary, cnvc_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load cnvc raster
#cnvc_algn_file <- file.path(cnvc_aligned_dir, "cnvc_algn_good.tif")
#cnvc_algn <- rast(cnvc_algn_file)
#plot(cnvc_algn)
#cnvc_algn

# Define the mapping from Broad_Cat -> native vs. tame (0 = not grass, 1 = tame, 2 = native)
bcat_to_nt <- function(bc) {
  if (is.na(bc)) return(NA_integer_)  # keep NA as NA
  
  switch(bc,
         "Other Native Upland Vegetation"           = 0L,
         "Wetland & Water"                          = 0L,
         "Ruderal Grassland"                        = 1L,
         "Alvar Native Grassland"                   = 2L,
         "Tallgrass Prairie Native Grassland"       = 2L,
         "Mixedgrass Prairie Native Grassland"      = 2L,
         "Rocky Mountain Montane Native Grassland"  = 2L,
         "Rocky Mountain Lowland Native Grassland"  = 2L,
         "Intermountain Steppe Native Grassland"    = 2L,
         "?"                                        = 0L,
         "Pacific Coastal Prairie Native Grassland" = 2L,
         "Atlantic Coastal Native Grassland"        = 2L,
         "Eastern Acidic Rocky Native Grassland"    = 2L,
         "Developed & Other Ag"                     = 0L,
         "Ag Pasture/Hay"                           = 1L,
         "Ruderal Vegetation"                       = 0L,
         # default if a Broad_Cat shows up that isn't listed:
         { warning(sprintf("Unmapped Broad_Cat encountered: '%s' -> set to NA", bc)); NA_integer_ }
  )
}

# Build the lookup table: Value -> native vs. tame
stopifnot(all(c("Value", "Broad_Cat") %in% names(cnvc_cat)))

# Normalize spaces just in case (optional but helpful)
cnvc_cat$Broad_Cat <- trimws(cnvc_cat$Broad_Cat)

# Vectorized mapping
nt_vec <- vapply(cnvc_cat$Broad_Cat, bcat_to_nt, integer(1))

lut_nt <- data.frame(
  from = cnvc_cat$Value,  # original raster codes
  to   = nt_vec          # 0/1/2/NA per Broad_Cat
)

# Use classify with a two-column data.frame for exact value mapping
# 'others=NA' means any raster value not present in 'from' becomes NA 
cnvc_algn_nattam <- classify(cnvc_algn, rcl = lut_nt, others = NA)

# Verify results and plot
cnvc_algn_nattam
plot(cnvc_algn_nattam, main = "Native vs. tame reclass by Broad_Cat (1 = tame, 2 = native)")

# Check unique output values (should be 0, 1, 2, and NA)
unique_vals_nt <- unique(values(cnvc_algn_nattam, mat = FALSE))
print(unique_vals_nt)

# Save the raster as 8-bit unsigned integer
cnvc_nt_out_path <- file.path(cnvc_aligned_dir, "cnvc_algn_natvtam.tif")
writeRaster(cnvc_algn_nattam, cnvc_nt_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################ isolate native

# Define the mapping from Broad_Cat -> native 
bcat_to_n <- function(bc) {
  if (is.na(bc)) return(NA_integer_)  # keep NA as NA
  
  switch(bc,
         "Other Native Upland Vegetation"           = 0L,
         "Wetland & Water"                          = 0L,
         "Ruderal Grassland"                        = 0L,
         "Alvar Native Grassland"                   = 1L,
         "Tallgrass Prairie Native Grassland"       = 1L,
         "Mixedgrass Prairie Native Grassland"      = 1L,
         "Rocky Mountain Montane Native Grassland"  = 1L,
         "Rocky Mountain Lowland Native Grassland"  = 1L,
         "Intermountain Steppe Native Grassland"    = 1L,
         "?"                                        = 0L,
         "Pacific Coastal Prairie Native Grassland" = 1L,
         "Atlantic Coastal Native Grassland"        = 1L,
         "Eastern Acidic Rocky Native Grassland"    = 1L,
         "Developed & Other Ag"                     = 0L,
         "Ag Pasture/Hay"                           = 0L,
         "Ruderal Vegetation"                       = 0L,
         # default if a Broad_Cat shows up that isn't listed:
         { warning(sprintf("Unmapped Broad_Cat encountered: '%s' -> set to NA", bc)); NA_integer_ }
  )
}

# Build the lookup table: Value -> native 
stopifnot(all(c("Value", "Broad_Cat") %in% names(cnvc_cat)))

# Normalize spaces just in case (optional but helpful)
cnvc_cat$Broad_Cat <- trimws(cnvc_cat$Broad_Cat)

# Vectorized mapping
n_vec <- vapply(cnvc_cat$Broad_Cat, bcat_to_n, integer(1))

lut_n <- data.frame(
  from = cnvc_cat$Value,  # original raster codes
  to   = n_vec          # 0/1/NA per Broad_Cat
)

# Use classify with a two-column data.frame for exact value mapping
# 'others=NA' means any raster value not present in 'from' becomes NA 
cnvc_algn_native <- classify(cnvc_algn, rcl = lut_n, others = NA)

# Verify results and plot
print(cnvc_algn_native)
plot(cnvc_algn_native, main = "Native grass")

# Check unique output values (should be 0, 1, and NA)
#unique_vals_n <- unique(values(cnvc_algn_native, mat = FALSE))
#print(unique_vals_n)

# Save the raster as 8-bit unsigned integer
cnvc_n_out_path <- file.path(cnvc_aligned_dir, "cnvc_algn_native.tif")
writeRaster(cnvc_algn_native, cnvc_n_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################ isolate ruderal

# Define the mapping from Broad_Cat -> ruderal
bcat_to_rud <- function(bc) {
  if (is.na(bc)) return(NA_integer_)  # keep NA as NA
  
  switch(bc,
         "Other Native Upland Vegetation"           = 0L,
         "Wetland & Water"                          = 0L,
         "Ruderal Grassland"                        = 1L,
         "Alvar Native Grassland"                   = 0L,
         "Tallgrass Prairie Native Grassland"       = 0L,
         "Mixedgrass Prairie Native Grassland"      = 0L,
         "Rocky Mountain Montane Native Grassland"  = 0L,
         "Rocky Mountain Lowland Native Grassland"  = 0L,
         "Intermountain Steppe Native Grassland"    = 0L,
         "?"                                        = 0L,
         "Pacific Coastal Prairie Native Grassland" = 0L,
         "Atlantic Coastal Native Grassland"        = 0L,
         "Eastern Acidic Rocky Native Grassland"    = 0L,
         "Developed & Other Ag"                     = 0L,
         "Ag Pasture/Hay"                           = 0L,
         "Ruderal Vegetation"                       = 0L,
         # default if a Broad_Cat shows up that isn't listed:
         { warning(sprintf("Unmapped Broad_Cat encountered: '%s' -> set to NA", bc)); NA_integer_ }
  )
}

# Build the lookup table: Value -> ruderal
stopifnot(all(c("Value", "Broad_Cat") %in% names(cnvc_cat)))

# Normalize spaces just in case (optional but helpful)
cnvc_cat$Broad_Cat <- trimws(cnvc_cat$Broad_Cat)

# Vectorized mapping
rud_vec <- vapply(cnvc_cat$Broad_Cat, bcat_to_rud, integer(1))

lut_rud <- data.frame(
  from = cnvc_cat$Value,  # original raster codes
  to   = rud_vec          # 0/1/NA per Broad_Cat
)

# Use classify with a two-column data.frame for exact value mapping
# 'others=NA' means any raster value not present in 'from' becomes NA 
cnvc_algn_ruderal <- classify(cnvc_algn, rcl = lut_rud, others = NA)

# Verify results and plot
print(cnvc_algn_ruderal)
plot(cnvc_algn_ruderal, main = "Ruderal grass")

# Check unique output values (should be 0, 1, and NA)
#unique_vals_rud <- unique(values(cnvc_algn_ruderal, mat = FALSE))
#print(unique_vals_rud)

# Save the raster as 8-bit unsigned integer
cnvc_rud_out_path <- file.path(cnvc_aligned_dir, "cnvc_algn_ruderal.tif")
writeRaster(cnvc_algn_ruderal, cnvc_rud_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################ isolate ag pasture/hay

# load spreadsheet with broad categories
cnvc_cat <- read_excel(cnvc_cat_file)
head(cnvc_cat)

# Inspect current levels (RAT)
lev <- levels(cnvc_algn)[[1]]   # data.frame or NULL if not categorical
if (is.null(lev)) {
  stop("cnvc_algn does not have a levels table (RAT). Use Option B to ratify first.")
}

# Peek at column names to find the correct join key in the RAT
print(names(lev))
# The join key in the RAT might be called "Value", "value", or similar.

# Harmonize a join key named "Value" in the RAT
if (!"Value" %in% names(lev)) {
  if ("value" %in% names(lev)) {
    names(lev)[names(lev) == "value"] <- "Value"
  } else if ("ID" %in% names(lev)) {
    # In some rasters, the RAT only has ID and no explicit Value column.
    # If your cnvc_cat maps ID -> Broad_Cat, you can join on ID instead.
    # Otherwise you need Option B (ratify) to ensure the original values are recorded.
    warning("No 'Value' column found in raster levels; found only 'ID'. If your cnvc_cat maps original pixel values, use Option B to ratify.")
  }
}

# Join: keep all raster categories; attach Broad_Cat by matching Value
lev_joined <- merge(
  lev,
  cnvc_cat[, c("Value", "Broad_Cat")],
  by = "Value",
  all.x = TRUE
)

# Write the joined table back to the raster
levels(cnvc_algn)[[1]] <- lev_joined

# check that new field is present
cats(cnvc_algn)

# Define the mapping from Broad_Cat -> ag pasture/hay
bcat_to_past <- function(bc) {
  if (is.na(bc)) return(NA_integer_)  # keep NA as NA
  
  switch(bc,
         "Other Native Upland Vegetation"           = 0L,
         "Wetland & Water"                          = 0L,
         "Ruderal Grassland"                        = 0L,
         "Alvar Native Grassland"                   = 0L,
         "Tallgrass Prairie Native Grassland"       = 0L,
         "Mixedgrass Prairie Native Grassland"      = 0L,
         "Rocky Mountain Montane Native Grassland"  = 0L,
         "Rocky Mountain Lowland Native Grassland"  = 0L,
         "Intermountain Steppe Native Grassland"    = 0L,
         "?"                                        = 0L,
         "Pacific Coastal Prairie Native Grassland" = 0L,
         "Atlantic Coastal Native Grassland"        = 0L,
         "Eastern Acidic Rocky Native Grassland"    = 0L,
         "Developed & Other Ag"                     = 0L,
         "Ag Pasture/Hay"                           = 1L,
         "Ruderal Vegetation"                       = 0L,
         # default if a Broad_Cat shows up that isn't listed:
         { warning(sprintf("Unmapped Broad_Cat encountered: '%s' -> set to NA", bc)); NA_integer_ }
  )
}

# Build the lookup table: Value -> ag pasture/hay
stopifnot(all(c("Value", "Broad_Cat") %in% names(cnvc_cat)))

# Normalize spaces just in case (optional but helpful)
cnvc_cat$Broad_Cat <- trimws(cnvc_cat$Broad_Cat)

# Vectorized mapping
past_vec <- vapply(cnvc_cat$Broad_Cat, bcat_to_past, integer(1))

lut_past <- data.frame(
  from = cnvc_cat$Value,  # original raster codes
  to   = past_vec          # 0/1/NA per Broad_Cat
)

# Use classify with a two-column data.frame for exact value mapping
# 'others=NA' means any raster value not present in 'from' becomes NA 
cnvc_algn_past <- classify(cnvc_algn, rcl = lut_past, others = NA)

# Verify results and plot
print(cnvc_algn_past)
plot(cnvc_algn_past, main = "Ag pasture/hay")

# Check unique output values (should be 0, 1, and NA)
#unique_vals_past <- unique(values(cnvc_algn_past, mat = FALSE))
#print(unique_vals_past)

# Save the raster as 8-bit unsigned integer
cnvc_past_out_path <- file.path(cnvc_aligned_dir, "cnvc_algn_agpasturehay.tif")
writeRaster(cnvc_algn_past, cnvc_past_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

####### for pilot area

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334

# crop to pilot area
crop_ext <- ext(pa)
cnvc_pa <- crop(r_cnvc_mask_trim, crop_ext)
cnvc_pa
plot (cnvc_pa)

# mask to pilot area
cnvc_masked <- mask(cnvc_pa, pa)
cnvc_masked
plot(cnvc_masked)

# print raster details 
origin(cnvc_masked) # 6.696121 -2.550334

#class       : SpatRaster 
#dimensions  : 1367, 1692, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -828263.3, -777503.3, 1111437, 1152447  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source(s)   : memory
#color table : 1 
#varname     : cnvc_crop_bbox_trim 
#categories  : ELCODE, Broad_Cat 
#name        : ELCODE 
#min value   :   G141 
#max value   :  G9901

# Save the raster as 8-bit unsigned integer
cnvc_pa_out_path <- file.path(cnvc_cropped_dir, "cnvc_crop_pa.tif")
writeRaster(cnvc_masked, cnvc_pa_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# crop binary raster to pilot area
cnvc_bin_pa <- crop(cnvc_trim_binary, crop_ext)
cnvc_bin_pa
plot (cnvc_bin_pa)

# mask to pilot area
cnvc_bin_masked <- mask(cnvc_bin_pa, pa)
cnvc_bin_masked
plot(cnvc_bin_masked)

#class       : SpatRaster 
#dimensions  : 1367, 1692, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -828263.3, -777503.3, 1111437, 1152447  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source(s)   : memory
#varname     : cnvc_crop_bbox_trim 
#name        : ELCODE 
#min value   :      0 
#max value   :      1 

# Save the raster as 8-bit unsigned integer
cnvc_bin_pa_out_path <- file.path(cnvc_pilot_dir, "cnvc_crop_pa_binary.tif")
writeRaster(cnvc_bin_masked, cnvc_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
cnvc_nt_pa <- crop(cnvc_algn_nattam, crop_ext)
cnvc_nt_pa
plot (cnvc_nt_pa)

# mask to pilot area
cnvc_nt_masked <- mask(cnvc_nt_pa, pa)
cnvc_nt_masked
plot(cnvc_nt_masked)

# Save the raster as 8-bit unsigned integer
cnvc_nt_pa_out_path <- file.path(cnvc_pilot_dir, "cnvc_algn_crop_pa_natvtam.tif")
writeRaster(cnvc_nt_masked, cnvc_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################

######### crop and mask to study area

# load raster
#cnvc_file <- file.path(cnvc_aligned_dir, "cnvc_algn_good.tif")
#cnvc <- rast(cnvc_file)

crop_ext <- ext(sa)

# crop raster to study area
cnvc_sa <- crop(cnvc_algn, crop_ext)
cnvc_sa
plot (cnvc_sa)

# mask to study area
cnvc_masked <- mask(cnvc_sa, sa)
cnvc_masked
plot(cnvc_masked)

cnvc_masked

# Save the raster as 8-bit unsigned integer
cnvc_sa_out_path <- file.path(cnvc_masked_dir, "cnvc_algn_masked.tif")
writeRaster(cnvc_masked, cnvc_sa_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

######### crop and mask to study area

# load binary raster
cnvc_binary_file <- file.path(cnvc_aligned_dir, "cnvc_algn_binary.tif")
cnvc_binary <- rast(cnvc_binary_file)

crop_ext <- ext(sa)

# crop raster to study area
cnvc_bin_sa <- crop(cnvc_binary, crop_ext)
cnvc_bin_sa
plot (cnvc_bin_sa)

# mask to study area
cnvc_bin_masked <- mask(cnvc_bin_sa, sa)
cnvc_bin_masked
plot(cnvc_bin_masked)

# Save the raster as 8-bit unsigned integer
cnvc_bin_sa_out_path <- file.path(cnvc_masked_dir, "cnvc_algn_binary_masked.tif")
writeRaster(cnvc_bin_masked, cnvc_bin_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load native vs tame raster
cnvc_nattam_file <- file.path(cnvc_aligned_dir, "cnvc_algn_natvtam.tif")
cnvc_nattam <- rast(cnvc_nattam_file)

crop_ext <- ext(sa)

# crop raster to study area
cnvc_nt_sa <- crop(cnvc_nattam, crop_ext)
cnvc_nt_sa
plot (cnvc_nt_sa)

# mask to study area
cnvc_nt_masked <- mask(cnvc_nt_sa, sa)
cnvc_nt_masked
plot(cnvc_nt_masked)

# Save the raster as 8-bit unsigned integer
cnvc_nt_sa_out_path <- file.path(cnvc_masked_dir, "cnvc_algn_natvtam_masked.tif")
writeRaster(cnvc_nt_masked, cnvc_nt_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load ag pasture/hay raster
cnvc_aph_file <- file.path(cnvc_aligned_dir, "cnvc_algn_agpasturehay.tif")
cnvc_aph <- rast(cnvc_aph_file)

crop_ext <- ext(sa)

# crop raster to study area
cnvc_aph_sa <- crop(cnvc_aph, crop_ext)
cnvc_aph_sa
plot (cnvc_aph_sa)

# mask to study area
cnvc_aph_masked <- mask(cnvc_aph_sa, sa)
cnvc_aph_masked
plot(cnvc_aph_masked)

# Save the raster as 8-bit unsigned integer
cnvc_aph_sa_out_path <- file.path(cnvc_masked_dir, "cnvc_algn_agpasturehay_masked.tif")
writeRaster(cnvc_aph_masked, cnvc_aph_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load native raster
cnvc_nat_file <- file.path(cnvc_aligned_dir, "cnvc_algn_native.tif")
cnvc_nat <- rast(cnvc_nat_file)

crop_ext <- ext(sa)

# crop raster to study area
cnvc_nat_sa <- crop(cnvc_nat, crop_ext)
cnvc_nat_sa
plot (cnvc_nat_sa)

# mask to study area
cnvc_nat_masked <- mask(cnvc_nat_sa, sa)
cnvc_nat_masked
plot(cnvc_nat_masked)

# Save the raster as 8-bit unsigned integer
cnvc_nat_sa_out_path <- file.path(cnvc_masked_dir, "cnvc_algn_native_masked.tif")
writeRaster(cnvc_nat_masked, cnvc_nat_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load ruderal raster
cnvc_rud_file <- file.path(cnvc_aligned_dir, "cnvc_algn_ruderal.tif")
cnvc_rud <- rast(cnvc_rud_file)

crop_ext <- ext(sa)

# crop raster to study area
cnvc_rud_sa <- crop(cnvc_rud, crop_ext)
cnvc_rud_sa
plot (cnvc_rud_sa)

# mask to study area
cnvc_rud_masked <- mask(cnvc_rud_sa, sa)
cnvc_rud_masked
plot(cnvc_rud_masked)

# Save the raster as 8-bit unsigned integer
cnvc_rud_sa_out_path <- file.path(cnvc_masked_dir, "cnvc_algn_ruderal_masked.tif")
writeRaster(cnvc_rud_masked, cnvc_rud_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

##################################################

# Crop to Prairies Ecozone for Don

# read into terra
pe <- vect(pe_file)
plot(pe)

# load native raster
cnvc_nat_file <- file.path(cnvc_aligned_dir, "cnvc_algn_native.tif")
cnvc_nat <- rast(cnvc_nat_file)
plot(cnvc_nat)

# --- Rasterize shapefile 'pe' to match raster 'sa' ---
# Reproject classification study area to match CNVC raster CRS
pe_vect_proj <- project(pe, crs(sa))

# field = 1 makes all polygon cells have value 1
pe_rast <- rasterize(
  pe_vect_proj,
  sa,
  field = 1,
  touches = TRUE     # include boundary cells
)

plot(pe_rast)
pe_rast

# --- Mask the CNVC raster using the new polygon raster ---
crop_ext <- ext(pe_rast)

# crop raster to PE
cnvc_nat_crop <- crop(cnvc_nat, crop_ext)
cnvc_nat_crop
plot (cnvc_nat_crop)

# mask to PE
cnvc_nat_pe <- mask(cnvc_nat_crop, pe_rast)
cnvc_nat_pe
plot(cnvc_nat_pe)

# --- Save result ---

# Save the raster as 8-bit unsigned integer
cnvc_nat_pe_out_path <- file.path(cnvc_masked_dir, "cnvc_algn_native_masked_pe.tif")
writeRaster(cnvc_nat_pe, cnvc_nat_pe_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Save the raster as 8-bit unsigned integer
pe_out_path <- file.path(pe_dir, "prairies_ecozone_algn_rast.tif")
writeRaster(pe_rast, pe_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

