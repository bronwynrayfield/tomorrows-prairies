# Reclassify AAFC LUTS rasters
# Reclass 2010 and 2020 to 4 classes: 1=Native; 3=Crop; 4=Settlement; 5=Other
# Reclass 2015 and 2020 to 5 classes: 1=Native; 2= Tame; 3=Crop; 4=Settlement; 5=Other

#--- Setup ----
source(file.path("scripts", "0-setup.R"))
source(file.path("scripts", "0-functions-general.R"))

#--- Load files ----
lulc2020 <- rast(tagged_path("./processeddatasets", "landcover/aafc/cropped", aafc2020_filename, "proj-naaea", "modelling-extent"))
lulc2015 <- rast(tagged_path("./processeddatasets", "landcover/aafc/cropped", aafc2015_filename, "proj-naaea", "modelling-extent"))
lulc2010 <- rast(tagged_path("./processeddatasets", "landcover/aafc/cropped", aafc2010_filename, "proj-naaea", "modelling-extent"))

#--- Reclassify to 4-class and 5-class LULC rasters----
# Prep shared Raster Attribute Table based on 2015 and 2020 AAFC LUTS rasters 
# 2015 includes class 85 (Newly-Detected Road <10 years)
# 2020 includes class 92 (Snow and Ice)
levs2015 <- levels(lulc2015)[[1]]
levs2020 <- levels(lulc2020)[[1]]
levs <- merge(levs2015, levs2020, by="Value", all=T)

# Coalesce (prefer 2020 label, fall back to 2015)
levs$class_name <- ifelse(is.na(levs$LU2020), levs$LU2015, levs$LU2020)
rat <- levs[, c("Value", "class_name")]
names(rat)[names(rat) == "Value"] <- "from"

# Define groups
# 5 groups
groups5 <- list(
  settlement = c(21,22,24,28,29,81,82,84,88,89),
  crop       = c(52,56),
  tame       = c(51,55),
  native     = c(63,64,65,66,67,68),
  other      = c(25,31,41,42,43,44,47,48,49,73,74,75,76,85,91,92)
)
# 4 groups
groups4 <- list(
  settlement = groups5$settlement,
  crop       = c(groups5$crop, groups5$tame),
  native     = groups5$native,
  other      = groups5$other
)

# --- Build 5-class table
to5 <- c(native = 1L, tame = 2L, crop = 3L, settlement = 4L, other = 5L)
tbl5 <- data.frame(from = unlist(groups5, use.names=FALSE),
                   group = rep(names(groups5), lengths(groups5)),
                   stringsAsFactors=FALSE)
tbl5$to <- unname(to5[tbl5$group])
qaqc5 <- merge(tbl5, rat, by="from", all.x=TRUE, sort=FALSE)
# qaqc5
# from      group to                                                class_name
# 1    21 settlement  4                                                Settlement
# 2    22 settlement  4                               High Reflectance Settlement
# 3    24 settlement  4                                         Settlement Forest
# 4    28 settlement  4                                      Vegetated Settlement
# 5    29 settlement  4                          Very High Reflectance Settlement
# 6    81 settlement  4                       Newly-Detected Settlement <10 years
# 7    82 settlement  4      Newly-Detected High Reflectance Settlement <10 years
# 8    84 settlement  4                Newly-Detected Settlement Forest <10 years
# 9    88 settlement  4             Newly-Detected Vegetated Settlement <10 years
# 10   89 settlement  4 Newly-Detected Very High Reflectance Settlement <10 years
# 11   52       crop  3                                           Annual Cropland
# 12   56       crop  3                         Land Converted to Annual Cropland
# 13   51       tame  2                                                  Cropland
# 14   55       tame  2                                Land Converted to Cropland
# 15   63     native  1                                         Prairie Grassland
# 16   64     native  1                                       Temperate Grassland
# 17   65     native  1                            Temperate Rainforest Grassland
# 18   66     native  1                                          Boreal Grassland
# 19   67     native  1                                         Montane Grassland
# 20   68     native  1                                   Alpine Tundra Grassland
# 21   25      other  5                                                     Roads
# 22   31      other  5                                                     Water
# 23   41      other  5                                                    Forest
# 24   42      other  5                                            Forest Wetland
# 25   43      other  5               Forest Regenerating after Harvest <20 years
# 26   44      other  5       Forest Wetland Regenerating after Harvest <20 years
# 27   47      other  5             Forest Regenerating after Harvest 20-29 years
# 28   48      other  5     Forest Wetland Regenerating after Harvest 20-29 years
# 29   49      other  5                  Forest Regenerating after Fire <20 years
# 30   73      other  5                                           Prairie Wetland
# 31   74      other  5                                         Temperate Wetland
# 32   75      other  5                              Temperate Rainforest Wetland
# 33   76      other  5                                            Boreal Wetland
# 34   85      other  5                             Newly-Detected Road <10 years
# 35   91      other  5                                                Other Land
# 36   92      other  5                                              Snow and Ice

rcl5 <- as.matrix(qaqc5[, c("from","to")])

# --- Build 4-class table
to4 <- c(native = 1L, crop = 3L, settlement = 4L, other = 5L)
tbl4 <- data.frame(from = unlist(groups4, use.names=FALSE),
                   group = rep(names(groups4), lengths(groups4)),
                   stringsAsFactors=FALSE)
tbl4$to <- unname(to4[tbl4$group])
qaqc4 <- merge(tbl4, rat, by="from", all.x=TRUE, sort=FALSE)
# qaqc4
# from      group to                                                class_name
# 1    21 settlement  4                                                Settlement
# 2    22 settlement  4                               High Reflectance Settlement
# 3    24 settlement  4                                         Settlement Forest
# 4    28 settlement  4                                      Vegetated Settlement
# 5    29 settlement  4                          Very High Reflectance Settlement
# 6    81 settlement  4                       Newly-Detected Settlement <10 years
# 7    82 settlement  4      Newly-Detected High Reflectance Settlement <10 years
# 8    84 settlement  4                Newly-Detected Settlement Forest <10 years
# 9    88 settlement  4             Newly-Detected Vegetated Settlement <10 years
# 10   89 settlement  4 Newly-Detected Very High Reflectance Settlement <10 years
# 11   52       crop  3                                           Annual Cropland
# 12   56       crop  3                         Land Converted to Annual Cropland
# 13   51       crop  3                                                  Cropland
# 14   55       crop  3                                Land Converted to Cropland
# 15   63     native  1                                         Prairie Grassland
# 16   64     native  1                                       Temperate Grassland
# 17   65     native  1                            Temperate Rainforest Grassland
# 18   66     native  1                                          Boreal Grassland
# 19   67     native  1                                         Montane Grassland
# 20   68     native  1                                   Alpine Tundra Grassland
# 21   25      other  5                                                     Roads
# 22   31      other  5                                                     Water
# 23   41      other  5                                                    Forest
# 24   42      other  5                                            Forest Wetland
# 25   43      other  5               Forest Regenerating after Harvest <20 years
# 26   44      other  5       Forest Wetland Regenerating after Harvest <20 years
# 27   47      other  5             Forest Regenerating after Harvest 20-29 years
# 28   48      other  5     Forest Wetland Regenerating after Harvest 20-29 years
# 29   49      other  5                  Forest Regenerating after Fire <20 years
# 30   73      other  5                                           Prairie Wetland
# 31   74      other  5                                         Temperate Wetland
# 32   75      other  5                              Temperate Rainforest Wetland
# 33   76      other  5                                            Boreal Wetland
# 34   85      other  5                             Newly-Detected Road <10 years
# 35   91      other  5                                                Other Land
# 36   92      other  5                                              Snow and Ice

rcl4 <- as.matrix(qaqc4[, c("from","to")])

# Reclassify
lulc2020_5 <- terra::classify(lulc2020, rcl5)
lulc2015_5 <- terra::classify(lulc2015, rcl5)
lulc2020_4 <- terra::classify(lulc2020, rcl4)
lulc2010_4 <- terra::classify(lulc2010, rcl4)

#--- QA/QC ----

# Original rasters
freqs_raw <- list(
  lulc2020 = freq(lulc2020, digits=0, bylayer=FALSE),
  lulc2015 = freq(lulc2015, digits=0, bylayer=FALSE),
  lulc2010 = freq(lulc2010, digits=0, bylayer=FALSE)
)

# Rename 'count' column to the list name, keep (value, count)
freqs_raw_named <- Map(function(df, nm) {
  df <- df[, c("value", "count")]
  names(df)[2] <- nm
  df
}, freqs_raw, names(freqs_raw))

# Merge into one table
freq_raw_all <- Reduce(function(x, y) merge(x, y, by="value", all=TRUE),
                       freqs_raw_named)

freq_raw_all[is.na(freq_raw_all)] <- 0
freq_raw_all <- freq_raw_all[order(freq_raw_all$value), ]

# Calculate totals
total_row <- data.frame(
  value = "TOTAL",
  t(colSums(freq_raw_all[, -1], na.rm = TRUE)),
  row.names = NULL
)

freq_raw_all <- rbind(freq_raw_all, total_row)
freq_raw_all
# value  lulc2020  lulc2015  lulc2010
# 1                                    Alpine Tundra Grassland    226050    226050    226056
# 2                                            Annual Cropland 280450646 280266203         0
# 3                                           Boreal Grassland    876617    877732    881451
# 4                                             Boreal Wetland  10178532  10186433  10227772
# 5                                                   Cropland 107199420 106535723 387697994
# 6                                                     Forest  75112242  75411321  76843978
# 7                   Forest Regenerating after Fire <20 years   1279581   1289447   1417370
# 8                Forest Regenerating after Harvest <20 years   5360976   5337535   5054091
# 9              Forest Regenerating after Harvest 20-29 years   1525623   1299281    414295
# 10                                            Forest Wetland  20850000  20791452  20922702
# 11       Forest Wetland Regenerating after Harvest <20 years    269245    280016    232618
# 12     Forest Wetland Regenerating after Harvest 20-29 years     59696    110010     73792
# 13                               High Reflectance Settlement   1501320   1208931   1181913
# 14                         Land Converted to Annual Cropland    913359   1159117         0
# 15                                Land Converted to Cropland   1077018   1468543   1183743
# 16                                         Montane Grassland   2019280   2019373   2021550
# 17      Newly-Detected High Reflectance Settlement <10 years    248663    339939    288565
# 18                             Newly-Detected Road <10 years         0     83946     85465
# 19                       Newly-Detected Settlement <10 years     69064    854191    895009
# 20                Newly-Detected Settlement Forest <10 years     41940    158983    121000
# 21             Newly-Detected Vegetated Settlement <10 years   1193447   2508041   1430102
# 22 Newly-Detected Very High Reflectance Settlement <10 years    118127    229751    211761
# 23                                                Other Land   2265200   2288533   2289068
# 24                                         Prairie Grassland  88861427  88941058  89801262
# 25                                           Prairie Wetland  31321229  31680786  32259391
# 26                                                     Roads  13433681  13349735  13348216
# 27                                                Settlement   3090585   2291372   2226878
# 28                                         Settlement Forest    727835    610388    599442
# 29                                              Snow and Ice     23326         0         0
# 30                                       Temperate Grassland     25436     25450     25461
# 31                            Temperate Rainforest Grassland      7334      7334      7358
# 32                              Temperate Rainforest Wetland      5285      5285      5291
# 33                                         Temperate Wetland    194729    194829    195134
# 34                                      Vegetated Settlement   5908615   4607356   4497227
# 35                          Very High Reflectance Settlement    733300    532160    510214
# 36                                                     Water  32559771  32552295  32552430
# 37                                                     TOTAL 689728599 689728599 689728599

# Reclassified rasters
freqs_recl <- list(
  lulc2020_5 = freq(lulc2020_5, digits=0, bylayer=FALSE),
  lulc2015_5 = freq(lulc2015_5, digits=0, bylayer=FALSE),
  lulc2020_4 = freq(lulc2020_4, digits=0, bylayer=FALSE),
  lulc2010_4 = freq(lulc2010_4, digits=0, bylayer=FALSE)
)

# Rename 'count' column to the list name, keep (value, count)
freqs_recl_named <- Map(function(df, nm) {
  df <- df[, c("value", "count")]
  names(df)[2] <- nm
  df
}, freqs_recl, names(freqs_recl))

# Merge into one table
freq_recl_all <- Reduce(function(x, y) merge(x, y, by="value", all=TRUE),
                        freqs_recl_named)

freq_recl_all[is.na(freq_recl_all)] <- 0
freq_recl_all <- freq_recl_all[order(freq_recl_all$value), ]

# Calculate totals
total_row <- data.frame(
  value = "TOTAL",
  t(colSums(freq_recl_all[, -1], na.rm = TRUE)),
  row.names = NULL
)

freq_recl_all <- rbind(freq_recl_all, total_row)
freq_recl_all
# value lulc2020_5 lulc2015_5 lulc2020_4 lulc2010_4
# 1     1   92016144   92096997   92016144   92963138
# 2     2  108276438  108004266          0          0
# 3     3  281364005  281425320  389640443  388881737
# 4     4   13632896   13341112   13632896   11962111
# 5     5  194439116  194860904  194439116  195921613
# 6 TOTAL  689728599  689728599  689728599  689728599


#--- Save to disk ----
write_rast_tagged(lulc2020_5, "./processeddatasets", "landcover/aafc/cropped", aafc2020_filename, "proj-naaea", "modelling-extent", "5-class")
write_rast_tagged(lulc2015_5, "./processeddatasets", "landcover/aafc/cropped", aafc2015_filename, "proj-naaea", "modelling-extent", "5-class")
write_rast_tagged(lulc2020_4, "./processeddatasets", "landcover/aafc/cropped", aafc2010_filename, "proj-naaea", "modelling-extent", "4-class")
write_rast_tagged(lulc2010_4, "./processeddatasets", "landcover/aafc/cropped", aafc2010_filename, "proj-naaea", "modelling-extent", "4-class")
