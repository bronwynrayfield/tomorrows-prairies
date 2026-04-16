
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
cnvc_root <- file.path(raw_grass_data_dir, "./Grassland/NatureServe_CNVC/updated version 23Oct2025/NS_CNVC_TIFF_COMPRESSED")
cnvc_aligned_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/aligned")
cnvc_cropped_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/cropped")
cnvc_masked_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/masked")
cnvc_pilot_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/pilot")
cnvc_cat_dir <- file.path(raw_grass_data_dir, "./Grassland/NatureServe_CNVC/updated version 23Oct2025")

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

#----- natureserve Plowprint ----#

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

# Save the raster as 8-bit unsigned integer
cnvc_crop_out_path <- file.path(cnvc_cropped_dir, "cnvc_crop_bbox.tif")
writeRaster(r_cnvc_mask, cnvc_crop_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# load cropped raster
cnvc_crop_file <- file.path(cnvc_cropped_dir, "cnvc_crop_bbox.tif")
cnvc_crop <- rast(cnvc_crop_file)
plot(cnvc_crop)

# reproject before resampling
#cnvc_reproj <- project(r_cnvc, crs(sa)) 

#freq_table2 <- freq(cnvc_reproj) 
#print(freq_table2) 

#origin(cnvc_reproj) # 1.8392815 0.4291859
#cnvc_reproj

# align to reference raster 
cnvc_algn <- resample(cnvc_crop, sa, method = "near")
plot(cnvc_algn)
cnvc_algn

# Get values and counts
#freq_table3 <- freq(cnvc_reproj_algn) 
#print(freq_table3) 

#origin(cnvc_reproj_algn) # 6.696121 -2.550334
#cnvc_reproj_algn

# Save the raster as 8-bit unsigned integer
cnvc_out_path <- file.path(cnvc_aligned_dir, "cnvc_algn_good.tif")
writeRaster(cnvc_algn, cnvc_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# Trim NA-only rows/columns
r_cnvc_mask_trim <- trim(r_cnvc_mask)

# visualize the trimmed raster
plot(r_cnvc_mask_trim)

# print raster details 
origin(r_cnvc_mask_trim) # 6.696121 -2.550334
r_cnvc_mask_trim

#class       : SpatRaster 
#dimensions  : 92747, 142291, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2242973, 2025757, -749792.6, 2032617  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_27d8e84469_10200.tif 
#color table : 1 
#varname     : NatureServe_CNVC_COMPRESSED_20251023 
#categories  : ELCODE, GNAME, GNAME2, Grank, RndGrank, ClassifCD, NATURAL_, WETLAND, NSX_LINK, NATIONS, SUBNATIONS, BiomeKey, BiomeName, SubBiKey, SubBiNam, EcoBiKey, EcoBiName, DivKey, DivName, MGKey, MGName, GroupKey, Group_, ClassName, HexCode, Red, Green, Blue 
#name        : ELCODE 
#min value   :   G007 
#max value   :  G9901

# Save the raster as 8-bit unsigned integer
cnvc_trim_out_path <- file.path(cnvc_cropped_dir, "cnvc_crop_bbox_trim.tif")
writeRaster(r_cnvc_mask_trim, cnvc_trim_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# load trimmed raster
cnvc_trim_file <- file.path(cnvc_cropped_dir, "cnvc_crop_bbox_trim.tif")
r_cnvc_mask_trim <- rast(cnvc_trim_file)

# what values are present outside of study area?
# Location in longitude/latitude
pt <- data.frame(lon = -122.4787, lat = 48.7491)  # example: Bellingham, WA
pt <- data.frame(lon = -116.1630295, lat = 52.4805298)  # example: just N of cnvc extent

# Make SpatVector point (declare CRS). Change to your point CRS if needed.
p <- vect(pt, geom = c("lon", "lat"), crs = "EPSG:4326")

# If raster is not in EPSG:4326, reproject the point to raster CRS
if (crs(r_cnvc_mask_trim) != crs(p)) {
  p <- project(p, crs(r_cnvc_mask_trim))
}

# Extract (nearest neighbor)
val_nearest <- terra::extract(r_cnvc_mask_trim, p)  # returns a data.frame
print(val_nearest)

# print frequency table
freq_table <- freq(r_cnvc_mask_trim) 
print(freq_table) 

    layer value      count
1       1  G007      56466
2       1  G009    3304371
3       1  G013     155396
4       1  G015   69324908
5       1  G016    1662678
6       1  G020   90364705
7       1  G021   16076243
8       1  G029     155055
9       1  G030   16809546
10      1  G031   55385986
11      1  G032   53107103
12      1  G033    5028483
13      1  G034   17166727
14      1  G036     138403
15      1  G037   10122771
16      1  G038    2193210
17      1  G039    7064084
18      1  G044       8573
19      1  G045    3273849
20      1  G046   47120251
21      1  G048   43480021
22      1  G059   53061308
23      1  G061      95433
24      1  G063     102423
25      1  G075    2979190
26      1  G089     115061
27      1  G101    3659652
28      1  G102     743067
29      1  G103   38286944
30      1  G105    1925923
31      1  G107     191407
32      1  G109        545
33      1  G113    2960524
34      1  G125    5861925
35      1  G130     520251
36      1  G133   98089571
37      1  G136     318450
38      1  G141   36156169
39      1  G144  130850410
40      1  G145    4399495
41      1  G146   10318706
42      1  G147   36888283
43      1  G154    4738452
44      1  G160    1427837
45      1  G161    2709301
46      1  G165   32187073
47      1  G166    5307512
48      1  G177          1
49      1  G178     119034
50      1  G179     387683
51      1  G180      36771
52      1  G183     242039
53      1  G184      44117
54      1  G186     145645
55      1  G187       1661
56      1  G190    1614274
57      1  G191     387948
58      1  G192    7617233
59      1  G195   20788357
60      1  G198     159930
61      1  G200    1562533
62      1  G201      65934
63      1  G202      21781
64      1  G203     104468
65      1  G205       1439
66      1  G206    2401001
67      1  G208   25522695
68      1  G209    4079382
69      1  G210  125801316
70      1  G211   17914993
71      1  G212    7158703
72      1  G213   32974636
73      1  G215   31535634
74      1  G216   23309762
75      1  G217   31942018
76      1  G218   82450436
77      1  G219   83313108
78      1  G220   51768958
79      1  G222   26947219
80      1  G223   51797423
81      1  G225   15349793
82      1  G226   13449442
83      1  G235    5738800
84      1  G237   14247228
85      1  G240   95220363
86      1  G241   32703860
87      1  G243    5031613
88      1  G248   16843255
89      1  G249    7390603
90      1  G252    1451427
91      1  G253   15124727
92      1  G257   23532860
93      1  G258      34572
94      1  G261    1031947
95      1  G264    2969968
96      1  G267    6435838
97      1  G268   12231208
98      1  G271   22617986
99      1  G272   20047679
100     1  G273   32621888
101     1  G276    7744795
102     1  G277   15656485
103     1  G281    2890376
104     1  G282    9452730
105     1  G285      17764
106     1  G292        145
107     1  G293     975366
108     1  G295   52466267
109     1  G296   58796429
110     1  G300   53925059
111     1  G301   14833036
112     1  G302  157454812
113     1  G303  249989359
114     1  G304  129929840
115     1  G305   16815640
116     1  G307    5557618
117     1  G308   34119652
118     1  G310   50937244
119     1  G311   21610076
120     1  G312   31047427
121     1  G314    4647785
122     1  G316    2879526
123     1  G317    8885012
124     1  G318    3541308
125     1  G319    6494469
126     1  G320      32972
127     1  G322      95544
128     1  G328   19750064
129     1  G329    2523978
130     1  G331  280778665
131     1  G332   10302944
132     1  G333    9037660
133     1  G334   32829806
134     1  G340        141
135     1  G342       1810
136     1  G344   61044688
137     1  G345     122067
138     1  <NA>      34485
139     1  G487      29445
140     1  G488    1627710
141     1  G490    1057044
142     1  G493     149884
143     1  G494      63285
144     1  G496    2460937
145     1  G497   36761514
146     1  G498     235482
147     1  G499    1337888
148     1  G505      38749
149     1  G506   12108147
150     1  G507    1651487
151     1  G510    2780209
152     1  G517    2923878
153     1  G520     159278
154     1  G521    6749000
155     1  G524    4171256
156     1  G526    5347883
157     1  G527    2415563
158     1  G531    4954213
159     1  G533      27800
160     1  G537   31861987
161     1  G538   28533870
162     1  G541      17474
163     1  G552     622081
164     1  G553    1556404
165     1  G556    3506315
166     1  G557     766796
167     1  G559    3630979
168     1  G563    1879343
169     1  G565   15458270
170     1  G566    7087237
171     1  G567    2083137
172     1  G569   12323101
173     1  G570   48147558
174     1  G571   21950162
175     1  G573    3151912
176     1  G582         72
177     1  G584        494
178     1  G599       8555
179     1  G600   75975945
180     1  G624   18823083
181     1  G632     856285
182     1  G638      22081
183     1  <NA>      72101
184     1  G647      13006
185     1  G648    8080597
186     1  G650  102025348
187     1  G652   13870276
188     1  G653    4957184
189     1  G654     431717
190     1  G655      48200
191     1  G657      20309
192     1  G658      19813
193     1  G660     145660
194     1  G661      84994
195     1  G662     585780
196     1  G663     156765
197     1  G671         92
198     1  G673   12688634
199     1  G675     637289
200     1  G677    2355174
201     1  G678    1564016
202     1  G679   41302454
203     1  G680    8852521
204     1  G681       1449
205     1  G742   70688608
206     1  G744   32878626
207     1  G749   11463987
208     1  G751   11554959
209     1  G753    1577598
210     1  G762      16308
211     1  G767      15541
212     1  G770    8341306
213     1  G775    2058132
214     1  G777        678
215     1  G784    2571528
216     1  G788      83642
217     1  G790       5904
218     1  G791     715143
219     1  G793       1863
220     1  G798     249496
221     1  G800    1218025
222     1  G804    1857881
223     1  G805     179860
224     1  G806   29756574
225     1  G807     211300
226     1  G809     125278
227     1  G831   15599777
228     1  G838      71280
229     1  G839    2132070
230     1  G841      22947
231     1  G842        543
232     1  <NA>       1905
233     1  G844       1256
234     1  G849   17649480
235     1  G851    8831669
236     1  G853     410058
237     1  G887   12798136
238     1  G888   74388540
239     1  G889   14032634
240     1  G891     734242
241     1  G893     781291
242     1  G895     985356
243     1  G899   49724821
244     1  G900   82177444
245     1  G902    3515835
246     1  G903     784990
247     1  G904   12867370
248     1  G905    8874121
249     1  G906    5957397
250     1  G907   26826936
251     1  G908    5294676
252     1  G909     157080
253     1  G910          1
254     1  G913     470317
255     1  G914     182015
256     1  G915      64089
257     1  G916      58973
258     1  G918    3079788
259     1  G919   17569442
260     1  G920   19521881
261     1  G921  127545194
262     1  G922   55870981
263     1  G951    9652010
264     1  G962       1028
265     1  G978       8164
266     1  G979    5232481
267     1  G982    1460121
268     1  G983    2088938
269     1  G984    1224900
270     1  G986     297797
271     1  G987   43673721
272     1  G988   55834273
273     1  G989   41713198
274     1  G990   57880071
275     1  G991   10573411
276     1  G999    1001288
277     1 G1171      85922
278     1 G1172    5937274
279     1 G1176     728020
280     1 G1177    6032553
281     1 G1178      15857
282     1 G1179       2752
283     1 G1180         84
284     1 G1181      74483
285     1 G1183    1077321
286     1 G1184    6912107
287     1 G1211   16162702
288     1 G1213    9289225
289     1 G1214    1377895
290     1  <NA>    7823152
291     1  <NA>     138128
292     1 G8801  363378738
293     1 G8802 1729468222
294     1 G8803   62876155
295     1 G8804  507693068
296     1 G8805  144250899
297     1 G8806    3034869
298     1 G8807  538155783
299     1 G8808        351
300     1 G8809   23155793
301     1 G9901     676282

# need broad categories to reclassify to binary

# load spreadsheet with broad categories
cnvc_cat <- read_excel(cnvc_cat_file)
head(cnvc_cat)

# Inspect current levels (RAT)
lev <- levels(r_cnvc_mask_trim)[[1]]   # data.frame or NULL if not categorical
if (is.null(lev)) {
  stop("r_cnvc_mask_trim does not have a levels table (RAT). Use Option B to ratify first.")
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
levels(r_cnvc_mask_trim)[[1]] <- lev_joined

# check that new field is present
cats(r_cnvc_mask_trim)

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
cnvc_trim_binary <- classify(r_cnvc_mask_trim, rcl = lut, others = NA)

# Verify results and plot
print(cnvc_trim_binary)
plot(cnvc_trim_binary, main = "Binary reclass by Broad_Cat (1 = grass)")

# Check unique output values (should be 0, 1, and NA)
unique_vals <- unique(values(cnvc_trim_binary, mat = FALSE))
print(unique_vals)

# Save the raster as 8-bit unsigned integer
cnvc_binary_out_path <- file.path(cnvc_cropped_dir, "cnvc_trim_binary.tif")
writeRaster(cnvc_trim_binary, cnvc_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

############### do for aligned raster 

# need broad categories to reclassify to binary

# load cropped raster
#cnvc_algn_file <- file.path(cnvc_aligned_dir, "cnvc_algn.tif")
#cnvc_algn <- rast(cnvc_algn_file)
plot(cnvc_algn)

# print frequency table
freq_table <- freq(cnvc_algn) 
print(freq_table) 

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
writeRaster(cnvc_algn_binary, cnvc_binary_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

# load cnvc raster
cnvc_algn_file <- file.path(cnvc_aligned_dir, "cnvc_algn_good.tif")
cnvc_algn <- rast(cnvc_algn_file)
plot(cnvc_algn)
cnvc_algn

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
print(cnvc_algn_nattam)
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
unique_vals_n <- unique(values(cnvc_algn_native, mat = FALSE))
print(unique_vals_n)

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
unique_vals_rud <- unique(values(cnvc_algn_ruderal, mat = FALSE))
print(unique_vals_rud)

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
unique_vals_past <- unique(values(cnvc_algn_past, mat = FALSE))
print(unique_vals_past)

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
cnvc_file <- file.path(cnvc_aligned_dir, "cnvc_algn_good.tif")
cnvc <- rast(cnvc_file)

crop_ext <- ext(sa)

# crop raster to study area
cnvc_sa <- crop(cnvc, crop_ext)
cnvc_sa
plot (cnvc_sa)

# mask to study area
cnvc_masked <- mask(cnvc_sa, sa)
cnvc_masked
plot(cnvc_masked)

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

