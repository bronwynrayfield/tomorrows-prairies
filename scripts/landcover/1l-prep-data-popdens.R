# =============================================================================
# Purpose: download census dissemination block boundaries and
#          Geographic Attribute Files (population density predictor)
#
# 2011: https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2011-eng.cfm
# Language=English, Format=ArcGIS(.shp) Boundary Files=Dissemination blocks, Cartographic Boundary File 
#
# 2021: https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21
# Language=English, Type=Cartographic Boundary Files
# Statistical boundaries=Dissemination blocks, Format=Shapefile
#
# Automatic downloads didn't work so I manually downloaded from the URLs
# =============================================================================
source(file.path("scripts/landcover", "0a-setup.R"))

options(timeout = 3600)

db_filenames  <- c("2011" = db2011_filename,  "2021" = db2021_filename)
gaf_filenames <- c("2011" = gaf2011_filename, "2021" = gaf2021_filename)
db_urls <- c(
  "2011" = "https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/gdb_000b11a_e.zip",
  "2021" = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/ldb_000b21a_e.zip"
)
gaf_urls <- c(
  "2011" = "https://www12.statcan.gc.ca/census-recensement/2011/geo/ref/files-fichiers/2011_92-151_XBB_txt.zip",
  "2021" = "https://www12.statcan.gc.ca/census-recensement/2021/geo/aip-pia/attribute-attribs/files-fichiers/2021_92-151_X.zip"
)

popdens_dir <- file.path("existingdatasets", capitalize_first(popdens_filepath))
dir.create(popdens_dir, recursive = TRUE, showWarnings = FALSE)

for (y in names(db_urls)) {
  f <- file.path(popdens_dir, paste0(db_filenames[y], ".zip"))
#  download.file(db_urls[y], f, mode = "wb")
  unzip(f, exdir = popdens_dir)
  
  f <- file.path(popdens_dir, paste0(gaf_filenames[y], ".zip"))
#  download.file(gaf_urls[y], f, mode = "wb")
  unzip(f, exdir = popdens_dir)
}
