https://www12.statcan.gc.ca/census-recensement/alternative_alternatif.cfm?l=eng&dispext=zip&teng=lpc_000a21a_e.zip&k=%20%20%20%20%204664&loc=//www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpc_000a21a_e.zip

https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?Year=21

https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpc_000a21a_e.zip
https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/gpc_000a11a_e.zip

source(file.path("scripts/landcover", "0a-setup.R"))


popctr_filenames <- c("2011" = popctr2011_filename, "2021" = popctr2021_filename)
pc_urls <- c(
  "2011" = "https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/gpc_000a11a_e.zip",
  "2021" = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpc_000a21a_e.zip"
)

pc_dir <- file.path("existingdatasets", capitalize_first(popctr_filepath))
dir.create(pc_dir, recursive = TRUE, showWarnings = FALSE)

for (y in names(pc_urls)) {
  f <- file.path(pc_dir, paste0(popctr_filenames[y], ".zip"))
  download.file(pc_urls[y], f, mode = "wb")
  unzip(f, exdir = pc_dir)
}
