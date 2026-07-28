
###################################################################################
##                 SETUP FOR GRASSLAND CONSENSUS MAPPING ANALYSIS                ##
###################################################################################

##----------------------- Alison Long - July 2026 -------------------------##

# step 8b: create redundancy penalty raster to assist creation of consensus grassland map 

# ===========================================
# SETUP & INPUTS
# ===========================================

# Temp dirs
tmp_root <- file.path(Sys.getenv("HOME"), "r_tmp")
#dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Load packages
library(dplyr)

# Force a fresh run instead of silently re-using old scratch outputs
FORCE_REBUILD <- TRUE

# Directories
out_dir    <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/results/consensus_map"

# Files
agreement_bin_global_file <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/??.csv"
agreement_nt_global_file  <- "s3://sagemaker-gst-stage.sharing/alison-long/tomorrows-prairies/processed/??.csv"
agreement_bin_ecodist_file
agreement_nt_ecodist_file

# Load
agreement_bin_global <- aws.s3::s3read_using(read_csv, object = agreement_bin_global_file, show_col_types = FALSE) %>%
  rename(file = file_s3)

agreement_nt_global <- aws.s3::s3read_using(read_csv, object = agreement_nt_global_file, show_col_types = FALSE) %>%
  rename(file = file_s3)

agreement_bin_ecodist_file
agreement_nt_ecodist_file

#--------------------------------------------------
# Convert agreement into a penalty
#
# agreement = 0  -> penalty = 1
# agreement = 1  -> penalty = floor
#
# floor prevents highly redundant rasters from
# becoming nearly useless
#--------------------------------------------------

agreement_to_penalty <- function(
    agreement,
    floor = 0.25
) {
  
  floor + (1 - floor) * (1 - agreement)
  
}

# ===========================================
# BINARY PENALTIES
# ===========================================

#--------- Global binary redundancy ---------
calc_global_bin_penalty <- function(
    agreement_bin_global,
    floor = 0.25
) {
  
  bind_rows(
    
    agreement_bin_global %>%
      transmute(
        raster = raster_i,
        agreement = mcc
      ),
    
    agreement_bin_global %>%
      transmute(
        raster = raster_j,
        agreement = mcc
      )
    
  ) %>%
    group_by(raster) %>%
    summarise(
      
      mean_agreement =
        mean(
          agreement,
          na.rm = TRUE
        ),
      
      n_pairs = n(),
      
      .groups = "drop"
      
    ) %>%
    mutate(
      
      penalty =
        floor +
        (1 - floor) *
        (1 - mean_agreement)
      
    )
  
}

#--------- Ecodistrict-level binary redundancy ---------
calc_ecodist_bin_penalty <- function(
    agreement_bin,
    floor = 0.25
) {
  
  bind_rows(
    
    agreement_bin_ecodist %>%
      transmute(
        Ecds_ID,
        raster = raster_i,
        agreement = mcc
      ),
    
    agreement_bin_ecodist %>%
      transmute(
        Ecds_ID,
        raster = raster_j,
        agreement = mcc
      )
    
  ) %>%
    group_by(
      Ecds_ID,
      raster
    ) %>%
    summarise(
      
      mean_agreement =
        mean(
          agreement,
          na.rm = TRUE
        ),
      
      .groups = "drop"
      
    ) %>%
    mutate(
      
      penalty =
        agreement_to_penalty(
          mean_agreement,
          floor = floor
        )
      
    )
  
}

#------------ Combined binary penalty -------------
calc_combined_bin_penalty <- function(
    agreement_bin_ecodist,
    agreement_bin_global,
    floor = 0.25
) {
  
  global_tbl <-
    calc_global_bin_penalty(
      agreement_bin_global,
      floor
    ) %>%
    rename(
      global_penalty = penalty
    )
  
  ecodist_tbl <-
    calc_ecodist_bin_penalty(
      agreement_bin_ecodist,
      floor
    ) %>%
    rename(
      ecodist_penalty = penalty
    )
  
  ecodist_tbl %>%
    left_join(
      global_tbl,
      by = "raster"
    ) %>%
    mutate(
      
      penalty =
        sqrt(
          global_penalty *
            ecodist_penalty
        )
      
    ) %>%
    select(
      Ecds_ID,
      raster,
      penalty
    )
  
}

# ===========================================
# TRICLASS PENALTIES
# ===========================================

#--------- Global triclass redundancy ----------
calc_global_nt_penalty <- function(
    agreement_nt_global,
    floor = 0.25
) {
  
  bind_rows(
    
    agreement_nt_global %>%
      transmute(
        raster = raster_i,
        agreement = mcc_multi
      ),
    
    agreement_nt_global %>%
      transmute(
        raster = raster_j,
        agreement = mcc_multi
      )
    
  ) %>%
    group_by(raster) %>%
    summarise(
      
      mean_agreement =
        mean(
          agreement,
          na.rm = TRUE
        ),
      
      n_pairs = n(),
      
      .groups = "drop"
      
    ) %>%
    mutate(
      
      penalty =
        floor +
        (1 - floor) *
        (1 - mean_agreement)
      
    )
  
}

#-------- Ecodistrict-level triclass redundancy ---------
calc_ecodist_nt_penalty <- function(
    agreement_nt_ecodist,
    floor = 0.25
) {
  
  bind_rows(
    
    agreement_nt_ecodist %>%
      transmute(
        Ecds_ID,
        raster = raster_i,
        agreement = mcc_multi
      ),
    
    agreement_nt_ecodist %>%
      transmute(
        Ecds_ID,
        raster = raster_j,
        agreement = mcc_multi
      )
    
  ) %>%
    group_by(
      Ecds_ID,
      raster
    ) %>%
    summarise(
      
      mean_agreement =
        mean(
          agreement,
          na.rm = TRUE
        ),
      
      .groups = "drop"
      
    ) %>%
    mutate(
      
      penalty =
        agreement_to_penalty(
          mean_agreement,
          floor = floor
        )
      
    )
  
}

#---------- Combined triclass redundancy ----------
calc_combined_nt_penalty <- function(
    agreement_nt_ecodist,
    agreement_nt_global,
    floor = 0.25
) {
  
  global_tbl <-
    calc_global_nt_penalty(
      agreement_nt_global,
      floor
    ) %>%
    rename(
      global_penalty = penalty
    )
  
  ecodist_tbl <-
    calc_ecodist_nt_penalty(
      agreement_nt_ecodist,
      floor
    ) %>%
    rename(
      ecodist_penalty = penalty
    )
  
  ecodist_tbl %>%
    left_join(
      global_tbl,
      by = "raster"
    ) %>%
    mutate(
      
      penalty =
        sqrt(
          global_penalty *
            ecodist_penalty
        )
      
    ) %>%
    select(
      Ecds_ID,
      raster,
      penalty
    )
  
}

# ===========================================
# BUILD TABLES
# ===========================================

# Option 1: combined
bin_penalty_tbl <-
  calc_combined_bin_penalty(
    agreement_bin
  )

nt_penalty_tbl <-
  calc_combined_nt_penalty(
    agreement_nt
  )

# Option 2: ecodistrict only
bin_penalty_tbl <-
  calc_ecodist_bin_penalty(
    agreement_bin
  )

nt_penalty_tbl <-
  calc_ecodist_nt_penalty(
    agreement_nt
  )

# Option 3: global only
global_bin <- calc_global_bin_penalty(
  agreement_bin_global
)

global_nt <- calc_global_nt_penalty(
  agreement_nt_global
)

# ===========================================
# BUILD RASTER STACKS
# ===========================================

#------- Ecodistrict-level or combined ------

# binary
penalty_bin_stack <- build_weight_stack(
  accuracy_table = bin_penalty_tbl %>%
    rename(
      raster_name = raster
    ),
  stack_names = names(bin_stack),
  value_col = "penalty",
  ecodist_r = ecodist_r,
  out_prefix = "penalty_bin"
)

# native vs tame
penalty_nt_stack <- build_weight_stack(
  accuracy_table = nt_penalty_tbl %>%
    rename(
      raster_name = raster
    ),
  stack_names = names(nt_stack),
  value_col = "penalty",
  ecodist_r = ecodist_r,
  out_prefix = "penalty_nt"
)

#------- Global only --------

# binary
penalty_bin_stack <- build_weight_stack(
  accuracy_table = global_bin %>%
    rename(
      raster_name = raster
    ),
  stack_names = names(bin_stack),
  value_col = "penalty",
  out_prefix = "penalty_bin"
)

# native vs tame
penalty_nt_stack <- build_weight_stack(
  accuracy_table = global_nt %>%
    rename(
      raster_name = raster
    ),
  stack_names = names(nt_stack),
  value_col = "penalty",
  out_prefix = "penalty_nt"
)

# ===========================================
# EXPORT
# ===========================================

# binary
pen_bin_out <- "/mnt/tmp/penalty_binary.tif"

writeRaster(
  penalty_bin_stack,
  pen_bin_out,
  overwrite = TRUE,
  wopt = list(
    datatype = "INT1U",
    gdal = c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = pen_bin_out,
  object = paste0(out_dir, "/penalty_binary.tif")
)

# triclass
pen_nt_out <- "/mnt/tmp/penalty_natvtam.tif"

writeRaster(
  penalty_nt_stack,
  pen_nt_out,
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = c("COMPRESS=ZSTD","TILED=YES","BIGTIFF=YES")
  )
)

aws.s3::put_object(
  file   = pen_nt_out,
  object = paste0(out_dir, "/penalty_natvtam.tif")
)


