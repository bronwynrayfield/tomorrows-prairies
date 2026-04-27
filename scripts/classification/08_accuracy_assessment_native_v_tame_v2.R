
#######################################################################################
##              PRELIMINARY ACCURACY ASSESSMENT: NATIVE AND TAME CLASSES             ##
#######################################################################################

##----------------------- Alison Long - April 2026 -------------------------

# step 8 (multiclass): assess accuracy of datasets mapping grasslands
# against ground-truthed point data (2015 to present), with 3 classes:
#   0 = not grass, 1 = tame grass, 2 = native grass
#
# THIS VERsION COMPARES:
#   1) native and tame points to native and tame raster values
#   2) native vs non-native (tame, crop, shrub, wetland) points to native vs non-native (tame, not grass) raster values
#      and tame vs non-tame (native, crop, shrub, wetland) points to tame vs non-tame (native, not grass) raster values
#
# Key changes from binary version:
#   - Uses LC_class3 (0/1/2) instead of LC_grass (0/1)
#   - Rasters expected as 0/1/2
#   - Time-aligned rules:
#         native (2): latest allowed_year <= Yer_Cll
#         tame   (1): allowed_year == Yer_Cll
#   - 3×3 confusion matrices and per-class metrics
#   - 2×2 confusion matrices and per-class metrics (native vs tame only)
#   - More accuracy metrics (not just recall)
