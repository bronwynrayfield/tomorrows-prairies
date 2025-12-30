# Tomorrow’s Prairies: Baselines and projections for grasslands in Canada

Purpose: To contribute to mitigation of climate change and biodiversity loss by researching historical baselines and projections in the extent and condition of Canadian grasslands.

Below is the folder structure for this project, showing the organization of data, scripts, and outputs.
```text
tomorrows-prairies/
├── existingdatasets/            # Source datasets from external providers
│   ├── biodiversity/            # Species, habitat, connectivity, ecosystem condition data
│   ├── boundaries/              # Jurisdictional and ecological boundaries
│   ├── carbon/                  # Parameterization, calibration, and validation carbon data
│   ├── climate/                 # Historical and projected climate data
│   ├── cropland/                # Cropland extent and properties
│   ├── grassland/               # Grassland classification datasets
│   ├── landcover/               # Land use land cover classification datasets
│   └── management/              # Protected area and other management data
├── processeddatasets/           # Datasets prepared for analysis and modelling
│   ├── studyarea/               # Spatial extents of analyses
│   ├── biodiversity/            
│   ├── boundaries/              
│   ├── carbon/                  
│   ├── climate/                 
│   ├── cropland/                
│   ├── grassland/               
│   ├── landcover/               
│   └── management/                                        
├── scripts/                     # Code for data processing and modelling
│   ├── classification/          # Scripts for grassland classification
│   └── landcover/               # Scripts for land cover processing
├── models/                      # Model configurations and outputs
├── results/                     # Final results and visualizations
│   ├── figures/                 # Plots, charts, and maps
│   └── tables/                  # Summary tables and statistics
├── docs/                        # Project documentation and notes
├── .gitignore                   # Git exclusions for files and folders
├── README.md                    # Project overview and instructions
└── LICENSE                      # Code licensing information
```
