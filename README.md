# Tomorrow’s Prairies: Baselines and projections for grasslands in Canada

Purpose: To contribute to mitigation of climate change and biodiversity loss by researching historical baselines and projections in the extent and condition of Canadian grasslands.

Below is the folder structure for this project, showing the organization of data, scripts, and outputs.
```text
tomorrows-prairies/
├── existingdatasets/            # Source datasets from external providers
│   ├── Climate/                 # Historical and projected climate data
│   ├── Cropland/                # Cropland extent and properties
│   ├── Landcover/               # Land use land cover classification datasets
│   ├── Management/              # Protected area and tenure data
│   ├── Socio-economic/          # Protected area and tenure data
│   ├── Soil/                    # Protected area and tenure data
│   └── Topography/              # Protected area and other management data
├── processeddatasets/           # Datasets prepared for analysis and modelling
│   ├── boundaries/              # Jurisdictional and ecological boundaries
│   ├── climate/                 
│   ├── cropland/                
│   ├── landcover/                  
│   ├── management/
│   ├── predictions/                  
│   ├── predictors/                  
│   ├── response/               
│   ├── socio-economic/               
│   ├── soil/               
│   ├── studyarea/               # Spatial extents of analyses
│   └── topography/                                        
├── scripts/                     # Code for data processing and modelling
│   ├── classification/          # Scripts for grassland classification
│   └── landcover/               # Scripts for land cover processing
├── output/                      # Model outputs and results visualizations
├── .gitignore                   # Git exclusions for files and folders
├── README.md                    # Project overview and instructions
└── LICENSE                      # Code licensing information
```
