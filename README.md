# Tomorrow’s Prairies: Baselines and scenario-based projections for grasslands in Canada

Tomorrow's Prairies is a collaborative research project mapping grasslands and projecting changes in biodiversity and carbon based on scenarios of grassland conversion.

Below is the folder structure for this project, showing the organization of data, scripts, and outputs.
```text
tomorrows-prairies/
├── existingdatasets/            # Source datasets from external providers
│   ├── Climate/                 # Historical and projected climate data
│   ├── Cropland/                # Cropland extent and properties
│   ├── Landcover/               # Land use land cover classification datasets
│   ├── Management/              # Protected area and tenure data
│   ├── Socio-economic/          # 
│   ├── Soil/                    # 
│   └── Topography/              # 
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
│   └── landcover/               # Scripts for land cover modelling
├── output/                      # Model outputs and results visualizations
│   ├── classification/          # Output for grassland classification
│   └── landcover/               # Output for land cover modelling
├── .gitignore                   # Git exclusions for files and folders
├── README.md                    # Project overview and instructions
└── LICENSE                      # Code licensing information
```
