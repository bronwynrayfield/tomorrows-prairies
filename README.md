# Grassland Futures in the Canadian Prairies
This project analyzes future potential grassland loss in the Canadian prairies under different socioeconomic and climate scenarios. The analysis explores the consequences for biodiversity and ecosystem carbon storage to inform conservation and land-use planning.

Below is the folder structure for this project, showing the organization of data, scripts, and outputs.
```text
tomorrows-prairies/
├── existingdatasets/            # Source datasets from external providers
│   ├── Biodiversity/            # Species and habitat data
│   ├── Boundaries/              # Jurisdictional and ecological boundaries
│   ├── Carbon/                  # Soil or biomass carbon data
│   ├── Climate/                 # Historical and projected climate data
│   ├── Cropland/                # Cropland extent and transitions
│   ├── Grassland/               # Native grassland distribution
│   ├── Landcover/               # Land use land cover classification datasets
│   └── Management/              # Protected area and other management data
├── processeddatasets/           # Datasets prepared for analysis and modeling
│   ├── Grassland/               # Cleaned and formatted grassland inputs
│   └── Studyarea/               # Clipped and harmonized spatial extent
├── scripts/                     # Code for data processing and modeling
│   ├── Classification/          # Scripts for land cover classification
│   └── Landcover/               # Scripts for land cover processing
├── models/                      # Model configurations and outputs
├── results/                     # Final results and visualizations
│   ├── Figures/                 # Plots, charts, and maps
│   └── Tables/                  # Summary tables and statistics
├── docs/                        # Project documentation and notes
├── .gitignore                   # Git exclusions for files and folders
├── README.md                    # Project overview and instructions
└── LICENSE                      # Code licensing information
```
