# Wetland Phenology Analysis - San Francisco Bay Area

Analysis of wetland phenology metrics for the San Francisco Bay area (Suisun region) using time series of Sentinel-2 satellite data.

Paper: Remotely sensed phenology reveals environmental and management controls on coastal wetland plant communities

## Overview

This project implements a comprehensive workflow for analyzing land surface phenology (LSP) metrics in wetland ecosystems. The analysis focuses on Plant Functional Types (PFTs) in the Suisun Bay area, extracting phenological patterns from Sentinel-2 time series data and applying machine learning techniques for classification and pattern recognition.

### Study Area
- **Location**: Suisun Bay, San Francisco Bay Area, California
- **Coordinate System**: UTM Zone 10N (EPSG:32610)
- **Time Period**: 2018-2020 (Sentinel-2 time series)

### Plant Functional Types (PFTs) Analyzed
The following PFT classes are included in the analysis:
- **Annual Grassland** (Class 1)
- **Perennial Grassland** (Class 4)
- **Pickleweed-Cordgrass-Saltbush** (Class 5)
- **Tule-Cattail** (Class 6)
- **Wet Meadows** (Class 9)
- **Other** (Class 3)

## Project Structure

```
Wetland_Phenology_Python/
├── 0_DataCube_getData.ipynb          # Extract Sentinel-2 data cube
├── 0_DEM.ipynb                        # Digital Elevation Model processing
├── 0_phenology_water_mask.ipynb      # Apply water mask to phenology data
├── 0_preparePFT.ipynb                 # Prepare Plant Functional Type data
├── 1_Autoencoder.ipynb                # Autoencoder for dimensionality reduction
├── 1_PCA.ipynb                        # Principal Component Analysis
├── 2_Get_Data.ipynb                   # Extract data at stratified sampling points
├── 3_KMeans.ipynb                     # K-Means clustering analysis
├── 4_UnifydataFrames.ipynb            # Unify and merge data frames
├── 5_clusterClasesAnalysis.ipynb      # Cluster and class analysis
├── 6_Variable_importance.ipynb       # Variable importance analysis
├── 7_diversity.ipynb                  # Diversity metrics analysis
├── get_random_points.ipynb            # Generate random sampling points
├── PhenoPy_example.ipynb              # PhenoPy library examples
├── stratified_points_in_polygon.py   # Stratified point sampling utility
├── data/                              # Input data files
│   ├── dem_clip.csv
│   ├── dem.csv
│   ├── KernelPCA_clip.csv
│   ├── lspAll.csv
│   ├── PFT_Classes_Info.txt
│   └── rmseAll.csv
└── outputs/                           # Analysis outputs
    ├── explained_variance_kernelPCA.txt
    ├── explained_variance_PCA.txt
    ├── mrpp.RData
    └── plsImp_veg.RData
```

## Workflow

The analysis follows a sequential workflow:

1. **Data Acquisition** (`0_DataCube_getData.ipynb`)
   - Extract Sentinel-2 time series data
   - Create data cube for analysis

2. **Preprocessing**
   - **DEM Processing** (`0_DEM.ipynb`): Process digital elevation model data
   - **Water Mask** (`0_phenology_water_mask.ipynb`): Apply water mask to phenology rasters
   - **PFT Preparation** (`0_preparePFT.ipynb`): Prepare Plant Functional Type classifications

3. **Dimensionality Reduction**
   - **PCA** (`1_PCA.ipynb`): Principal Component Analysis
   - **Autoencoder** (`1_Autoencoder.ipynb`): Deep learning-based dimensionality reduction

4. **Data Extraction** (`2_Get_Data.ipynb`)
   - Create stratified sampling points within PFT classes
   - Extract phenology variables and DEM at sampling points

5. **Clustering** (`3_KMeans.ipynb`)
   - K-Means clustering analysis on phenology data

6. **Data Integration** (`4_UnifydataFrames.ipynb`)
   - Merge and unify data frames from different sources

7. **Analysis**
   - **Cluster Analysis** (`5_clusterClasesAnalysis.ipynb`): Analyze cluster patterns and classes
   - **Variable Importance** (`6_Variable_importance.ipynb`): Assess importance of phenology variables
   - **Diversity Analysis** (`7_diversity.ipynb`): Calculate diversity metrics

## Phenology Metrics

The analysis extracts various land surface phenology (LSP) metrics including:

- **SOS**: Start of Season (DOY)
- **POS**: Peak of Season (DOY)
- **EOS**: End of Season (DOY)
- **LOS**: Length of Season
- **MSP**: Mid Spring (DOY)
- **MAU**: Mid Autumn (DOY)
- **AOS**: Amplitude of Season
- **IOS**: Integral of Season [SOS-EOS]
- **ROG**: Rate of Greening [slope SOS-POS]
- **ROS**: Rate of Senescence [slope POS-EOS]
- **SW**: Skewness of Growing Season [SOS-EOS]

## Dependencies

Key Python libraries used:
- `rasterio` - Geospatial raster I/O
- `geopandas` - Geospatial data operations
- `xarray` / `rioxarray` - Multi-dimensional arrays and geospatial data
- `pandas` / `numpy` - Data manipulation and numerical computing
- `scikit-learn` - Machine learning (PCA, KMeans, etc.)
- `seaborn` / `matplotlib` - Data visualization
- `rasterstats` - Zonal statistics
- `tqdm` - Progress bars

## Usage

1. **Setup**: Ensure all required dependencies are installed
2. **Data Preparation**: Place your Sentinel-2 data and PFT shapefiles in the appropriate directories
3. **Run Workflow**: Execute notebooks sequentially (0 → 1 → 2 → ... → 7)
4. **Results**: Check the `outputs/` directory for analysis results

## Author

**Javier Lopatin**  
Date: 2025-11-10  
Version: 1.0

## License

[Add license information if applicable]

## Citation

If you use this code in your research, please cite appropriately.
