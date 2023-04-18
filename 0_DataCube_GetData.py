
# For testing inside the Datacube Platform!
import datacube
import numpy as np
import sys, os
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from datacube.utils.masking import make_mask, mask_invalid_data, describe_variable_flags
from datacube.utils.rio import configure_s3_access
from dask.distributed import Client

sys.path.append('/home/javierlopatin/Documents/GitHub/PhenoPY')
from phenoxr.phenoXr import Pheno

# client = Client()
# configure_s3_access(aws_unsigned=False, requester_pays=True, client=client)
configure_s3_access(aws_unsigned=False, requester_pays=True)

dc = datacube.Datacube(app="Pheno_test")
dc.list_products().name

veg_proxy = 'NDVI'
dates = ('2019-01-01', '2020-12-31')
inProduct = 'landsat8_c2l2_sr'
resolution = 30

# load shapefile 
shp = gpd.read_file('../data/shapefiles/shapefile.shp')

# get shapefile bounding box
shp_bouning_box = shp.total_bounds

# get current working directory
cwd = os.getcwd()

# get coordinates of the shapefile bounding box
X_min = shp_bouning_box[0]
X_max = shp_bouning_box[2]
Y_min = shp_bouning_box[1]
Y_max = shp_bouning_box[3]

# range of X and Y coordinates
X_range = X_max - X_min
Y_range = Y_max - Y_min

# list of the X and Y coordinates with a grid of 4 blocks
X_list = np.linspace(X_min, X_max, 4)
Y_list = np.linspace(Y_min, Y_max, 4)

# query to load data using the list of X and Y coordinates
query = {
    "x": (X_list[0], X_list[1]),
    "y": (Y_list[0], Y_list[1]),
    "time": ('2018-01-01', '2020-12-31'),
    "output_crs": "EPSG:32610",
    "resolution": (-10, 10),
    "dask_chunks": {"time": 1},
    "group_by":"solar_day",
    'product': "s2_l2a",
    'skip_broken_datasets': True,
    'crs': 'EPSG:32610',

    'measurements': ['B04', 'B08'],
}

ds_s2 = dc.load(
    **query,
)

# describe the SCL variable
from datacube.utils import masking
masking.describe_variable_flags(ds_s2.SCL).values

# Multiple flags are combined as logial OR using the | symbol
cloud_free_mask = (
    masking.make_mask(ds_s2.SCL, qa="vegetation") | 
    masking.make_mask(ds_s2.SCL, qa="bare soils") |
    masking.make_mask(ds_s2.SCL, qa="water") |
    masking.make_mask(ds_s2.SCL, qa="snow or ice")
)

cloud_free_mask

# Calculate proportion of good pixels
valid_pixel_proportion = cloud_free_mask.sum(dim=("x", "y"))/(cloud_free_mask.shape[1] * cloud_free_mask.shape[2])
valid_threshold = 0.7
observations_to_keep = (valid_pixel_proportion >= valid_threshold)

# Mask the data
ds_s2_valid = ds_s2.where(cloud_free_mask)

# only keep observations above the good pixel proportion threshold
# The .compute() step means the values will be loaded into memory. This step may take some time
ds_s2_keep = ds_s2_valid.sel(time=observations_to_keep).compute()
ds_s2_keep

# plot example band of one observation
ds_s2_keep.B04.isel(time=0).plot()

# order ds_s2_keep by day of year
ds_s2_keep = ds_s2_keep.sortby(ds_s2_keep.time.dt.dayofyear)

# calculate EVI Index for each observation using the formula
ds_s2_keep['EVI'] = (2.5 * (ds_s2_keep.B08 - ds_s2_keep.B04)) / (ds_s2_keep.B08 + (2.4 * ds_s2_keep.B04) + 1)

# plot example EVI Index of one observation
ds_s2_keep.EVI.isel(time=0).plot()

# time series plot of EVI Index for the central pixel of the image, using dots to represent the observations
ds_s2_keep.EVI.isel(x=ds_s2_keep.x.shape[0]//2, y=ds_s2_keep.y.shape[0]//2).plot(marker='o')

# PhenoShape
ans = ds_s2_keep.EVI.pheno.PhenoShape()
ans

# plot the example date
ds_s2_keep.EVI.isel(time=0).plot()

# unstck a dictionary of xarrays to a single xarray
ans = xr.Dataset(ans)    





# list of files
files = glob.glob(folder + pattern)

# load and merge list of path files using rioxarray
from rioxarray.merge import merge_arrays

pheno1 = rio.open_rasterio(files[0])
pheno2 = rio.open_rasterio(files[1])
pheno3 = rio.open_rasterio(files[2])
pheno4 = rio.open_rasterio(files[3])
pheno5 = rio.open_rasterio(files[4])

phenoAll = merge_arrays([pheno1, pheno2, pheno3, pheno4, pheno5])

    

      