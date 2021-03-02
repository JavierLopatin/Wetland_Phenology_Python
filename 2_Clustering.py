#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:46:05 2020

@author: javier
"""

#%%
import os
import pandas as pd
import geopandas as gpd
import rasterio
import rasterio.mask
from rasterio.merge import merge
from rasterstats import zonal_stats
import numpy as np
from tqdm import tqdm       
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# set workind directory
os.chdir('/home/javier/Documents/SF_delta/Sentinel/Siusun')

# -----------------------------------------------------------------------------
#%%
def stratified_points_in_polygon(shp, number=100, dist=10, maxiter=1000):
    '''
    Generate stratified points inside polygons
    '''    
    from shapely.geometry import Point
    import geopandas as gpd
    import random
    import multiprocessing
    from joblib import Parallel, delayed
    from tqdm import tqdm

    # read shapefile and get labels
    shp = gpd.read_file(shp)

    def my_function(j,shp):
    # ancillary funciton to be passed to the parallel processing
        polygon = shp.geometry[j]
        points = []
        min_x, min_y, max_x, max_y = polygon.bounds
        i = 0
        while len(points) < number:
            point = Point(random.uniform(min_x, max_x), random.uniform(min_y, max_y))
            if polygon.contains(point):
                if i == 0:
                    points.append(point)
                if i != 0:
                    try:
                        points.append(point)
                        dist = [point.distance(x) for x in points][:-1]
                        # delete last point if distance to other points > dist
                        if np.any(np.array(dist) < dist):
                            points.pop() 
                    except i > maxiter:
                        sys.exit(1)                   
            i += 1       
        df = gpd.GeoDataFrame(crs=shp.crs , geometry=points)
        df['class'] = shp['NVCSName'][j]
        return df
        
    # list all geometries of the shapefile
    #myList = list(shp.geometry)
    
    # parallel processing through the list
    num_cores = multiprocessing.cpu_count()
    outlist = Parallel(n_jobs=num_cores)(delayed(my_function)(j,shp) for j in tqdm(range(len(shp))))
    
    return pd.concat(outlist) # concatenate all data
    
def ExtractPointValues(raster, shp, bandNames, labels, ID='FID'):
    from rasterstats import point_query
    import shapefile
    import multiprocessing
    from joblib import Parallel, delayed
    """ Extract raster values by a shapefile point mask.
    """

    # Shapefile management
    shape = shapefile.Reader(shp)
    records = pd.DataFrame(shape.records())
    n = pd.DataFrame(shape.fields)[0].values.tolist().index(ID)
    _id = records[n-1]

    # empty matrix to store results
    matrix = np.empty((len(records), len(bandNames)+1), dtype=object)
    matrix[:,0] = _id

    # Extract values
    for i in range(len(bandNames)):
        # get values with parallel processing
        num_cores = multiprocessing.cpu_count()
    # ancillary function to be passed to parallel processing
    def _funtion(i,shp,raster):
        stats = point_query(shp, raster, band=i+1)
        return pd.DataFrame(stats)
    # parallel processing
    stats = Parallel(n_jobs=num_cores)(delayed(_funtion)(i,shp,raster) 
                                       for i in tqdm(range(len(bandNames))))
    x = pd.concat(stats, axis=1) # concatenate all dataframes into one
    x.columns = bandNames # add colum names
    x.index = records[0]  # and shapefile ID to index
    # save data to .CSV file
    name = os.path.basename(raster)
    x.to_csv(name[:-4] + ".csv", index=False, header=True, na_rep='NA')

    return x  
    
# -----------------------------------------------------------------------------
#%% 
# load shapefile
# load Suisun area
shp = 'shapefiles/SuisunMarsh_NVCSName_diss.shp'
suisun = gpd.read_file(shp)
labels = suisun['NVCSName']; labels

# generate 50 stratified points within the first county polygon in geodata
samples = stratified_points_in_polygon(shp, number=50)
samples.to_file('samplepoints.shp')
#%%
# Load and merge rasters
import fiona

with fiona.open("shapefiles/SuisunMarsh_diss.shp", "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]


with rasterio.open('TSS_LSP_all.tif') as src:
    out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
    out_meta = src.meta
    out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})

with rasterio.open('nRMSE_all_fill_clip.tif') as src2:
    out_image2, out_transform2 = rasterio.mask.mask(src2, shapes, crop=True)

 
out_image.shape
out_image2.shape

merged = np.concatenate([out_image, out_image2], axis=0).squeeze()
merged.shape

out_meta.update(count = out_image.shape[0]+out_image2.shape[0])

with rasterio.open('LSP_nRMSE.tif', 'w', **out_meta) as dst:
       dst.write(merged)

    
#%%
# -----------------------------------------------------------------------------
# extract point values from rasters
# samples using the LSP and %RMSE metrics
LSPbands = ['SOS', 'POS', 'EOS', 'vSOS', 'vPOS', 'vEOS', 'LOS', 'MSP', 'MAU',
            'vMSP', 'vMAU', 'AOS', 'IOS', 'ROG', 'ROS', 'SW', 'nRMSE_all',
            'nRMSE_beginning','nRMSE_middle','nRMSE_end']

df_lsp1 = ExtractPointValues(raster='LSP_nRMSE.tif', shp='samplepoints.shp', 
                             bandNames=LSPbands, labels=labels, ID='class')
df_lsp1
#df_lsp2 = ExtractPointValues(raster='TSS_LSP_nRMSE_over03.tif', shp='samplepoints.shp', 
#                             bandNames=LSPbands, labels=labels, ID='class')

# samples using the phenoshape fits
xnew = np.linspace(5, 360, 46, dtype='int16')
bandNames = []
for i in range(46):
    a = "DOY_" + str(xnew[i])
    bandNames.append(a)

df_phen1 = ExtractPointValues(raster='TSS_phen_all.tif', shp='samplepoints.shp', 
                             bandNames=bandNames, labels=labels, ID='class')
df_phen2 = ExtractPointValues(raster='TSS_phen_over0.3.tif', shp='samplepoints.shp', 
                             bandNames=bandNames, labels=labels, ID='class')


# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#%%
#### clustering and PCA

# Parameterization
# ----------------------
inData = df_lsp1.dropna()
variables = 'lsp_nRMSE_all'
inRaster = 'LSP_nRMSE.tif'
# ----------------------

#%%

# correlation
from sklearn.preprocessing import StandardScaler

corrmat = inData.corr()
corrmat[corrmat == 1] = np.nan

plt.figure()
sns.heatmap(corrmat)
plt.title('Correlation among variables')
plt.tight_layout()
plt.savefig('corrplot_'+variables+'.png', res=300)

#%%
# preprocessing
# standardize variables
ss = StandardScaler()
X_std = ss.fit_transform(inData)

# Calculate the Eigenvalues and Eigenvectors from the Covariance Matrix
cov_mat = np.cov(X_std.T)
eigen_vals, eigen_vecs = np.linalg.eig(cov_mat)
'''
The eigen_vecs variable I’ve defined here represent the principal componants, 
or direction of maximum variance, whereas the eigen_vals is simply a scalar 
that defines their magnitude.
'''
tot = sum(eigen_vals)
# var_exp ratio is fraction of eigen_val to total sum
var_exp = [(i / tot) for i in sorted(eigen_vals, reverse=True)]
# calculate the cumulative sum of explained variances
cum_var_exp = np.cumsum(var_exp)
'''
The cum_var_exp variable is just the cumulative sum of the explained variance 
and the var_exp is the ratio of the eigenvalue to the total sum of eigenvalues. 
I plotted both of these values below in order to see what percentage of the 
total variance is explained by each principal component. Since the eigenvalues 
are sorted by decreasing order we can see the impact of of adding an additional 
principal component.
'''

## plot accumulative variance
plt.figure()
plt.bar(range(1, len(var_exp)+1), var_exp, alpha=0.75, align='center',
        label='individual explained variance')
plt.step(range(1, len(var_exp)+1), cum_var_exp, where='mid',
         label='cumulative explained variance')
plt.ylim(0, 1.1)
plt.xlabel('Principal components')
plt.ylabel('Explained variance ratio')
plt.legend(loc='best')
plt.title('Correlation among variables')
plt.tight_layout()
plt.savefig('cumsunPCA_'+variables+'.png', res=300)


# -----------------------------------------------------------------------------
#%%
### Select number of components, and process again
from sklearn.decomposition import PCA

# number of PCA components to use
n_components = 6

# PCA with 2 primary components
pca = PCA(n_components=n_components) #PCA with 3 primary components

# fit and transform both PCA models
X_pca = pca.fit_transform(X_std)

variance = pca.explained_variance_ratio_

# -----------------------------------------------------------------------------
#%%
### Using K-means++ to Cluster the Principal Components
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

n = list(range(3,10)) + list(range(10,70,5)) # sequence of clusters to run
distortions = []  # sum of squared error within the each cluster
silhouette = []   # The best value is 1 and the worst value is -1. Values near 0 indicate overlapping clusters.
for i in tqdm(n):
    km = KMeans(n_clusters=i,
               init='k-means++',
               n_init=10,
               max_iter=300,
               random_state=0,
               n_jobs=6)
    km.fit(X_pca)
    cluster_labels = km.fit_predict(X_pca)
    distortions.append(km.inertia_)
    silhouette.append(silhouette_score(inData, cluster_labels))
clust = pd.concat([pd.DataFrame(n), pd.DataFrame(distortions), pd.DataFrame(silhouette)], axis=1)
clust.columns = ['clusters', 'distortions', 'silhouette']
clust.to_csv('clusters'+variables+'.csv')
print(clust['silhouette'])

# plot
fig, ax1 = plt.subplots()
color = 'tab:red'
plt.plot(n, distortions, marker='o', alpha=0.75, color=color)
plt.title('Clustering')
plt.ylabel('Distortions', color=color)
plt.xlabel('Number of clusters')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Silhouette index', color=color)  # we already handled the x-label with ax1
ax2.plot(n, silhouette, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  
plt.savefig('clustering_'+variables+'.png', res=300)

### select number of clusters
#%%
import matplotlib.colors as pltc
n_clust = 5

km = KMeans(n_clusters=n_clust,
           init='k-means++',
           n_init=10,
           max_iter=300,
           tol=1e-04,
           random_state=0)

y_km = km.fit_predict(X_pca)
centers = km.cluster_centers_

# color labels
cmapp5 = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e']
cmapp = pltc.LinearSegmentedColormap.from_list(y_km, cmapp5)

# plot figure
fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(X_pca.T[0], X_pca.T[1], alpha=0.3, c=y_km, cmap=cmapp)
plt.title('PCA ('+ str(np.round(np.sum(variance)*100,2))+' %) - Clustering')
plt.xlabel('PC1 ('+str(np.round(variance[0]*100,2))+' %)')
plt.ylabel('PC2 ('+str(np.round(variance[1]*100,2))+' %)')
plt.scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5)
plt.tight_layout()
plt.savefig('clusteringColor_'+variables+'.png', res=400)

# -----------------------------------------------------------------------------
#%% Predict Clusters to map

## load raster 
with rasterio.open(inRaster) as src:
    arr = src.read() # as array
    meta = src.meta.copy() # save metadata info

# transpose 
arr2 = np.transpose(arr, [1,2,0]) # get to (row, column, band) shape
# reshape to 2D
h, w, numBands = arr2.shape
arr2 = np.reshape(arr2, (w*h, numBands))
# get location of NaN values
isnan = np.isnan(arr2[:,0]) 
# predict classes
arr_std = ss.fit_transform(arr2) # normalization
# Nan Values to zero
arr_pca = pca.fit_transform(np.nan_to_num(arr_std, 0)) # PCA decomposition
arr_cluster = km.fit_predict(arr_pca)
arr_cluster2 = arr_cluster.astype('float32')
arr_cluster2[isnan] = np.nan
# reshaped to 3D
arr_cluster2 = np.reshape(arr_cluster2, (h, w))

# save predict
meta.update(count = 1)
with rasterio.open('cluster_'+variables+'.tif', 'w', **meta) as dst:
       dst.write(arr_cluster2, 1)

# -----------------------------------------------------------------------------
#%% Extract cluster values in 

cluster_stats = zonal_stats('shapefiles/SuisunMarsh_NVCSName_diss.shp', 'cluster_'+variables+'.tif', stats='majority', band=1)
cluster_stats2 = pd.DataFrame(cluster_stats)
cluster_stats2.index = labels
cluster_stats2 = cluster_stats2.sort_values('majority')
cluster_stats2.to_csv('clusterClases_PFTs_masked'+variables+'.csv')
