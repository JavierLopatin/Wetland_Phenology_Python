#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 17:19:43 2020

@author: javier
"""

from shapely.geometry import Point
import numpy as np
import pandas as pd
import rtree
from multiprocessing import Manager, Process

def improved_function(shp, number, maxiter, cassNames):
    # Create a spatial index to quickly check containment
    index = rtree.index.Index()
    for i, polygon in enumerate(shp.geometry):
        index.insert(i, polygon.bounds)

    # Create a shared queue to store the results
    manager = Manager()
    outlist = manager.list()

    def worker(j, outlist):
        # get the j-th polygon
        polygon = shp.geometry[j]
        # get the bounds of the polygon
        min_x, min_y, max_x, max_y = polygon.bounds

        # Generate a batch of points
        batch_size = 100
        points = np.array([[random.uniform(min_x, max_x), random.uniform(min_y, max_y)] for _ in range(batch_size)])
        inside_points = []

        for point in points:
            point = Point(point)
            # check if point is inside the polygon
            if polygon.contains(point):
                # Use the spatial index to quickly find the nearest points
                nearest = index.nearest((point.x, point.y, point.x, point.y), objects=True)
                for n in nearest:
                    # check the distance of the point with the nearest points
                    if point.distance(shp.geometry[n.object]) < number:
                        inside_points.append(point)
                        break

        # create a dataframe for the points inside the polygon
        df = pd.DataFrame({"geometry": inside_points})
        df = df.assign(class=shp.class[j])
        outlist.append(df)

    # Use all CPU cores minus one for parallel processing
    num_cores = multiprocessing.cpu_count() - 1
    with multiprocessing.Pool(num_cores) as pool:
        pool.starmap(worker, [(j, outlist) for j in range(len(shp))])

    return pd.concat(outlist)


if __name__ == "__main__":
    
    # Read a shapefile
    shp = gpd.read_file("Suisun_PFT_diss.shp")

    # Define the number of points to generate
    number = 1000

    # Define the maximum number of iterations
    maxiter = 100

    # Define the column name of the class
    cassNames = "OBJECTID"

    # Generate points inside the polygons
    df = improved_function(shp, number, maxiter, cassNames)

    df.to_file('samplepoints.shp')
