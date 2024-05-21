import sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
import geopandas as gpd
import scanpy as sc
import anndata
import polars as pl
import pickle

from tifffile import imread, imwrite
from csbdeep.utils import normalize
from stardist.models import StarDist2D
from shapely.geometry import Polygon, Point
from scipy import sparse
from matplotlib.colors import ListedColormap
import pyarrow

from csbdeep.io import save_tiff_imagej_compatible
from stardist import export_imagej_rois


def main(args):

    # Read full image
    img = imread(args.image)

    # crop and save if coordinates exist
    if args.rmax is not None:
        cropped_img = img[args.rmin : args.rmax, args.cmin : args.cmax]
        save_tiff_imagej_compatible(args.out+'/img_Stardist.tif', cropped_img, axes='YXC')

    # Load the pretrained model
    model = StarDist2D.from_pretrained('2D_versatile_he')

    # Percentile normalization of the image
    min_percentile = 5
    max_percentile = 95
    img = normalize(img, min_percentile, max_percentile)

    # Run segmentation for the full image (we will select and crop later)
    labels, polys = model.predict_instances_big(img, axes='YXC', block_size=4096, prob_thresh=0.01,nms_thresh=0.001, min_overlap=128, context=128, normalizer=None, n_tiles=(4,4,1))

     # Save the results to save time if required later.
    f1 = open(args.out+'/labels_FullSection.pckl', 'wb')
    pickle.dump(labels, f1)
    f1.close()

    f2 = open(args.out+'/polys_FullSection.pckl', 'wb')
    pickle.dump(polys, f2)
    f2.close()

    if args.rmax is not None:
        polys2 = polys
        index = []
        for nuclei in range(len(polys2['coord'])):
            if np.all(polys2['coord'][nuclei][0] > args.rmin) and np.all(polys2['coord'][nuclei][0] < args.rmax) and np.all(polys2['coord'][nuclei][1] > args.cmin) and np.all(polys2['coord'][nuclei][1] < args.cmax) :
                polys2['coord'][nuclei][0] = polys2['coord'][nuclei][0] - args.rmin
                polys2['coord'][nuclei][1] = polys2['coord'][nuclei][1] - args.cmin
                index.append(nuclei)
        
        # Save masks for Subset
        export_imagej_rois(args.out+'/img_rois_Stardist_Subset.zip', polys2['coord'][index])


    # Creating a list to store Polygon geometries
    geometries = []

    # Iterating through each nuclei in the 'polys' DataFrame 
    for nuclei in range(len(polys['coord'])):
        coords = [(y, x) for x, y in zip(polys['coord'][nuclei][0], polys['coord'][nuclei][1])]

        # Creating a Polygon geometry from the coordinates
        geometries.append(Polygon(coords))

    # Creating a GeoDataFrame using the Polygon geometries
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf["id"] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]

    # Load Visium HD data
    raw_h5_file = args.xenadir + "/filtered_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(raw_h5_file)

    # Load the Spatial Coordinates
    tissue_position_file = args.xenadir + "/spatial/tissue_positions.parquet"
    df_tissue_positions = pd.read_parquet(tissue_position_file)
    # Set the index of the dataframe to the barcodes
    df_tissue_positions = df_tissue_positions.set_index("barcode")

    # Create an index in the dataframe to check joins
    df_tissue_positions["index"] = df_tissue_positions.index

    # Adding the tissue positions to the meta data
    adata.obs = pd.merge(adata.obs, df_tissue_positions, left_index=True, right_index=True)

    # Create a GeoDataFrame from the DataFrame of coordinates
    geometry = [
        Point(xy)
        for xy in zip(
            df_tissue_positions["pxl_col_in_fullres"],
            df_tissue_positions["pxl_row_in_fullres"],
        )
    ]

    gdf_coordinates = gpd.GeoDataFrame(df_tissue_positions, geometry=geometry)

    # Identify BCs within a nuclei mask and overlaps as well
    result_spatial_join = gpd.sjoin(gdf_coordinates, gdf, how="left", predicate="within")
    result_spatial_join["is_within_polygon"] = ~result_spatial_join["index_right"].isna()
    barcodes_in_overlaping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=["index"])]["index"])
    result_spatial_join["is_overlap"] = result_spatial_join["index"].isin(barcodes_in_overlaping_polygons)

    Result = result_spatial_join.loc[result_spatial_join['is_within_polygon']]

    Result[['index','id','is_within_polygon','is_overlap']].to_csv(args.out+"/Nuclei_Barcode_Map.csv")



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-i','--image', type=str, help='Path to .btf image')
    parser.add_argument('-r1','--rmin', type=int, help='row min')
    parser.add_argument('-r2','--rmax', type=int, help='row max')
    parser.add_argument('-c1','--cmin', type=int, help='column min')
    parser.add_argument('-c2','--cmax', type=int, help='column max')
    parser.add_argument('-x','--xenadir', type=str, help='Path to spaceranger outs')
    parser.add_argument('-o',"--out", type=str,help="Path to output directory")
    args = parser.parse_args()
    main(args)