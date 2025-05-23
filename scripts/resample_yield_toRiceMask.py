"""
@author: Enzo
"""
from osgeo import gdal
import os
import glob

def get_all_tif_files(root_folder):
    # Use glob to recursively find all .tif files
    tif_files = glob.glob(os.path.join(root_folder, '**', '*rice_Topt5.tif'), recursive=True)
    return tif_files

# Open the target dataset to get projection, geotransform, and extent
tgt_ds = gdal.Open('Data/Asia mask/asia-rice-extent.tif')
if tgt_ds is None:
    raise FileNotFoundError("Target dataset not found.")

tgt_proj = tgt_ds.GetProjection()
tgt_geotransform = tgt_ds.GetGeoTransform()

origin_x, origin_y = tgt_geotransform[0], tgt_geotransform[3]  # Origin (minimum x, minimum y)
pixel_size_x, pixel_size_y = tgt_geotransform[1], tgt_geotransform[5]  # Pixel sizes (delta x, delta y)
print(f"Pixel size: rx={pixel_size_x}, ry={pixel_size_y}")

# Calculate bounding box of the target dataset
minx = origin_x
maxx = origin_x + tgt_ds.RasterXSize * pixel_size_x
miny = origin_y + tgt_ds.RasterYSize * pixel_size_y
maxy = origin_y

bbox = (minx, miny, maxx, maxy)
print(f"BBox: {bbox}")

folder_path = "Data/Model results/Crop_yield"
files = os.listdir(folder_path)
# Filter GeoTIFF files
tif_files = get_all_tif_files(folder_path)
# for file in tif_files:
#     print(file)

for f in tif_files:
    src_filename = f
    print(src_filename)
    src_ds = gdal.Open(src_filename)

    output_file = src_filename.replace("rice", "rice_rs")
    # output_file = output_file.replace("_yield/", "_yield/rs/")
    # Perform the resampling
    resample_options = gdal.WarpOptions(
        dstSRS=tgt_proj,
        xRes=pixel_size_x,
        yRes=pixel_size_y,
        resampleAlg='near',
        outputBounds=bbox  # Ensuring the output bounds match the target dataset's extent
    )
    gdal.Warp(output_file, src_ds, options=resample_options)
    # Cleanup
    src_ds = None        
print("Resampling completed.")
    
# Cleanup   
tgt_ds = None