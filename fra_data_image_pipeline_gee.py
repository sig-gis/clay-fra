import os
import geopandas as gpd
from shapely import Point
import shapely
import ee
import sys
import argparse
import rasterio
from rasterio.transform import from_bounds
import tqdm

sys.path.append("../")


ee.Initialize(project='pc530-fao-fra-rss',
              opt_url='https://earthengine-highvolume.googleapis.com')

classes ={
        'Stable Non Forest': 0,
        'Stable Forest': 1,
        'Forest Loss': 2,
        'Forest Gain': 3,
    }
# read feature geoJSON (these are per GEZ)
def transform(gdf):

    
    list_of_features = []
    for i,row in gdf.iterrows():
        f_dict = {
            'PLOTID': row.PL_PLOTID,
            'geometry': row.geometry,
            'top-left': Point(
                [list(row.geometry.bounds)[0],
                 list(row.geometry.bounds)[3]]),

            'example1': 
            {'label': row.CHANGE0010.replace(' ',''),
            'label_int': classes[row.CHANGE0010],
            't1start':'2000-01-01',
            't1end':'2000-12-31', 
            't2start': '2010-01-01',
            't2end': '2010-12-31',
            },

            'example2': 
            {'label': row.CHANGE1018.replace(' ',''),
            'label_int': classes[row.CHANGE1018],
            't1start':'2010-01-01',
            't1end':'2010-12-31', 
            't2start': '2018-01-01',
            't2end': '2018-12-31',
            }
            
            
            }
        
        list_of_features.append(f_dict)
    
    return list_of_features

def get_landsat_composite(region:shapely.Polygon,
                          start:str,
                          end:str):
    # Define the region of interest as a bounding box.
    region = ee.Geometry.Polygon(list(region.exterior.coords))
    band_mapper = {
        'l5':{
            'bands': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
            'band_names': ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']},
        'l7':{
            'bands': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
            'band_names': ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']},
        'l8':{
            'bands': ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
            'band_names': ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']},
        'l9':{
            'bands': ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
            'band_names': ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
        }
    }
    
    # Applies scaling factors.
    def apply_scale_factorsl5l7(image):
       optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
       thermal_bands = image.select('ST_B6').multiply(0.00341802).add(149.0)
       return image.addBands(optical_bands, None, True).addBands(
           thermal_bands, None, True
           )
    def apply_scale_factorsl8l9(image):
       optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
       thermal_bands = image.select('ST_B10').multiply(0.00341802).add(149.0)
       return image.addBands(optical_bands, None, True).addBands(
           thermal_bands, None, True
           )
    
    def qal5l7(image):
        """Custom QA masking method for Landsat9 surface reflectance dataset"""
        qa_band = image.select("QA_PIXEL")
        qa_flag = int('111111',2)
        sat_mask = image.select('QA_RADSAT').eq(0);
        mask = qa_band.bitwiseAnd(qa_flag).eq(0).And(sat_mask)
        return apply_scale_factorsl5l7(image).updateMask(mask)
    
    def qal8l9(image):
        """Custom QA masking method for Landsat9 surface reflectance dataset"""
        qa_band = image.select("QA_PIXEL")
        qa_flag = int('111111',2)
        sat_mask = image.select('QA_RADSAT').eq(0);
        mask = qa_band.bitwiseAnd(qa_flag).eq(0).And(sat_mask)
        return apply_scale_factorsl8l9(image).updateMask(mask)
    
    # Create a Landsat image collection for the specified date range.
    l5 = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
          .filterBounds(region)
          .filterDate(start, end)
          .map(qal5l7)
          .select(band_mapper['l5']['bands'], band_mapper['l5']['band_names'])
          )
    l7 = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
            .filterBounds(region)
            .filterDate(start, end)
            .map(qal5l7)
            .select(band_mapper['l7']['bands'], band_mapper['l7']['band_names'])
            )
    l8 = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
            .filterBounds(region)
            .filterDate(start, end)
            .map(qal8l9)
            .select(band_mapper['l8']['bands'], band_mapper['l8']['band_names'])
            )
    l9 = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
            .filterBounds(region)
            .filterDate(start, end)
            .map(qal8l9)
            .select(band_mapper['l9']['bands'], band_mapper['l9']['band_names'])
            )
    
    collection = l5.merge(l7).merge(l8).merge(l9)
    image = collection.median()
    return image

def download_tif_from_geom_centr(image:ee.Image,
                 geom:shapely.Geometry,
                 size:int,
                 gsd:int,
                 out_dir:str,
                 fileprefix:str):
    
    def get_innermost_H_W(array, height=32, width=32):
        # Get the shape of the input array
        rows, cols = array.shape
        
        # Calculate the starting indices
        start_row = (rows - height) // 2
        start_col = (cols - width) // 2
        
        # Slice the array to get the innermost 32x32 values
        inner_array = array[start_row:start_row + height, start_col:start_col + width]
        
        return inner_array
    
    if os.path.exists(f'{out_dir}/{fileprefix}.tif'):
        print(f'{fileprefix}.tif already exists')
        return
    hex = ee.Geometry.Polygon(list(geom.exterior.coords))
    # bounds = hex.bounds() # gives us little bit less than a 32x32px img 
    centroid_buffer = hex.centroid().buffer(gsd*size/2)
    image = image.reproject("EPSG:4326", scale=30).clipToBoundsAndScale(centroid_buffer, scale=30)

    # request image as numpy array, do some reshaping to 32,32
    request = {'expression': image,
               'fileFormat': 'NUMPY_NDARRAY'}
    data = ee.data.computePixels(request)
    # Assuming data is a structured array with multiple bands
    bands = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
    data_innermost = {band: get_innermost_H_W(data[band],size,size) for band in bands}
    
    # Get the CRS and transform info from the ee.Image
    crs = image.projection().crs().getInfo()
    coords = centroid_buffer.bounds().coordinates().getInfo()[0]
    west = coords[0][0]
    east = coords[2][0]
    south = coords[1][1]
    north = coords[3][1]
    transform = from_bounds(west, south, east, north, 32, 32)
    
    # Write the (32,32) array to a GeoTIFF
    out_path = f'{out_dir}/{fileprefix}.tif'
    with rasterio.open(
        out_path,
        'w',
        driver='GTiff',
        height=size,
        width=size,
        count=len(bands),
        dtype=data_innermost[bands[0]].dtype,
        crs=crs,
        transform=transform
    ) as dst:
        for i, band in enumerate(bands, start=1):
            dst.write(data_innermost[band], i)
    
    print(f"Downloaded {out_path}")
    return 

#### MAIN ####
def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-file', 
                        type=str, 
                        required=True, 
                        help='Path to the geojson file')
    parser.add_argument('-sample-count', 
                        type=int, 
                        required=False, 
                        help='Number of samples to process')

    args = parser.parse_args()

    features = gpd.read_file(args.file)
    if args.sample_count:
        features = features.sample(n=args.sample_count) # seed needs to be set for reproducibility
    print(f"Processing {len(features)} features")

    basename = str(os.path.basename(args.file).split('.')[0])
    outdir = os.path.join(os.getcwd(), 'data', 'classify-fao-fra', f'{basename}_tifs')
    subfolders = [i.replace(' ', '') for i in classes.keys()]
    subdirs = [os.path.join(outdir, i) for i in subfolders]
    for dir in subdirs:
        os.makedirs(dir, exist_ok=True)

    size = 32
    gsd = 30

    for i in tqdm.tqdm(transform(features)):
        # each sample from GeoJSON will result in 4 downloaded .tifs 
        # (2 transition samples x 2 .tifs per transition)
        plotid = i['PLOTID']
        geometry = i['geometry']

        for ex in ['example1', 'example2']:
            example = i[ex]
            ex_folder = f"{example['label']}"
            fileprefix = os.path.join(ex_folder, f"PL{plotid}_{example['label']}_{example['t1start'][:4]}_{example['t2start'][:4]}")

            image1 = get_landsat_composite(region=geometry, start=example['t1start'], end=example['t1end'])
            download_tif_from_geom_centr(image=image1, geom=geometry, gsd=gsd, size=size, out_dir=outdir, fileprefix=f"{fileprefix}_t1")
            
            image2 = get_landsat_composite(region=geometry, start=example['t2start'], end=example['t2end'])
            download_tif_from_geom_centr(image=image2, geom=geometry, gsd=gsd, size=size, out_dir=outdir, fileprefix=f"{fileprefix}_t2")
        
if __name__ == '__main__':
    main()
