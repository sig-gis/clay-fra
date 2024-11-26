import os
import geopandas as gpd
from shapely import Point
import shapely
import ee
import sys

sys.path.append("../")


ee.Initialize(project='pc530-fao-fra-rss',
              opt_url='https://earthengine-highvolume.googleapis.com')


# read feature geoJSON (these are per GEZ)
def transform(gdf):

    classes ={
        'Stable Non Forest': 0,
        'Stable Forest': 1,
        'Forest Loss': 2,
        'Forest Gain': 3,
    }
    list_of_features = []
    for i,row in gdf.iterrows():
        
        f_dict = {
            'PLOTID': row.PL_PLOTID,
            'geometry': row.geometry,

            'example1': 
            {'label': row.CHANGE0010,
            'label_int': classes[row.CHANGE0010],
            't1start':'2000-01-01',
            't1end':'2000-12-31', 
            't2start': '2010-01-01',
            't2end': '2010-12-31',
            },

            'example2': 
            {'label': row.CHANGE1018,
            'label_int': classes[row.CHANGE1018],
            't1start':'2010-01-01',
            't1end':'2010-12-31', 
            't2start': '2018-01-01',
            't2end': '2018-12-31',
            }
            
            
            }
        
        list_of_features.append(f_dict)
    
    return list_of_features

def get_landsat_composite(region:shapely.Point,
                          start:str,
                          end:str):
    # Define the region of interest as a bounding box.
    region = ee.Geometry.Point([region.x, region.y])
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
    def apply_scale_factors(image):
       optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
       thermal_bands = image.select('ST_B6').multiply(0.00341802).add(149.0)
       return image.addBands(optical_bands, None, True).addBands(
           thermal_bands, None, True
           )
    def qa(image):
        """Custom QA masking method for Landsat9 surface reflectance dataset"""
        qa_band = image.select("QA_PIXEL")
        qa_flag = int('111111',2)
        sat_mask = image.select('QA_RADSAT').eq(0);
        mask = qa_band.bitwiseAnd(qa_flag).eq(0).And(sat_mask)
        return apply_scale_factors(image).updateMask(mask)
    
    # Create a Landsat image collection for the specified date range.
    l5 = (ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
          .filterBounds(region)
          .filterDate(start, end)
          .map(qa)
          .select(band_mapper['l5']['bands'], band_mapper['l5']['band_names'])
          )
    l7 = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
            .filterBounds(region)
            .filterDate(start, end)
            .map(qa)
            .select(band_mapper['l7']['bands'], band_mapper['l7']['band_names'])
            )
    l8 = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
            .filterBounds(region)
            .filterDate(start, end)
            .map(qa)
            .select(band_mapper['l8']['bands'], band_mapper['l8']['band_names'])
            )
    l9 = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
            .filterBounds(region)
            .filterDate(start, end)
            .map(qa)
            .select(band_mapper['l9']['bands'], band_mapper['l9']['band_names'])
            )
    
    collection = l5.merge(l7).merge(l8).merge(l9)
    image = collection.median()
    return image#.reproject("EPSG:3857", scale=30).clipToBoundsAndScale(region, scale=30)
    

def download_tif(image:ee.Image,
                 region:shapely.Point,
                 size:int,
                 gsd:int,
                 out_dir:str,
                 fileprefix:str):
    proj = ee.Projection('EPSG:4326').atScale(gsd).getInfo()

    # region is Point(x,y)
    coord_x = region.x
    coord_y = region.y
    
    # Get scales out of the transform.
    scale_x = proj['transform'][0]
    scale_y = -proj['transform'][4]
    
    # Retrieve the pixel values for the bounding box.
    request = {'expression': image,
               'fileFormat': 'GEO_TIFF',#'NUMPY_NDARRAY',
               # Make a projection to discover the scale in degrees.
               'grid': {
                   'dimensions': {
                       'width': size,
                       'height': size
                       },
                       'affineTransform': {
                           'scaleX': scale_x,
                           'shearX': 0,
                           'translateX': coord_x,
                           'shearY': 0,
                           'scaleY': scale_y,
                           'translateY': coord_y
                           },
                           'crsCode': proj['crs'],
                        }
                }
    data = ee.data.computePixels(request)   
    with open(f'{out_dir}/{fileprefix}.tif', 'wb') as f:
            f.write(data)
    return 

#### MAIN ####
file = "data/classify-fao-fra/features/gez15_centroids.geojson"
features = gpd.read_file(file)
basename = str(os.path.basename(file).split('.')[0])
outdir=os.path.join(os.getcwd(),'data','classify-fao-fra',f'{basename}_tifs')
os.makedirs(outdir, exist_ok=True)

size = 32
gsd = 30


for i in transform(features)[:10]:
    plotid = i['PLOTID']
    geometry = i['geometry']

    # TODO: loop through [t1start/t1end, t2start/t2end] tuples, 
    # we have to download 2 images per example, 2 examples per plot, so 4 tifs per plot
    # make sure each folder and .tif is named descriptively enough 
    # so that when we feed t1.tif and t2.tif in to model to get pairwise embeddings, 
    # we know the label is 'Stable Forest', 'Forest Loss', etc.
    example1 = i['example1']
    image1 = get_landsat_composite(region=geometry,
                                   start=example1['t1start'],
                                   end=example1['t1end'])
    
    tif1 = download_tif(image=image1,
                        region=geometry,
                        size=size,
                        gsd=gsd,
                        out_dir=outdir,
                        fileprefix=f"PL{plotid}_example1_{example1['label']}")
    
    example2 = i['example2']
    image2 = get_landsat_composite(region=geometry,
                                   start=example2['t1start'],
                                   end=example2['t1end'])
    tif2 = download_tif(image=image2,
                        region=geometry,
                        size=size,
                        gsd=gsd,
                        out_dir=outdir,
                        fileprefix=f"PL{plotid}_example2_{example2['label']}")
    # break
print()
