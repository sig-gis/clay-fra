import stackstac
import pystac_client
from rasterio.enums import Resampling
import geopandas as gpd
import pandas as pd
from shapely import Point

# Point over Monchique Portugal
lat, lon = 37.30939, -8.57207

# Dates of a large forest fire
start = "2018-07-01"
end = "2018-09-01"

# THIS WORKS
STAC_API = "https://earth-search.aws.element84.com/v1"
COLLECTION = "sentinel-2-l2a" 

# BUT, swap to a landsat collection from the same STAC API
# COLLECTION = "landsat-c2-l2"
# we get the following error:
# RuntimeError: Error opening 's3://usgs-landsat/collection02/level-2/standard/oli-tirs/2018/203/034/LC08_L2SP_203034_20180830_20200831_02_T1/LC08_L2SP_203034_20180830_20200831_02_T1_SR_B7.TIF': 
# RasterioIOError('AWS_SECRET_ACCESS_KEY and AWS_NO_SIGN_REQUEST configuration options not defined, and /home/kyle/.aws/credentials not filled')

# use a totally different STAC API and COLLECTION 
STAC_API = "https://landsatlook.usgs.gov/stac-server"
COLLECTION = "landsat-c2l2-sr"
# we get a differnt error
# RuntimeError: Error opening 'https://landsatlook.usgs.gov/data/collection02/level-2/standard/oli-tirs/2018/203/034/LC08_L2SP_203034_20180830_20200831_02_T1/LC08_L2SP_203034_20180830_20200831_02_T1_SR_B7.TIF': 
# RasterioIOError("'/vsicurl/https://landsatlook.usgs.gov/data/collection02/level-2/standard/oli-tirs/2018/203/034/LC08_L2SP_203034_20180830_20200831_02_T1/LC08_L2SP_203034_20180830_20200831_02_T1_SR_B7.TIF' 
# not recognized as being in a supported file format.")

# Search the catalogue
catalog = pystac_client.Client.open(STAC_API)
search = catalog.search(
    collections=[COLLECTION],
    datetime=f"{start}/{end}",
    bbox=(lon - 1e-5, lat - 1e-5, lon + 1e-5, lat + 1e-5),
    max_items=10,
    query={"eo:cloud_cover": {"lt": 80}},
)

all_items = search.get_all_items()

# Reduce to one per date (there might be some duplicates
# based on the location)
items = []
dates = []
for item in all_items:
    if item.datetime.date() not in dates:
        items.append(item)
        dates.append(item.datetime.date())

print(f"Found {len(items)} items")

# Extract coordinate system from first item
epsg = items[0].properties["proj:epsg"]

# Convert point of interest into the image projection
# (assumes all images are in the same projection)
poidf = gpd.GeoDataFrame(
    pd.DataFrame(),
    crs="EPSG:4326",
    geometry=[Point(lon, lat)],
).to_crs(epsg)

coords = poidf.iloc[0].geometry.coords[0]

# Create bounds in projection
size = 32
gsd = 30
bounds = (
    coords[0] - (size * gsd) // 2,
    coords[1] - (size * gsd) // 2,
    coords[0] + (size * gsd) // 2,
    coords[1] + (size * gsd) // 2,
)

# Retrieve the pixel values, for the bounding box in
# the target projection. In this example we use only
# the RGB and NIR bands.
stack = stackstac.stack(
    items,
    bounds=bounds,
    snap_bounds=False,
    epsg=epsg,
    resolution=gsd,
    dtype="float64",
    rescale=False,
    fill_value=0,
    assets=["blue", "green", "red", "nir08","swir16","swir22"],
    resampling=Resampling.nearest,
)
print(stack)

stack = stack.compute()
print(stack.data.shape)