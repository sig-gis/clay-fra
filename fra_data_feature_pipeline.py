import ee
import google.auth
import shapely.geometry
from pyproj import Transformer

"""
Step 1 do all the FC processing and export a GeoJSON containing PLOTID, Change1 label, Change2 label, GEZ
this way we do the part of the sampling requireing EE requests only once, do the subsequent subsampling in geopandas befroe sending computePixels requests..
"""
PROJECT = "pc530-fao-fra-rss"  # change to your cloud project name

## INIT WITH HIGH VOLUME ENDPOINT
credentials, _ = google.auth.default()
ee.Initialize(
    credentials,
    project=PROJECT,
    opt_url="https://earthengine-highvolume.googleapis.com",
)
    
# USE HEXAGONS TO MAKE PATCH BOUNDS ######################################################
gez_zones = range(11,17)
n = 100
for i in gez_zones:
    
    hex_zone = (ee.FeatureCollection("projects/pc530-fao-fra-rss/assets/reference/hexWCenPropertiesTropics")
    .filterMetadata("GEZ", "equals", i)
    )

    hex_zone_centroid = hex_zone.map(lambda h: ee.Feature(h.geometry().centroid()).copyProperties(h))
    gdf = (ee.data.computeFeatures({'expression': hex_zone_centroid,
                                          'fileFormat': 'GEOPANDAS_GEODATAFRAME'})
                                                .set_crs("EPSG:4326")
                                                # .to_crs('EPSG:3857')
                                                )
    gdf_subset = gdf[['PL_PLOTID','CHANGE0010','CHANGE1018','GEZ','geometry']]
    
    
    gdf_subset.to_file(f"data/classify-fao-fra/features/gez{i}_centroids.geojson", driver='GeoJSON')