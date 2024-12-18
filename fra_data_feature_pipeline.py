import ee
import google.auth
import shapely.geometry
from pyproj import Transformer
import argparse

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
    
def main():
    parser = argparse.ArgumentParser(description="Export FRA features from GEE to GeoJSON.")
    parser.add_argument('-gez', 
                        nargs='+', 
                        type=int, 
                        required=True, 
                        help='GEZ zones to process')
    parser.add_argument('-samples', 
                        nargs='+', 
                        type=int, 
                        required=False, 
                        help='Number of samples to limit per class in order of Stable Non Forest, Stable Forest, Forest Loss, Forest Gain')

    args = parser.parse_args()
    gez_zones = args.gez
    n = args.samples

    for i in gez_zones:
        hex_zone = (ee.FeatureCollection("projects/pc530-fao-fra-rss/assets/reference/hexWCenPropertiesTropics")
                    .filterMetadata("GEZ", "equals", i)
                    )
        if n:
            if not isinstance(n,list) or len(n) != 4:
                raise ValueError("Please provide a list of samples for each of the 4 classes (list length 4)")
            
            hex_zone_snf = (hex_zone.filter(ee.Filter.Or(
                ee.Filter.eq("CHANGE0010", "Stable Non Forest"),
                ee.Filter.eq("CHANGE1018", "Stable Non Forest")))
                .limit(n[0]))
            hex_zone_sf = (hex_zone.filter(ee.Filter.Or(
                ee.Filter.eq("CHANGE0010", "Stable Forest"),
                ee.Filter.eq("CHANGE1018", "Stable Forest")))
                .limit(n[1]))
            hex_zone_fl = (hex_zone.filter(ee.Filter.Or(
                ee.Filter.eq("CHANGE0010", "Forest Loss"),
                ee.Filter.eq("CHANGE1018", "Forest Loss")))
                .limit(n[2]))
            hex_zone_fg = (hex_zone.filter(ee.Filter.Or(
                ee.Filter.eq("CHANGE0010", "Forest Gain"),
                ee.Filter.eq("CHANGE1018", "Forest Gain")))
                .limit(n[3]))
            hex_zone = ee.FeatureCollection(hex_zone_snf.merge(hex_zone_sf).merge(hex_zone_fl).merge(hex_zone_fg))
        
        gdf = (ee.data.computeFeatures({'expression': hex_zone,
                                        'fileFormat': 'GEOPANDAS_GEODATAFRAME'})
                .set_crs("EPSG:4326"))

        gdf_subset = gdf[['PL_PLOTID', 'CHANGE0010', 'CHANGE1018', 'GEZ', 'geometry']]
        gdf_subset.to_file(f"data/classify-fao-fra/features/gez{i}_hex_subsamples_{n[0]}_{n[1]}_{n[2]}_{n[3]}.geojson", driver='GeoJSON')

if __name__ == "__main__":
    main()