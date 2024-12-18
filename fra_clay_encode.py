import math

# import geopandas as gpd
import numpy as np
# import pandas as pd
# import pystac_client
# import stackstac
import torch
import yaml
from box import Box
from matplotlib import pyplot as plt
from rasterio.enums import Resampling
from shapely import Point
from sklearn import decomposition, svm
from torchvision.transforms import v2

from src.model import ClayMAEModule

import rasterio
import os
import datetime

device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
torch.set_default_device(device)

def read_tif(file_path):
    with rasterio.open(file_path) as src:
        array = src.read()
    return array

def tif_center(file_path):
    with rasterio.open(file_path) as src:
        bounds = src.bounds

    # Assuming you have the latitude and longitude of the center of the footprint
    lat_center = (bounds.top + bounds.bottom) / 2
    lon_center = (bounds.left + bounds.right) / 2
    return lat_center, lon_center

def rearrange_bands(array):
    """reorders bands to be in the order expected by the model
        R,G,B,NIR,SWIR1,SWIR2
    """
    return array[[2,1,0,3,4,5],:,:]

def get_transform(platform:str="landsat-c2l1",
                  metadata_yml:str="configs/metadata.yaml",
                  bands:list=["red","green","blue","nir08","swir16","swir22"]
                  ):
                  
    metadata = Box(yaml.safe_load(open(metadata_yml)))
    mean = []
    std = []
    waves = []

    # Use the band names to get the correct values in the correct order.
    for band in bands:
        mean.append(metadata[platform].bands.mean[band])
        std.append(metadata[platform].bands.std[band])
        waves.append(metadata[platform].bands.wavelength[band])

    # Prepare the normalization transform function using the mean and std values.
    transform = v2.Compose(
        [
            v2.Normalize(mean=mean, std=std),
        ]
    )

    print(mean)
    print(std)
    print(waves)
    print(transform)
    return transform, waves

def normalize_timestamp(date):
    """Prep datetimes embedding using a normalization function from the model code"""
    week = date.isocalendar().week * 2 * np.pi / 52
    hour = date.hour * 2 * np.pi / 24

    return (math.sin(week), math.cos(week)), (math.sin(hour), math.cos(hour))

def normalize_latlon(lat, lon):
    """# Prep lat/lon embedding using the lat/lon of sample"""
    lat = lat * np.pi / 180
    lon = lon * np.pi / 180

    return (math.sin(lat), math.cos(lat)), (math.sin(lon), math.cos(lon))

def make_datacube(platform:str,
                  pixels:np.ndarray,
                  week_norm:tuple,
                  hour_norm:tuple,
                  lat_norm:tuple,
                  lon_norm:tuple,
                  gsd:int,
                  waves:list,
                  device=None):
    """Create a datacube from the input data
    Can be one or multiple samples
    """
    datacube = {
        "platform": platform,
        "time": torch.tensor(
            np.hstack((week_norm, hour_norm)),
            dtype=torch.float32,
            device=device,
        ),
        "latlon": torch.tensor(
            np.hstack((lat_norm, lon_norm)), dtype=torch.float32, device=device
        ),
        "pixels": pixels.to(device),
        "gsd": torch.tensor(gsd, device=device),
        "waves": torch.tensor(waves, device=device),
    }
    return datacube

def load_model(ckpt_file:str,
               metadata_path:str):
    
    ckpt = ckpt_file
    model = ClayMAEModule.load_from_checkpoint(
        ckpt, 
        metadata_path=metadata_path,#"../../configs/metadata.yaml", 
        shuffle=False, 
        mask_ratio=0
    )
    model.eval()

    model = model.to(device)
    return model

def encode(model,
           datacube):
    with torch.no_grad():
        unmsk_patch, unmsk_idx, msk_idx, msk_matrix = model.model.encoder(datacube)

    # The first embedding is the class token, which is the
    # overall single embedding. We extract that for PCA below.
    embeddings = unmsk_patch[:, 0, :].cpu().numpy()
    return embeddings

###### TEsting
dir = "/home/kyle/code_repos/clay-fra/data/classify-fao-fra/gez16_hex_subsamples_50_50_200_200_tifs/StableForest"
files = os.listdir(dir)[:2]


samples = [read_tif(os.path.join(dir,file)) for file in files]
samples = [rearrange_bands(sample) for sample in samples]

transform, waves = get_transform()
pixels = torch.from_numpy(np.stack(samples,axis=0).astype(np.float32))
pixels = transform(pixels)

# normalize time
dates = []
for file in files:
    if 't1' in file:
        date_str = os.path.basename(file).split('_')[-3]
    else:
        date_str = os.path.basename(file).split('_')[-2]

    date = datetime.datetime(int(date_str), 1, 1,0, 0, 0)
    dates.append(date)


datetimes = [datetime.datetime(int(date_str), 1, 1,0, 0, 0) for date in dates]
times = [normalize_timestamp(dat) for dat in datetimes]
week_norm = [dat[0] for dat in times]
hour_norm = [dat[1] for dat in times]

# normalize location
lat_norm=[]
lon_norm=[]
for file in files:
    lat,lon = tif_center(os.path.join(dir,file))
    norm = normalize_latlon(lat, lon)
    lat_norm.append(norm[0])
    lon_norm.append(norm[1])


datacube = make_datacube(platform="landsat-c2l1",
                         pixels=pixels,
                         week_norm=week_norm,
                         hour_norm=hour_norm,
                         lat_norm=lat_norm,
                         lon_norm=lon_norm,
                         waves=waves,
                         gsd=30,
                         device=device)

# load model and run sample thru it
model = load_model("checkpoints/clay-v1-base.ckpt",
                   "configs/metadata.yaml")

embeddings = encode(model,datacube)
print(embeddings.shape)
print(embeddings)
print('done')
