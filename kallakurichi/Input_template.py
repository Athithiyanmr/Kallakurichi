# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import geopandas as gpd
import pandas as pd
import numpy as np

from osgeo import ogr, gdal
from osgeo import gdal_array
from osgeo import gdalconst
import geopandas as gpd
from geopandas import overlay
import pandas as pd
from matplotlib import pyplot
import matplotlib.pyplot as plt

from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.font_manager as fm
from pygc import great_distance
import matplotlib.patches as mpatches
import rasterio
import rasterio.plot
from rasterio.plot import show
import rasterio.mask

import fiona
import matplotlib.ticker as ticker
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors

import folium
from folium import plugins,Vega
from folium import raster_layers

import branca.colormap as cm
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
from folium.plugins import MarkerCluster,FloatImage, Draw, MeasureControl
import geopandas as gpd
import os
from matplotlib import colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import openpyxl
from openpyxl import Workbook,load_workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.chart import LineChart, Reference
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

from matplotlib.patches import ConnectionPatch

from matplotlib.transforms import ScaledTranslation

from matplotlib.patches import FancyArrow
from shapely.geometry import shape


# %%
# !pip freeze > requirements.txt

# %% [markdown]
# ## kallakurichi

# %%
def get_rooted(stem):
    return "D:/Lila_Kallakuruchi/" + stem
def read_df_UT(stem):
    return gpd.read_file(get_rooted(stem)).to_crs(epsg = 4326)
def get_in_workdir(stem):
        return get_rooted('workdir/' + stem)
def get_in_output(stem):
    return get_rooted("output/"+ stem)
def read_rastergdal_UT(stem):
    return gdal.Open(get_rooted(stem))  
def read_raster_UT(stem):
    return rasterio.open(get_rooted(stem))


# %%
def plot_cities(fig, ax):
    shp_cities = read_df_UT("workdir/extra_inputs/shp_major_cities/Kallakurichi_major_towns-point.shp")
    shp_cities['coords'] = shp_cities['geometry'].apply(lambda x: x.representative_point().coords[:])
    shp_cities['coords'] = [coords[0] for coords in shp_cities['coords']]
    shp_cities["coords"].tolist()
    shp_cities[['lat', 'lon', 'zero']] = gpd.GeoDataFrame(shp_cities['coords'].tolist(), index=shp_cities.index)

    x = shp_cities["lat"]
    y = shp_cities["lon"]    
    
    labels =shp_cities["Name"]

    for i in range(0,len(shp_cities)):
        plt.text(x[i]+0.008,y[i]+0.008,labels[i],fontsize=4,color = '#c3c3c3', ha = 'center',zorder=10)

    shp_cities.plot(ax=ax, markersize=3, color='grey',zorder=11)


# %%
def plot_common_features(fig, ax):
    plt.rcParams['font.family'] = 'Helvetica';
    plt.grid(color="grey",linestyle = '--', linewidth = 0.1)


    scalebar = AnchoredSizeBar(ax.transData,
                           0.1835161781605187,
                           '20 km',  
                           loc='lower left',
                           frameon=False,
                           size_vertical=0.0125,
                           color='lightgrey',
                           bbox_to_anchor=(0.1, 0.0125),
                           bbox_transform=ax.transAxes,
                           sep = 1.5,
                           fontproperties={'family': 'Helvetica', 'size': 5.5})

    ax.add_artist(scalebar)
    
    ax.tick_params(axis='x', colors='grey', labelsize=3, labeltop = 'True', labelrotation = 270)
    ax.tick_params(axis='y', colors='grey', labelsize=3, labelright = 'True')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlim(78.50, 79.50)
    ax.set_ylim(11.40, 12.30)

    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none') 
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')

    x, y, arrow_length = 0.5, 0.99, 0.1
    ax.annotate('N',color= "lightgrey", xy=(x, y), xytext=(x, y-arrow_length), arrowprops=dict(facecolor='none',edgecolor="lightgrey", width=4, headwidth=12), ha='center', va='center', fontsize=10, xycoords=ax.transAxes)


# %%
def area_acres(df):
    crs_utm = 32644 
    df = df.to_crs(crs_utm)
    df["area_acres"] = (((df.geometry.area)/10**6)*247.105)
    a = df["area_acres"].max()
    def area_class(df):
        if df["area_acres"] < 5:
            return "D"
        elif 5 <= df["area_acres"] <= 20:
            return "A"
        elif 20 < df["area_acres"] <= 100:
            return "B"
        else:
            return "C"
    df["area_class"] =df.apply(area_class, axis=1)
    df = df.to_crs(4326)
    return df


# %%
def add_area_and_class(gdf):
    if not isinstance(gdf.geometry, gpd.GeoSeries):
        raise ValueError("Input dataframe does not contain a valid geometry column.")

    gdf_utm = gdf.to_crs(32644)
    gdf_utm["area_acres"] = (((gdf_utm.geometry.area) / 10 ** 6) * 247.105)

    conditions = [
        gdf_utm["area_acres"] < 5,
        (gdf_utm["area_acres"] >= 5) & (gdf_utm["area_acres"] <= 20),
        (gdf_utm["area_acres"] > 20) & (gdf_utm["area_acres"] <= 100),
        gdf_utm["area_acres"] > 100,
    ]
    choices = ["D", "A", "B", "C"]
    gdf_utm["area_class"] = np.select(conditions, choices)

    gdf_wgs84 = gdf_utm.to_crs(4326)
    return gdf_wgs84


# %%
def add_area_and_class_forest(gdf):
    if not isinstance(gdf.geometry, gpd.GeoSeries):
        raise ValueError("Input dataframe does not contain a valid geometry column.")

    gdf_utm = gdf.to_crs(32644)
    gdf_utm["area_acres"] = (((gdf_utm.geometry.area) / 10 ** 6) * 247.105)

    conditions = [
        gdf_utm["area_acres"] < 2.47,
        (gdf_utm["area_acres"] >= 2.47) & (gdf_utm["area_acres"] <= 20),
        (gdf_utm["area_acres"] > 20) & (gdf_utm["area_acres"] <= 100),
        gdf_utm["area_acres"] > 100,
    ]
    choices = ["D", "A", "B", "C"]
    gdf_utm["area_class"] = np.select(conditions, choices)

    gdf_wgs84 = gdf_utm.to_crs(4326)
    return gdf_wgs84


# %%


def add_area_and_class_agri(gdf):
    if not isinstance(gdf.geometry, gpd.GeoSeries):
        raise ValueError("Input dataframe does not contain a valid geometry column.")

    gdf_utm = gdf.to_crs(32644)
    gdf_utm["area_acres"] = (((gdf_utm.geometry.area) / 10 ** 6) * 247.105)

    conditions = [
        gdf_utm["area_acres"] < 1,
        (gdf_utm["area_acres"] >= 1) & (gdf_utm["area_acres"] <= 5),
        (gdf_utm["area_acres"] > 5) & (gdf_utm["area_acres"] <= 10),
        gdf_utm["area_acres"] > 10,
    ]
    choices = ["D", "A", "B", "C"]
    gdf_utm["area_class"] = np.select(conditions, choices)

    gdf_wgs84 = gdf_utm.to_crs(4326)
    return gdf_wgs84


# %%
def add_area_and_class_housing(gdf):

    if not isinstance(gdf.geometry, gpd.GeoSeries):
        raise ValueError("Input dataframe does not contain a valid geometry column.")

    gdf_utm = gdf.to_crs(32644)
    gdf_utm["area_acres"] = (((gdf_utm.geometry.area) / 10 ** 6) * 247.105)

    conditions = [
        gdf_utm["area_acres"] < 0.25,
        (gdf_utm["area_acres"] >= 0.25) & (gdf_utm["area_acres"] <= 2),
        (gdf_utm["area_acres"] > 2) & (gdf_utm["area_acres"] <= 5),
        gdf_utm["area_acres"] > 5,
    ]
    choices = ["D", "A", "B", "C"]
    gdf_utm["area_class"] = np.select(conditions, choices)

    gdf_wgs84 = gdf_utm.to_crs(4326)
    return gdf_wgs84


# %%
def add_area_and_class_industry(gdf):

    if not isinstance(gdf.geometry, gpd.GeoSeries):
        raise ValueError("Input dataframe does not contain a valid geometry column.")

    gdf_utm = gdf.to_crs(32644)
    gdf_utm["area_acres"] = (((gdf_utm.geometry.area) / 10 ** 6) * 247.105)

    conditions = [
        gdf_utm["area_acres"] < 0.3,
        (gdf_utm["area_acres"] >= 0.3) & (gdf_utm["area_acres"] <= 20),
        (gdf_utm["area_acres"] > 20) & (gdf_utm["area_acres"] <= 100),
        gdf_utm["area_acres"] > 100,
    ]
    choices = ["D", "A", "B", "C"]
    gdf_utm["area_class"] = np.select(conditions, choices)

    gdf_wgs84 = gdf_utm.to_crs(4326)
    return gdf_wgs84


# %% [markdown]
# ### define func for creating Area_hect and area_class

# %%
def area_hect(df):
    crs_utm = 32644 
    df = df.to_crs(crs_utm)
    df["area_hect"] = ((df.geometry.area)/10**4)
    a = df["area_hect"].max()
    def area_class(df):
        if df["area_acres"] < 5:
            return "D"
        elif 5 <= df["area_acres"] < 20:
            return "A"
        elif 20 <= df["area_acres"] < 100:
            return "B"
        else:
            return "C"
    df["area_class"] =df.apply(area_class, axis=1)
    df = df.to_crs(4326)
    return df


# %% [markdown]
# ### define func for the difference and intersection between shape files

# %%
def intersection(df,df1,dist):
    df = gpd.overlay(dist,df,how ="intersection")
    df1 = gpd.overlay(dist,df1,how ="intersection")
    df2 = gpd.overlay(df,df1,how ="intersection")
    return df2   


# %% [markdown]
# ### def fun for Calculating overlap area

# %%
def find_overlap_area(df,tag,fdf2):
    crs_utm = 32644    
    df = df.to_crs(crs_utm)
    df1 = pd.DataFrame(columns = ['olap%'+tag,'olaparea'+tag,'count'+tag])
    df1['olap%'+tag]=df1['olap%'+tag].astype('object')
    df1['olaparea'+tag]=df1['olaparea'+tag].astype('object')
    df1['count'+tag] = df1['count'+tag].astype('object')
    fdf2=fdf2.to_crs(crs_utm)
    #set spatial index for data for faster processing
    sindex = fdf2.sindex
    for i in range(len(df)):
        geometry = df.iloc[i]['geometry']
        fids = list(sindex.intersection(geometry.bounds))
        if fids:
            olaparea = ((fdf2.iloc[fids]['geometry'].intersection(geometry)).area).sum()
            count = (fdf2.iloc[fids]['geometry'].intersection(geometry)).count()
            olap_perc = olaparea*100/geometry.area
            olaparea = (olaparea/10**6)*247.1               
        else:
            olaparea = 0
            olap_perc = 0
            count = 0
        df1.at[i,'olap%'+tag] =  olap_perc      
        df1.at[i,'olaparea'+tag] = olaparea
        df1.at[i,'count'+tag] = count
    return pd.concat([df,df1], axis= 1)


# %%
def find_overlap(df,tag,fdf2):
    crs_utm = 32644
    df = df.to_crs(crs_utm)
    df1 = pd.DataFrame(columns = ['op%'+tag,'opar'+tag,'cnt'+tag])
    df1['op%'+tag]=df1['op%'+tag].astype('object')
    df1['opar'+tag]=df1['opar'+tag].astype('object')
    df1['cnt'+tag] = df1['cnt'+ tag].astype('object')
    fdf2=fdf2.to_crs(crs_utm)

    sindex = fdf2.sindex

    for i in range(len(df)):
        geometry = df.iloc[i]['geometry']
        intersection = fdf2[fdf2.intersects(geometry)]
        count = len(intersection)
        if count:
            olaparea = intersection['geometry'].intersection(geometry).area.sum()
            olap_perc = olaparea*100/geometry.area
            olaparea = (olaparea/10**6)*247.1               
        else:
            olaparea = 0
            olap_perc = 0
        df1.at[i,'op%'+tag] =  olap_perc 
        df1.at[i,'cnt'+tag] = count
        df1.at[i,'opar'+tag] = olaparea
    return pd.concat([df,df1], axis= 1)


# %%
def find_overlap_count(df,tag,fdf2):
    crs_utm = 32644
    df = df.to_crs(crs_utm)
    df1 = pd.DataFrame(columns = ['op%'+tag,'opar'+tag,'cnt'+tag])
    df1['op%'+tag]=df1['op%'+tag].astype('object')
    df1['opar'+tag]=df1['opar'+tag].astype('object')
    df1['cnt'+tag] = df1['cnt'+ tag].astype('object')
    fdf2=fdf2.to_crs(crs_utm)
    
    #set spatial index for data for faster processing
    sindex = fdf2.sindex

    for i in range(len(df)):
        geometry = df.iloc[i]['geometry']
        intersection = fdf2[fdf2.intersects(geometry)]
        count = len(intersection)
        df1.at[i,'cnt'+tag] = count
    return pd.concat([df,df1], axis= 1)


# %%
def competing_acres(df,df1):
    df = df["geometry"]
    df1 =df1["geometry"]
    df = gpd.GeoDataFrame(df)
    df1 = gpd.GeoDataFrame(df1)
    df2 = gpd.overlay(df,df1,how ="intersection",keep_geom_type=True)
    df2 = gpd.GeoDataFrame(df2)
    df2 = df2.to_crs(4326)
    return df2


# %%
def calculate_slope(df, input_raster):
    df.geometry = df.geometry.buffer(0)
    outputdf = pd.DataFrame()

    for i in range(len(df)):
        input_shp =  get_rooted('workdir/temp.shp')
        selection = df['geometry'][i:i+1]
        selection.to_file(input_shp)

        output_raster = get_rooted('workdir/temp.tif')

        ds = gdal.Warp(output_raster,
                      input_raster,
                      format = 'GTiff',
                      cutlineDSName = input_shp,
                      cropToCutline=True,
                      )
        ds = None

        raster = gdal.Open(output_raster, gdal.GA_ReadOnly)
        rasterarr = raster.ReadAsArray()
        rasterarr = rasterarr[rasterarr!= -9999]

        if (np.size(rasterarr)==0):
            outputdf.at[i, 'min']=np.nan
            outputdf.at[i , 'max']=np.nan
            outputdf.at[i , 'mean']=np.nan
            outputdf.at[i , '25percentile']=np.nan
            outputdf.at[i , '75percentile']=np.nan

        else:
            outputdf.at[i, 'min']=rasterarr.min()
            outputdf.at[i , 'max']=rasterarr.max()
            outputdf.at[i , 'mean']=rasterarr.mean()
            outputdf.at[i , '25percentile']=np.percentile(rasterarr,25)
            outputdf.at[i , '75percentile']=np.percentile(rasterarr,75)
    df = pd.concat([df,outputdf], axis= 1)
    return df
