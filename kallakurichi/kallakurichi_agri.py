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

# %% [markdown]
# # Agriculture Report

# %%
import Input_template
from Input_template import *

# %%
## boundary files
district_boundary = "workdir/extra_inputs/shp_district/shp_district.shp"
village = "workdir/extra_inputs/shp_settlement/shp_settlement.shp"
taluk = "workdir/extra_inputs/shp_taluk/shp_taluk.shp"
## module files
theoritical = "workdir/agri/shp_agri_theo/shp_agri_theo.shp"
technical ="workdir/agri/shp_agri_tech/shp_agri_tech.shp"
high = "workdir/agri/shp_agri_high/shp_agri_high.shp"
med  = "workdir/agri/shp_agri_med/shp_agri_med.shp"
## other module high files
water_high = "workdir/water/shp_water_high/shp_water_high.shp"
forest_high = "workdir/forest/shp_forest_high/shp_forest_high.shp"
housing_high = "workdir/housing/shp_housing_high/shp_housing_high.shp"
solar_high = "workdir/solar/shp_solar_high/shp_solar_high.shp"
industry_high = "workdir/industry/shp_industry_high/shp_industry_high.shp"
Unused_land = "workdir/extra_inputs/Unused/Unused.shp"
## landcover shape files
Unused_land = "workdir/extra_inputs/Unused/Unused.shp" 
Sparseveg = "workdir/extra_inputs/shp_sparseveg/shp_sparseveg.shp"
Cropland = "workdir/extra_inputs/shp_cropland/shp_cropland.shp"
Forest = "workdir/extra_inputs/shp_forest/shp_forest.shp"
waterbodies = "workdir/extra_inputs/shp_water_bodies/shp_water_bodies.shp"
urban ="workdir/extra_inputs/shp_builtup/shp_builtup.shp"
access ="workdir/extra_inputs/shp_access/shp_access.shp"
## dist by size
S1 = "workdir/agri/S1/S1.shp"
S2 = "workdir/agri/S2/S2.shp"
S3 ="workdir/agri/S3/S3.shp"
## extra inputs
slope =get_rooted("workdir/raster/slope.tif")
shp_soil_erosion = read_df_UT("workdir/agri/soil_erosion/soil_erosion.shp")

# %%
## boundary files
shp_district = read_df_UT(district_boundary)
shp_village = read_df_UT(village)
shp_taluk = read_df_UT(taluk)
## module files
shp_agri_theo = read_df_UT(theoritical)
shp_agri_tech = read_df_UT(technical)
shp_agri_high =read_df_UT(high)
shp_agri_med = read_df_UT(med)
## other module high files
shp_water_high = read_df_UT(water_high)
shp_forest_high = read_df_UT(forest_high)
shp_housing_high = read_df_UT(housing_high)
shp_solar_high = read_df_UT(solar_high)
shp_indus_high =read_df_UT(industry_high)
## landcover shape files
shp_unused = read_df_UT(Unused_land)
## dist by size
S1 = read_df_UT(S1)
S2 = read_df_UT(S2)
S3 = read_df_UT(S3)

# %%
shp_unused = add_area_and_class_agri(shp_unused)


# %% [markdown]
# ## Report Tables

# %% [markdown]
# ### Stats

# %%
def summarize_data(df, name):
    summary = df.groupby("area_class")["area_acres"].agg(["sum", "count"]).reset_index()
    summary.columns = ["area_class", f"{name}_sum", f"{name}_count"]
    return summary

data_frames = [shp_unused, shp_agri_theo, shp_agri_tech, shp_agri_high, shp_agri_med]
data_names = ['unused', 'theo', 'tech', 'high', 'med']

summary_dict = {}

for df, name in zip(data_frames, data_names):
    summary_dict[name] = summarize_data(df, name)

merged_summary = summary_dict['unused']
for name in data_names[1:]:
    merged_summary = pd.merge(merged_summary, summary_dict[name], on="area_class", how="outer")

merged_summary = merged_summary.fillna(0)
merged_summary

# %% [markdown]
# ### Competing use (need to improve)

# %%
shp_agri_water= find_overlap(shp_agri_tech,"water",shp_water_high)
shp_agri_forest = find_overlap(shp_agri_tech,"forest",shp_forest_high)
shp_agri_housing =find_overlap(shp_agri_tech,"hsing",shp_housing_high)
shp_agri_indus =find_overlap(shp_agri_tech,"indus",shp_indus_high)
shp_agri_solar =find_overlap(shp_agri_tech,"solar",shp_solar_high)
combined = pd.concat([shp_indus_high,shp_forest_high,shp_water_high,shp_solar_high,shp_housing_high])
combined.reset_index(inplace =True,drop =True)
inter = gpd.overlay(shp_agri_tech,combined,how ="intersection",keep_geom_type=True)
inter = inter.dissolve()
inter = add_area_and_class(inter)
data = {
    'agri_forest': shp_agri_forest.oparforest.sum(),
    'agri_water': shp_agri_water.oparwater.sum(),
    'agri-housing':shp_agri_housing.oparhsing.sum(),
    'agri-indus': shp_agri_indus.oparindus.sum(),
    'agri-solar': shp_agri_solar.oparsolar.sum(),
    'Competing use':inter.area_acres.sum(),
}
competing_use = pd.DataFrame(data, index=[0])

# %%
competing_use

# %% [markdown]
# ### Distribution by type

# %%
S_list = [S1, S2, S3]

overlap_data = {}

land_use_categories = ['forest', 'water','hsing', 'indus','solar']
high_dataframes = [shp_forest_high, shp_water_high,shp_housing_high, shp_indus_high,shp_solar_high]


for idx, S in enumerate(S_list):
    S_op_dict = {}
    for category, high_df in zip(land_use_categories, high_dataframes):
        S_op_dict[category] = find_overlap(S, category, high_df)
    overlap_data[f"S{idx+1}"] = S_op_dict

data = []
for idx, S in enumerate(S_list):
    row_data = [
        overlap_data[f"S{idx+1}"]["forest"].oparforest.sum(),
        overlap_data[f"S{idx+1}"]["water"].oparwater.sum(),
        overlap_data[f"S{idx+1}"]["hsing"].oparhsing.sum(),
        overlap_data[f"S{idx+1}"]["indus"].oparindus.sum(),
        overlap_data[f"S{idx+1}"]["solar"].oparsolar.sum()
    ]
    data.append(row_data)

Dist_by_type = pd.DataFrame(data, columns=['Forest', 'Water', 'housing', 'industry','solar'], index=['S1', 'S2', 'S3'])

# %%
Dist_by_type

# %%
largest_plot = shp_agri_tech.area_acres.max()

# %%
largest_plot

# %% [markdown]
# ## Soil erosion

# %% [markdown]
# ### the correct results are in the bottom 

# %%
shp_soil_erosion.erosion_ca.value_counts()

# %%
shp_moderately = shp_soil_erosion[shp_soil_erosion["erosion_ca"] =="moderately"]
shp_slightly = shp_soil_erosion[shp_soil_erosion["erosion_ca"] =="slightly"]
shp_severely = shp_soil_erosion[shp_soil_erosion["erosion_ca"] =="severely"]

# %%
shp_soil_erosion = read_df_UT("workdir/agri/soil_erosion/soil_erosion.shp")

shp_agri_tech = read_df_UT(technical)

shp_soil_erosion = shp_soil_erosion[["erosion_ca","geometry"]]

shp_agri_tech = shp_agri_tech[["area_acres","geometry"]]

categories = shp_soil_erosion
lands = shp_agri_tech
lands["Category"] =""
# Loop through the joined data and assign categories
for index, row in lands.iterrows():
    for cat_index, cat_row in categories.iterrows():
        if row['geometry'].within(cat_row['geometry']):
            lands.at[index, 'Category'] = cat_row['erosion_ca']
            break
    else:
        lands.at[index, 'Category'] = 'Mixed'

lands.Category.value_counts().sum()

df = lands.groupby("Category")["area_acres"].agg(["sum", "count"])

df

# %%
template = load_workbook(get_rooted('workdir/agri/agri_temp.xlsx'))
ws = template['Sheet1']

# Define the ranges where to insert the dataframes
ranges = ['A2:K5','A10:F10','B13:F15','B20:C23']

# Insert the dataframes into the workbook
for x, r in zip([merged_summary,competing_use,Dist_by_type,df], ranges):
    start_col, start_row, end_col, end_row = openpyxl.utils.cell.range_boundaries(r)
    for i, row in enumerate(x.values):
        for j, value in enumerate(row):
            ws.cell(row=start_row + i, column=start_col + j, value=value)
ws.cell(row=17, column=2, value=largest_plot)           
template.save(get_in_output('agri/agri_temp_filled.xlsx'))

# %% [markdown]
# ## Top15

# %%
shp_agri_high_sorted = shp_agri_high.sort_values(by=["area_acres"],ascending = False)
shp_agri_top15  = shp_agri_high_sorted[:15]
shp_agri_top15.reset_index(inplace =True)

# %%
df = shp_agri_top15
slope =get_rooted("workdir/raster/slope.tif")
top15_agri = calculate_slope(df, slope)

# %%
top15_agri['coords'] = top15_agri['geometry'].apply(lambda x: x.representative_point().coords[:])
top15_agri['coords'] = [coords[0] for coords in top15_agri['coords']]
top15_agri["coords"].tolist()
top15_agri[['lat', 'lon']] = gpd.GeoDataFrame(top15_agri['coords'].tolist(), index=top15_agri.index) 

# %%
overlap = find_overlap(top15_agri,"solar",shp_solar_high)
overlap = find_overlap(overlap,"forest",shp_forest_high)
overlap = find_overlap(overlap,"water",shp_water_high)
overlap = find_overlap(overlap,"hsing",shp_housing_high)
overlap = find_overlap(overlap,"indus",shp_indus_high)

# %%
overlap.columns

# %%
overlap.info()

# %%
overlap.iloc[:,34:49] = overlap.iloc[:,34:49].astype(float)

# %%
overlap_top15 = overlap[["lat","lon","area_acres","wtmindist","wtdist","erosion_cl","min","max",'op%solar', 'oparsolar',
       'cntsolar', 'op%forest', 'oparforest', 'cntforest', 'op%water',
       'oparwater', 'cntwater', 'op%hsing', 'oparhsing', 'cnthsing', 'op%indus',
       'oparindus', 'cntindus','geometry']]

# %%
overlap_top15.to_file(get_in_output("agri/Top15_final"))

# %%
overlap_top15_excel = overlap[["lat","lon","area_acres","wtmindist","wtdist","erosion_cl","min","max",'op%solar', 'oparsolar',
       'cntsolar', 'op%forest', 'oparforest', 'cntforest', 'op%water',
       'oparwater', 'cntwater', 'op%hsing', 'oparhsing', 'cnthsing', 'op%indus',
       'oparindus', 'cntindus']]

# %%
overlap_top15_excel.to_excel(get_in_output("agri/Top15_final.xlsx"))

# %% [markdown]
# # Visuals

# %% [markdown]
# ## Technical Suitability - Technical, theoretical, and no potential lands

# %%
fig8, ax8 = plt.subplots(figsize=(5, 5))

plot_common_features(fig8, ax8)
plot_cities(fig8, ax8)


shp_district.plot(figsize=(5,5),color="none", ax=ax8, linewidth = 0.5, zorder=5)

shp_unused.plot(color="#424242",ax =ax8, label='No Potential')
shp_agri_theo.plot(color="#a09a22",ax =ax8, label='Theoretical Potential')
shp_agri_tech.plot(color="#f3ff8b",ax =ax8, label='Technical Potential')


No_P = mpatches.Patch(color='#424242', label='No potential')
Theo_P = mpatches.Patch(color='#a09a22', label='Theoretical potential')
Tech_P = mpatches.Patch(color='#f3ff8b', label='Technical potential')
    
plt.legend(handles = [No_P, Theo_P, Tech_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)


print(plt.rcParams['font.family'])


plt.savefig(get_in_output("images/agri/Technical_suitability.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ## Distribution by type

# %%
fig9, ax9 = plt.subplots(figsize=(5, 5))

plot_common_features(fig9, ax9)
plot_cities(fig9, ax9)


shp_district.plot(figsize=(5,5),color="none", ax=ax9, linewidth = 0.5)

S1.plot(color="#393928",ax =ax9, label='>1 to 5 acres')
S2.plot(color="#789821",ax =ax9, label='>5 to 10 acres')
S3.plot(color="#d9f389",ax =ax9, label='>10 acres')


S1 = mpatches.Patch(color='#393928', label='>1 to 5 acres')
S2 = mpatches.Patch(color='#789821', label='>5 to 10 acres')
S3 = mpatches.Patch(color='#d9f389', label='>10 acres')
    
plt.legend(handles = [S1, S2, S3], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

print(plt.rcParams['font.family']);


plt.savefig(get_in_output("images/agri/Distribution by size.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ## High Potential

# %%
fig10, ax10 = plt.subplots(figsize=(5, 5))

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)

shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

shp_agri_tech.plot(color="#4f4d0e",ax =ax10, label='Low Potential')
shp_agri_med.plot(color="#a6c12d",ax =ax10, label='Medium Potential')
shp_agri_high.plot(color="#dbe7a4",ax =ax10, label='High Potential')


Low_P = mpatches.Patch(color='#4f4d0e', label='Low potential')
Med_P = mpatches.Patch(color='#a6c12d', label='Medium potential')
High_P = mpatches.Patch(color='#dbe7a4', label='High potential')
    
plt.legend(handles = [Low_P, Med_P, High_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

print(plt.rcParams['font.family'])


plt.savefig(get_in_output("images/agri/High Potential_H_M_L.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# # Settlement analysis

# %%
popdf = pd.DataFrame()
for j in range(len(shp_village)):
    input_shp =  get_rooted('workdir/temp_agri.shp')

    selection = shp_village.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted("workdir/raster/population_10_lon_70_general-v1.5.tif")

        output_raster = get_rooted('workdir/temp_agri.tif')
        ds = gdal.Warp(output_raster,
                      input_raster,
                      format = 'GTiff',
                      cutlineDSName = input_shp,
                      cropToCutline=True,
                      )
        ds = None

        raster = gdal.Open(output_raster, gdal.GA_ReadOnly)
        rasterarr = raster.ReadAsArray()
        #Set -9999 as no data values
        rasterarr = np.where(rasterarr==-9999, np.nan,rasterarr)
        #remove nodata values
        rasterarr = rasterarr[~np.isnan(rasterarr)]


    if (np.size(rasterarr)==0):
        popdf.at[j, "totpop"]=0


    else:    

        pop_sum = rasterarr.sum()

        popdf.at[j, "totpop"]=pop_sum

shp_village_final = pd.concat([shp_village, popdf], axis = 1)

# %%
shp_village_final.columns

# %%
shp_village_final.totpop.sum()

# %%
shp_village_final['coords'] = shp_village_final['geometry'].apply(lambda x: x.representative_point().coords[:])
shp_village_final['coords'] = [coords[0] for coords in shp_village_final['coords']]
shp_village_final["coords"].tolist()
shp_village_final[['lat', 'lon']] = gpd.GeoDataFrame(shp_village_final['coords'].tolist(), index=shp_village_final.index) 

# %%
shp_village_final = shp_village_final.to_crs(32644)
shp_village_final["TGA(acres)"] = ((shp_village_final.geometry.area)/10**6)*247.105
shp_village_final = shp_village_final.to_crs(4326)

# %%
shp_village_tech = find_overlap(shp_village_final,"tech",shp_agri_tech)

# %%
# shp_village_tech =  find_overlap(shp_village_tech,"unused",shp_unused)

# %%
shp_village_tech = shp_village_tech.to_crs(4326)

# %%
shp_village_tech.info()

# %%
shp_village_tech = shp_village_tech.iloc[:, [14,15,16,17] + list(range(44, 55))]

# %%
shp_village_tech.info()

# %%
shp_village_tech.iloc[:,12:15] = shp_village_tech.iloc[:,12:15].astype(float)

# %%
shp_village_tech= shp_village_tech[['p_name', 'd_name', 'b_name', 'p_name_rd', 'area_acres', 'area_class'
       , 'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech',
       'opartech', 'cnttech','geometry']]

# %%
shp_village_tech.to_file(get_in_output("agri/settlement"))

# %%
shp_village_tech.columns

# %%
shp_village = shp_village_tech[['p_name', 'd_name', 'b_name', 'p_name_rd', 'area_acres', 'area_class',
       'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech', 'opartech', 'cnttech']]

# %%
shp_village.to_excel(get_in_output("agri/230517_settlement.xlsx"))

# %%
shp_village.sort_values("p_name",inplace =True)

# %%
shp_village.reset_index(drop =True,inplace =True)

# %%
shp_village_tech.drop("geometry",axis =1,inplace =True)

# %%
shp_village.to_excel(get_in_output("agri/230517_settlement_ordered.xlsx"))

# %% [markdown]
# # Settlement visuals

# %%
settlement = read_df_UT("output/agri/settlement/settlement.shp")

# %%
S3 = LinearSegmentedColormap.from_list('testCmap1', colors=["#fbfada", "#def399", "#c4d85e", "#7a8737", "#4f541e"], N=256)

# %%

fig, ax = plt.subplots(figsize=(5, 5))

plot_common_features(fig, ax)  # assuming this function is defined elsewhere

shp_district.plot(figsize=(5,5), color="none", ax=ax, linewidth=0.5, zorder=5)
settlement.plot(column='op%tech', cmap=S3, ax=ax)
a = settlement["op%tech"].min()
b = settlement["op%tech"].max()
sm = plt.cm.ScalarMappable(cmap=S3)
cbaxes = fig.add_axes([0.7, 0.18, 0.2, 0.02]) 
cbar = plt.colorbar(sm, orientation = 'horizontal', cax=cbaxes, shrink = 0.2)
cbar.mappable.set_clim(vmin = a, vmax = b)
cbar.ax.tick_params(labelsize=3, color = 'grey')
cbar.outline.set_visible(False)
cbar.ax.set_title('Technical potential \n (% of TGA)', fontsize=5)

print(a,b)

plt.savefig(get_in_output("images/agri/Technical potential area(%).jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ### Taluk analysis

# %%
from shapely.geometry import shape

# %%
shp_taluk = "workdir/extra_inputs/shp_taluk/shp_taluk.shp"

# %%
shp_taluk = read_df_UT(shp_taluk)

# %%
shp_taluk = gpd.overlay(shp_district,shp_taluk,how="intersection")

# %%
# shp_taluk.to_file(get_in_workdir("extra_inputs/shp_taluk"))

# %%
popdf = pd.DataFrame()
for j in range(len(shp_taluk)):
    input_shp =  get_rooted('workdir/temp_t_h.shp')

    selection = shp_taluk.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted("workdir/raster/population_10_lon_70_general-v1.5.tif")

        output_raster = get_rooted('workdir/temp_t_h.tif')
        ds = gdal.Warp(output_raster,
                      input_raster,
                      format = 'GTiff',
                      cutlineDSName = input_shp,
                      cropToCutline=True,
                      )
        ds = None

        raster = gdal.Open(output_raster, gdal.GA_ReadOnly)
        rasterarr = raster.ReadAsArray()
        #Set -9999 as no data values
        rasterarr = np.where(rasterarr==-9999, np.nan,rasterarr)
        #remove nodata values
        rasterarr = rasterarr[~np.isnan(rasterarr)]


    if (np.size(rasterarr)==0):
        popdf.at[j, "totpop"]=0


    else:    

        pop_sum = rasterarr.sum()

        popdf.at[j, "totpop"]=pop_sum

shp_taluk_final = pd.concat([shp_taluk, popdf], axis = 1)

# %%
shp_taluk_final

# %%
shp_taluk_final['coords'] = shp_taluk_final['geometry'].apply(lambda x: x.representative_point().coords[:])
shp_taluk_final['coords'] = [coords[0] for coords in shp_taluk_final['coords']]
shp_taluk_final["coords"].tolist()
shp_taluk_final[['lat', 'lon']] = gpd.GeoDataFrame(shp_taluk_final['coords'].tolist(), index=shp_taluk_final.index) 

# %%
shp_taluk_final = shp_taluk_final.to_crs(32644)
shp_taluk_final["TGA(acres)"] = ((shp_taluk_final.geometry.area)/10**6)*247.105
shp_taluk_final = shp_taluk_final.to_crs(4326)

# %%
shp_taluk_final = find_overlap(shp_taluk_final,"tech",shp_agri_tech)

# %%
# shp_taluk_final['geometry'] = shp_taluk_final['geometry'].apply(lambda x: shape(x).buffer(0).buffer(0.0000000000001))
# shp_taluk_finalshp_taluk_final = gpd.GeoDataFrame(shp_taluk_final, geometry='geometry')
# shp_taluk_final = shp_taluk_final.loc[shp_taluk_final.is_valid]

# %%
# shp_taluk_final =  find_overlap(shp_taluk_final,"unused",shp_unused)

# %%
shp_taluk_final = shp_taluk_final.to_crs(4326)

# %%
shp_taluk_final= shp_taluk_final[['Taluk_name', 'area_acres', 'area_class'
       , 'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech',
       'opartech', 'cnttech','geometry']]

# %%
shp_taluk_final.totpop.sum()

# %%
shp_taluk_final.to_file(get_in_output("agri/shp_taluk"))

# %%
shp_taluk_final.drop("geometry",axis =1,inplace =True)

# %%
shp_taluk_final.to_excel(get_in_output("agri/taluk_analysis.xlsx"))

# %%
shp_taluk_final

# %%
# Create a new column that maps categories to specific colors
shp_soil_erosion['color'] = np.where(shp_soil_erosion['erosion_ca'] == 'slightly', '#fff4aa',
                       np.where(shp_soil_erosion['erosion_ca'] == 'moderately', '#ffd65a',
                                np.where(shp_soil_erosion['erosion_ca'] == 'severely', '#cb8b54', 'white')))


# %%

fig, ax = plt.subplots(figsize=(5, 5))

plot_common_features(fig, ax) 


shp_agri_tech.plot(figsize=(5,5),facecolor ="none",edgecolor="#4f4d0e" ,ax=ax, linewidth=0.1, zorder=7)
shp_agri_med.plot(figsize=(5,5),facecolor ="none",edgecolor="#a6c12d" ,ax=ax, linewidth=0.1, zorder=8)
shp_agri_high.plot(figsize=(5,5),facecolor ="none",edgecolor="#dbe7a4" ,ax=ax, linewidth=0.1, zorder=9)

shp_district.plot(figsize=(5,5), color="none", ax=ax, linewidth=0.5, zorder=5)
shp_soil_erosion.plot(color=shp_soil_erosion['color'], ax=ax)
Severely = mpatches.Patch(color="#cb8b54", label = 'Severely eroded')
Moderately = mpatches.Patch(color="#ffd65a", label = 'Moderately eroded')
Slightly= mpatches.Patch(color="#fff4aa", label = 'Slightly eroded')

Low_P = mpatches.Patch(color='#4f4d0e', label='Low potential')
Med_P = mpatches.Patch(color='#a6c12d', label='Medium potential')
High_P = mpatches.Patch(color='#dbe7a4', label='High potential')

plt.legend(handles = [Severely, Moderately, Slightly,Low_P,Med_P,High_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.27), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)
plt.savefig(get_in_output("images/agri/soil_erosion.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ### Soil erosion stats

# %%
shp_soil_erosion = read_df_UT("workdir/agri/soil_erosion/soil_erosion.shp")

shp_agri_tech = read_df_UT(technical)

shp_soil_erosion = shp_soil_erosion[["erosion_ca","geometry"]]

shp_agri_tech = shp_agri_tech[["area_acres","geometry"]]

shp_soil_erosion["erosion_ca"].fillna("Others",inplace =True)

categories = shp_soil_erosion
lands = shp_agri_tech

lands["Category"] =""
# Loop through the joined data and assign categories
for index, row in lands.iterrows():
    for cat_index, cat_row in categories.iterrows():
        if row['geometry'].within(cat_row['geometry']):
            lands.at[index, 'Category'] = cat_row['erosion_ca']
            break
    else:
        lands.at[index, 'Category'] = 'Mixed'



df = lands.groupby("Category")["area_acres"].agg(["sum", "count"])

df

# %%
df

# %%
shp_soil_erosion = read_df_UT("workdir/agri/soil_erosion/soil_erosion.shp")

shp_agri_tech = read_df_UT(technical)

shp_soil_erosion = shp_soil_erosion[["erosion_ca","geometry"]]
shp_agri_tech = shp_agri_tech[["area_acres","geometry"]]

shp_soil_erosion["erosion_ca"].fillna("other", inplace=True)

# %%
shp_soil_erosion["erosion_ca"].value_counts()

# %%
shp_slightly = shp_soil_erosion[shp_soil_erosion["erosion_ca"] =="slightly"]
shp_slightly.reset_index(drop =True,inplace =True)
shp_moderately = shp_soil_erosion[shp_soil_erosion["erosion_ca"] =="moderately"]
shp_moderately.reset_index(drop =True,inplace =True)
shp_severely = shp_soil_erosion[shp_soil_erosion["erosion_ca"] =="severely"]
shp_severely.reset_index(drop =True,inplace =True)
shp_other = shp_soil_erosion[shp_soil_erosion["erosion_ca"] =="other"]
shp_other.reset_index(drop =True,inplace =True)

# %%
soil_eros_res = find_overlap(shp_agri_tech,"sl",shp_slightly)

# %%
soil_eros_res = find_overlap(soil_eros_res,"mo",shp_moderately)

# %%
soil_eros_res = find_overlap(soil_eros_res,"sv",shp_severely)

# %%
soil_eros_res = find_overlap(soil_eros_res,"ot",shp_other)


# %%
def get_category(soil_eros_res):
    if soil_eros_res['op%ot'] > 0:
        return 'Others'
    elif soil_eros_res['op%sl'] > 99:
        return 'slightly'
    elif soil_eros_res['op%mo'] > 99:
        return 'moderately'
    elif soil_eros_res['op%sv'] > 99:
        return 'severely'
    else:
        return "Mixed"
    
# apply the function to create the new column
soil_eros_res['category'] = soil_eros_res.apply(get_category, axis=1)

# %%
soil_eros_res.groupby("category")["area_acres"].agg(["sum", "count"])

# %%
soil_eros_res.category.value_counts()

# %%
shp_slightly = soil_eros_res[soil_eros_res["category"] =="slightly"]
shp_moderately = soil_eros_res[soil_eros_res["category"] =="moderately"]
shp_severely = soil_eros_res[soil_eros_res["category"] =="severely"]
shp_others = soil_eros_res[soil_eros_res["category"] =="Others"]
shp_mixed = soil_eros_res[soil_eros_res["category"] =="Mixed"]

# %%
shp_slightly = shp_slightly.to_crs(4326)
shp_moderately = shp_moderately.to_crs(4326)
shp_severely = shp_severely.to_crs(4326)
shp_others = shp_others.to_crs(4326)
shp_mixed = shp_mixed.to_crs(4326)

# %%
shp_slightly.reset_index(inplace=True,drop =True)
shp_moderately.reset_index(inplace =True,drop =True)
shp_severely.reset_index(inplace =True,drop =True)
shp_others.reset_index(inplace =True,drop =True)
shp_mixed.reset_index(inplace =True,drop =True)

# %%
# shp_slightly.to_file(get_in_output("agri/shp_slightly"))
# shp_moderately.to_file(get_in_output("agri/shp_moderately"))
# shp_severely.to_file(get_in_output("agri/shp_severely"))
# shp_others.to_file(get_in_output("agri/shp_others"))
# shp_mixed.to_file(get_in_output("agri/shp_mixed"))

# %%

# %%

# %%
nan_count = lands['Category'].isna().sum()

# Calculate area of land with NaN category values
nan_area = lands.loc[lands['Category'].isna(), 'area_acres'].sum()

# %%
print(nan_area,nan_count)

# %%
shp_agri_tech.area_acres.sum(),shp_agri_tech.shape

# %%
df.sum()

# %%
79869.70836315978-79653.110594

# %%
10727-10652

# %%
# lands.to_file(get_in_output("agri/soil_erosion_final"))
