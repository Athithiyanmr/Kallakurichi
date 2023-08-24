#!/usr/bin/env python
# coding: utf-8
# %% [markdown]
# # Kallakurichi Forest report

# %%
import Input_template
from Input_template import *

# %% [markdown]
# ## Defining variables with input paths

# %%
## boundary files
district_boundary = "workdir/extra_inputs/shp_district/shp_district.shp"
village = "workdir/extra_inputs/shp_settlement/shp_settlement.shp"
taluk = "workdir/extra_inputs/shp_taluk/shp_taluk.shp"
## module files
theoritical = "workdir/forest/shp_forest_theo/shp_forest_theo.shp"
technical ="workdir/forest/shp_forest_tech/shp_forest_tech.shp"
high = "workdir/forest/shp_forest_high/shp_forest_high.shp"
med  = "workdir/forest/shp_forest_med/shp_forest_med.shp"
## other module high files
water_high = "workdir/water/shp_water_high/shp_water_high.shp"
solar_high = "workdir/solar/shp_solar_high/shp_solar_high.shp"
agri_high = "workdir/agri/shp_agri_high/shp_agri_high.shp"
housing_high = "workdir/housing/shp_housing_high/shp_housing_high.shp"
industry_high = "workdir/industry/shp_industry_high/shp_industry_high.shp"
## landcover shape files
Unused_land = "workdir/extra_inputs/Unused/Unused.shp"
Sparseveg = "workdir/extra_inputs/shp_sparseveg/shp_sparseveg.shp"
Cropland = "workdir/extra_inputs/shp_cropland/shp_cropland.shp"
Forest = "workdir/extra_inputs/shp_forest/shp_forest.shp"
Waterbodies = "workdir/extra_inputs/shp_water_bodies/shp_water_bodies.shp"
builtup ="workdir/extra_inputs/shp_builtup/shp_builtup.shp"
## dist by size
S1 = "workdir/forest/S1/S1.shp"
S2 = "workdir/forest/S2/S2.shp"
S3 ="workdir/forest/S3/S3.shp"
## extra inputs
slope = get_rooted("workdir/raster/slope.tif")
roads_primary  =  "workdir/extra_inputs/shp_roads_primary/shp_roads_primary.shp"
roads_secondary = "workdir/extra_inputs/shp_roads_secondary/shp_roads_secondary.shp"
roads_tertiary  =  "workdir/extra_inputs/shp_roads_tertiary/shp_roads_tertiary.shp"
railways ="workdir/extra_inputs/shp_railways/shp_railways.shp"
powerlines= "workdir/extra_inputs/shp_powerlines/shp_powerlines.shp"
substation = "workdir/extra_inputs/shp_substations/shp_substations.shp"
population_raster ="workdir/raster/population_10_lon_70_general-v1.5.tif"


# %%
## boundary files
shp_district = read_df_UT(district_boundary)
shp_village = read_df_UT(village)
shp_taluk = read_df_UT(taluk)
## module files
shp_forest_theo = read_df_UT(theoritical)
shp_forest_tech = read_df_UT(technical)
shp_forest_high =read_df_UT(high)
shp_forest_med = read_df_UT(med)
## other module high files
shp_water_high = read_df_UT(water_high)
shp_solar_high = read_df_UT(solar_high)
shp_agri_high = read_df_UT(agri_high)
shp_housing_high = read_df_UT(housing_high)
shp_industry_high = read_df_UT(industry_high)
## landcover shape files
shp_unused = read_df_UT(Unused_land)
shp_sparseveg = read_df_UT(Sparseveg)
shp_cropland = read_df_UT(Cropland)
shp_Forest = read_df_UT(Forest)
shp_waterbodies = read_df_UT(Waterbodies)
shp_builtup = read_df_UT(builtup)
## dist by size
S1 = read_df_UT(S1)
S2 = read_df_UT(S2)
S3 = read_df_UT(S3)
## extra inputs
shp_roads_primary  =  read_df_UT(roads_primary)
shp_roads_secondary = read_df_UT(roads_secondary)
shp_railways = read_df_UT(railways)
shp_powerlines= read_df_UT(powerlines)
shp_substation = read_df_UT(substation)


# %% [markdown]
# ### Adding area_acres and area_class

# %%
shp_forest_theo = add_area_and_class_forest(shp_forest_theo)
shp_forest_tech = add_area_and_class_forest(shp_forest_tech)
shp_forest_high =add_area_and_class_forest(shp_forest_high)
shp_forest_med = add_area_and_class_forest(shp_forest_med)
shp_unused = add_area_and_class_forest(shp_unused)


# %% [markdown]
# ## Report Tables

# %%
def summarize_data(df, name):
    summary = df.groupby("area_class")["area_acres"].agg(["sum", "count"]).reset_index()
    summary.columns = ["area_class", f"{name}_sum", f"{name}_count"]
    return summary

data_frames = [shp_unused, shp_forest_theo, shp_forest_tech, shp_forest_high, shp_forest_med]
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
# #### Competing Use

# %%
shp_forest_water= find_overlap(shp_forest_tech,"water",shp_water_high)
shp_forest_solar = find_overlap(shp_forest_tech,"solar",shp_solar_high)
shp_forest_agri = find_overlap(shp_forest_tech,"agri",shp_agri_high)
shp_forest_housing = find_overlap(shp_forest_tech,"hsing",shp_housing_high)
shp_forest_industry = find_overlap(shp_forest_tech,"indus",shp_industry_high)
combined = pd.concat([shp_agri_high,shp_industry_high,shp_water_high,shp_solar_high,shp_housing_high])
combined.reset_index(inplace =True,drop =True)
inter = gpd.overlay(shp_forest_tech,combined,how ="intersection",keep_geom_type=True)
inter = inter.dissolve()
inter = add_area_and_class(inter)
data = {
    'forest_water': shp_forest_water.oparwater.sum(),
    'forest_solar': shp_forest_solar.oparsolar.sum(),
    'forest_agri': shp_forest_agri.oparagri.sum(),
    'forest_housing': shp_forest_housing.oparhsing.sum(),
    'forest_industry': shp_forest_industry.oparindus.sum(),
    'Competing use':inter.area_acres.sum(),
}
competing_use = pd.DataFrame(data, index=[0])
competing_use


# %%
S_list = [S1, S2, S3]

overlap_data = {}

land_use_categories = ['solar', 'water', 'agri', 'hsing', 'indus']
high_dataframes = [shp_solar_high, shp_water_high, shp_agri_high, shp_housing_high, shp_indus_high]

for idx, S in enumerate(S_list):
    S_op_dict = {}
    for category, high_df in zip(land_use_categories, high_dataframes):
        S_op_dict[category] = find_overlap(S, category, high_df)
    overlap_data[f"S{idx+1}"] = S_op_dict

data = []
for idx, S in enumerate(S_list):
    row_data = [
        overlap_data[f"S{idx+1}"]["solar"].oparsolar.sum(),
        overlap_data[f"S{idx+1}"]["water"].oparwater.sum(),
        overlap_data[f"S{idx+1}"]["agri"].oparagri.sum(),
        overlap_data[f"S{idx+1}"]["hsing"].oparhsing.sum(),
        overlap_data[f"S{idx+1}"]["indus"].oparindus.sum()
    ]
    data.append(row_data)

Dist_by_type = pd.DataFrame(data, columns=['solar', 'Water', 'agri', 'housing', 'industry'], index=['S1', 'S2', 'S3'])

# %%
Dist_by_type


# %%
largest_plot = shp_forest_tech.area_acres.max()

# %%
largest_plot

# %% [markdown]
# #### for reference excel file

# %%
template = load_workbook(get_rooted('workdir/forest/Forest_temp.xlsx'))
ws = template['Sheet1']

# Define the ranges where to insert the dataframes
ranges = ['A2:K5','A10:E10','B13:F15']

# Insert the dataframes into the workbook
for x, r in zip([merged_summary,competing_use,Dist_by_type], ranges):
    start_col, start_row, end_col, end_row = openpyxl.utils.cell.range_boundaries(r)
    for i, row in enumerate(x.values):
        for j, value in enumerate(row):
            ws.cell(row=start_row + i, column=start_col + j, value=value)
ws.cell(row=17, column=2, value=largest_plot)           
template.save(get_in_output('forest/Forest_temp_filled.xlsx'))

# %% [markdown]
# #### Top 15 
#
# ###### Assuming high file has > 15 lands
# ###### Assuming centroid is there

# %%
shp_forest_high = shp_forest_high.sort_values(by=["area_acres"],ascending = False)
shp_forest_high = shp_forest_high[:15]
shp_forest_high.reset_index(inplace =True,drop =True)
overlap_water = find_overlap(shp_forest_high,"water",shp_water_high)
overlap_solar = find_overlap(overlap_water,"solar",shp_solar_high)
overlap_agri = find_overlap(overlap_solar,"agri",shp_agri_high)
overlap_housing = find_overlap(overlap_agri,"hsing",shp_housing_high)
overlap_industry = find_overlap(overlap_housing,"indus",shp_industry_high)
overlap_industry = overlap_industry.to_crs(4326)

df = overlap_industry
slope = get_rooted("workdir/raster/slope.tif")
outputdf = calculate_slope(df, slope)


# %%
outputdf.columns

# %%
outputdf['coords'] = outputdf['geometry'].apply(lambda x: x.representative_point().coords[:])
outputdf['coords'] = [coords[0] for coords in outputdf['coords']]
outputdf["coords"].tolist()
outputdf[['lat', 'lon']] = gpd.GeoDataFrame(outputdf['coords'].tolist(), index=outputdf.index) 

# %%
outputdf_excel = outputdf[["lat","lon",'area_acres', 'area_class', 'op%water', 'oparwater',
       'cntwater', 'op%solar', 'oparsolar', 'cntsolar', 'op%agri', 'oparagri',
       'cntagri', 'op%hsing', 'oparhsing', 'cnthsing', 'op%indus', 'oparindus',
       'cntindus', 'min', 'max', 'mean']]

# outputdf_excel.to_excel(get_rooted("output/forest/Top15_forest.xlsx"))

# %% [markdown]
# ## Visuals

# %%
fig8, ax8 = plt.subplots(figsize=(5, 5))

plot_common_features(fig8, ax8)
plot_cities(fig8, ax8)


shp_district.plot(figsize=(5,5),color="none", ax=ax8, linewidth = 0.5, zorder=5)

shp_unused.plot(color="#424242",ax =ax8, label='No Potential')
shp_forest_theo.plot(color="#15915C",ax =ax8, label='Theoretical Potential')
shp_forest_tech.plot(color="#99CC66",ax =ax8, label='Technical Potential')


No_P = mpatches.Patch(color='#424242', label='No potential')
Theo_P = mpatches.Patch(color='#15915C', label='Theoretical potential')
Tech_P = mpatches.Patch(color='#99CC66', label='Technical potential')
    
plt.legend(handles = [No_P, Theo_P, Tech_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)


print(plt.rcParams['font.family'])


plt.savefig(get_in_output("forest/Technical_suitability.jpg"),dpi =1500)
plt.show()

# %%
fig9, ax9 = plt.subplots(figsize=(5, 5))

plot_common_features(fig9, ax9)
plot_cities(fig9, ax9)


shp_district.plot(figsize=(5,5),color="none", ax=ax9, linewidth = 0.5)

S1.plot(color="#455555",ax =ax9, label='>2.47 to 20 acres')
S2.plot(color="#54AD64",ax =ax9, label='>20 to 100 acres')
S3.plot(color="#99CC66",ax =ax9, label='>100 acres')


S1 = mpatches.Patch(color='#455555', label='>2.47 to 20 acres')
S2 = mpatches.Patch(color='#54AD64', label='>20 to 100 acres')
S3 = mpatches.Patch(color='#99CC66', label='>100 acres')
    
plt.legend(handles = [S1, S2, S3], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

print(plt.rcParams['font.family']);


plt.savefig(get_in_output("forest/Distribution by size.jpg"),dpi =1500)
plt.show()


# %%
fig10, ax10 = plt.subplots(figsize=(5, 5))

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)

shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

shp_forest_tech.plot(color="#21583B",ax =ax10, label='Low Potential')
shp_forest_med.plot(color="#009541",ax =ax10, label='Medium Potential')
shp_forest_high.plot(color="#99CC66",ax =ax10, label='High Potential')


Low_P = mpatches.Patch(color='#21583B', label='Low potential')
Med_P = mpatches.Patch(color='#009541', label='Medium potential')
High_P = mpatches.Patch(color='#99CC66', label='High potential')
    
plt.legend(handles = [Low_P, Med_P, High_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

print(plt.rcParams['font.family'])


plt.savefig(get_in_output("forest/High Potential_H_M_L.jpg"),dpi =1500)
plt.show()


# %% [markdown]
# ## Settlement Analysis

# %%
popdf = pd.DataFrame()
for j in range(len(shp_village)):
    input_shp =  get_rooted('temp_f.shp')

    selection = shp_village.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted(population_raster)

        output_raster =  get_rooted('temp_f.tif')
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
        popdf.at[j, "totpop_ras"]=0


    else:    

        pop_sum = rasterarr.sum()

        popdf.at[j, "totpop_ras"]=pop_sum

shp_village_final = pd.concat([shp_village, popdf], axis = 1)

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
shp_village_tech = find_overlap(shp_village_final,"tech",shp_forest_tech)
shp_village_tech =  find_overlap(shp_village_tech,"unused",shp_unused)

# %% [markdown]
# #### add / change inputs 

# %%
shp_village_final = shp_village_tech[['p_name','d_name', 'b_name', 'p_name_rd',
       'pancha_id', 'block_id_o', 'dist_id', 'district_r', 'gp_code', "totpop",
       'lat', 'lon', 'totpix', 'f_totarkm2', 'totarkm2', 'fcover%', 'fper1000', 'totpop_ras',
       'TGA(acres)', 'op%tech', 'opartech', 'cnttech', 'op%unused',
       'oparunused', 'cntunused', 'geometry']]

# %%
# shp_village_final.to_file(get_in_output("forest/settlement"))

# %%
shp_village = shp_village_tech[['p_name', 'd_name', 'b_name', 'p_name_rd', 'pancha_id', 'block_id_o',
       'dist_id', 'district_r', 'gp_code', 'totpop', 'lat', 'lon', 'totpix',
       'f_totarkm2', 'totarkm2', 'fcover%', 'fper1000', 'totpop_ras',
       'TGA(acres)', 'op%tech', 'opartech', 'cnttech', 'op%unused',
       'oparunused', 'cntunused', 'geometry']]

# %%
shp_village.sort_values("p_name",inplace =True)

# %%
shp_village.reset_index(drop =True,inplace =True)

# %%
shp_village.drop("geometry",axis =1,inplace =True)

# %%
# shp_village.to_excel(get_in_output("forest/230517_settlement_ordered.xlsx"))

# %% [markdown]
# #### the below code is for clipping 

# %%
# shp_village_tech = shp_village_tech[shp_village_tech["TGA(acres)"] > 100]

# %% [markdown]
# ### Visuals

# %%
settlement_final = read_df_UT("output/forest/settlement/settlement.shp")

# %%
settlement_final.columns

# %%
C1 = LinearSegmentedColormap.from_list('testCmap1', colors=["#F5FAF7", "#A4EDC1", "#46C779", "#149169", "#306151"], N=256)

# %%
fig, ax = plt.subplots(figsize=(5, 5))

plot_common_features(fig, ax)  # assuming this function is defined elsewhere

shp_district.plot(figsize=(5,5), color="none", ax=ax, linewidth=0.5, zorder=5)
settlement_final.plot(column='op%tech', cmap=C1, ax=ax)
a = settlement_final["op%tech"].min()
b = settlement_final["op%tech"].max()
sm = plt.cm.ScalarMappable(cmap=C1)
cbaxes = fig.add_axes([0.7, 0.18, 0.2, 0.02]) 
cbar = plt.colorbar(sm, orientation = 'horizontal', cax=cbaxes, shrink = 0.2)
cbar.mappable.set_clim(vmin = a, vmax = b)
cbar.ax.tick_params(labelsize=3, color = 'grey')
cbar.outline.set_visible(False)
cbar.ax.set_title('Technical potential area(%)\n per settlement', fontsize=5)

plt.savefig(get_in_output("forest/Technical potential area(%).jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ## Taluk

# %% [markdown]
# same as settlement except with taluk boundaries

# %%
shp_taluk = read_df_UT('workdir/extra_inputs/shp_taluk/shp_taluk.shp')

# %%
shp_taluk = gpd.overlay(shp_district, shp_taluk, how ='intersection')

# %%
popdf = pd.DataFrame()
for j in range(len(shp_taluk)):
    input_shp =  get_rooted('temp5.shp')

    selection = shp_taluk.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted(population_raster)

        output_raster =  get_rooted('temp5.tiff')
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
        popdf.at[j, "totpop_ras"]=0


    else:    

        pop_sum = rasterarr.sum()

        popdf.at[j, "totpop_ras"]=pop_sum

shp_taluk_final = pd.concat([shp_taluk, popdf], axis = 1)

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
shp_taluk_tech = find_overlap(shp_taluk_final,"tech",shp_forest_tech)


# %%
shp_taluk_tech = shp_taluk_tech.to_crs(4326)

# %%
shp_taluk_tech= shp_taluk_tech[['Taluk_name','lat', 'lon','totpop_ras',
       'TGA(acres)', 'op%tech','opartech', 'cnttech', 'geometry']]

# %%
shp_taluk_tech = shp_taluk_tech.to_crs(4326)

# %%
shp_taluk_tech.to_file(get_in_output("forest/shp_taluk_tech"))

# %%
shp_taluk_tech.columns

# %%
shp_taluk_tech.drop("geometry", axis = 1)
# shp_taluk_tech.to_excel(get_in_output("forest/Taluk_forest.xlsx"))
