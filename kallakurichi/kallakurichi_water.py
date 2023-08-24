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
# # Water Report

# %%
import Input_template
from Input_template import *

# %%
## boundary files
district_boundary = "workdir/extra_inputs/shp_district/shp_district.shp"
village = "workdir/extra_inputs/shp_settlement/shp_settlement.shp"
taluk = "workdir/extra_inputs/shp_taluk/shp_taluk.shp"
## module files
technical ="workdir/water/shp_water_tech/shp_water_tech.shp"
high = "workdir/water/shp_water_high/shp_water_high.shp"
med  = "workdir/water/shp_water_med/shp_water_med.shp"
low ="workdir/water/shp_water_low/shp_water_low.shp"
## other module high files
agri_high = "workdir/agri/shp_agri_high/shp_agri_high.shp"
forest_high = "workdir/forest/shp_forest_high/shp_forest_high.shp"
housing_high = "workdir/housing/shp_housing_high/shp_housing_high.shp"
solar_high = "workdir/solar/shp_solar_high/shp_solar_high.shp"
industry_high = "workdir/industry/shp_industry_high/shp_industry_high.shp"
## landcover shape files
Unused_land = "workdir/extra_inputs/Unused/Unused.shp"
waterbodies = "workdir/extra_inputs/shp_water_bodies/shp_water_bodies.shp"
## extra inputs
slope =get_rooted("workdir/raster/slope.tif")

# %%
## boundary files
shp_district = read_df_UT(district_boundary)
shp_village = read_df_UT(village)
shp_taluk = read_df_UT(taluk)
## module files
shp_water_tech = read_df_UT(technical)
shp_water_high =read_df_UT(high)
shp_water_med = read_df_UT(med)
shp_water_low = read_df_UT(low)
## other module high files
shp_agri_high = read_df_UT(agri_high)
shp_forest_high = read_df_UT(forest_high)
shp_housing_high = read_df_UT(housing_high)
shp_solar_high = read_df_UT(solar_high)
shp_indus_high =read_df_UT(industry_high)
## landcover shape files
shp_unused = read_df_UT(Unused_land)
shp_waterbodies = read_df_UT(waterbodies)


# %%
shp_water_low = add_area_and_class(shp_water_low)
shp_unused = add_area_and_class(shp_unused)
shp_waterbodies = add_area_and_class(shp_waterbodies)

# %% [markdown]
# ## exporting this files for HTML

# %%
shp_water_tech[['WD_Rating','RunBL_Rat',
       'Fin_Rating']].head(25)

# %%
# shp_WD_high = shp_water_tech[shp_water_tech["WD_Rating"] == "H"]
# shp_WD_med = shp_water_tech[shp_water_tech["WD_Rating"] == "M"]
# shp_WD_low = shp_water_tech[shp_water_tech["WD_Rating"] == "L"]
# shp_RO_high = shp_water_tech[shp_water_tech["RunBL_Rat"] == "H"]
# shp_RO_med = shp_water_tech[shp_water_tech["RunBL_Rat"] == "M"]
# shp_RO_low = shp_water_tech[shp_water_tech["RunBL_Rat"] == "L"]

# %%
# shp_WD_high.to_file(get_in_workdir("water/shp_WD_high"))
# shp_WD_med.to_file(get_in_workdir("water/shp_WD_med"))
# shp_WD_low.to_file(get_in_workdir("water/shp_WD_low"))
# shp_RO_high.to_file(get_in_workdir("water/shp_RO_high"))
# shp_RO_med.to_file(get_in_workdir("water/shp_RO_med"))
# shp_RO_low.to_file(get_in_workdir("water/shp_RO_low"))

# %% [markdown]
# ## surface and ground water files

# %%
shp_sw_tech = read_df_UT("workdir/water/sw_tech/sw_tech.shp")
shp_sw_high = read_df_UT("workdir/water/sw_high/sw_high.shp")
shp_sw_med = read_df_UT("workdir/water/sw_med/sw_med.shp")
shp_sw_low = read_df_UT("workdir/water/sw_low/sw_low.shp")

# %%
shp_gw_tech = read_df_UT("workdir/water/gw_tech/gw_tech.shp")
shp_gw_high = read_df_UT("workdir/water/gw_high/gw_high.shp")
shp_gw_med = read_df_UT("workdir/water/gw_med/gw_med.shp")
shp_gw_low = read_df_UT("workdir/water/gw_low/gw_low.shp")

# %%
shp_sw_tech = add_area_and_class(shp_sw_tech)
shp_sw_high = add_area_and_class(shp_sw_high)
shp_sw_med = add_area_and_class(shp_sw_med)
shp_sw_low = add_area_and_class(shp_sw_low)
shp_gw_tech = add_area_and_class(shp_gw_tech)
shp_gw_high = add_area_and_class(shp_gw_high)
shp_gw_med = add_area_and_class(shp_gw_med)
shp_gw_low = add_area_and_class(shp_gw_low)

# %%
shp_sw_tech = shp_sw_tech[["area_acres","geometry"]]
shp_sw_high = shp_sw_high[["area_acres","geometry"]]
shp_sw_med = shp_sw_med[["area_acres","geometry"]]
shp_sw_low = shp_sw_low[["area_acres","geometry"]]
shp_gw_tech = shp_gw_tech[["area_acres","geometry"]]
shp_gw_high = shp_gw_high[["area_acres","geometry"]]
shp_gw_med = shp_gw_med[["area_acres","geometry"]]
shp_gw_low = shp_gw_low[["area_acres","geometry"]]

# %%
shp_sw_tech.area_acres.sum(),shp_sw_tech.shape

# %%
shp_water_tech.area_acres.sum(),shp_water_tech.shape


# %%
# check the code

# %%
# shp_sw_tech = gpd.overlay(shp_sw_tech,shp_water_tech,how= "intersection")
# shp_sw_tech.reset_index(inplace =True,drop =True)
# shp_sw_tech = add_area_and_class(shp_sw_tech)
# shp_sw_tech = shp_sw_tech[shp_sw_tech["area_acres"] > 0.01]

# %%
# shp_sw_low = gpd.overlay(shp_sw_low,shp_water_tech,how= "intersection")
# shp_sw_low.reset_index(inplace =True,drop =True)
# shp_sw_low = add_area_and_class(shp_sw_low)
# shp_sw_low = shp_sw_low[shp_sw_low["area_acres"] > 0.01]

# %% [markdown]
# ## Report Tables

# %% [markdown]
# ### Stats

# %%
def summarize_data(df, name):
    summary = df.groupby("area_class")["area_acres"].agg(["sum", "count"]).reset_index()
    summary.columns = ["area_class", f"{name}_sum", f"{name}_count"]
    return summary

data_frames = [shp_unused,shp_water_tech, shp_water_high, shp_water_med, shp_water_low]
data_names = ['unused', 'tech','high', 'med'.'low']

summary_dict = {}

for df, name in zip(data_frames, data_names):
    summary_dict[name] = summarize_data(df, name)

merged_summary = summary_dict['unused']
for name in data_names[1:]:
    merged_summary = pd.merge(merged_summary, summary_dict[name], on="area_class", how="outer")

merged_summary = merged_summary.fillna(0)
merged_summary

# %% [markdown]
# ## Additional stats (Run off)

# %%
shp_water_tech.area_acres.sum(),shp_water_tech.Run_Tot.sum(),shp_water_tech.shape

# %%
shp_water_high.area_acres.sum(),shp_water_high.Run_Tot.sum(),shp_water_high.shape

# %%
shp_water_med.area_acres.sum(),shp_water_med.Run_Tot.sum(),shp_water_med.shape

# %%
shp_water_low.area_acres.sum(),shp_water_low.Run_Tot.sum(),shp_water_low.shape

# %% [markdown]
# ## Distribution by type
#

# %%
overlap= find_overlap(shp_water_tech,"sw",shp_sw_tech)
overlap = find_overlap(overlap,"gw",shp_gw_tech)
conditions = [
    (overlap['oparsw'] > 0) & (overlap['opargw'] == 0),
    (overlap['opargw'] > 0) & (overlap['oparsw'] == 0),
    (overlap['oparsw'] > 0) & (overlap['opargw'] > 0),
    (overlap['oparsw'] == 0) & (overlap['opargw'] == 0)
]
labels = ['SW', 'GW', 'both', 'nothing']
# use numpy.select to create new column
overlap['Type'] = np.select(conditions, labels, default='unknown')
df = overlap.groupby("Type")["area_acres"].agg(["sum", "count"])
print("AREA :",df)
df_run = overlap.groupby("Type")["Run_Tot"].agg(["sum", "count"])
print("Run off :",df_run)

# %% [markdown]
# ## exporting the files

# %%
overlap = overlap[['oparsw','opargw',"Type","geometry"]]
overlap = overlap.to_crs(4326)
sw = overlap[overlap["Type"] == "SW"]
gw = overlap[overlap["Type"] == "GW"]
both = overlap[overlap["Type"] == "both"]
sw.reset_index(drop=True,inplace =True)
gw.reset_index(drop=True,inplace =True)
both.reset_index(drop=True,inplace =True)

# %%
# sw.to_file(get_in_output("water/sw"))
# gw.to_file(get_in_output("water/gw"))
# both.to_file(get_in_output("water/both"))

# %%
overlap_low= find_overlap(shp_water_low,"sw",shp_sw_low)
overlap_low = find_overlap(overlap_low,"gw",shp_gw_low)
conditions = [
    (overlap_low['oparsw'] > 0) & (overlap_low['opargw'] == 0),
    (overlap_low['opargw'] > 0) & (overlap_low['oparsw'] == 0),
    (overlap_low['oparsw'] > 0) & (overlap_low['opargw'] > 0),
    (overlap_low['oparsw'] == 0) & (overlap_low['opargw'] == 0)
]
labels = ['SW', 'GW', 'both', 'nothing']

# use numpy.select to create new column
overlap_low['Type'] = np.select(conditions, labels, default='unknown')
df_low = overlap_low.groupby("Type")["area_acres"].agg(["sum", "count"])
print("AREA :",df_low)
df_low_run = overlap_low.groupby("Type")["Run_Tot"].agg(["sum", "count"])
print("Run off :",df_low_run)

# %%
overlap_med= find_overlap(shp_water_med,"sw",shp_sw_med)
overlap_med = find_overlap(overlap_med,"gw",shp_gw_med)
conditions = [
    (overlap_med['oparsw'] > 0) & (overlap_med['opargw'] == 0),
    (overlap_med['opargw'] > 0) & (overlap_med['oparsw'] == 0),
    (overlap_med['oparsw'] > 0) & (overlap_med['opargw'] > 0),
    (overlap_med['oparsw'] == 0) & (overlap_med['opargw'] == 0)
]
labels = ['SW', 'GW', 'both', 'nothing']

# use numpy.select to create new column
overlap_med['Type'] = np.select(conditions, labels, default='unknown')
df_med = overlap_med.groupby("Type")["area_acres"].agg(["sum", "count"])
print("AREA :",df_med)
df_med_run = overlap_med.groupby("Type")["Run_Tot"].agg(["sum", "count"])
print("Run off :",df_med_run)

# %%
overlap_high= find_overlap(shp_water_high,"sw",shp_sw_high)
overlap_high = find_overlap(overlap_high,"gw",shp_gw_high)
conditions = [
    (overlap_high['oparsw'] > 0) & (overlap_high['opargw'] == 0),
    (overlap_high['opargw'] > 0) & (overlap_high['oparsw'] == 0),
    (overlap_high['oparsw'] > 0) & (overlap_high['opargw'] > 0),
    (overlap_high['oparsw'] == 0) & (overlap_high['opargw'] == 0)
]
labels = ['SW', 'GW', 'both', 'nothing']

overlap_high['Type'] = np.select(conditions, labels, default='unknown')
df_high = overlap_high.groupby("Type")["area_acres"].agg(["sum", "count"])
print("AREA :",df_high)
df_high_run = overlap_high.groupby("Type")["Run_Tot"].agg(["sum", "count"])
print("Run off :",df_high_run)

# %% [markdown]
# ### Competing use (need to improve)

# %%
shp_water_agri= find_overlap(shp_water_tech,"agri",shp_agri_high)
shp_water_forest = find_overlap(shp_water_tech,"forest",shp_forest_high)
shp_water_housing =find_overlap(shp_water_tech,"hsing",shp_housing_high)
shp_water_indus =find_overlap(shp_water_tech,"indus",shp_indus_high)
shp_water_solar =find_overlap(shp_water_tech,"solar",shp_solar_high)
data = {
    'water_forest': shp_water_forest.oparforest.sum(),
    'water_agri': shp_water_agri.oparagri.sum(),
    'water-housing':shp_water_housing.oparhsing.sum(),
    'water-indus': shp_water_indus.oparindus.sum(),
    'water-solar': shp_water_solar.oparsolar.sum()
}
competing_use = pd.DataFrame(data, index=[0])
competing_use

# %%
largest_plot = shp_water_tech.area_acres.max()
largest_plot

# %%
template = load_workbook(get_rooted('workdir/water/water_temp.xlsx'))
ws = template['Sheet1']

# Define the ranges where to insert the dataframes
ranges = ['A2:K5','A10:E10','B21:C23','B27:C9']

# Insert the dataframes into the workbook
for x, r in zip([merged_summary,competing_use,df,df_low], ranges):
    start_col, start_row, end_col, end_row = openpyxl.utils.cell.range_boundaries(r)
    for i, row in enumerate(x.values):
        for j, value in enumerate(row):
            ws.cell(row=start_row + i, column=start_col + j, value=value)
ws.cell(row=17, column=2, value=largest_plot)           
template.save(get_in_output('water/water_temp_filled.xlsx'))

# %% [markdown]
# ## Top15

# %%
shp_water_high_sorted = shp_water_high.sort_values(by=["area_acres"],ascending = False)

# %%
shp_water_top15 = shp_water_high_sorted[:15]
shp_water_top15.reset_index(inplace =True,drop =True)

# %%
df = shp_water_top15
slope =get_rooted("workdir/raster/slope.tif")
top15_water = calculate_slope(df, slope)

# %%
top15_water['coords'] = top15_water['geometry'].apply(lambda x: x.representative_point().coords[:])
top15_water['coords'] = [coords[0] for coords in top15_water['coords']]
top15_water["coords"].tolist()
top15_water[['lat', 'lon']] = gpd.GeoDataFrame(top15_water['coords'].tolist(), index=top15_water.index)

# %%
overlap_final = find_overlap(top15_water,"solar",shp_solar_high)
overlap_final = find_overlap(overlap_final,"forest",shp_forest_high)
overlap_final = find_overlap(overlap_final,"agri",shp_agri_high)
overlap_final = find_overlap(overlap_final,"hsing",shp_housing_high)
overlap_final = find_overlap(overlap_final,"indus",shp_indus_high)

# %%
overlap_final.iloc[:,56:71] = overlap_final.iloc[:,56:71].astype(float)

# %%
overlap_final_top15 = overlap_final[["lat","lon","area_acres","Run_Tot","HiETAr%","min","max",'op%solar', 'oparsolar',
       'cntsolar', 'op%forest', 'oparforest', 'cntforest', 'op%agri',
       'oparagri', 'cntagri', 'op%hsing', 'oparhsing', 'cnthsing', 'op%indus',
       'oparindus', 'cntindus','geometry']]

# %%
overlap_final_top15.to_file(get_in_output("water/Top15_final"))

# %%
overlap_final_excel = overlap_final[["lat","lon","area_acres","Run_Tot","HiETAr%","min","max",'op%solar', 'oparsolar',
       'cntsolar', 'op%forest', 'oparforest', 'cntforest', 'op%agri',
       'oparagri', 'cntagri', 'op%hsing', 'oparhsing', 'cnthsing', 'op%indus',
       'oparindus', 'cntindus']]

# %%
# overlap_final_excel.to_excel(get_in_output("water/Top15_final.xlsx"))

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
shp_water_tech.plot(color="#07c4ff",ax =ax8, label='Technical Potential')



No_P = mpatches.Patch(color='#424242', label='No potential')
Tech_P = mpatches.Patch(color='#07c4ff', label='Technical potential')
    
plt.legend(handles = [No_P, Tech_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)


print(plt.rcParams['font.family'])


plt.savefig(get_in_output("images/water/Technical_suitability.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ## High Potential

# %%
fig10, ax10 = plt.subplots(figsize=(5, 5))

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)

shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

shp_water_low.plot(color="#0558a6",ax =ax10, label='Low Potential')
shp_water_med.plot(color="#a8a2a2",ax =ax10, label='Medium Potential')
shp_water_high.plot(color="#07c4ff",ax =ax10, label='High Potential')


Low_P = mpatches.Patch(color='#0558a6', label='Low potential')
Med_P = mpatches.Patch(color='#a8a2a2', label='Medium potential')
High_P = mpatches.Patch(color='#07c4ff', label='High potential')
    
plt.legend(handles = [Low_P, Med_P, High_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

print(plt.rcParams['font.family'])


plt.savefig(get_in_output("images/water/High Potential_H_M_L.jpg"),dpi =1500)
plt.show()

# %%
fig10, ax10 = plt.subplots(figsize=(5, 5))

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)

shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

shp_WD_low.plot(color="#335342",ax =ax10, label='Low water demand')
shp_WD_med.plot(color="#77b574",ax =ax10, label='Medium water demand')
shp_WD_high.plot(color="#c5d3ae",ax =ax10, label='High water demand')


Low_P = mpatches.Patch(color='#335342', label='Low water demand')
Med_P = mpatches.Patch(color='#77b574', label='Medium water demand')
High_P = mpatches.Patch(color='#c5d3ae', label='High water demand')
    
plt.legend(handles = [Low_P, Med_P, High_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

print(plt.rcParams['font.family'])
plt.savefig(get_in_output("images/water/water_demand_H_M_L.jpg"),dpi =1500)
plt.show()

# %%
fig10, ax10 = plt.subplots(figsize=(5, 5))

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)

shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

shp_RO_low.plot(color="#95daff",ax =ax10, label='Low run off')
shp_RO_med.plot(color="#185fbb",ax =ax10, label='Medium run off')
shp_RO_high.plot(color="#123e71",ax =ax10, label='High run off')


Low_P = mpatches.Patch(color='#123e71', label='Low run off')
Med_P = mpatches.Patch(color='#185fbb', label='Medium run off')
High_P = mpatches.Patch(color='#95daff', label='High run off')
    
plt.legend(handles = [Low_P, Med_P, High_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)
print(plt.rcParams['font.family'])
plt.savefig(get_in_output("images/water/run_off_H_M_L.jpg"),dpi =1500)
plt.show()

# %%
fig10, ax10 = plt.subplots(figsize=(5, 5))

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)

shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

sw.plot(color="#56aeff",ax =ax10, label='Surface water')
gw.plot(color="#326e80",ax =ax10, label='Ground water')
both.plot(color="#6ce6f8",ax =ax10, label='Surface & ground water')


sw_P = mpatches.Patch(color='#56aeff', label='Surface water')
gw_P = mpatches.Patch(color='#326e80', label='Ground water')
both_P = mpatches.Patch(color='#6ce6f8', label='Surface & ground water')
    
plt.legend(handles = [sw_P, gw_P, both_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

print(plt.rcParams['font.family'])


plt.savefig(get_in_output("images/water/SW_GW_both.jpg"),dpi =1500)
plt.show()

# %%
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

from matplotlib.patches import ConnectionPatch

fig10, ax10 = plt.subplots(figsize=(5, 5))

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)

shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

shp_water_low.plot(color="#0558a6",ax =ax10, label='Low Potential')
shp_water_med.plot(color="#a8a2a2",ax =ax10, label='Medium Potential')
shp_water_high.plot(color="#07c4ff",ax =ax10, label='High Potential')

Low_P = mpatches.Patch(color='#0558a6', label='Low potential')
Med_P = mpatches.Patch(color='#a8a2a2', label='Medium potential')
High_P = mpatches.Patch(color='#07c4ff', label='High potential')

plt.legend(handles=[Low_P, Med_P, High_P], loc='upper left', bbox_to_anchor=(1.2, 0.2), title='Legend\n', fontsize=5.5, markerscale=2, title_fontsize=5.5, framealpha=0, borderpad=0.3, handletextpad=0.5, handlelength=1.0)


axins1 = zoomed_inset_axes(ax10, zoom=4)  # create new set of axes
shp_water_low.plot(color="#0558a6",ax =axins1, label='Low Potential')

shp_water_med.plot(color="#a8a2a2",ax =axins1, label='Medium Potential')
shp_water_high.plot(color="#07c4ff",ax =axins1, label='High Potential') # plot data on new axes
   
axins1.tick_params(axis='x', colors='grey', labelsize=3, labelrotation = 270)
axins1.tick_params(axis='y', colors='grey', labelsize=3) #reducing the size of the axis values
axins1.xaxis.set_major_formatter(FormatStrFormatter('%.2f')) #axis value formatting for both axis
axins1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins1.set_xlim(79.38, 79.44)  # set limits of new axes
axins1.set_ylim(11.68, 11.74)

axins1.spines['bottom'].set_color('grey')
axins1.spines['top'].set_color('grey') 
axins1.spines['right'].set_color('grey')
axins1.spines['left'].set_color('grey')


plt.savefig(get_in_output("images/water/test.jpg"),dpi =1500)
plt.show()

# %%
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

fig10, ax10 = plt.subplots()

# plot common features and cities
plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)
fig10.subplots_adjust(left=0.05, right=0.8)
# plot water potentials
shp_district.plot(figsize=(5,5), color="none", ax=ax10, linewidth=0.5, zorder=1)
shp_water_low.plot(color="#0558a6", ax=ax10, label='Low Potential')
shp_water_med.plot(color="#a8a2a2", ax=ax10, label='Medium Potential')
shp_water_high.plot(color="#07c4ff", ax=ax10, label='High Potential')

# create legend
Low_P = mpatches.Patch(color='#0558a6', label='Low potential')
Med_P = mpatches.Patch(color='#a8a2a2', label='Medium potential')
High_P = mpatches.Patch(color='#07c4ff', label='High potential')
plt.legend(handles=[Low_P, Med_P, High_P], loc='upper left', bbox_to_anchor=(0.8, 0.2), title='Legend\n', fontsize=5.5, markerscale=2, title_fontsize=5.5, framealpha=0, borderpad=0.3, handletextpad=0.5, handlelength=1.0)

# create zoomed inset axes
axins1 = zoomed_inset_axes(ax10, zoom=4, bbox_to_anchor=(1,0.75 ), bbox_transform=ax10.figure.transFigure)

# plot water potentials on zoomed inset axes
shp_water_low.plot(color="#0558a6", ax=axins1, label='Low Potential')
shp_water_med.plot(color="#a8a2a2", ax=axins1, label='Medium Potential')
shp_water_high.plot(color="#07c4ff", ax=axins1, label='High Potential')

# adjust tick parameters and limits of zoomed inset axes
axins1.tick_params(axis='both', colors='grey', labelsize=3)
axins1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins1.set_xlim(79.38, 79.44)
axins1.set_ylim(11.68, 11.74)

# Modify the edge color of the inset axes
axins1.spines['top'].set_edgecolor('grey')
axins1.spines['right'].set_edgecolor('grey')
axins1.spines['bottom'].set_edgecolor('grey')
axins1.spines['left'].set_edgecolor('grey')
# Modify the line thickness of the inset axes spines
axins1.spines['top'].set_linewidth(0.5)
axins1.spines['right'].set_linewidth(0.5)
axins1.spines['bottom'].set_linewidth(0.5)
axins1.spines['left'].set_linewidth(0.5)

# create zoomed inset axes
axins2 = zoomed_inset_axes(ax10, zoom=4, bbox_to_anchor=(1,0.5), bbox_transform=ax10.figure.transFigure)

# plot water potentials on zoomed inset axes
shp_water_low.plot(color="#0558a6", ax=axins2, label='Low Potential')
shp_water_med.plot(color="#a8a2a2", ax=axins2, label='Medium Potential')
shp_water_high.plot(color="#07c4ff", ax=axins2, label='High Potential')

# Modify the edge color of the inset axes
axins2.spines['top'].set_edgecolor('grey')
axins2.spines['right'].set_edgecolor('grey')
axins2.spines['bottom'].set_edgecolor('grey')
axins2.spines['left'].set_edgecolor('grey')
# Modify the line thickness of the inset axes spines
axins2.spines['top'].set_linewidth(0.5)
axins2.spines['right'].set_linewidth(0.5)
axins2.spines['bottom'].set_linewidth(0.5)
axins2.spines['left'].set_linewidth(0.5)


# adjust tick parameters and limits of zoomed inset axes
axins2.tick_params(axis='both', colors='grey', labelsize=3)
axins2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins2.set_xlim(79.15, 79.21)
axins2.set_ylim(11.85, 11.91)

plt.savefig(get_in_output("images/water/test.jpg"), dpi=1500)
plt.show()

# %% [markdown]
# # Settlement analysis

# %%
popdf = pd.DataFrame()
for j in range(len(shp_village)):
    input_shp =  get_rooted('workdir/temp_water.shp')

    selection = shp_village.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted("workdir/raster/population_10_lon_70_general-v1.5.tif")

        output_raster = get_rooted('workdir/temp_water.tif')
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
shp_village_final['coords'] = shp_village_final['geometry'].apply(lambda x: x.representative_point().coords[:])
shp_village_final['coords'] = [coords[0] for coords in shp_village_final['coords']]
shp_village_final["coords"].tolist()
shp_village_final[['lat', 'lon']] = gpd.GeoDataFrame(shp_village_final['coords'].tolist(), index=shp_village_final.index) 

# %%
shp_village_final = shp_village_final.to_crs(32644)
shp_village_final["TGA(acres)"] = ((shp_village_final.geometry.area)/10**6)*247.105
shp_village_final = shp_village_final.to_crs(4326)

# %%
shp_village_tech = find_overlap(shp_village_final,"tech",shp_water_tech)
shp_village_tech = find_overlap(shp_village_tech,"sw",sw)
shp_village_tech = find_overlap(shp_village_tech,"gw",gw)
shp_village_tech = find_overlap(shp_village_tech,"both",both)
shp_village_tech = find_overlap(shp_village_tech,"WB",shp_waterbodies)

# %%
# shp_village_tech =  find_overlap(shp_village_tech,"unused",shp_unused) ## calculating this in other modules

# %%
shp_village_tech = shp_village_tech.to_crs(4326)

# %%
shp_village_tech = shp_village_tech.iloc[:, [14,15,16,17] + list(range(44, 67))]

# %%
shp_village_tech.columns

# %%
shp_village_tech.iloc[:,12:27] = shp_village_tech.iloc[:,12:27].astype(float)

# %%
shp_village_tech= shp_village_tech[['p_name', 'd_name', 'b_name', 'p_name_rd', 'area_acres', 'area_class'
       , 'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech',
       'opartech', 'cnttech', 'op%sw', 'oparsw', 'cntsw', 'op%gw', 'opargw',
       'cntgw','op%both', 'oparboth', 'cntboth', 'op%WB', 'oparWB', 'cntWB','geometry']]

# %%
# shp_village_tech.to_file(get_in_output("water/settlement"))

# %%
shp_village_tech.columns

# %%
shp_village = shp_village_tech[['p_name', 'd_name', 'b_name', 'p_name_rd', 'area_acres', 'area_class',
       'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech', 'opartech', 'cnttech','op%sw', 'oparsw', 'cntsw', 'op%gw', 'opargw',
       'cntgw','op%both', 'oparboth', 'cntboth', 'op%WB', 'oparWB', 'cntWB']]

# %%
# shp_village.to_excel(get_in_output("water/settlement.xlsx"))

# %%
shp_village.sort_values("p_name",inplace =True)

# %%
shp_village.reset_index(drop =True,inplace =True)

# %%
shp_village_tech.drop("geometry",axis =1,inplace =True)

# %%
# shp_village.to_excel(get_in_output("water/230517_settlement_ordered.xlsx"))

# %% [markdown]
# # Settlement visuals

# %%
settlement = read_df_UT("output/water/settlement/settlement.shp")

# %%
S3 = LinearSegmentedColormap.from_list('testCmap1', colors=["#ecfaff", "#c3efff", "#7abddc", "#5292b6", "#3b5573"], N=256)

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
cbar.ax.set_title('Technical potential area(%)\n per settlement', fontsize=5)

print(a,b)

plt.savefig(get_in_output("images/water/Technical potential area(%).jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ### Taluk analysis

# %%
from shapely.geometry import shape


# %%
shp_taluk = "workdir/extra_inputs/shp_taluk/shp_taluk.shp"

# %%
popdf = pd.DataFrame()
for j in range(len(shp_taluk)):
    input_shp =  get_rooted('workdir/temp_water_new.shp')

    selection = shp_taluk.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted("workdir/raster/population_10_lon_70_general-v1.5.tif")

        output_raster = get_rooted('workdir/temp_water_new.tif')
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
shp_taluk_final['coords'] = shp_taluk_final['geometry'].apply(lambda x: x.representative_point().coords[:])
shp_taluk_final['coords'] = [coords[0] for coords in shp_taluk_final['coords']]
shp_taluk_final["coords"].tolist()
shp_taluk_final[['lat', 'lon']] = gpd.GeoDataFrame(shp_taluk_final['coords'].tolist(), index=shp_taluk_final.index) 

# %%
shp_taluk_final = shp_taluk_final.to_crs(32644)
shp_taluk_final["TGA(acres)"] = ((shp_taluk_final.geometry.area)/10**6)*247.105
shp_taluk_final = shp_taluk_final.to_crs(4326)

# %%
shp_taluk_final = find_overlap(shp_taluk_final,"tech",shp_water_tech)

# %%
# shp_taluk_final['geometry'] = shp_taluk_final['geometry'].apply(lambda x: shape(x).buffer(0).buffer(0.0000000000001))
# shp_taluk_finalshp_taluk_final = gpd.GeoDataFrame(shp_taluk_final, geometry='geometry')
# shp_taluk_final = shp_taluk_final.loc[shp_taluk_final.is_valid]

# %%
# shp_waterbodies = read_df_UT("workdir/extra_inputs/shp_waterbodies/shp_waterbodies.shp")
shp_sw =  read_df_UT("output/water/sw/sw.shp")
shp_gw = read_df_UT("output/water/gw/gw.shp")
shp_both = read_df_UT("output/water/both/both.shp")

# %%
shp_taluk_final =  find_overlap(shp_taluk_final,"WB",shp_waterbodies)

# %%
shp_taluk_final =  find_overlap(shp_taluk_final,"sw",shp_sw)
shp_taluk_final =  find_overlap(shp_taluk_final,"gw",shp_gw)
shp_taluk_final =  find_overlap(shp_taluk_final,"both",shp_both)

# %%
shp_taluk_final = shp_taluk_final.to_crs(4326)

# %%
shp_taluk_final= shp_taluk_final[['Taluk_name', 'area_acres', 'area_class'
       , 'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech',
       'opartech', 'cnttech','op%WB','oparWB','cntWB','cntsw',"cntgw","cntboth",'geometry']]

# %%
shp_taluk_final.totpop.sum()

# %%
# shp_taluk_final.to_file(get_in_output("water/shp_taluk"))

# %%
shp_taluk_final.drop("geometry",axis =1,inplace =True)

# %%
# shp_taluk_final.to_excel(get_in_output("water/taluk_analysis.xlsx"))
