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
# # Industry Report

# %%
import Input_template
from Input_template import *

# %%
## boundary files
district_boundary = "workdir/extra_inputs/shp_district/shp_district.shp"
village = "workdir/extra_inputs/shp_settlement/shp_settlement.shp"
taluk = "workdir/extra_inputs/shp_taluk/shp_taluk.shp"
## module files
theoritical = "workdir/industry/shp_industry_theo/shp_industry_theo.shp"
technical ="workdir/industry/shp_industry_tech/shp_industry_tech.shp"
high = "workdir/industry/shp_industry_high/shp_industry_high.shp"
med  = "workdir/industry/shp_industry_med/shp_industry_med.shp"
## other module high files
water_high = "workdir/water/shp_water_high/shp_water_high.shp"
forest_high = "workdir/forest/shp_forest_high/shp_forest_high.shp"
agri_high = "workdir/agri/shp_agri_high/shp_agri_high.shp"
solar_high = "workdir/solar/shp_solar_high/shp_solar_high.shp"
housing_high = "workdir/housing/shp_housing_high/shp_housing_high.shp"
## landcover shape files
Unused_land = "workdir/extra_inputs/Unused/Unused.shp"
Sparseveg = "workdir/extra_inputs/shp_sparseveg/shp_sparseveg.shp"
Cropland = "workdir/extra_inputs/shp_cropland/shp_cropland.shp"
Forest = "workdir/extra_inputs/shp_forest/shp_forest.shp"
waterbodies = "workdir/extra_inputs/shp_water_bodies/shp_water_bodies.shp"
urban ="workdir/extra_inputs/shp_builtup/shp_builtup.shp"
access ="workdir/extra_inputs/shp_access/shp_access.shp"
## dist by size
S1 = "workdir/industry/S1/S1.shp"
S2 = "workdir/industry/S2/S2.shp"
S3 ="workdir/industry/S3/S3.shp"
## extra inputs
slope =get_rooted("workdir/raster/slope.tif")
slope_4326 = get_rooted("workdir/raster/slope_4326.tif")
GHI = get_rooted("workdir/raster/GHI_kallakurichi_cut.tif")
roads_primary  =  "workdir/extra_inputs/shp_roads_primary/shp_roads_primary.shp"
roads_secondary = "workdir/extra_inputs/shp_roads_secondary/shp_roads_secondary.shp"
railways ="workdir/extra_inputs/shp_railways/shp_railways.shp"
powerlines= "workdir/extra_inputs/shp_powerlines/shp_powerlines.shp"
substation = "workdir/extra_inputs/shp_substations/shp_substations.shp"
population ="workdir/raster/population_10_lon_70_general-v1.5.tif"
slope =get_rooted("workdir/raster/slope.tif")

# %%
## boundary files
shp_district = read_df_UT(district_boundary)
shp_village = read_df_UT(village)
shp_taluk = read_df_UT(taluk)
## module files
shp_indus_theo = read_df_UT(theoritical)
shp_indus_tech = read_df_UT(technical)
shp_indus_high =read_df_UT(high)
shp_indus_med = read_df_UT(med)
## other module high files
shp_water_high = read_df_UT(water_high)
shp_forest_high = read_df_UT(forest_high)
shp_agri_high = read_df_UT(agri_high)
shp_solar_high = read_df_UT(solar_high)
shp_housing_high =read_df_UT(housing_high)
## landcover shape files
shp_unused = read_df_UT(Unused_land)
shp_access = read_df_UT(access)
## dist by size
S1 = read_df_UT(S1)
S2 = read_df_UT(S2)
S3 = read_df_UT(S3)


# %%
shp_unused = add_area_and_class_industry(shp_unused)


# %% [markdown]
# ## Report Tables

# %% [markdown]
# ### Stats

# %%
def summarize_data(df, name):
    summary = df.groupby("area_class")["area_acres"].agg(["sum", "count"]).reset_index()
    summary.columns = ["area_class", f"{name}_sum", f"{name}_count"]
    return summary

data_frames = [shp_unused, shp_indus_theo, shp_indus_tech, shp_indus_high, shp_indus_med]
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
shp_indus_water= find_overlap(shp_indus_tech,"water",shp_water_high)
shp_indus_forest = find_overlap(shp_indus_tech,"forest",shp_forest_high)
shp_indus_agri =find_overlap(shp_indus_tech,"agri",shp_agri_high)
shp_indus_housing =find_overlap(shp_indus_tech,"hsing",shp_housing_high)
shp_indus_solar =find_overlap(shp_indus_tech,"solar",shp_solar_high)
combined = pd.concat([shp_agri_high,shp_forest_high,shp_water_high,shp_solar_high,shp_housing_high])
combined.reset_index(inplace =True,drop =True)
inter = gpd.overlay(shp_indus_tech,combined,how ="intersection",keep_geom_type=True)
inter = inter.dissolve()
inter = add_area_and_class(inter)
data = {
    'indus_forest': shp_indus_forest.oparforest.sum(),
    'indus_water': shp_indus_water.oparwater.sum(),
    'indus-agri':shp_indus_agri.oparagri.sum(),
    'indus-housing': shp_indus_housing.oparhsing.sum(),
    'indus-solar': shp_indus_solar.oparsolar.sum(),
    'Competing use':inter.area_acres.sum(),
}
competing_use = pd.DataFrame(data, index=[0])
competing_use

# %% [markdown]
# ### Distribution by type

# %%
S_list = [S1, S2, S3]

overlap_data = {}

land_use_categories = ['forest', 'water', 'agri', 'hsing', 'solar']
high_dataframes = [shp_forest_high, shp_water_high, shp_agri_high, shp_housing_high, shp_solar_high]


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
        overlap_data[f"S{idx+1}"]["agri"].oparagri.sum(),
        overlap_data[f"S{idx+1}"]["hsing"].oparhsing.sum(),
        overlap_data[f"S{idx+1}"]["solar"].oparsolar.sum()
    ]
    data.append(row_data)

Dist_by_type = pd.DataFrame(data, columns=['Forest', 'Water', 'agri', 'housing', 'solar'], index=['S1', 'S2', 'S3'])
Dist_by_type

# %%
largest_plot = shp_indus_tech.area_acres.max()

# %%
largest_plot

# %%
template = load_workbook(get_rooted('workdir/industry/indus_temp.xlsx'))
ws = template['Sheet1']

# Define the ranges where to insert the dataframes
ranges = ['A2:K5','A10:F10','B13:F15']

# Insert the dataframes into the workbook
for x, r in zip([merged_summary,competing_use,Dist_by_type], ranges):
    start_col, start_row, end_col, end_row = openpyxl.utils.cell.range_boundaries(r)
    for i, row in enumerate(x.values):
        for j, value in enumerate(row):
            ws.cell(row=start_row + i, column=start_col + j, value=value)
ws.cell(row=17, column=2, value=largest_plot)           
template.save(get_in_output('industry/indus_temp_filled.xlsx'))

# %% [markdown]
# ## Top15

# %%
shp_indus_high = shp_indus_high.sort_values(by=["area_acres"],ascending = False)
shp_indus_top15  = shp_indus_high[:15]
shp_indus_top15.reset_index(inplace =True)

# %%
shp_indus_top15

# %%
df = shp_indus_top15
slope =get_rooted("workdir/raster/slope.tif")
top15_indus = calculate_slope(df, slope)

# %%
top15_indus['coords'] = top15_indus['geometry'].apply(lambda x: x.representative_point().coords[:])
top15_indus['coords'] = [coords[0] for coords in top15_indus['coords']]
top15_indus["coords"].tolist()
top15_indus[['lat', 'lon']] = gpd.GeoDataFrame(top15_indus['coords'].tolist(), index=top15_indus.index) 

# %%
overlap = find_overlap(top15_indus,"solar",shp_solar_high)
overlap = find_overlap(overlap,"forest",shp_forest_high)
overlap = find_overlap(overlap,"water",shp_water_high)
overlap = find_overlap(overlap,"agri",shp_agri_high)
overlap = find_overlap(overlap,"hsing",shp_housing_high)

# %%
overlap.columns

# %%
overlap.info()

# %%
overlap.iloc[:,39:54] = overlap.iloc[:,39:54].astype(float)

# %%
overlap.info()

# %%
overlap_top15 = overlap[["lat","lon","area_acres","minDistRd","rdmindist","rdmintype","rdtype","min","max",'op%solar', 'oparsolar',
       'cntsolar', 'op%forest', 'oparforest', 'cntforest', 'op%water',
       'oparwater', 'cntwater', 'op%agri', 'oparagri', 'cntagri', 'op%hsing',
       'oparhsing', 'cnthsing','geometry']]

# %%
overlap_top15.to_file(get_in_output("industry/Top15_final"))

# %%
overlap_top15_excel = overlap[["lat","lon","area_acres","minDistRd","rdmindist","rdmintype","rdtype","min","max",'op%solar', 'oparsolar',
       'cntsolar', 'op%forest', 'oparforest', 'cntforest', 'op%water',
       'oparwater', 'cntwater', 'op%agri', 'oparagri', 'cntagri', 'op%hsing',
       'oparhsing', 'cnthsing']]

# %%
overlap_top15_excel.to_excel(get_in_output("industry/Top15_final.xlsx"))

# %% [markdown]
# # Visuals

# %% [markdown]
# ## Technical Suitability - Technical, theoretical, and no potential lands

# %%
fig8, ax8 = plt.subplots(figsize=(5, 5))

plot_common_features(fig8, ax8)
plot_cities(fig8, ax8)


shp_district.plot(figsize=(5,5),color="none", ax=ax8, linewidth = 0.5, zorder=5)

shp_unused.plot(color="#424242",ax =ax8, label='No Potential',zorder=1)
shp_indus_theo.plot(color="#a95c38",ax =ax8, label='Theoretical Potential',zorder=2)
shp_indus_tech.plot(color="#ffaa75",ax =ax8, label='Technical Potential',zorder=3)


No_P = mpatches.Patch(color='#424242', label='No potential')
Theo_P = mpatches.Patch(color='#a95c38', label='Theoretical potential')
Tech_P = mpatches.Patch(color='#ffaa75', label='Technical potential')
    
plt.legend(handles = [No_P, Theo_P, Tech_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)


print(plt.rcParams['font.family'])


plt.savefig(get_in_output("images/industry/Technical_suitability.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ## Distribution by type

# %%
fig9, ax9 = plt.subplots()

plot_common_features(fig9, ax9)
plot_cities(fig9, ax9)
fig9.subplots_adjust(left=0.05, right=0.8)

shp_district.plot(figsize=(5,5),color="none", ax=ax9, linewidth = 0.5)

S1.plot(color="#4a372d",ax =ax9, label='>0.3 to 20 acres')
S2.plot(color="#d87357",ax =ax9, label='>20 to 100 acres')
S3.plot(color="#ffb9a5",ax =ax9, label='>100 acres')


S1_p = mpatches.Patch(color='#4a372d', label='>0.3 to 20 acres')
S2_p = mpatches.Patch(color='#d87357', label='>20 to 100 acres')
S3_p = mpatches.Patch(color='#ffb9a5', label='>100 acres')
    
plt.legend(handles = [S1_p, S2_p, S3_p], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

# create zoomed inset axes
axins1 = zoomed_inset_axes(ax9, zoom=4, bbox_to_anchor=(1,0.9 ), bbox_transform=ax9.figure.transFigure)

S1.plot(color="#4a372d",ax =axins1, label='>0.3 to 20 acres')
S2.plot(color="#d87357",ax =axins1, label='>20 to 100 acres')
S3.plot(color="#ffb9a5",ax =axins1, label='>100 acres')

# adjust tick parameters and limits of zoomed inset axes
axins1.tick_params(axis='both', colors='grey', labelsize=3)
axins1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins1.set_xlim(79.16, 79.22)
axins1.set_ylim(11.90, 11.96)

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

zoom_patch1 = ConnectionPatch(xyA=(79.19, 11.93), coordsA=ax9.transData,
                             xyB=(79.15, 11.93), coordsB=axins1.transData,
                             arrowstyle="->", shrinkA=5, shrinkB=5,
                             mutation_scale=15, fc="none", ec="grey")
zoom_patch1.set_linestyle('dashed')
zoom_patch1.set_linewidth(0.2)
ax9.add_artist(zoom_patch1)

# create zoomed inset axes
axins2 = zoomed_inset_axes(ax9, zoom=4, bbox_to_anchor=(1,0.65), bbox_transform=ax9.figure.transFigure)

# plot water potentials on zoomed inset axes
S1.plot(color="#4a372d",ax =axins2, label='>0.3 to 20 acres')
S2.plot(color="#d87357",ax =axins2, label='>20 to 100 acres')
S3.plot(color="#ffb9a5",ax =axins2, label='>100 acres')

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
axins2.set_xlim(78.88, 78.94)
axins2.set_ylim(11.80, 11.86)
zoom_patch2 = ConnectionPatch(xyA=(78.91, 11.83), coordsA=ax9.transData,
                             xyB=(78.87, 11.83), coordsB=axins2.transData,
                             arrowstyle="->", shrinkA=5, shrinkB=5,
                             mutation_scale=15, fc="none", ec="grey")
zoom_patch2.set_linestyle('dashed')
zoom_patch2.set_linewidth(0.2)
ax9.add_artist(zoom_patch2)


axins3 = zoomed_inset_axes(ax9, zoom=4, bbox_to_anchor=(1,0.4), bbox_transform=ax9.figure.transFigure)

# plot water potentials on zoomed inset axes
S1.plot(color="#4a373d",ax =axins3, label='>0.3 to 20 acres')
S2.plot(color="#d87357",ax =axins3, label='>20 to 100 acres')
S3.plot(color="#ffb9a5",ax =axins3, label='>100 acres')

# Modify the edge color of the inset axes
axins3.spines['top'].set_edgecolor('grey')
axins3.spines['right'].set_edgecolor('grey')
axins3.spines['bottom'].set_edgecolor('grey')
axins3.spines['left'].set_edgecolor('grey')
# Modify the line thickness of the inset axes spines
axins3.spines['top'].set_linewidth(0.5)
axins3.spines['right'].set_linewidth(0.5)
axins3.spines['bottom'].set_linewidth(0.5)
axins3.spines['left'].set_linewidth(0.5)

# adjust tick parameters and limits of zoomed inset axes
axins3.tick_params(axis='both', colors='grey', labelsize=3)
axins3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins3.set_xlim(78.90, 78.96)
axins3.set_ylim(11.68, 11.74)
zoom_patch3 = ConnectionPatch(xyA=(78.91, 11.71), coordsA=ax9.transData,
                             xyB=(78.89, 11.71), coordsB=axins3.transData,
                             arrowstyle="->", shrinkA=5, shrinkB=5,
                             mutation_scale=15, fc="none", ec="grey")
zoom_patch3.set_linestyle('dashed')
zoom_patch3.set_linewidth(0.2)
ax9.add_artist(zoom_patch3)


plt.savefig(get_in_output("images/industry/Distribution by size.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ## High Potential

# %%
fig10, ax10 = plt.subplots()

plot_common_features(fig10, ax10)
plot_cities(fig10, ax10)
fig10.subplots_adjust(left=0.05, right=0.8)
shp_district.plot(figsize=(5,5),color="none", ax=ax10, linewidth = 0.5, zorder= 1)

shp_indus_tech.plot(color="#1e1e1e",ax =ax10, label='Low Potential')
shp_indus_med.plot(color="#bb584d",ax =ax10, label='Medium Potential')
shp_indus_high.plot(color="#ffb485",ax =ax10, label='High Potential')


Low_P = mpatches.Patch(color='#1e1e1e', label='Low potential')
Med_P = mpatches.Patch(color='#bb584d', label='Medium potential')
High_P = mpatches.Patch(color='#ffb485', label='High potential')
    
plt.legend(handles = [Low_P, Med_P, High_P], loc = 'upper left', bbox_to_anchor=(0.8, 0.2), title = 'Legend\n', fontsize = 5.5, markerscale = 2, title_fontsize = 5.5, framealpha= 0, borderpad = 0.3, handletextpad = 0.5, handlelength = 1.0)

# create zoomed inset axes
# create zoomed inset axes
axins1 = zoomed_inset_axes(ax10, zoom=4, bbox_to_anchor=(1,0.9 ), bbox_transform=ax10.figure.transFigure)


shp_indus_tech.plot(color="#1e1e1e",ax =axins1, label='Low Potential')
shp_indus_med.plot(color="#bb584d",ax =axins1, label='Medium Potential')
shp_indus_high.plot(color="#ffb485",ax =axins1, label='High Potential')

# adjust tick parameters and limits of zoomed inset axes
axins1.tick_params(axis='both', colors='grey', labelsize=3)
axins1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins1.set_xlim(79.16, 79.22)
axins1.set_ylim(11.90, 11.96)

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

zoom_patch1 = ConnectionPatch(xyA=(79.19, 11.93), coordsA=ax10.transData,
                             xyB=(79.15, 11.93), coordsB=axins1.transData,
                             arrowstyle="->", shrinkA=5, shrinkB=5,
                             mutation_scale=15, fc="none", ec="grey")
zoom_patch1.set_linestyle('dashed')
zoom_patch1.set_linewidth(0.2)
ax10.add_artist(zoom_patch1)

# create zoomed inset axes
axins2 = zoomed_inset_axes(ax10, zoom=4, bbox_to_anchor=(1,0.65), bbox_transform=ax10.figure.transFigure)

shp_indus_tech.plot(color="#1e1e1e",ax =axins2, label='Low Potential')
shp_indus_med.plot(color="#bb584d",ax =axins2, label='Medium Potential')
shp_indus_high.plot(color="#ffb485",ax =axins2, label='High Potential')
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
axins2.set_xlim(78.88, 78.94)
axins2.set_ylim(11.80, 11.86)
zoom_patch2 = ConnectionPatch(xyA=(78.91, 11.83), coordsA=ax10.transData,
                             xyB=(78.87, 11.83), coordsB=axins2.transData,
                             arrowstyle="->", shrinkA=5, shrinkB=5,
                             mutation_scale=15, fc="none", ec="grey")
zoom_patch2.set_linestyle('dashed')
zoom_patch2.set_linewidth(0.2)
ax10.add_artist(zoom_patch2)


axins3 = zoomed_inset_axes(ax10, zoom=4, bbox_to_anchor=(1,0.4), bbox_transform=ax10.figure.transFigure)

shp_indus_tech.plot(color="#1e1e1e",ax =axins3, label='Low Potential')
shp_indus_med.plot(color="#bb584d",ax =axins3, label='Medium Potential')
shp_indus_high.plot(color="#ffb485",ax =axins3, label='High Potential')

# Modify the edge color of the inset axes
axins3.spines['top'].set_edgecolor('grey')
axins3.spines['right'].set_edgecolor('grey')
axins3.spines['bottom'].set_edgecolor('grey')
axins3.spines['left'].set_edgecolor('grey')
# Modify the line thickness of the inset axes spines
axins3.spines['top'].set_linewidth(0.5)
axins3.spines['right'].set_linewidth(0.5)
axins3.spines['bottom'].set_linewidth(0.5)
axins3.spines['left'].set_linewidth(0.5)

# adjust tick parameters and limits of zoomed inset axes
axins3.tick_params(axis='both', colors='grey', labelsize=3)
axins3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axins3.set_xlim(78.90, 78.96)
axins3.set_ylim(11.68, 11.74)
zoom_patch3 = ConnectionPatch(xyA=(78.91, 11.71), coordsA=ax10.transData,
                             xyB=(78.89, 11.71), coordsB=axins3.transData,
                             arrowstyle="->", shrinkA=5, shrinkB=5,
                             mutation_scale=15, fc="none", ec="grey")
zoom_patch3.set_linestyle('dashed')
zoom_patch3.set_linewidth(0.2)
ax10.add_artist(zoom_patch3)



plt.savefig(get_in_output("images/industry/High Potential_H_M_L.jpg"),dpi =1500)
plt.show()

# %% [markdown]
# # Settlement analysis

# %%
shp_village

# %%
popdf = pd.DataFrame()
for j in range(len(shp_village)):
    input_shp =  get_rooted('workdir/temp_ind.shp')

    selection = shp_village.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted("workdir/raster/population_10_lon_70_general-v1.5.tif")

        output_raster = get_rooted('workdir/temp_ind.tif')
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
shp_village_tech = find_overlap(shp_village_final,"tech",shp_indus_tech)

# %%
shp_village_tech =  find_overlap(shp_village_tech,"unused",shp_unused)

# %%
shp_village_tech = shp_village_tech.to_crs(4326)

# %%
shp_village_tech.info()

# %%
shp_village_tech = shp_village_tech.iloc[:, [14,15,16,17] + list(range(44, 58))]

# %%
shp_village_tech.info()

# %%
shp_village_tech.iloc[:,12:18] = shp_village_tech.iloc[:,12:18].astype(float)

# %%
shp_village_tech= shp_village_tech[['p_name', 'd_name', 'b_name', 'p_name_rd', 'area_acres', 'area_class'
       , 'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech',
       'opartech', 'cnttech', 'op%unused', 'oparunused', 'cntunused','geometry']]

# %%
shp_village_tech.to_file(get_in_output("industry/settlement"))

# %%
shp_village_tech.columns

# %%
shp_village = shp_village_tech[['p_name', 'd_name', 'b_name', 'p_name_rd', 'area_acres', 'area_class',
       'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech', 'opartech', 'cnttech',
       'op%unused', 'oparunused', 'cntunused']]

# %%
shp_village.to_excel(get_in_output("industry/settlement.xlsx"))

# %%
shp_village_tech.sort_values("p_name",inplace =True)

# %%
shp_village_tech.reset_index(drop =True,inplace =True)

# %%
shp_village_tech.drop("geometry",axis =1,inplace =True)

# %%
shp_village_tech.to_excel(get_in_output("industry/230405_settlement_ordered.xlsx"))

# %% [markdown]
# # Settlement visuals

# %%
settlement = read_df_UT("output/industry/settlement/settlement.shp")


# %%
S3 = LinearSegmentedColormap.from_list('testCmap1', colors=["#f8e2dd", "#e1b5ad", "#cb8f84", "#90534a", "#4c2a26"], N=256)

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



plt.savefig(get_in_output("images/industry/Technical potential area(%).jpg"),dpi =1500)
plt.show()

# %% [markdown]
# ### Taluk analysis

# %%
shp_taluk = gpd.overlay(shp_district,shp_taluk,how="intersection")

# %%
popdf = pd.DataFrame()
for j in range(len(shp_taluk)):
    input_shp =  get_rooted('workdir/temp_ind.shp')

    selection = shp_taluk.geometry[j:j+1]
    if selection.geometry.is_empty.bool():
        rasterarr = []
    else:
        selection.to_file(input_shp)

        input_raster= get_rooted("workdir/raster/population_10_lon_70_general-v1.5.tif")

        output_raster = get_rooted('workdir/temp1.tif')
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
shp_taluk_final = find_overlap(shp_taluk_final,"tech",shp_indus_tech)

# %%
shp_taluk_final['geometry'] = shp_taluk_final['geometry'].apply(lambda x: shape(x).buffer(0).buffer(0.0000000000001))
shp_taluk_finalshp_taluk_final = gpd.GeoDataFrame(shp_taluk_final, geometry='geometry')
shp_taluk_final = shp_taluk_final.loc[shp_taluk_final.is_valid]

# %%
shp_taluk_final =  find_overlap(shp_taluk_final,"unused",shp_unused)

# %%
shp_taluk_final = shp_taluk_final.to_crs(4326)

# %%
shp_taluk_final= shp_taluk_final[['Taluk_name', 'area_acres', 'area_class'
       , 'totpop', 'lat', 'lon', 'TGA(acres)', 'op%tech',
       'opartech', 'cnttech', 'op%unused', 'oparunused', 'cntunused','geometry']]

# %%
shp_taluk_final.totpop.sum()

# %%
# shp_taluk_final.to_file(get_in_output("industry/shp_taluk"))

# %%
shp_taluk_final.drop("geometry",axis =1,inplace =True)

# %%
shp_taluk_final.to_excel(get_in_output("industry/taluk_analysis.xlsx"))

# %%
