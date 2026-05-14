# 🗺️ Kallakurichi Land Suitability Assessment

> **Multi-criteria land suitability mapping for Kallakurichi district, Tamil Nadu — identifying optimal land use opportunities across forestation, agriculture, water harvesting, housing, industry, and solar energy.**

[![Python](https://img.shields.io/badge/Python-3776AB?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![Jupyter](https://img.shields.io/badge/Jupyter-F37626?style=flat-square&logo=jupyter&logoColor=white)](https://jupyter.org)
[![Sentinel-2](https://img.shields.io/badge/Sentinel--2-003DA5?style=flat-square)](https://sentinel.esa.int)
[![GeoPandas](https://img.shields.io/badge/GeoPandas-139C5A?style=flat-square)](https://geopandas.org)
[![MIT License](https://img.shields.io/badge/License-MIT-green?style=flat-square)](LICENSE)

---

## 📌 What Is This?

Unused and fallow lands represent significant untapped potential for sustainable development — but their optimal use varies by terrain, ecology, infrastructure, and climate. This project maps and scores unused land across Kallakurichi district for **six competing land use categories** using satellite-derived LULC classification, change detection, and multi-criteria spatial analysis.

The study was commissioned by the **Tamil Nadu State Planning Commission** (under Auroville Consulting) and outputs were delivered to state and district government departments for evidence-based spatial planning.

---

## 🌟 Key Results

| Metric | Value |
|---|---|
| **Study area** | Kallakurichi district, Tamil Nadu |
| **Analysis period** | 2019 – 2023 (5-year LULC change detection) |
| **Satellite imagery** | Sentinel-2 + Landsat time series |
| **Land use categories assessed** | 6 (forestation, agriculture, water harvesting, housing, industry, solar) |
| **Ground truthing** | Field validation across district |
| **Deliverables** | Interactive online map + detailed district report |
| **Client** | Tamil Nadu State Planning Commission (Government of Tamil Nadu) |

**Policy impact:** Results were presented to relevant state and local departments and informed district-level spatial planning decisions. An interactive online map was created for ongoing use by government stakeholders.

---

## 🎯 Project Objectives

1. **Identify unused and fallow lands** across Kallakurichi district using satellite-based LULC classification
2. **Detect land use change** (2019–2023) to understand trends in land availability
3. **Assess multi-criteria suitability** of identified lands for 6 development and conservation categories
4. **Ground truth** findings through field validation
5. **Disseminate results** to government departments via interactive maps and detailed reports

---

## 🌱 Six Land Use Categories

| Category | Key Criteria |
|---|---|
| **Forestation** | Ecological connectivity, slope, soil type, proximity to existing forest |
| **Agriculture** | Soil quality, water availability, accessibility, existing agricultural patterns |
| **Water Harvesting** | Watershed position, DEM (drainage accumulation), soil permeability |
| **Housing** | Proximity to roads and settlements, slope, flood risk, infrastructure |
| **Industrial Development** | Road/rail access, flat terrain, distance from sensitive ecological zones |
| **Ground-mounted Solar** | Solar irradiance, slope (<3%), land availability, grid proximity |

---

## 🔄 Methodology

```
1. Sentinel-2 + Landsat imagery acquisition (2019 & 2023)
       ↓
2. LULC classification (supervised, multi-class)
       ↓
3. LULC change detection (2019 → 2023)
       ↓
4. Unused / fallow land identification
       ↓
5. Multi-criteria suitability scoring per land use category
   (weighted overlay using spatial layers: DEM, soil, hydrology,
    road networks, ecological constraints)
       ↓
6. Competing use analysis (overlapping suitability zones)
       ↓
7. Ground truthing and field validation
       ↓
8. Interactive map + district-level report generation
```

---

## 🗂️ Repository Structure

```
Kallakurichi/
├── kallakurichi/
│   ├── Kallakuruchi_analysis.ipynb       # Master LULC analysis notebook
│   ├── Kallakuruchi_analysis.py
│   ├── Kallakurichi_competinguse.ipynb   # Competing land use overlay analysis
│   ├── Kallakurichi_competinguse.py
│   ├── Kallakuruchi_forest.ipynb         # Forestation suitability
│   ├── kallakurichi_agri.ipynb           # Agriculture suitability
│   ├── kallakurichi_housing.ipynb        # Housing suitability
│   ├── kallakurichi_industry.ipynb       # Industrial suitability
│   ├── kallakurichi_solar.ipynb          # Solar energy suitability
│   ├── kallakurichi_water.ipynb          # Water harvesting suitability
│   ├── kallakuruchi_input.ipynb          # Input data preparation
│   └── Input_template.ipynb              # Reusable analysis template
├── requirements.txt
└── README.md
```

Each suitability category has a **dedicated notebook** with its own criteria, weights, and spatial outputs — enabling modular updates when policy priorities change.

---

## 🛠️ Setup

```bash
# Create environment
conda create --name lila python=3.10
conda activate lila

# Install dependencies
conda install jupyter nbconvert
conda install --file requirements.txt -c conda-forge
pip install dask-geopandas geopandas datashader
```

---

## ▶️ Usage

Open the relevant notebook for the land use category you want to analyse:

```bash
jupyter notebook kallakurichi/Kallakuruchi_analysis.ipynb   # Master analysis
jupyter notebook kallakurichi/kallakurichi_solar.ipynb      # Solar suitability
jupyter notebook kallakurichi/Kallakurichi_competinguse.ipynb  # Competing use overlay
```

Inputs required per notebook:
- District boundary shapefile (EPSG:32644)
- Sentinel-2 / Landsat imagery stack (2019, 2023)
- Ancillary spatial layers: DEM, road network, soil map, hydrological data

---

## 📤 Key Outputs

| Output | Description |
|---|---|
| LULC classification maps | Land cover for 2019 and 2023 |
| Change detection map | Transitions between LULC classes (2019→2023) |
| Unused land inventory | Mapped fallow/unused parcels with area statistics |
| Suitability maps (x6) | Scored land parcels per development category |
| Competing use map | Spatial overlay showing conflicting suitability zones |
| Interactive online map | Web-deliverable for government stakeholders |
| District report | Detailed PDF report with maps, charts, and policy recommendations |

---

## 🌍 Applications

- District-level spatial planning and land governance
- Solar energy potential assessment for distributed generation
- Forestation and ecological restoration planning
- Agricultural land identification and rural development
- Evidence-based policy for climate adaptation at district scale

---

## 👤 Author

**Athithiyan M R** — Geospatial Data Scientist | Remote Sensing | Climate Analytics

This work was carried out at **Auroville Consulting** in collaboration with the **Tamil Nadu State Planning Commission**.

[![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-0A66C2?style=flat-square&logo=linkedin)](https://www.linkedin.com/in/athithiyan-m-r-/)
[![GitHub](https://img.shields.io/badge/GitHub-Athithiyanmr-181717?style=flat-square&logo=github)](https://github.com/Athithiyanmr)

---

## 📜 License

MIT License © 2026 Athithiyan M R
