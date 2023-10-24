# Codes to reconstruct Sea Surface Salinity observations made by the Soil Moisture Active Passive (SMAP) satellite in the northwest Atlantic Ocean

This repository contains the codes used to conduct the study of Barceló-Llull, B., Drushka, K., &amp; Gaube, P. (2021). Lagrangian reconstruction to extract small-scale salinity variability from SMAP observations. Journal of Geophysical Research: Oceans, 126, e2020JC016477. https://doi.org/10.1029/2020JC016477

## Description of the repository:

- `simulation_notebooks/`:  contains the codes to perform the backward simulations.
  - `sim_Oleander_backward_altimetry_Bigger_Domain_dt.ipynb`: Oleander simulations.
  - `sim_weekly_backward_altimetry_Bigger_Domain_nrt_dates_corrected.ipynb`: weekly simulations for the dates in which only near-real-time altimetry data was available.
  - `sim_weekly_backward_altimetry_Bigger_Domain_dt.ipynb`: weekly simulations for the dates in which delayed-time altimetry data was available.