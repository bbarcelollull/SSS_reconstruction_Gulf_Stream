# Codes to reconstruct Sea Surface Salinity observations made by the Soil Moisture Active Passive (SMAP) satellite in the northwest Atlantic Ocean

This repository contains the codes used to conduct the study of Barceló-Llull, B., Drushka, K., &amp; Gaube, P. (2021). Lagrangian reconstruction to extract small-scale salinity variability from SMAP observations. Journal of Geophysical Research: Oceans, 126, e2020JC016477. https://doi.org/10.1029/2020JC016477

## Description of the repository:

- `simulation_notebooks/`:  contains the codes to perform the Lagrangian simulations.
  
  - `sim_Oleander_backward_altimetry_Bigger_Domain_dt.ipynb`: Oleander simulations.
  - `sim_weekly_backward_altimetry_Bigger_Domain_nrt_dates_corrected.ipynb`: weekly simulations for the dates in which only near-real-time altimetry data was available.
  - `sim_weekly_backward_altimetry_Bigger_Domain_dt.ipynb`: weekly simulations for the dates in which delayed-time altimetry data was available.
  - `sim_Oleander_backward_Oscar_check_bigger_domain_1day.ipynb`: test of a simulation using Oscar currents.
    
- `analysis_Oleander_simulations/`:  codes for the analysis of the Oleander set of simulations.
  
  - `sim_Oleander_BiggerDomain_all_step1_tag_SSS_to_particles_final.py`: Tag SSS from SMAP to each particle at T0 of the backward simulation.
  - `sim_Oleander_BiggerDomain_all_step2_comp_with_TSG_save_and_plot.py`: Save advected and SMAP data interpolated onto 
TSG data points for each simulation to detect fronts.
  - `sim_Oleander_BiggerDomain_all_step4_method_detect_fronts_save.py`: Code to detect and quantify fronts. 
  - `sim_Oleander_BiggerDomain_all_step5_plot_fronts_figures_examples_figurePAPER_Fig3.py`: create Figure 3.
  - `sim_Oleander_BiggerDomain_all_step6_plot_fronts_check_intensity.py`: check relative intensity
  - `sim_Oleander_BiggerDomain_all_step6_plot_fronts_statistics_abs_int_figurePAPER_Fig4.py`: create Figure 4.
  - `sim_Oleander_BiggerDomain_all_step6_plot_fronts_statistics_abs_int_figurePAPER_rev_Fig5.py`: create Figure 5.

- `analysis_weekly_simulations/`:  codes for the analysis of the weekly set of simulations.
  
  - ``:
  - ``:
  - ``:
  - ``:
  - ``:
  - ``:

- `analysis_test_Oscar_currents/`:  codes for the analysis of the simulation done with Oscar currents (it is a test).

  - `sim_Oleander_BiggerDomain_1day_Osc_tag_SSS_to_particles.py`: Tag SSS from SMAP to each particle at T0 of the backward simulation.
  - `sim_Oleander_BiggerDomain_1day_Osc_comp_with_TSG.py`: Code to open the simulation file with SSS tagged and compare with TSG data. 

- `generate_additional_figures/`:  codes to create figures 1 and 2 of the paper.
  
  - `Fig_PAPER_1_sim_Oleander_BiggerDomain_1day_alt_comp_with_TSG_resampling_rev1_save.py`: resample and smooth S_adv and save data in a file
  - `Fig_PAPER_1_sim_Oleander_BiggerDomain_1day_alt_comp_with_TSG_resampling_rev1_plot.py`: make Figure 1
  - `Fig_PAPER_2_Oleander_number_trans_STATISTICS_Bigger_Domain_QCed.py`: make Figure 2
 
In each folder are included the toolboxes used for the respective codes: `sim_toolbox.py`, `sss_toolbox.py` and `Tools_OI.py`.

