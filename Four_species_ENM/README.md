## Spatiotemporal dynamics of questing activity by four tick species

This repository serves to store R scripts to reproduce analyses and figures presented in the manuscript **"Spatiotemporal dynamic models of questing activity by four tick species in the Central Great Plains."**

Authors: Marlon E. Cobos, Joanna L. Corimanya, Abigail C. Perkins, Zenia Ruiz-Utrilla, Claudia Nuñez-Penichet, Eric Ng’eno, A. Townsend Peterson*, Kathryn T. Duncan

Contact: town@ku.edu

<br>

### Citation

In process...

<br>

### Description

Note: The environmental raster layers used in data preparation analyses, was downloaded and prepared using Google Earth Engine. You can find the GEE script at this <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/General/GEE_Daymet_data_prep.txt" target="_blank">link</a> 

This repository sub folder contains the following:

- <a href="https://github.com/marlonecobos/Tick_KSOK/tree/main/Four_species_ENM/Data" target="_blank">Data</a> (created using <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/General/KS_OK_tickdata.csv" target="_blank">all data</a> and <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/A_americanum_ENM/Scripts/Data_preparation_Aa.R" target="_blank">data preparation script)</a>
  - Processed tick record data and associated environmental variable values for each tick species (a CSV file)
    - *Amblyomma maculatum* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/am_presence_absence_env.csv" target="_blank">Link</a>
    - *Dermacentor albipictus* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/da_presence_absence_env.csv" target="_blank">Link</a>
    - *Dermacentor variabilis* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/dv_presence_absence_env.csv" target="_blank">Link</a>
    - *Ixodes scapularis* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/is_presence_absence_env.csv" target="_blank">Link</a>
  - Processed tick record data and associated principal component variable values for each tick species (a CSV file)
    - *A. maculatum* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/am_presence_absence_pcs.csv" target="_blank">Link</a>
    - *D. albipictus* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/da_presence_absence_pcs.csv" target="_blank">Link</a>
    - *D. variabilis* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/dv_presence_absence_pcs.csv" target="_blank">Link</a>
    - *I. scapularis* <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Data/is_presence_absence_pcs.csv" target="_blank">Link</a>
- <a href="https://github.com/marlonecobos/Tick_KSOK/tree/main/Four_species_ENM/Scripts" target="_blank">Scripts</a> 
  - Code to run models and obtain predictions using all *A. maculatum* data (an R script) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Scripts/Risk_modeling_Am.R" target="_blank">Link</a>
  - Code to run models and obtain predictions using all *D. albipictus* data (an R script) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Scripts/Risk_modeling_Da.R" target="_blank">Link</a>
  - Code to run models and obtain predictions using all *D. variabilis* data (an R script) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Scripts/Risk_modeling_Dv.R" target="_blank">Link</a>
  - Code to run models and obtain predictions using all *I. scapularis* data (an R script) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/Four_species_ENM/Scripts/Risk_modeling_Is.R" target="_blank">Link</a>
