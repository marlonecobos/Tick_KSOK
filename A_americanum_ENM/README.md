## Spatiotemporal dynamics of A. americanum questing activity

This repository serves to store R scripts to reproduce analyses and figures presented in the manuscript **"Modeling spatiotemporal dynamics of *Amblyomma americanum* questing activity in the Central Great Plains."**

Authors: Marlon E. Cobos*, Taylor Winters,Ismari Martinez,Yuan Yao,Xiangming Xiao,Anuradha Ghosh,Kellee Sundstrom,Kathryn Duncan,Robert E. Brennan,Susan E. Little,A. Townsend Peterson

Contact: manubio13@gmail.com

<br>

### Citation

Cobos ME, Winters T, Martinez I, Yao Y, Xiao X, Ghosh A, Sundstrom K, Duncan K, Brennan RE, Little SE, Peterson AT (2024) Modeling spatiotemporal dynamics of *Amblyomma americanum* questing activity in the central Great Plains. PLoS ONE 19(10): e0304427. https://doi.org/10.1371/journal.pone.0304427

<br>

### Description

Note: The environmental raster layers used in data preparation analyses, was downloaded and prepared using Google Earth Engine. You can find the GEE script in this <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/General/GEE_Daymet_data_prep.txt" target="_blank">Link</a> 

This repository sub folder contains the following:

- <a href="https://github.com/marlonecobos/Tick_KSOK/tree/main/A_americanum_ENM/Data" target="_blank">Data</a>  
  - Tick records collected in the entire project (focused in *A. americanum*; a CSV file) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/A_americanum_ENM/Data/KS_OK_tickdata_Aa.csv" target="_blank">Link</a> 
- <a href="https://github.com/marlonecobos/Tick_KSOK/tree/main/A_americanum_ENM/Scripts" target="_blank">Scripts</a> 
  - Code to process data before running analyses (an R script) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/A_americanum_ENM/Scripts/Data_preparation_Aa.R" target="_blank">Link</a>
  - Code to run models and obtain predictions using all *A. americanum* data (an R script) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/A_americanum_ENM/Scripts/Risk_modeling_Aa.R" target="_blank">Link</a>
  - Code to run models and obtain predictions using data separated by *A. americanum* life stage (an R script) <a href="https://github.com/marlonecobos/Tick_KSOK/blob/main/A_americanum_ENM/Scripts/Risk_modeling_Aa_life_stages.R" target="_blank">Link</a>
