# Repository_For_Data_Scripts_Code_Used_In_The_TEMPEST-D_Data_Assimilation_Workflow_Manuscript

## Description:
This repository contains data and scripts used for the production of the Wu et al. (2020) JTECH paper.

In addition, access to online code repositories as well as publicly avaialbe dataset that were used in the production of the Wu et al. (2020) is also documented in this README. 
 

## Data:

  - data/tempestd_retrievals/
    One sampe file of the TEMPEST-D retrieval from Schulte et al. (2020) (netCDF)

Data that is publicly available, therefore, not included in this repository:

  - TEMPEST-D data in HDF5 format: 
    https://tempest.colostate.edu/data

  - GDAS analysis output (operational production: nemsio)
    NOAA NCEI GDAS-Daily Tar Files: https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-data-assimilation-system-gdas

  - MHS radiance dataset (BUFR)
    NOAA NCEI GDAS-Daily BUFR Files: https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-data-assimilation-system-gdas

  - GMI dataset in HDF5 format:
    NASA GES DISC: https://disc.gsfc.nasa.gov/datasets/GPM_1CGPMGMI_05/summary 

  - FV3GFS output files (various formats):
    NOAA NCEP NOMADS Data Archive: https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod
    (gfs.tNNz.atmanl.nemsio, gfs.tNNz.atmfXXX.nemsio, and gfs.tNNz.XXX.tm00.bufr_d were used)
 
## Scripts: 

  - gsi_rad_diag/ 
    read GSI radiance diagnostic file (run_read_diag_rad.sh)
    and produce O-B/O-A plots (readgsi_diagrad_multiple.py)

  - read_fv3gfs_output/
    extract FV3GFS variables from nemsio output file: forecast, analysis, and production analysis 
    (run_read_forecast.sh, run_read_analysis.sh, and run_read_production.sh) 
    and compute analysis and forecast bias relative to the operational GDAS analysis 
    (analysis_bias_comparison.py and forecast_bias_comparison.py) 
 
  - anomaly_correlation/
    compute anomaly correlation coefficient (compute_anomaly_correlation.py and compile_acc.py) 
    and the script that read and write NCEP/NCAR Reanalysis dataset (read_clima.py)

  - create_synthetic_gmi/
    to create synthetic GMI brightness temperatures (Tbs) using the output from FV3GFS and the CRTM

  - tempestd_gmi_analysis/
    to perform comparison between observed and synthetic GMI Tbs generated from TEMPEST-D assimilation experiments
 
 
## Code:

  - TempestD_converter/ 
    the NCEP-BUFR encoder to create TEMPEST-D BUFR data files

  
  - py-ncepbufr-docker/
    to create a Docker image that has the environment to run py-ncepbufr script for TempestD_converter
    
  - tempestd_preproc/
    the script that performs quality assurance on the TEMPEST-D raw data to avoid unphysical info flowing into data assimilation


Code maintained by other version control manager, therefore, not included in this repository:

  - FV3GFS global-workflow on GitHub Repository: 
    https://github.com/NOAA-EMC/global-workflow

  - NCEP/EMC GFS forecast verification package on Fanglin Yang’s GitHub repository: 
    https://github.com/yangfanglin/gfs_verif

  - GSI tempestd_dev branch on NOAA Virtual Lab Git Repository: 
    https://vlab.ncep.noaa.gov/redmine/projects/comgsi/repository?utf8=%E2%9C%93&rev=tempestd_dev
