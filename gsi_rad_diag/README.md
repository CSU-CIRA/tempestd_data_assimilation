# Read_GSI_Radiance_Diagnostic_File

## Purpose: 
Read GSI generated radiance diagnostic file (ascii)
Perform statistics calculation of observation minus background (O-B) and observation minus analysis (O-A) 
Make plots to show fit-to-observation.

## Description:
The use of this script assumes the following:
1. pre-existing GSI generated radiance diagnostic file in ascii
2. If a binary diagnostic file is available, a fortran-written converter is required.
   A fortran-written diagnostic converter is available in GSI source code (EMC maintained GitHub repository):
   gsi/util/Analysis_Utilities/read_diag/read_diag_rad.f90

## Note:
Currently, the script is written for TEMPEST-D and MHS because the diatnostic file is sensor specific. 
Expanding the script to work fow other sensors is straightforward as long as one is familiar with the GSI source code.

