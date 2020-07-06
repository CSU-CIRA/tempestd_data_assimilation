## Overview

Create synthetic GPM GMI brightness temperatures using the Community
Radiative Transfer Model (CRTM).

## Run instructions
### Environment setup
The main driver script is `crtm_scripts/crtm_gmi_control.sh`. Within the
script, several variables and paths need to be changed prior to execution.
`count` and `count_max` should be the first and last forecast hour to
process. `atm_prefix` and `sfc_prefix` should reflect the nemsio output
for the forecast run. `begin_yyyymmdd` and `begin_hour` refer to the
cycle start day and time. The `top_directory` should be the one that
houses each of the subdirectories contained in this repository and the
`gmi_orbit` directory should point to the location of the binary-format
geolocation files for the GMI orbit that most closely corresponds with
the forecast day and time. `nemsio_source` is the directory that contains
output from the cycled experiment. The `case` is a subdirectory that must
exist within the `crtm_interface`, `nemsio_source`, and `nemsio_directory`
directories. The directory described by `nemsio_directory` must also
be defined in `gfdl_microphysics/gfdl_diams.f90` and
`crtm_interface/Get_CRTM.f90`. In addition to the path to the
`nemsio_directory`, will need to have paths pointing to the CRTM folder
where the coefficient file is located and needs to match the path defined
in `crtm_interface/Makefile`.

Updated executables need to be built in `gfdl_microphysics/` and
`crtm_interface/` using `gfdl_microphysics/Make.gfdl_diams` and
`crtm_interface/Makefile`, respectively, prior to running
`crtm_scripts/crtm_gmi_control.sh`

`crtm_scripts/crtm_gmi_control.sh` can be invoked from the command line
with no arguments or by updating and using `crtm_scripts/submit.sh` to
submit jobs on a Slurm cluster.
