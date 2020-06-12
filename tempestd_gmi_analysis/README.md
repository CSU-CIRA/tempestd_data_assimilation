## Overview

Analyze brightness temperatures from TEMPEST-D assimilation
experiments and compare to synthetic GMI brightness temperatures
created by the CRTM.

## Code Components

**tempestd_assim_vis.py**: Plots histograms (line plots), difference
maps, and mean global brightness temperatures for each forecast
hour with standard deviation

**histogram_generation.py**: Creates histograms of brightness
temperature bin counts for each experiment or observations for
a given case

## Run instructions
### Environment setup
All components of this code are intended to run in an environment
described by the Dockerfile in this repository. To build the
container image use `docker run -it --name tempestd_postproc`
in the top level directory. To create and run a container from
the image, use `docker run -it --name tempestd_postproc_cont
tempestd_postproc`, adding -v flags before the image name as
necessary to provide any data mounts needed for the code
to be able to access the necessary input data as well as the
code itself.

If the container is stopped (you can check by running `docker ps`), 
start it by running `docker start tempestd_postproc_cont` and
then enter it using `docker exec -it tempestd_postproc_cont /bin/bash`. 
This will land you in the `/` directory.

The code will be available in the container at the location
you mount.

### Creating histograms
- Get help: 
```
python histogram_generation.py -h
```
- Create histogram:
1. Update paths, case, and options as necessary
at the top-level definitions in histogram_generation.py
2. 
```
python histogram_generation.py 

```
- Create plots:
- Get help: 
```
python tempestd_assim_vis.py -h
```
1. Update paths, case, and options as necessary
at the top-level definitions in tempestd_assim_vis.py
2.
```
python tempestd_assim_vis.py 
    -f function_to_run
```
