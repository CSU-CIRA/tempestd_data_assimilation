# py-ncepbufr-docker

Main program file: buildDocker

Other project files: Dockerfile, README.md, runDocker, buildf

CIRA project dependencies: None.

### Purpose:
To create a Docker image that has the GitHub JCSDA/py-ncepbufr project built
and installed, and that can be used to create a container that will run
Python/py-ncepbufr scripts in a host directory.

### Description:
The buildDocker script builds the image using the Dockerfile and tags it
pyncepbufr:latest. The build downloads the GitHub JCSDA/py-ncepbufr project and
its dependencies, runs the project's setup.py script to build the NCEP BUFR
library from the code included in the project and gets it installed under the
miniconda that was installed as one of the dependencies.

This results in a tagged Docker image that has py-ncepbufr available in its
Python.

It also includes a gfortran build environment and a libbufr.a, which can be
used to build and test example NCEPBUFR FORTRAN programs.

### Installation:
Install the pyncepbufr image in the host computer's Docker environment.

For those at CIRA this is best done by pulling the image from CIRA's GitLab
Docker registry on bear.

For those not at CIRA build it using this project. Using the buildDocker script
is recommended.

The included runDocker script will create a container that can be used to test
the build. It will run bash in the container and attach it to the invoking
terminal.

Optionally cd to /root/py-ncepbufr-master/test and run the test python
scripts there.

Here are the pull instructions for those at CIRA:

```
docker pull bear.cira.colostate.edu:4567/small-sat/py-ncepbufr-docker/pyncepbufr:latest
```

Then adding the original tag is recommended - the tag the image would get if it
were built it on the local host machine.

```
docker tag bear.cira.colostate.edu:4567/small-sat/py-ncepbufr-docker/pyncepbufr:latest pyncepbufr:latest
```

### Usage:
Create a Docker container that mounts a local Python/py-ncepbufr project
directory using this image. An example runDocker script for this would be:

```
#!/bin/bash

# Get the CWD this script is run from
currWD=$(pwd)
target=/$(basename $currWD)

docker run -t -i --rm -w $target \
  --network=host \
    --mount type=bind,source=$currWD,target=$target \
      pyncepbufr:latest
```

Run this in the project directory.

The invoking terminal will now be running bash in the container as root. Cd to
/<basename project_dir\> and run the project's Python/py-ncepbufr scripts.

Optionally build and run NCEP BUFR f90 programs using the buildf script. This
script is already in PATH. It can only handle single source file programs. Look
at it in /usr/local/bin/buildf as a starting point for more complicated builds.

Give it the f90 program name without the .f90 extension. Example:
```
buildf bufr_encode_sample
```
Builds bufr\_encode\_sample.exe from bufr\_encode\_sample.f90.

### JCSDA/py-ncepbufr Documentation:

Clone or ZIP download the JCSDA/py-ncepbufr project from
https://github.com/JCSDA/py-ncepbufr.

In a browser load:
<install_parent_dir\>/py-ncepbufr/docs/ncepbufr/index.html.
