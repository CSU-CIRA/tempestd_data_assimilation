# TempestD_converter

## Purpose:
Converts Preprocessed Tempest-D HDF5 files to BUFR formatted files.

## Description:
The converter consists of a containerized set of processes that do the
conversions. Everything runs in a single Docker container. Python
Multiprocessing is used to convert multiple files in parallel.

This depends on the SmallSat py-ncepbufr-docker project in that the base Docker
image used in the Dockerfile is the py-ncepbufr-docker image.

## Installation:
Install the tempd-bufr-converter:latest image in the host computer's Docker
environment.

For those at CIRA this is best done by pulling the image from CIRA's GitLab
Docker registry on bear. Those not at CIRA will need to build the image from
this project. Using the buildDocker script is recommended.

## Usage:
### Development:
The runDocker script, which encapsulates `docker-compose -f docker-compose.dev.yml up`,
will build the tempd-bufr-converter:latest image if needed and then start up a
tempestd\_convert\_1 container based on that image to run the converter code.
The code will start running when the container is created. It will convert
multiple HDF5 files in parallel.

### Operations:
On an operational machine (ssdate at CIRA), running the converter with just
`docker-compose up` is recommended. The docker-compose.yml file is configured
for the operational environment. Again, it will start up a tempestd\_convert\_1
container, which will immediately start the converter code.

### Configuration:
Environment variables are used for the configuration settings, and are set in
the docker-compose yaml files. Edit the docker-compose.dev.yml or
docker-compose.yml file to change the settings. If the input or output
directories need to change, then the "volumes:" setting may also need to
change.

### Workflow:
The converter requires that three directories be associated with the
preprocessed Tempest-D HDF5 files. These are specified in the docker-compose
configuration using the three environment variables: `H5_DATA_DIR`,
`H5_PROCESSED_DIR` and `H5_BAD_DIR`. It is recommended that they be put under
the same parent directory. A `BUFR_OUTPUT_DIR` must also be specified. Here is
the example from the docker-compose.yml file:
```
- H5_DATA_DIR=/data/preprocessor_output/output
- H5_PROCESSED_DIR=/data/preprocessor_output/used
- H5_BAD_DIR=/data/preprocessor_output/failed
- BUFR_OUTPUT_DIR=/data/bufr_converter_output/output
```

It is recommended that all of the H5 directories be under the same parent
directory as shown, and that the `BUFR_OUTPUT_DIR` be part of an equivalent
directory structure with `used/` and `failed/` subdirectories under
`bufr_converter_output/` in this example. Subsequent processing that uses the
BUFR files can then use these directories in the same way the converter uses
the `preprocessor_output/` subdirectories.

The converter reads the HDF5 files in `preprocessor_output/output/` and writes
the BUFR files to `bufr_converter_output/output/`. When a conversion is
complete the HDF5 file is moved to the `used/` subdirectory. If a conversion
fails the HDF5 file is moved to `failed/`. This keeps track of which files have
been processed and which have not.

The converter continually reads `preprocessor_output/output/` looking for new
HDF5 files. It only exits when the tempestd\_convert\_1 container is manually
terminated. It expects any new HDF5 file names to end with an extra .part
extension while they are being written so that they will only be read after
writing is complete. The converter also adds .part to the end of the BUFR file
names as they are being written for the same reason.
