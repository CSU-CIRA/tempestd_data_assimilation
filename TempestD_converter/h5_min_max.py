# coding: utf-8

# Stock modules
import argparse
import numpy as np
import h5py

'''
"python h5_min_max.py <h5_file_path> <h5_data_set_name>"
Reports the minimum and maximum values for an HDF5 dataset.
'''

#
# Main Program
#
parser = argparse.ArgumentParser()
parser.add_argument(
  'h5_file_path',
  help = 'The path of the HDF5 file'
)
parser.add_argument(
  'h5_data_set_name',
  help = 'The name of the dataset to get the min and max values of'
)
parser_args = parser.parse_args()

h5_path = parser_args.h5_file_path

print('opening', h5_path)
hdf5 = h5py.File(h5_path, 'r')

h5_dataset = hdf5[parser_args.h5_data_set_name][:]

# Insure all missing values are nan
# This is somewhat Tempest-D specific
h5_dataset[h5_dataset < -990] = np.nan

themin = np.nanmin(h5_dataset)
themax = np.nanmax(h5_dataset)

print(f'min: {themin}, max: {themax}')
