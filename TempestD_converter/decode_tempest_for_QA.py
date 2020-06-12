# coding: utf-8

# Stock modules
import os
import sys
import math
import argparse
import glob
import re
import numpy as np
import numpy.ma as ma
import datetime
import pytz
import h5py
import logging

# Downloaded specialty modules
import ncepbufr

# Local modules
from tempest_h5_to_bufr import TD_record

'''
"python decode_tempest_for_QA.py bufr_file h5_data_dir"
Decodes the Tempest-D BUFR file and compares the pixel values in every subset
to the corresponding pixel values in the HDF5 file.
'''

# For converting the iwp & lwp values
KG_TO_G = 1000

class TD_HDF5_File:
    '''
    Encapsulates an HDF5 file. Specifically keys it to the time range given in
    the file name.
    '''

    # Free the file data if this many match failures have been seen
    MAX_MATCH_FAILURES = 10

    # TEMPEST_L1_LG_pub_20181208T030000_20181208T090000_v1.41.h5
    REGEX = 'TEMPEST.*\_(\d{8}T\d{6})\_(\d{8}T\d{6})\_.*\.h5'
    FNAME_MATCHER = re.compile(REGEX)

    def __init__(self, filepath):
        self.filepath = filepath

        # Init the file object
        log.info('opening %s', filepath)
        self.file_obj = h5py.File(filepath, 'r')

        # Read the data in only when we need it
        self.h5_data = None

        # Determine the file's time range by parsing its name
        fname_match = self.FNAME_MATCHER.search(filepath)

        dt = datetime.datetime.strptime(
          fname_match.group(1),
          '%Y%m%dT%H%M%S'
        )
        # The pixel records we want to compare these to use timezone aware
        # datetimes so these need to be timezone aware too
        self.dt_start = pytz.utc.localize(dt)

        dt = datetime.datetime.strptime(
          fname_match.group(2),
          '%Y%m%dT%H%M%S'
        )
        self.dt_end = pytz.utc.localize(dt)

        # Count how many times the time range match has failed
        self.matchFailures = 0

    def read_file(self, record):
        '''
        Read the h5 datasets into memory using the dataset keys in the record
        data array
        '''
        self.h5_data = {}
        for dset_name in record.data:
            self.h5_data[dset_name] = self.file_obj[dset_name][:]

    def free_data(self):
        '''
        At least try to free the file data from memory
        '''
        for dset_name, h5_d in self.h5_data.items():
            del h5_d

        # If the above loop doesn't do it, this should
        self.h5_data = None

    def in_range(self, record):
        '''
        Returns True if the record was in file's time range, and False if not.
        '''
        if record.datetime() >= self.dt_start and record.datetime() <= self.dt_end:
            return True
        else:
            self.matchFailures += 1
            if self.matchFailures >= self.MAX_MATCH_FAILURES and self.h5_data:
                # Data in the BUFR files should be pretty contiguous by
                # date/time, so we probably don't need this hdf5 data anymore
                self.free_data()

                # But just in case we do - start over
                self.matchFailures = 0
            return False

    def check_pixel(self, record):
        '''
        We should know that the hdf5 file time range matches the record time,
        so the record's data should be in here and it should match up properly
        '''
        if not self.h5_data:
            self.read_file(record)

        row = record.sangle_idx
        col = record.sline_idx
        for dset_name, h5_d in self.h5_data.items():
            pixel_val = record.data[dset_name]
            h5_val = h5_d[row, col]
            #log.debug('pixel_val: %s, h5_val: %s', pixel_val, h5_val)
            if pixel_val == 10E10:
                # h5 val should be nan
                if np.isnan(h5_val) or h5_val < -990:
                    result = True
                else:
                    result = False
            elif h5_val < -300:
                    result = False
            elif np.issubdtype(h5_d.dtype, np.integer):
                # For h5 integer datasets check for equality
                result = (int(pixel_val) == h5_val)
            elif dset_name == '/iwp' or dset_name == '/lwp':
                # For the water path datasets check for closeness within 1000th
                result = math.isclose(pixel_val, h5_val,
                                      abs_tol=0.001)
            else:
                # For the other float datasets check for closeness within
                # 10000th
                result = math.isclose(pixel_val, h5_val,
                                      abs_tol=0.0001)

            if not result:
                log.warning('BUFR & HDF5 values don\'t match:')
                log.warning('HDF5 file: %s, row: %s, col: %s', self.filepath, row, col)
                log.warning('dataset: %s', dset_name)
                log.warning('bufr val: %s, h5 val: %s\n', pixel_val, h5_val)
                #sys.exit(1)

class PixelChecker:
    '''
    Checks the accuracy of the pixel data from a BUFR file against one or more
    input HDF5 files, to QA the conversion of the HDF5 files to BUFR
    '''
    def __init__(self, h5_dir):
        # Open and initialize all the HDF5 files
        h5_glob = os.path.join(h5_dir, '*.h5')
        h5_paths = glob.glob(h5_glob)

        self.filelist = []
        for h5_path in h5_paths:
            self.filelist.append(TD_HDF5_File(h5_path))

        # Sort by starting datetime
        self.filelist.sort(key = lambda x: x.dt_start)

    def check(self, pixel):
        '''
        Makes sure the original hdf5 file pixel data is correctly matched by
        the BUFR pixel data contained in the pixel record
        '''
        was_range_matched = False
        was_verified = False
        for h5_file in self.filelist:
            if h5_file.in_range(pixel):
                was_range_matched = True
                if h5_file.check_pixel(pixel):
                    was_verified = True

        if not was_range_matched:
            literal = (
              'No HDF5 file\'s time range matched this pixel\'s datetime:'
            )
            log.warning('%s\n%s', literal, str(pixel.__dict__))

        return was_verified

#
# Main Program
#
# Setup logging. Note that the output directly from BUFRLIB goes to stdout,
# so the logging output has to do the same for both to go to the same
# place.
log = logging.getLogger(__name__)
logging.basicConfig(
  format='%(asctime)s %(levelname)-8s%(name)s: %(message)s',
  level=logging.INFO
)

parser = argparse.ArgumentParser()
parser.add_argument(
  'bufr_file',
  help='The BUFR file to check'
)
parser.add_argument(
  'input_hdf5_dir',
  help='The directory that contains the Tempest-D HDF5 files to check against'
)
#parser.add_argument(
#  'increment',
#  help='The size of the steps through the BUFR file pixels'
#)
parser_args = parser.parse_args()

# Open the HDF5 files and get them set up for pixel data checking
checker = PixelChecker(parser_args.input_hdf5_dir)

bufr = ncepbufr.open(parser_args.bufr_file)
while bufr.advance() == 0:
    bufr_header = '{:10d}{:6d}{:^10}'.format(
      bufr.msg_date, bufr.msg_counter, bufr.msg_type)
    log.info(bufr_header)

    while bufr.load_subset() == 0:
        pixel = TD_record()
        pixel.data = {}

        #scalarstr1 = 'SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH CHSQ'
        # Squeeze out any extra singleton dimensions and fill in any masked
        # values with the fill value - 10E10
        hdr = bufr.read_subset(pixel.scalarstr1).squeeze().filled()
        #log.debug('hdr: %s', hdr)

        # Put the time info in both the pixel time members and the data dict so
        # that both getting the object's datetime and looking up the data in
        # the hdf5 files is easy
        pixel.year = pixel.data['/year'] = int(hdr[1])
        pixel.month = pixel.data['/month'] = int(hdr[2])
        pixel.day = pixel.data['/day'] = int(hdr[3])
        pixel.hour = pixel.data['/hour'] = int(hdr[4])
        pixel.min = pixel.data['/minute'] = int(hdr[5])
        pixel.sec = pixel.data['/second'] = int(hdr[6])
        #log.debug('hdr[7]: %s, %s', hdr[7], type(hdr[7]))
        pixel.data['/pixel latitude'] = hdr[7]
        pixel.data['/pixel longitude'] = hdr[8]
        pixel.data['/chi'] = hdr[9]

        #scalarstr2 = 'CLAVR SAZA BEARAZ SOZA SOLAZI SANG FOVN SLNM'
        hdr = bufr.read_subset(pixel.scalarstr2).squeeze().filled()
        pixel.data['/converge'] = hdr[0]
        pixel.data['/zenith_angle'] = hdr[1]
        pixel.data['/scan_angle'] = hdr[5]

        pixel.sangle_idx = int(hdr[6])
        pixel.sline_idx = int(hdr[7])

        #ilwpstr = 'COLN ILWP'
        obs = bufr.read_subset(pixel.ilwpstr, rep=True).squeeze().filled()
        # Convert back to g m^-1
        obs[1, :][obs[1, :] < 9E9] *= KG_TO_G
        #log.debug('ilwp obs: %s', obs)
        pixel.data['/iwp'] = obs[1, 0]
        pixel.data['/lwp'] = obs[1, 1]

        #tmbrstr = 'CHNM TMBR'
        obs = bufr.read_subset(pixel.tmbrstr, rep=True).squeeze().filled()
        #log.debug('tmbr obs: %s', obs)
        pixel.data['/Tb 89 GHz'] = obs[1, 0]
        pixel.data['/Tb 165 GHz'] = obs[1, 1]
        pixel.data['/Tb 176 GHz'] = obs[1, 2]
        pixel.data['/Tb 180 GHz'] = obs[1, 3]
        pixel.data['/Tb 182 GHz'] = obs[1, 4]

        checker.check(pixel)

bufr.close()
