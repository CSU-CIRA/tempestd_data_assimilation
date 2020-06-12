# coding: utf-8

# Stock modules
import os
import sys
import re
import datetime
import numpy as np
import h5py
import logging
import shutil
import traceback

# Downloaded specialty modules
import ncepbufr

'''
Call function:
"h5_to_bufr(h5_path)" from another Python program after importing in this file.
Converts one of Louie's preprocessed HDF5 Tempest-D data files to BUFR.
'''

# Concerns:
# Cannot set the message receipt_time, but GSI doesn't require the receipt_time
#
# The BUFR message size limit of 10000 bytes cannot be changed with py-netcdf,
# but again, this is not an issue for GSI.

# For converting the iwp & lwp values
G_TO_KG = 10**-3

# For catching bufr errors
class BufrError(Exception):
    def __init__(self, bufr_path_part):
        self.bufr_path_part = bufr_path_part
        Exception.__init__(self)

class TD_record:
    '''
    Encapsulates the input record data. The record variables will get filled
    in for the objects when the hdf5 files are read.
    '''
    scalarstr1 = 'SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH CHSQ'
    scalarstr2 = 'CLAVR SAZA BEARAZ SOZA SOLAZI SANG FOVN SLNM'

    # "ILWPSEQ"2 
    ilwpstr = 'COLN ILWP'
    ilwptup = (1, 2)

    # "BRIT"5
    tmbrstr = 'CHNM TMBR'
    tmbrtup = (1, 2, 3, 4, 5)

    subset = 'NC021026'
    said = 1024
    bearaz = 25
    soza = 25
    solazi = 25

    def __init__(self):
        self.dt = None

    def datetime(self):
        '''
        Return the datatime object for this record. It may need to be created
        first. This assumes that the individual data time fields are defined.
        '''
        if not self.dt:
            self.dt = datetime.datetime(
              self.year, self.month, self.day,
              self.hour, self.min, self.sec
            )

        return self.dt

class TD_BUFR_File:
    '''
    Encapsulates the BUFR file for a given 6 hour time range. A new object is
    created when the first record of a new time range is added.
    '''

    # Set the file time range upper limits in seconds
    HOUR3 = 3*3600
    HOUR9 = 9*3600
    HOUR15 = 15*3600
    HOUR21 = 21*3600
    HOUR24 = 24*3600
    SEARCH_HOURS = (HOUR3, HOUR9, HOUR15, HOUR21, HOUR24)

    def __init__(self, log, version_str, h5_fname, record):
        self.log = log
        self.h5_fname = h5_fname

        # Determine the bufr file time range that would include this record
        # Use it to create the file name
        # Go through the upper limits of the time ranges in seconds
        sec_of_day = record.sec + record.min*60 + record.hour*3600
        for test_sec in self.SEARCH_HOURS:
            if sec_of_day <= test_sec:
                break
        if test_sec == self.HOUR24:
            end_hour = 3
        else:
            end_hour = test_sec // 3600

        self.dt_end = datetime.datetime(record.year, record.month, record.day,
                                        end_hour, 0, 0)
        if test_sec == self.HOUR24:
            # end_hour is 3 on the next day
            self.dt_end = self.dt_end + datetime.timedelta(days=1)

        self.dt_start = self.dt_end - datetime.timedelta(hours=6)

        # Create a time stamp for the end of the gfs file
        self.dt_gfs_end = self.dt_end - datetime.timedelta(minutes=45)

        dt_center = self.dt_end - datetime.timedelta(hours=3)

        # Use the dt_center fields to create the gfs bufr file name
        # The gfs file will be written first, so the next few object variables
        # are set accordingly. They will all be changed when it is time to
        # write the gdas file.
        self.bufr_file = (
          'gfs.t{:02d}z.tempd.{:d}{:02d}{:02d}{}.bufr_d'.format(
            dt_center.hour, dt_center.year, dt_center.month, dt_center.day,
            version_str
          )
        )

        # Determine the relative lower level directory path this file should
        # end up in.
        self.bufr_subdirs = 'gfs.{:d}{:02d}{:02d}/{:02d}'.format(
          dt_center.year,
          dt_center.month,
          dt_center.day,
          dt_center.hour
        )

        # But, initially leave out the subdirs and give the file name a
        # .<hdf5_file_name>.part extension. This will be the file's path while
        # it is being written. The <hdf5_file_name> is there to keep the file
        # name from colliding with any bufr file of the same name that might be
        # created by a different converter
        self.bufr_path_part = os.path.join(
          os.environ['BUFR_OUTPUT_DIR'],
          self.bufr_file + '.' + self.h5_fname + '.part'
        )

        # See if the .part file already exists
        if os.path.isfile(self.bufr_path_part):
            # Probably a BUFR error on the previous run
            raise BufrError(self.bufr_path_part)

        # Open the bufr file
        self.log.info('Opening: %s', self.bufr_path_part)
        self.bufr = ncepbufr.open(self.bufr_path_part,'w',table='tempest-D.table')

        # The file data-hour is needed when the first bufr message is created
        self.datehour = int('{:d}{:02d}{:02d}{:02d}'.format(
          dt_center.year, dt_center.month, dt_center.day, dt_center.hour
        ))

        # Open the first bufr message
        # Again, the message length limit will cause messages to be closed and
        # new ones opened automatically, all with the same date-hour.
        self.bufr.open_message(TD_record.subset, self.datehour)
    
        self.subsetcnt = 0

        # Write this first record to the file
        self.write_a_record(record)

    def write_a_record(self, record):
        '''
        Writes the record to a bufr subset
        '''

        # SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH CHSQ
        scalars1 = np.array((
          record.said, record.year, record.month, record.day, record.hour,
          record.min, record.sec, record.lat, record.lon, record.chi
        ), dtype = np.float64)

        # CLAVR SAZA BEARAZ SOZA SOLAZI SANG FOVN SLNM
        scalars2 = np.array((
          record.cld_mask, record.zenith, record.bearaz, record.soza,
          record.solazi, record.scan, record.sangle_idx, record.sline_idx
        ), dtype = np.float64)

        # Write the header values to BUFR
        scalarstr_tup = (record.scalarstr1, record.scalarstr2)
        scalars_tup   = (scalars1, scalars2)
        for scalarstr, scalars in zip(scalarstr_tup, scalars_tup):
            # "bufrfy" the missing values
            scalars[scalars < -998] = 10E10

            for name, value in zip(scalarstr.split(' '), scalars):
                self.log.debug('%-10s%s', name, '{:.4f}'.format(value))

            self.bufr.write_subset(scalars, scalarstr)

        # "Replicate" the water paths and channels
        # ILWP TMBR
        for replstr, repltup, repl in zip((record.ilwpstr, record.tmbrstr),
                                          (record.ilwptup, record.tmbrtup),
                                          (record.ilwp, record.Tb)):
            # "bufrfy" the missing values
            repl[repl < -998] = 10E10

            repl = np.array((repltup, repl), dtype = np.float64)
            self.log.debug('%s:\n%s', replstr, repl)

            subs_end = (replstr == 'CHNM TMBR')
            self.log.debug('subs_end: %s', subs_end)
            self.bufr.write_subset(repl, replstr, rep=True, end=subs_end)

        self.subsetcnt += 1

        if self.subsetcnt % 500000 == 0:
            self.log.info('%d subsets written to %s', self.subsetcnt,
                          os.path.basename(self.bufr_path_part))

    def write_record(self, record):
        # First, the date/time should be within this bufr file's time range
        if not (record.datetime() > self.dt_start and
                record.datetime() <= self.dt_end):
            raise ValueError(
              (
                'It is incorrect for this record: %s, angle_idx: %d, line_idx:'
                ' %d, to be written to this BUFR file: %s'
              ) % (record.datetime().isoformat(), record.sangle_idx,
                   record.sline_idx, os.path.basename(self.bufr_path))
            )

        if record.datetime() <= self.dt_gfs_end:
            # The file name better still start with 'gfs'
            if 'gfs' != os.path.basename(self.bufr_path_part)[0:3]:
                raise ValueError(
                  (
                    'This record: %s, angle_idx: %d, line_idx:'
                    ' %d, should be written to the gfs file instead of: %s'
                  ) % (record.datetime().isoformat(), record.sangle_idx,
                       record.sline_idx, os.path.basename(self.bufr_path_part))
                )

        elif 'gfs' == os.path.basename(self.bufr_path_part)[0:3]:
            # It's time to switch to the gdas file

            self.bufr.close_message()
            self.bufr.close()

            self.log.info('All %d subsets written to %s', self.subsetcnt,
                          os.path.basename(self.bufr_path_part))

            # Move this to the appropriate subdirs without the .part extensions
            bufr_dir = os.path.join(
              os.environ['BUFR_OUTPUT_DIR'],
              self.bufr_subdirs
            )
            os.makedirs(bufr_dir, exist_ok=True)

            gfs_path = os.path.join(bufr_dir, self.bufr_file)
            os.rename(self.bufr_path_part, gfs_path)
            
            # Switch to the gdas file
            self.bufr_subdirs = self.bufr_subdirs.replace('gfs.', 'gdas.', 1)
            self.bufr_file = self.bufr_file.replace('gfs.', 'gdas.', 1)

            self.bufr_path_part = os.path.join(
              os.environ['BUFR_OUTPUT_DIR'],
              self.bufr_file + '.' + self.h5_fname + '.part'
            )

            # Make sure it doesn't already exist
            if os.path.isfile(self.bufr_path_part):
                raise ValueError(
                  (
                    'This record: %s, angle_idx: %s, line_idx:'
                    ' %s, should be the first record written to the gdas file,'
                    ' but the gdas file "%s" already exists'
                  ) % record.datetime().isoformat(), record.sangle_idx,
                      record.sline_idx, os.path.basename(self.bufr_path_part)
                )

            shutil.copy(gfs_path, self.bufr_path_part)

            # Open the gdas bufr file for appending
            self.log.info('Opening: %s for appending', self.bufr_path_part)
            self.bufr = ncepbufr.open(self.bufr_path_part, 'a')

            # Open the first bufr message
            self.bufr.open_message(TD_record.subset, self.datehour)

        # Write the record to the current bufr file
        self.write_a_record(record)

    def close(self):
        self.bufr.close_message()
        self.bufr.close()

        self.log.info('All %d subsets written to %s', self.subsetcnt,
                      os.path.basename(self.bufr_path_part))

        # Move this to the appropriate subdirs without the .part extensions
        bufr_dir = os.path.join(
          os.environ['BUFR_OUTPUT_DIR'],
          self.bufr_subdirs
        )
        os.makedirs(bufr_dir, exist_ok=True)

        gdas_path = os.path.join(bufr_dir, self.bufr_file)
        os.rename(self.bufr_path_part, gdas_path)

def read_in_dataset(dataset):
    '''
    Read the Tempest-D HDF5 dataset into a numpy array and then convert all the
    missing values in the array to -999.
    '''
    ds_array = dataset[:]

    ds_array[np.isnan(ds_array)] = -999
    ds_array[ds_array < -990] = -999

    return ds_array

def h5_to_bufr(h5_path):
    # Setup logging. With the HDF5 file name as the log message prefix. Note
    # that the output directly from BUFRLIB goes to stdout, so stdout has to be
    # redirected so that logging output for both BUFRLIB and python go to the
    # same place.
    log = logging.getLogger(os.path.basename(h5_path))

    try:

        # Get the version string from the hdf5 file name if there is one
        h5_version = ''
        h5_version_match = re.search(os.environ['H5_VERSION_REGEX'], h5_path)
        if h5_version_match:
            pieces = h5_version_match.groupdict()
            h5_version = pieces.get('version', '')

        log.info('Opening: %s', h5_path)
        hdf5 = h5py.File(h5_path, 'r')

        # Read all the hdf5 datasets into memory - it's just too slow otherwise
        # Use the year as the key - if a pixel's year is missing, skip it
        h5_year = read_in_dataset(hdf5['/year'])

        log.debug('year[:, 1]: %s', h5_year[:, 1])
        log.debug('year[:, 1] dtype, itemsize, shape: %s, %s, %s',
                  h5_year[:, 1].dtype, h5_year[:, 1].itemsize,
                  h5_year[:, 1].shape)
        log.debug('h5_year.shape[0]: %d', h5_year.shape[0])

        if h5_year.shape[0] == 400:
            # Data includes the full scan - reduce it to the middle 160 scan
            # line indexes. Full scan middle index = 199
            # scan_e is exactly the last index wanted
            scan_s = 199 - 79
            scan_e = 199 + 80
        elif h5_year.shape[0] == 160:
            scan_s = 0
            scan_e = 159
        else:
            raise ValueError(
              'Bad number of scan angles: %d' % h5_year.shape[0]
            )
        
        h5_month = read_in_dataset(hdf5['/month'])
        h5_day = read_in_dataset(hdf5['/day'])
        h5_hour = read_in_dataset(hdf5['/hour'])
        h5_min = read_in_dataset(hdf5['/minute'])
        h5_sec = read_in_dataset(hdf5['/second'])

        h5_lat = read_in_dataset(hdf5['/pixel latitude'])
        h5_lon = read_in_dataset(hdf5['/pixel longitude'])

        h5_chi = read_in_dataset(hdf5['/chi'])

        # Convert the water path's g m^-2 to kg m^-2 after they are read in
        h5_iwp = read_in_dataset(hdf5['/iwp'])
        h5_iwp[h5_iwp > -998] *= G_TO_KG

        h5_lwp = read_in_dataset(hdf5['/lwp'])
        h5_lwp[h5_lwp > -998] *= G_TO_KG

        h5_converge = read_in_dataset(hdf5['/converge'])

        h5_zenith = read_in_dataset(hdf5['/zenith_angle'])
        h5_scan = read_in_dataset(hdf5['/scan_angle'])

        h5_89 = read_in_dataset(hdf5['/Tb 89 GHz'])
        h5_165 = read_in_dataset(hdf5['/Tb 165 GHz'])
        h5_176 = read_in_dataset(hdf5['/Tb 176 GHz'])
        h5_180 = read_in_dataset(hdf5['/Tb 180 GHz'])
        h5_182 = read_in_dataset(hdf5['/Tb 182 GHz'])

        bfile = None
        recordcnt = 1
        # Loop through the scan lines
        sline_end = 0
        for sline_idx in range(h5_year.shape[1]):
            # Restrict records written for testing
            #if recordcnt > 50:
            #    break
            year_scan = h5_year[:, sline_idx]

            # Skip all scan lines where the year values are all invalid
            if np.all(year_scan == -999):
                continue

            if not sline_idx == sline_end:
                log.info('valid scan line gap: %s:%s', sline_end, sline_idx - 1)
            # Set it to next expected scan line index
            sline_end = sline_idx + 1

            log.debug('year_scan: %s', year_scan)

            for sangle_idx in range(h5_year.shape[0]):
                # Skip any pixels outside the valid scan angle index range
                if sangle_idx < scan_s or sangle_idx > scan_e:
                    continue

                # Also skip any pixels with missing lat or lon values as well
                # as those with a missing year
                lat_val = h5_lat[sangle_idx, sline_idx]
                lon_val = h5_lon[sangle_idx, sline_idx]
                if year_scan[sangle_idx] < -998 or \
                  lat_val < -998 or lon_val < -998:
                    continue

                # Restrict records written for testing
                #if recordcnt > 50:
                #    break

                # Create the record obj
                record = TD_record()

                #if recordcnt < 200000:
                #    # Indent 'record.year =' line to cause python error
                record.year = year_scan[sangle_idx]
                record.month = h5_month[sangle_idx, sline_idx]
                record.day = h5_day[sangle_idx, sline_idx]
                record.hour = h5_hour[sangle_idx, sline_idx]
                record.min = h5_min[sangle_idx, sline_idx]
                record.sec = h5_sec[sangle_idx, sline_idx]

                record.lat = lat_val
                record.lon = lon_val

                record.chi = h5_chi[sangle_idx, sline_idx]

                record.ilwp = np.array((
                  h5_iwp[sangle_idx, sline_idx],
                  h5_lwp[sangle_idx, sline_idx]
                ), dtype = np.float64)

                record.cld_mask = h5_converge[sangle_idx, sline_idx]

                record.zenith = h5_zenith[sangle_idx, sline_idx]
                record.scan = h5_scan[sangle_idx, sline_idx]

                # Keep the scan angle index to 0-159 in the BUFR file
                record.sangle_idx = sangle_idx - scan_s
                record.sline_idx = sline_idx

                record.Tb = np.array((
                  h5_89[sangle_idx, sline_idx],
                  h5_165[sangle_idx, sline_idx],
                  h5_176[sangle_idx, sline_idx],
                  h5_180[sangle_idx, sline_idx],
                  h5_182[sangle_idx, sline_idx]
                ), dtype = np.float64)

                log.debug(
                  ('record lat, lon, year, month, day, hour, min, sec, scan,'
                   ' zenith, Tb: %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s'),
                  record.lat, record.lon, record.year, record.month,
                  record.day, record.hour, record.min, record.sec, record.scan,
                  record.zenith, record.Tb
                )

                #sys.exit(0)

                if not bfile:
                    # The bufr file object does not exist - create it
                    bfile = TD_BUFR_File(log, h5_version,
                                         os.path.basename(h5_path), record)
                else:
                    bfile.write_record(record)

                recordcnt += 1

                #if recordcnt > 200000:
                #    # Cause BUFR error
                #    TD_record.scalarstr1 = 'SAID YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH CHSQ CLOWN'

        if sline_idx - sline_end > 1:
            log.info('valid scan line gap: %s:%s', sline_end, sline_idx)

        bfile.close()

        hdf5.close()
            
        # All the above seems to have succeeded so move the hdf5 file to the
        # processed dir
        shutil.move(h5_path, os.environ['H5_PROCESSED_DIR'])

    except BufrError as err:
        # The current run is a retry of a failed run. Abort
        # The bufr_path_part for the file will be in the exception
        log.error(
          'This conversion has already failed once, probably due to a BUFR'
          ' error. Aborting!'
        )
        log.error(
          'Moving %s to %s, and deleting %s',
          os.path.basename(h5_path), os.environ['H5_BAD_DIR'],
          err.bufr_path_part
        )
        shutil.move(h5_path, os.environ['H5_BAD_DIR'])
        os.remove(err.bufr_path_part)
        raise

    except (Exception, KeyboardInterrupt) as err:
        log.error('Exception occured:')
        exc_info = sys.exc_info()
        traceback.print_exception(*exc_info)
        if isinstance(err, KeyboardInterrupt):
            log.error('Termination was deliberate. Leaving %s where it is',
                      h5_path)
        else:
            log.error('moving %s to: %s',
                      os.path.basename(h5_path), os.environ['H5_BAD_DIR'])
            shutil.move(h5_path, os.environ['H5_BAD_DIR'])

        if os.path.isfile(bfile.bufr_path_part):
            log.error('deleting %s', bfile.bufr_path_part)
            os.remove(bfile.bufr_path_part)
        else:
            log.error(
              'The \'.part\' bufr file was never created -'
              ' no need to clean it up'
            )
        raise
