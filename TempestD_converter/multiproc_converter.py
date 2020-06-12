# coding: utf-8

# Stock modules
import os
import sys
import glob
import time
import re
import logging
import multiprocessing as mproc
import signal

# Local modules
import tempest_h5_to_bufr

'''
"python multiproc_converter.py"
Converts Louie's preprocessed HDF5 Tempest-D data files to BUFR.
Required env variables:
'LOGLEVEL' - The logging level (optional: defaults to 'info')
'H5_DATA_DIR' - The directory that contains the Tempest-D HDF5 files
'H5_PROCESSED_DIR' - The directory the Tempest-D HDF5 files will be moved to
  after they area processed successfully
'H5_BAD_DIR' - The directory the Tempest-D HDF5 files will be moved to if the
  conversion fails
'BUFR_OUTPUT_DIR' - The directory that the BUFR files will be written to
'SLEEP_SECONDS' - Wait time between checks for input files or for cpus to
  become available (optional: defaults to 60 sec)
'CPU_FRACTION' - The fraction of the node's total cpus that can be used
  (optional: defaults to 0.5)
'H5_VERSION_REGEX' - The HDF5 file name version extraction regex
  (optional: defaults to '(?P<version>\_v(\d|\.)*\d)'
'''

def signal_handler(sig, frame):
    print('__main__: Signal caught:', str(sig))
    print('__main__: children:', mproc.active_children())
    for proc in mproc.active_children():
        os.kill(proc.pid, signal.SIGINT)
        # Try to separate the log message blocks
        time.sleep(.001)
    os.kill(os.getpid(), signal.SIGINT)

def set_environ(log):
    for req_env in ('H5_DATA_DIR', 'H5_PROCESSED_DIR', 'H5_BAD_DIR',
                    'BUFR_OUTPUT_DIR'):
        if req_env not in os.environ:
            log.error('%s must be defined', req_env)
            sys.exit(1)
    log.debug('H5_DATA_DIR: %s', os.environ['H5_DATA_DIR'])

    if not os.environ.get('SLEEP_SECONDS'):
        os.environ['SLEEP_SECONDS'] = '60'

    if not os.environ.get('CPU_FRACTION'):
        os.environ['CPU_FRACTION'] = '.5'

    if not os.environ.get('H5_VERSION_REGEX'):
        os.environ['H5_VERSION_REGEX'] = '(?P<version>\_v(\d|\.)*\d)'

def spawn_converters():
    # Setup logging.
    logging.basicConfig(
      format='%(asctime)s %(levelname)-8s%(name)s: %(message)s',
      level=os.environ.get('LOGLEVEL', 'INFO').upper()
    )
    log = logging.getLogger(__name__)

    set_environ(log)

    # Init the glob pattern for reading the h5 file names from the input dir
    h5_glob = os.path.join(os.environ['H5_DATA_DIR'], '*.h5')

    # Wait SLEEP_SECONDS amount of time for processes to complete
    sleep_seconds = int(os.environ['SLEEP_SECONDS'])

    # Set the number of processes limit to half the number of cpus
    log.info('cpu_fraction: %s', os.environ['CPU_FRACTION'])
    proc_limit = int(mproc.cpu_count() *
                     float(os.environ['CPU_FRACTION']))
    if proc_limit <= 0:
        log.error(
          'proc_limit calculates to 0 CPUs. CPU_FRACTION must be increased'
        )
        return 1
    log.info('proc_limit: %s', proc_limit)

    while True:
        # Loop through the input files
        h5_paths = glob.glob(h5_glob)

        # Get the list of active converters if any
        proclist = mproc.active_children()

        # Insure these hdf5 files are not already being processed
        # The proc names will be set to the hdf5 paths later on
        proc_names = [proc.name for proc in proclist]
        h5_paths[:] = [path for path in h5_paths if path not in proc_names]
        if len(h5_paths) > 0:
            log.info('h5_paths: %s', h5_paths)

        did_sleep = False
        for h5_path in h5_paths:
    
            # See if we should wait until we are below the process limit
            proclist = mproc.active_children()
            while len(proclist) >= proc_limit:
                log.info('Doing proc_limit sleep')
                time.sleep(sleep_seconds)
                did_sleep = True
                proclist = mproc.active_children()

            # Start the conversion sub-process for this file
            log.info('Starting converter for %s', h5_path)
            proc = mproc.Process(
              target = tempest_h5_to_bufr.h5_to_bufr,
              args = (h5_path,),
              name = h5_path
            )
            proc.start()

        # Insure a sleep occurs before every recheck of the input dir
        if not did_sleep:
            log.debug('Doing input dir check sleep')
            time.sleep(sleep_seconds)

if __name__ == '__main__':
    # Convert SIGTERM to SIGINT for graceful termination
    signal.signal(signal.SIGTERM, signal_handler)
    spawn_converters()
