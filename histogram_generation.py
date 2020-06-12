from glob import glob
import re
import os
import sys
import numpy as np
import yaml
import h5py
import pandas as pd

gmi_data_dir = "/mnt/ssdatenas/gmi_native"
case = "may2019"
crtm_gmi_data_dir = "/mnt/ssdatenas/crtm_output"
total_fcst_hours = 240
fcst_step = 3

#Brightness temperature range for histograms
#arange is inclusive of the first number but exclusive of the last
TB_RANGE = np.arange(180, 311)

#CRTM landcover masking options; True means that data will be filtered out
CRTM_MASK_DICT = {"water": False, "land": True, "ice": True}
upper_case_keys = [k.capitalize() for k,v in CRTM_MASK_DICT.items() if v]
mask_words_for_plot_name = "".join(upper_case_keys)
if not mask_words_for_plot_name:
     mask_words_for_plot_name= "None"

gmi_regex = re.compile("(?P<preamble>[1C.GPM.GMI.XCAL2016-C.](.*))(?P<yyyymmdd>[0-9]{8})-S(?P<start_hhmmss>[0-9]{6})-E(?P<end_hhmmss>[0-9]{6}).(?P<granule>[0-9]{6}).(?P<version>[A-Z0-9]{4}).HDF5$")
relevant_hours = ["00", "03", "06", "09", "12", "15", "18", "21"]
forecast_sec_dict = {"00": 0*60*60, "03": 3*60*60, "06": 6*60*60, "09": 9*60*60, "12": 12*60*60, "15": 15*60*60, "18": 18*60*60, "21": 21*60*60}
    
#FV3 grid: [1536,768]
fv3_lat_step = 180./768.
fv3_lon_step = 360./1536.

#South to North
LAT_GRID_SOUTH = np.arange(-90, 90 + fv3_lat_step, fv3_lat_step, dtype=float)[:-1]
LAT_GRID_NORTH = np.arange(-90, 90 + fv3_lat_step, fv3_lat_step, dtype=float)[1:]

#West to East
LON_GRID_WEST = np.arange(0, 360 + fv3_lon_step, fv3_lon_step, dtype=float)[:-1]
LON_GRID_EAST = np.arange(0, 360 + fv3_lon_step, fv3_lon_step, dtype=float)[1:]

#South to North, starting half a grid box North of the southern most bin line and ending half a grid box North of the northern most bin line
LAT_CENTERS = np.arange(LAT_GRID_SOUTH[0] + fv3_lat_step / 2., LAT_GRID_NORTH[-1], fv3_lat_step, dtype=float)
#West to East, starting half a grid box East of the western most bin line and ending half a grid box east of the easternmost bin line
LON_CENTERS = np.arange(LON_GRID_WEST[0] + fv3_lon_step / 2., LON_GRID_EAST[-1], fv3_lon_step, dtype=float)
     
def create_crtm_landcover_mask(landcover_data):

    """
    Create mask for array elements based on criteria defined by
    CRTM_MASK_DICT

    Parameters
    ----------
    landcover_data : ndarray
        The landcover information from the CRTM
        0=water (ocean + inland), 1=land, 2=ice

    Returns
    -------
    ndarray
        Boolean array indicating which elements are masked/unmasked by
        the criteria defined by CRTM_MASK_DICT

    """

    #Three possible values: 0=water (ocean + inland), 1=land, 2=ice
    possible_filters = {"water":0, "land":1, "ice":2}
    filter_vals = []
    vals_to_mask = {filter_vals.append(v) for k, v in possible_filters.items() if CRTM_MASK_DICT[k]}
    mask = np.logical_or.reduce([landcover_data == v for v in filter_vals])
    return mask
    
def read_crtm_data(crtm_file, fill_val=None):
    """
    Read in synthetic imagery data from the CRTM-produced binary file

    Parameters
    ----------
    crtm_file : str
        Full path to the CRTM file
    fill_val : int or float
        Optional parameter to supply the array fill value for masking

    Returns
    -------
    ndarray
        CRTM synthetic data masked according to the fill_val

    """

    dat = np.fromfile(crtm_file, dtype=np.float32)
    dat = dat.reshape([768,1536])
    dat = np.ma.masked_where(dat == fill_val, dat)

    return dat
    
def get_hdf5_data(var, hdf5_file):
    """
    Extract given variable data from an HDF-5 file.

    Parameters
    ----------
    var : str
        Full path to data field within an HDF-5 file
    hdf5_file : str
        Full path to an HDF-5 file

    Returns
    -------
    ndarray
        Data contained in the HDF-5 field

    """

    f = h5py.File(hdf5_file, "r")
    try: 
        data = f[var][:]
    except:
        print("Problem retrieving " + var + " from " + hdf5_file)
        f.close()
        return False
    f.close()
    return data
    
def fv3_gridify(native_data, native_lat, native_lon):
    """
    Regrid data to the FV3-GFS grid defined by 
    LAT_CENTERS and LON_CENTERS

    Parameters
    ----------
    native_data : ndarray
        Data array to regrid to the FV3-GFS grid
    native_lat : ndarray
        Latitude array corresponding to data to regrid
    native_lon : ndarray
        Longitude array corresponding to data to regrid

    Returns
    -------
    ndarray
        Data regridded to the FV3-GFS grid

    """
    
    grid = np.empty([len(LON_CENTERS), len(LAT_CENTERS)], dtype=np.object)
    
    native_latlondat = zip(list(native_lat.flatten()), list(native_lon.flatten()), list(native_data.flatten()))
    
    for ll in native_latlondat:
        fv3_yidx = np.where(np.logical_and(ll[0] >= LAT_GRID_SOUTH, ll[0] <= LAT_GRID_NORTH))[0][0]
        fv3_xidx = np.where(np.logical_and(ll[1] >= LON_GRID_WEST, ll[1] <= LON_GRID_EAST))[0][0]

        if grid[fv3_xidx,fv3_yidx] is None:
            grid[fv3_xidx,fv3_yidx] = [ll[2]]
        else:
            grid[fv3_xidx,fv3_yidx].append(ll[2])
    
    x_action, y_action = np.nonzero(grid)

    for x, y in zip(x_action, y_action):
        if grid[x,y] is not None:
            grid[x,y] = np.mean(grid[x,y])
    
    return grid.T.astype(np.float32)
    
def create_gmi_fcst_time_lut():
    """
    Create the look up table of GMI files that correspond to each
    forecast time

    Parameters
    ----------
    none

    Returns
    -------
    dict
        GMI files for each cycle start and forecast step
    yaml file
        dictionary above written to a file for future ease of use 

    """
    dates = list(set([gmi_regex.search(os.path.basename(f)).groupdict()["yyyymmdd"] for f in glob(os.path.join(gmi_data_dir, "*", "*"))]))
    data_dict = {d+h:"" for h in relevant_hours for d in dates}
    for f in sorted(glob(os.path.join(gmi_data_dir, "*", "*"))):
        print(os.path.basename(f))
        try:
            filename_dict = gmi_regex.search(os.path.basename(f)).groupdict()
            #print(filename_dict)
        except:
            print(f + " does not appear to be a GMI file.")
            print("Regex: " + gmi_regex)
            print("Exiting.")
            sys.exit()            
        filename_start_hour = filename_dict["start_hhmmss"][0:2]
        if filename_start_hour in relevant_hours:
            print("Relevant hour")
            print("Inserting " + f + " into " + filename_dict["yyyymmdd"] + ":" + filename_start_hour)
            data_dict[filename_dict["yyyymmdd"]+filename_start_hour] = f
        elif int(filename_start_hour) >= 1 and int(filename_start_hour) < 2:
            if data_dict[filename_dict["yyyymmdd"]+"00"]:
                print("File already found for " + filename_dict["yyyymmdd"] + ":00")
                #print(data_dict[filename_dict["yyyymmdd"]]["00"])
                continue
            else:
                print("Start more than 1 but less than 3")
                day_before = int(filename_dict["yyyymmdd"])-1
                file_to_insert = sorted(glob(os.path.join(gmi_data_dir, str(day_before)[2:], "*"+str(day_before)+"*")))[-1]
                print("Inserting " + file_to_insert + " into " + filename_dict["yyyymmdd"] + ":00")
                data_dict[filename_dict["yyyymmdd"]+"00"] = file_to_insert
        else: 
            start_sec = int(filename_dict["start_hhmmss"][0:2])*60*60 + int(filename_dict["start_hhmmss"][2:4])*60 + int(filename_dict["start_hhmmss"][4:6])
            end_sec = int(filename_dict["end_hhmmss"][0:2])*60*60 + int(filename_dict["end_hhmmss"][2:4])*60 + int(filename_dict["end_hhmmss"][4:6])
            test = list({k for k, v in forecast_sec_dict.items() if v > start_sec and v < end_sec})
            print(test)
            if test:
                if data_dict[filename_dict["yyyymmdd"]+test[0]]:
                    print("File already found for " + filename_dict["yyyymmdd"] + ":" + test[0])
                    print(data_dict[filename_dict["yyyymmdd"]+test[0]])
                    continue
                else:
                    print("Inserting " + f + " into " + filename_dict["yyyymmdd"] + ":" + test[0])
                    data_dict[filename_dict["yyyymmdd"]+test[0]] = f
            
    with open("corresponding_file_info.yml" , "w") as outfile:
        yaml.dump(data_dict, outfile)
    return data_dict

def gather_forecast_results_from_cycles():
    """
    Create the look up table of experiment files for each cycle
    for each forecast time for a given case

    Parameters
    ----------
    none

    Returns
    -------
    dict
        Synthetic imagery files for each experiment for each
        forecast hour and each cycle start for a given case
    yaml file
        dictionary above written to a file for future ease of use 

    """
    if glob("corresponding_file_info.yml"):
        with open("corresponding_file_info.yml" , "r") as infile:        
            gmi_file_dict = yaml.full_load(infile)
    else:
        gmi_file_dict = create_gmi_fcst_time_lut()
    
    outfile_dict = {exp:{"{0:0=3d}".format(h):[] for h in range(0, total_fcst_hours+1, fcst_step)} for exp in os.listdir(os.path.join(crtm_gmi_data_dir, case))}
    outfile_dict["Observed"] = {"{0:0=3d}".format(h):[] for h in range(0, total_fcst_hours+1, fcst_step)}
    for exp, hr_dict in outfile_dict.items():
        for fcst_hr in hr_dict.keys():
            outfile_dict[exp][fcst_hr] = glob(os.path.join(crtm_gmi_data_dir, case, exp, "*", "tb_gmi_89VGhz_*_" + fcst_hr))
    #key off Control to get GMI files
    for fcst_hr, file_list in outfile_dict["Control"].items():
        for f in file_list:
            obs_time = re.split("_", os.path.basename(f))[3]
            outfile_dict["Observed"][fcst_hr].append(gmi_file_dict[obs_time])
    
    with open("forecast_hour_file_groups_" + case + ".yml" , "w") as outfile:
        yaml.dump(outfile_dict, outfile)
    return outfile_dict
    
def build_hist_df(exp, exp_file_dict):
    """
    Creates a histogram of brightness temperature bin counts
    for a given experiment or observations for a given case

    Parameters
    ----------
    exp : str
        Experiment name
    exp_file_dict: dict
        Synthetic imagery files for each experiment for each
        forecast hour and each cycle start for a given case
        (output of gather_forecast_results_from_cycles function)
        

    Returns
    -------
    boolean
        Returns true if the function completes
        (assumes the file has been written)

    """
        
    hist_df = pd.DataFrame(index=TB_RANGE[:-1], columns=exp_file_dict.keys())
    
    for fcst_hr, file_list in exp_file_dict.items():
        temp_data_array = np.ma.empty([])
        for i, f in enumerate(file_list):
            #print(f)
            if exp == "Observed":
                obs_tb_native = get_hdf5_data("/S1/Tc", f)
                obs_tb89GhzV_native = obs_tb_native[:, :, 7]

                obs_lat = get_hdf5_data("/S1/Latitude", f)
                obs_lon = get_hdf5_data("/S1/Longitude", f)
                #Convert -180 to 180 -> 0 to 360
                obs_lon[np.where(obs_lon < 0.)] = obs_lon[np.where(obs_lon < 0.)] + 360.

                #Regrid
                data_temp = fv3_gridify(obs_tb89GhzV_native, obs_lat, obs_lon)
            else:    
                data_temp = read_crtm_data(f, fill_val=0)
            temp_data_array = np.ma.append(temp_data_array, np.ma.masked_where(crtm_landcover_mask, data_temp))
        hist_df[fcst_hr] = np.histogram(temp_data_array.compressed(), bins=TB_RANGE)[0]
    
    hist_df.to_csv("Tb_frequency_" + case + "_" + exp + "_Mask" + mask_words_for_plot_name + ".csv")
    
    return True

CRTM_LANDCOVER_FILE = "LAND000000"
crtm_landcover_data = read_crtm_data(CRTM_LANDCOVER_FILE)
crtm_landcover_mask = create_crtm_landcover_mask(crtm_landcover_data)

if glob("forecast_hour_file_groups_" + case + ".yml"):
    with open("forecast_hour_file_groups_" + case + ".yml" , "r") as infile:        
            file_dict = yaml.full_load(infile)
else:
    file_dict = gather_forecast_results_from_cycles()

for exp, exp_dict in file_dict.items():
    build_hist_df(exp, exp_dict)

