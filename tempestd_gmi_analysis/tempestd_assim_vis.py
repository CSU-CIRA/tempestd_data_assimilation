from glob import glob
import re
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
import yaml
import h5py
import pandas as pd
import datetime
import argparse

gmi_data_dir = "/mnt/ssdatenas/gmi_native"
case = "dec2018"
crtm_gmi_data_dir = "/mnt/ssdatenas/crtm_output"
total_fcst_hours = 240
fcst_step = 3
cycle_starts = {"dec2018": ["2018120812", "2018120900", "2018120912", "2018121000",
               "2018121012", "2018121100", "2018121112", "2018121200", "2018121212"],
               "may2019": ["2019051212", "2019051300", "2019051312", "2019051400",
               "2019051412", "2019051500", "2019051512", "2019051600", "2019051612",
               "2019051700", "2019051712", "2019051800", "2019051812", "2019051900",
               "2019051912", "2019052000", "2019052012", "2019052100", "2019052112",
               "2019052200", "2019052212"]}
base_hist_plot_dir = "/visualization_tools/src/pub_testing"#"/mnt/ssdatenas/histogram_plots"
base_map_plot_dir = "/visualization_tools/src/pub_testing"#"/mnt/ssdatenas/diff_maps"
base_timeseries_plot_dir = "/visualization_tools/src/pub_testing"#"/mnt/ssdatenas/tempestd_gmi_global_diff_timeseries"

LINE_PLOT_DICT = {"Observed": {"color": "black"},
                 "Control": {"color": "blue"},
                 "AddMHS": {"color": "green"},
                 "AddTEMPESTD": {"color": "orange"}}

CRTM_LANDCOVER_FILE = "LAND000000"

#Brightness temperature range for histograms
#arange is inclusive of the first number but exclusive of the last
TB_RANGE = np.arange(180, 311)

FIGSIZE = (12,6)
#CRTM landcover masking options; True means that data will be filtered out
CRTM_MASK_DICT = {"water": False, "land": True, "ice": True}
upper_case_keys = [k.capitalize() for k,v in CRTM_MASK_DICT.items() if v]
mask_words_for_plot_name = "".join(upper_case_keys)
if not mask_words_for_plot_name:
     mask_words_for_plot_name= "None"

gmi_regex = re.compile("(?P<preamble>[1C.GPM.GMI.XCAL2016-C.](.*))(?P<yyyymmdd>[0-9]{8})-S(?P<start_hhmmss>[0-9]{6})-E(?P<end_hhmmss>[0-9]{6}).(?P<granule>[0-9]{6}).(?P<version>[A-Z0-9]{4}).HDF5$")
relevant_hours = ["00", "03", "06", "09", "12", "15", "18", "21"]
forecast_sec_dict = {"00": 0*60*60, "03": 3*60*60, "06": 6*60*60, "09": 9*60*60, "12": 12*60*60, "15": 15*60*60, "18": 18*60*60, "21": 21*60*60}

orbit_window_regex = re.compile("gmi_orbit_window_(?P<yyyymmddhh>[0-9]{10})")
tb_regex = re.compile("tb_gmi_89VGhz_(?P<yyyymmddhh>[0-9]{10})_(?P<fcst_step>[0-9]{3})")

dates = list(set([gmi_regex.search(os.path.basename(f)).groupdict()["yyyymmdd"] for f in glob(os.path.join(gmi_data_dir, "*", "*"))]))
data_dict = {d+h:"" for h in relevant_hours for d in dates}

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

def plot_data_on_map(data, outfile, central_lon, vmax, vmin, cmap, cbar_label, figsize=FIGSIZE):
    """
    Creates map plot of the given data

    Parameters
    ----------
    data_dict : ndarray
        The data to be plotted on the map
    outfile : str
        Full path to the output plot name
    central_lon : float
        The central longitude for the map
    vmax : float
        The maximum value for the colorbar/data
    vmin : float
        The minimum value for the colorbar/data
    cmap : matplotlib colormap
        The colormap to use for the plot
    cbar_label : str
        The label for the colorbar
    figsize: tuple
        Optional parameter defining the figure size (X, Y)
        

    Returns
    -------
    boolean
        Returns true if the function completes
        (assumes the plot has been created)
    """
    
    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
    
    im = ax.imshow(data, origin='lower', transform=ccrs.PlateCarree(central_longitude=central_lon), cmap=cmap, vmin=vmin, vmax=vmax)
    ax.coastlines(resolution='110m')
    
    cbar = fig.colorbar(im, ax=ax, pad=0.01, shrink=0.85)
    cbar.ax.set_xlabel(cbar_label)
    fig.savefig(outfile, bbox_inches='tight')
    plt.close()
    
    return True
    
def plot_histogram(hist_data_dict, figsize=FIGSIZE):

    """
    Creates a plot of the frequency of brightness temperature
    bin counts for each given experiment/observations for a
    given case

    Parameters
    ----------
    hist_data_dict : dict
        The histogram dictionary, which contains subdictionaries of
        brightness temperature bin counts for all experiments/observations
        for the given case
    figsize: tuple
        Optional parameter defining the figure size (X, Y)
        

    Returns
    -------
    boolean
        Returns true if the function completes
        (assumes the plot has been created)

    """

    for fcst_hr in range(0, total_fcst_hours+1, fcst_step):
        hr = "{0:0=3d}".format(fcst_hr)
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        for exp in hist_data_dict.keys():
            #bin edges have one less element than the counts
            #print(hist_data_dict[exp][hr][0])
            plt.plot(TB_RANGE[:-1], hist_data_dict[exp][hr], LINE_PLOT_DICT[exp]["color"], linewidth=2)
        plt.legend(hist_data_dict.keys(), loc='best', fontsize=15)
        ax.tick_params(axis='both', labelsize=15)
        ax.set_title("Tb (t=" + hr + ")", fontsize=16)
        ax.set_xlabel('Tb [K]',fontsize=15)
        ax.set_ylabel('Count [#]',fontsize=15)
        if case == "dec2018":
            ax.set_ylim([0, 22500])
        if case == "may2019":
            ax.set_ylim([0, 55000])
        ax.grid(which='both')
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)
    
        fig.savefig(os.path.join(base_hist_plot_dir, case, "TBhistogram_Mask" + mask_words_for_plot_name + "_" + hr + ".png"))
        plt.close()
    
    return

def plot_mean_diff_along_forecast(data_dict, figsize=FIGSIZE):

    """
    Creates line plots of the global mean and standard deviation
    for each experiment for each forecast hour

    Parameters
    ----------
    data_dict : dict
        The data dictionary, which contains the global mean
        and standard deviation at each forecast hour for all
        cycles for the given case
    figsize: tuple
        Optional parameter defining the figure size (X, Y)
        

    Returns
    -------
    boolean
        Returns true if the function completes
        (assumes the plot has been created)

    """

    case_title_dict = {"dec2018": "December 2018", "may2019": "May 2019"}
    ncycles = len(cycle_starts[case])
    hours = ["{0:0=3d}".format(h) for h in range(0, 120+1, fcst_step)]

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    for exp, mean_stddev_dict in data_dict.items():
        plt.plot(np.arange(0, 41), mean_stddev_dict["mean"][0:41], LINE_PLOT_DICT[exp]["color"], linewidth=2, label=exp)
        plt.plot(np.arange(0, 41), mean_stddev_dict["stddev"][0:41], LINE_PLOT_DICT[exp]["color"], linestyle="--")
        
    plt.legend(loc='best', fontsize=15)
    ax.tick_params(axis='both', labelsize=15)
    ax.set_title("Global Mean Difference from GMI " + str(ncycles) + " Cycles\n" + case_title_dict[case], fontsize=16)
    ax.set_xlabel('Forecast Hour',fontsize=15)
    ax.set_xticks(np.arange(0, 41))
    ax.set_xticklabels(hours)
    plt.xticks(rotation=30, fontsize=7)
    ax.set_ylabel('Mean Difference(Tb)',fontsize=15)
    ax.grid(which='both')
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)

    fig.savefig(os.path.join(base_timeseries_plot_dir, "AllCycles_ExpMinusGMI_Global_Mean_"+case+".png"))
    plt.close()
    
    return True

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
    
def plot_daily_composites(data_dict, outfile, central_lon, vmax, vmin, cmap, cbar_label, figsize=FIGSIZE):
    """
    Creates paneled map plots of the given data(one per experiment,
    vertically, as indicated by the number of dictionary keys
    in the data_dict)

    Parameters
    ----------
    data_dict : dict
        The data dictionary, which contains subdictionaries of
        the data to be plotted for each experiment
    outfile : str
        Full path to the output plot name
    central_lon : float
        The central longitude for the map
    vmax : float
        The maximum value for the colorbar/data
    vmin : float
        The minimum value for the colorbar/data
    cmap : matplotlib colormap
        The colormap to use for the plot
    cbar_label : str
        The label for the colorbar
    figsize: tuple
        Optional parameter defining the figure size (X, Y)
        

    Returns
    -------
    boolean
        Returns true if the function completes
        (assumes the plot has been created)

    """
    n_experiments = len(data_dict.keys())
    fig, ax = plt.subplots(n_experiments, 1, figsize=(12,n_experiments*5), subplot_kw={'projection': ccrs.PlateCarree()})
    
    for i, (exp, data) in enumerate(data_dict.items()):
        for layer in range(0, data.shape[2]):
            im = ax[i].imshow(data[:,:,layer], origin='lower', transform=ccrs.PlateCarree(central_longitude=central_lon), cmap=cmap, vmin=vmin, vmax=vmax)
        ax[i].coastlines(resolution='110m')
        ax[i].set_title(exp + " - GMI 89V GHz Tb", fontsize=12)

        cbar = fig.colorbar(im, ax=ax[i], pad=0.03, shrink=0.95)
        cbar.ax.set_xlabel(cbar_label)
    fig.savefig(outfile, bbox_inches='tight')
    plt.close()
    
    return

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

def dt_range(start, end, delta):
    """
    Creates paneled map plots of the given data(one per experiment,
    vertically, as indicated by the number of dictionary keys
    in the data_dict)

    Parameters
    ----------
    start : datetime
        Starting datetime
    end : datetime
        Ending datetime
    delta : timedelta
        Timestep
        

    Returns
    -------
    iterator
        Datetimes between the start and end at the given step
    
    """
    
    curr = start
    while curr < end:
        yield curr
        curr += delta
        
def create_global_mean_stddev_for_forecast_hour():
    """
    Calculates the global mean and standard deviation at each
    forecast hour for all cycles for the given case

    Parameters
    ----------
    none
    

    Returns
    -------
    boolean
        Returns true if the function completes
        (assumes the file has been written)

    """
    
    with open("forecast_hour_file_groups_" + case + ".yml" , "r") as infile:        
                file_dict = yaml.full_load(infile)

    mean_stddev_dfs_dict = {exp: pd.DataFrame(columns=["mean", "stddev"]) for exp in ["Control", "AddMHS", "AddTEMPESTD"]}

    for h, fcst_hr in enumerate(file_dict["Control"].keys()):
    #for h, fcst_hr in enumerate(["000", "003"]):
        print(fcst_hr)
        exp_data = {"Control": sorted(file_dict["Control"][fcst_hr]),
                   "AddMHS": sorted(file_dict["AddMHS"][fcst_hr]),
                   "AddTEMPESTD": sorted(file_dict["AddTEMPESTD"][fcst_hr]),
                   "Observed": sorted(file_dict["Observed"][fcst_hr])}
        temp_dict = {"Control": [],
                   "AddMHS": [],
                   "AddTEMPESTD": []}
        for i, f in enumerate(exp_data["Observed"]):
            if os.path.basename(f) == "1C.GPM.GMI.XCAL2016-C.20190521-S143353-E160627.029697.V05A.HDF5":
                #This file has corrupt data that hasn't been masked
                continue
            #print(f)
            #print(exp_data["Control"][i])
            #print(exp_data["AddMHS"][i])
            #print(exp_data["AddTEMPESTD"][i])
            #print()
            obs_tb_native = get_hdf5_data("/S1/Tc", f)
            obs_tb89GhzV_native = obs_tb_native[:, :, 7]

            obs_lat = get_hdf5_data("/S1/Latitude", f)
            obs_lon = get_hdf5_data("/S1/Longitude", f)
            #Convert -180 to 180 -> 0 to 360
            obs_lon[np.where(obs_lon < 0.)] = obs_lon[np.where(obs_lon < 0.)] + 360.

            #Regrid
            gmi_data = fv3_gridify(obs_tb89GhzV_native, obs_lat, obs_lon)

            #Mask
            gmi_data = np.ma.masked_where(crtm_landcover_mask, gmi_data)
            gmi_data = np.ma.masked_where(gmi_data == -9999.9, gmi_data)

            #plot_data_on_map(gmi_data, "Obs_CycleStart" + "2019052112" + "_FcstStep"+ fcst_hr + ".png", central_lon=180., vmax=300, vmin=200, cmap="viridis", cbar_label="Tb (K)")
            #sys.exit()

            for exp in ["Control", "AddMHS", "AddTEMPESTD"]:
                print(exp)
                #cycle_start = re.split("_", os.path.basename(exp_data[exp][i]))[3]
                fcst_data = read_crtm_data(exp_data[exp][i], fill_val=0)
                fcst_data = np.ma.masked_where(crtm_landcover_mask, fcst_data)

                diff = fcst_data - gmi_data
                #print("Mean:", np.nanmean(diff))
                #print("Stddev:", np.nanstd(diff))

                #if i == 0:
                if len(temp_dict[exp]) == 0:
                    temp_dict[exp] = diff
                else:
                    temp_dict[exp]  = np.ma.dstack((temp_dict[exp], diff))
        #print()
        for exp, all_cycle_diffs in temp_dict.items():
            #print(exp)
            #print(np.nanmean(all_cycle_diffs))
            mean_stddev_dfs_dict[exp].at[h, "mean"] = np.nanmean(all_cycle_diffs)
            #print(np.nanstd(all_cycle_diffs))
            mean_stddev_dfs_dict[exp].at[h, "stddev"] = np.nanstd(all_cycle_diffs)
            #mean_stddev_dict[exp]["mean"].append(np.nanmean(all_cycle_diffs))
            #mean_stddev_dict[exp]["stddev"].append(np.nanstd(all_cycle_diffs))
    for exp, df in mean_stddev_dfs_dict.items():
        df.to_csv(exp + "_global_diff_mean_stddev_" + case + ".csv")

    return True

if __name__ == "__main__":

    """
    Run the given visualization function and produce output plots

    """

    parser = argparse.ArgumentParser(
        description="TEMPEST-D Data Assimilation Plotting",
        prefix_chars="-")
    parser.add_argument(
        "-f", "--function",
        help="Function to run. Options:\
             plot_histogram,\
             plot_mean_diff_along_forecast,\
             plot_daily_composites",
        required=True)
    
    args = parser.parse_args()
    fn = args.function
    
    crtm_landcover_data = read_crtm_data(CRTM_LANDCOVER_FILE)
    crtm_landcover_mask = create_crtm_landcover_mask(crtm_landcover_data)
    
    cycle_start_regex = re.compile("(?P<yyyy>[0-9]{4})(?P<mm>[0-9]{2})(?P<dd>[0-9]{2})(?P<hh>[0-9]{2})")
    cycle_datetime_fmt = "%Y%m%d%H"

    
    if fn == "plot_histogram":

        if (not glob("Tb_frequency_" + case + "_Control_Mask" + mask_words_for_plot_name + ".csv") or
          not glob("Tb_frequency_" + case + "_AddMHS_Mask" + mask_words_for_plot_name + ".csv") or
          not glob("Tb_frequency_" + case + "_AddTEMPESTD_Mask" + mask_words_for_plot_name + ".csv") or
          not glob("Tb_frequency_" + case + "_Observed_Mask" + mask_words_for_plot_name + ".csv")):
            print("One or more histograms do not exist. Please create with \
                histogram_generation.py")
        
        control_df = pd.read_csv("Tb_frequency_" + case + "_Control_Mask" + mask_words_for_plot_name + ".csv", index_col=0)
        mhs_df = pd.read_csv("Tb_frequency_" + case + "_AddMHS_Mask" + mask_words_for_plot_name + ".csv", index_col=0)
        td_ncd_c1_df = pd.read_csv("Tb_frequency_" + case + "_AddTEMPESTD_Mask" + mask_words_for_plot_name + ".csv", index_col=0)
        obs_df = pd.read_csv("Tb_frequency_" + case + "_Observed_Mask" + mask_words_for_plot_name + ".csv", index_col=0)

        hist_dict = {"Control": control_df.to_dict(orient="list"),
                    "AddMHS": mhs_df.to_dict(orient="list"),
                    "AddTEMPESTD": td_ncd_c1_df.to_dict(orient="list"),
                    "Observed": obs_df.to_dict(orient="list")}

        plot_histogram(hist_dict)
        
        
    elif fn == "plot_mean_diff_along_forecast":

        if (not glob("Control_global_diff_mean_stddev_" + case + ".csv") or
          not glob("AddMHS_global_diff_mean_stddev_" + case + ".csv") or
          not glob("AddTEMPESTD_global_diff_mean_stddev_" + case + ".csv")):
            create_global_mean_stddev_for_forecast_hour()

        control_df = pd.read_csv("Control_global_diff_mean_stddev_" + case + ".csv", index_col=0)
        mhs_df = pd.read_csv("AddMHS_global_diff_mean_stddev_" + case + ".csv", index_col=0)
        td_ncd_c1_df = pd.read_csv("AddTEMPESTD_global_diff_mean_stddev_" + case + ".csv", index_col=0)

        
        mean_stddev_dict = {"Control": control_df.to_dict(orient="list"),
                           "AddMHS": mhs_df.to_dict(orient="list"),
                           "AddTEMPESTD": td_ncd_c1_df.to_dict(orient="list")}
        plot_mean_diff_along_forecast(mean_stddev_dict) 
         

    elif fn == "plot_daily_composites":
    
        vmin = -20.0
        vmax = 20.0
        binsize = 2
        ncolors = ((vmax - vmin) / binsize)
        discrete_cmap = cm.get_cmap("RdBu_r", ncolors)
        
        if glob("corresponding_file_info.yml"):
            with open("corresponding_file_info.yml" , "r") as infile:        
                file_dict = yaml.full_load(infile)
        else:
            file_dict = create_gmi_fcst_time_lut()

        for cs in cycle_starts[case]:
            daily_data_composite_dict = {}
            substr_dict = cycle_start_regex.search(cs).groupdict()
            start_datetime = datetime.datetime(int(substr_dict["yyyy"]), int(substr_dict["mm"]), int(substr_dict["dd"]), int(substr_dict["hh"]), 0, 0)
            end_datetime = start_datetime + datetime.timedelta(days=1)
            #print(start_datetime.strftime(cycle_datetime_fmt))
            #print(end_datetime.strftime(cycle_datetime_fmt))
            print("Cycle Start:", cs)
            cycle_day_range = [dt.strftime(cycle_datetime_fmt) for dt in dt_range(start_datetime, end_datetime, datetime.timedelta(hours=3))]
            #print(cycle_day_range)
            for exp in ["Control", "AddMHS", "AddTEMPESTD"]:
                print("Experiment:", exp)
                for i, fcst_time in enumerate(cycle_day_range):
                    print("Forecast time:", fcst_time)
                    cycle_dir = os.path.join(crtm_gmi_data_dir, case, exp, cs)
                    fcst_file = glob(os.path.join(cycle_dir, "tb_gmi_89VGhz_" + fcst_time + "_*"))[0]
                    gmi_file = file_dict[fcst_time]
                    #print(gmi_file)
                    fcst_data = read_crtm_data(fcst_file, fill_val=0)
                    obs_tb_native = get_hdf5_data("/S1/Tc", gmi_file)
                    obs_tb89GhzV_native = obs_tb_native[:, :, 7]

                    obs_lat = get_hdf5_data("/S1/Latitude", gmi_file)
                    obs_lon = get_hdf5_data("/S1/Longitude", gmi_file)
                    #Convert -180 to 180 -> 0 to 360
                    obs_lon[np.where(obs_lon < 0.)] = obs_lon[np.where(obs_lon < 0.)] + 360.

                    #Regrid
                    gmi_data = fv3_gridify(obs_tb89GhzV_native, obs_lat, obs_lon)

                    fcst_data = np.ma.masked_where(crtm_landcover_mask, fcst_data)
                    gmi_data = np.ma.masked_where(crtm_landcover_mask, gmi_data)
                    gmi_data = np.ma.masked_where(gmi_data == -9999.9, gmi_data)


                    diff = fcst_data - gmi_data

                    if i == 0:
                        daily_data_composite = diff
                    else:
                        daily_data_composite = np.ma.dstack((daily_data_composite, diff))

                daily_data_composite_dict[exp] = daily_data_composite
   
            plot_daily_composites(daily_data_composite_dict, os.path.join(base_map_plot_dir, case, "ExperimentsMinusObs_DailyComposite_CycleStart" + cs + ".png"), central_lon=180., vmax=vmax, vmin=vmin, cmap=discrete_cmap, cbar_label="Tb (K)")

    
    else:
        print("Function " + fn + "is not supported. Exiting.")
        sys.exit()


