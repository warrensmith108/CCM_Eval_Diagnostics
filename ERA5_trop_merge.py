

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from read_STRATOCLIM_merge import get_merge_array
import glob
import netCDF4 as nc4
import jdcal
from scipy.interpolate import RegularGridInterpolator as RGI
from read_STRATOCLIM_WAS import load_WAS
import datetime as dt

PLOTDIR = '/home/wsmith/STRATOCLIM/plots/'
SAVEDIR = '/home/wsmith/STRATOCLIM/vars/'
WAS_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/WAS_merge/'
DATAFILES = sorted(glob.glob('/UTLS/Field_Data/STRATOCLIM_Data/merge_data/stratoclim*.nc'))
TROPFILES = sorted(glob.glob('/data/wsmith/ERA5_tropopause/e*'))
trop_choice = 'wmo_1st_z'
#trop_choice = 'clp_z'

# Read an example merge array to get the StratoClim flight track information
var_data, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_merge_array('FOZAN:O3', DATAFILES)

# Read the tropopause files
file0_str = TROPFILES[0]

curr_year  = int(file0_str[-13:-9])
curr_month = int(file0_str[-8:-6])
curr_day   = int(file0_str[-5:-3])
jd_array = list(np.sum(jdcal.gcal2jd(curr_year,curr_month,curr_day)) + (np.arange(24)/24.))


file0 = nc4.Dataset(file0_str)  # First file

lon = file0.variables['lon'][:]
lat = file0.variables['lat'][:]
trop = file0.variables[trop_choice][:]

print('Reading the tropopause...')
for f in range(len(TROPFILES)-1):

    curr_file = TROPFILES[f+1]

    curr_year  = int(curr_file[-13:-9])
    curr_month = int(curr_file[-8:-6])
    curr_day   = int(curr_file[-5:-3])
    curr_jd_array = list(np.sum(jdcal.gcal2jd(curr_year,curr_month,curr_day)) + (np.arange(24)/24.))
    
    jd_array.extend(curr_jd_array)
    
    file = nc4.Dataset(curr_file)
    
    curr_trop = file.variables[trop_choice][:]
    
    trop = np.append(trop, curr_trop, axis=0)
    
# Create the interpolator function
f_tropz = RGI((jd_array,lat,lon), trop, bounds_error=False, fill_value=-999, method='linear')

trop_array = []
jd_array = []

# Now loop over all the aircraft points and do the interpolation
for pt in range(len(AVIONIK_LAT)):
    curr_lon = AVIONIK_LON[pt]
    curr_lat = AVIONIK_LAT[pt]
    
    curr_datestr = str(dateint[pt])
    curr_year = int(curr_datestr[0:4])
    curr_month = int(curr_datestr[4:6])
    curr_day = int(curr_datestr[6:8])

    curr_jd = np.sum(jdcal.gcal2jd(curr_year, curr_month, curr_day)) + int(curr_datestr[8:10])/24. + int(curr_datestr[10:12])/1440. + int(curr_datestr[12:])/86400.
    
    curr_trop = f_tropz((curr_jd, curr_lat, curr_lon))
    
    trop_array.append(curr_trop)
    jd_array.append(curr_jd)

print(np.min(trop_array), np.max(trop_array))

print(len(AVIONIK_LAT))
print(len(trop_array))

# Save the merged tropopause   
with open(SAVEDIR + trop_choice + '_merge.pkl','wb') as f:
    pkl.dump((trop_array),f)
    



######################  WE ALSO NEED A SEPARATE MERGE FOR THE WAS POINTS



all_data, var_headers = load_WAS(r'/home/wsmith/STRATOCLIM/STRATOCLIM_WAS_Times.xlsx')

merge_WAS_arr = []  # Create a blank array to be filled with the merged array
nWAS = (np.shape(all_data)[1]-1)/2  # The number of canisters
#print(nWAS)

for msr in range(int(nWAS)):
    
    print('************', msr)
    
    curr_ind = msr*2+1    # 1, 3, 5, etc
    
    curr_open_date    = all_data[1][curr_ind]
    curr_open_time    = all_data[2][curr_ind]            
    curr_open_dt      = curr_open_date + dt.timedelta(hours = curr_open_time.hour, minutes = curr_open_time.minute, seconds = curr_open_time.second)
    curr_open_jd      = sum(jdcal.gcal2jd(curr_open_dt.year, curr_open_dt.month, curr_open_dt.day)) + curr_open_dt.hour/24. + curr_open_dt.minute/1440. + curr_open_dt.second/86400.
    #curr_open_dateint = int(str(curr_open_dt.year) + str(curr_open_dt.month).zfill(2) + str(curr_open_dt.day).zfill(2) + \
    #                   str(curr_open_dt.hour).zfill(2) + str(curr_open_dt.minute).zfill(2) + str(curr_open_dt.second).zfill(2))

    curr_close_date    = all_data[1][curr_ind+1]
    curr_close_time    = all_data[2][curr_ind+1]
    curr_close_dt      = curr_close_date + dt.timedelta(hours = curr_close_time.hour, minutes = curr_close_time.minute, seconds = curr_close_time.second)
    curr_close_jd      = sum(jdcal.gcal2jd(curr_close_dt.year, curr_close_dt.month, curr_close_dt.day)) + curr_close_dt.hour/24. + curr_close_dt.minute/1440. + curr_close_dt.second/86400.    
    #curr_close_dateint = int(str(curr_close_dt.year) + str(curr_close_dt.month).zfill(2) + str(curr_close_dt.day).zfill(2) + \
    #                     str(curr_close_dt.hour).zfill(2) + str(curr_close_dt.minute).zfill(2) + str(curr_close_dt.second).zfill(2))                    
    
    curr_WAS_arr = []  # This will be filled with all the merge data falling within each open/close, to be averaged at the end        
    #curr_prs_arr = []
    for merge_t in range(len(dateint)):
        curr_jd = jd_array[merge_t]
        
        #print(curr_open_jd, curr_close_jd, curr_jd)
        
        if curr_jd >= curr_open_jd and curr_jd <= curr_close_jd:
        
    #        curr_prs_arr.append(AVIONIK_PRS[merge_t])
            curr_WAS_arr.append(trop_array[merge_t])
            
    
    
    print(np.nanmean(curr_WAS_arr))
    
    #merge_prs_arr.append(np.nanmean(curr_prs_arr))
    merge_WAS_arr.append(np.nanmean(curr_WAS_arr))
    

savefile = WAS_SAVEDIR + 'ERA5_' + trop_choice + '_merge.pkl'    
#print(len(merge_WAS_arr))
with open(savefile, 'wb') as f:
    pkl.dump(([], [], merge_WAS_arr, []), f)
            
