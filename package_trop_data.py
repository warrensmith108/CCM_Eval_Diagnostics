
import numpy as np
import pickle as pkl
import glob
import netCDF4 as nc4
import datetime as dt
import jdcal

OUTROOT = '/home/wsmith/STRATOCLIM/vars/'
TROPFILES = sorted(glob.glob('/data/wsmith/ERA5_tropopause/era5*'))
time_ref = dt.datetime(2000,1,1)

OUTFILE_z  = OUTROOT + 'trop_alt_file.pkl'
OUTFILE_p  = OUTROOT + 'trop_prs_file.pkl'
OUTFILE_t  = OUTROOT + 'trop_temp_file.pkl'
OUTFILE_th = OUTROOT + 'trop_theta_file.pkl'

# Read the first file and get the dimensions you need
TROPFILE0 = TROPFILES[0]
Tfile = nc4.Dataset(TROPFILE0)
trop_lon = Tfile.variables['lon'][:]
trop_lat = Tfile.variables['lat'][:]

# Create blank arrays
wmo_1st_z  = np.zeros((1,len(trop_lat),len(trop_lon)))-999
wmo_1st_p  = np.zeros((1,len(trop_lat),len(trop_lon)))-999
wmo_1st_t  = np.zeros((1,len(trop_lat),len(trop_lon)))-999
wmo_1st_th = np.zeros((1,len(trop_lat),len(trop_lon)))-999

trop_jd = []

for TROPFILE in TROPFILES:

    print('***************')
    print(TROPFILE)

    Tfile = nc4.Dataset(TROPFILE)
    curr_seconds = np.asarray(Tfile.variables['time'][::3])
    curr_trop_z  = np.asarray(Tfile.variables['wmo_1st_z'][::3,:,:])
    curr_trop_p  = np.asarray(Tfile.variables['wmo_1st_p'][::3,:,:])
    curr_trop_t  = np.asarray(Tfile.variables['wmo_1st_t'][::3,:,:])
    curr_trop_th = curr_trop_t*((1000./curr_trop_p)**(287./1004.))
    
    for t in range(len(curr_seconds)):
        curr_dt = time_ref+dt.timedelta(seconds=curr_seconds[t])
        curr_jd = np.sum(jdcal.gcal2jd(curr_dt.year,curr_dt.month,curr_dt.day))+(curr_dt.hour/24.)
        trop_jd.append(curr_jd)
   
    wmo_1st_z  = np.append(wmo_1st_z,  curr_trop_z,  axis=0)
    wmo_1st_p  = np.append(wmo_1st_p,  curr_trop_p,  axis=0)
    wmo_1st_t  = np.append(wmo_1st_t,  curr_trop_t,  axis=0)
    wmo_1st_th = np.append(wmo_1st_th, curr_trop_th, axis=0)
    
# Remove the first zeroed time slice
wmo_1st_z  = wmo_1st_z[1:,:,:]
wmo_1st_p  = wmo_1st_p[1:,:,:]
wmo_1st_t  = wmo_1st_t[1:,:,:]
wmo_1st_th = wmo_1st_th[1:,:,:]

# Save the arrays
with open(OUTFILE_z, 'wb') as f:
    pkl.dump((trop_lon, trop_lat, trop_jd, wmo_1st_z),f)
with open(OUTFILE_p, 'wb') as f:
    pkl.dump((trop_lon, trop_lat, trop_jd, wmo_1st_p),f)
with open(OUTFILE_t, 'wb') as f:
    pkl.dump((trop_lon, trop_lat, trop_jd, wmo_1st_t),f)
with open(OUTFILE_th, 'wb') as f:
    pkl.dump((trop_lon, trop_lat, trop_jd, wmo_1st_th),f)
    
# For general interest, print the ranges
print(np.nanmin(wmo_1st_z), np.nanmax(wmo_1st_z))
print(np.nanmin(wmo_1st_p), np.nanmax(wmo_1st_p))
print(np.nanmin(wmo_1st_t), np.nanmax(wmo_1st_t))
print(np.nanmin(wmo_1st_th), np.nanmax(wmo_1st_th))
    
    

