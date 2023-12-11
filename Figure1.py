
############################################
#
# This script is intended to create a map plot of the STRATOCLIM observation
# locations with mean UT GPH overplotted
#
############################################

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import glob
import os
from read_STRATOCLIM_WAS import load_WAS
from read_STRATOCLIM_merge import get_merge_array
import cartopy as cp
import cartopy.crs as ccrs
import pygrib
import matplotlib
import scipy
import scipy.ndimage
import scale_hgt as sh

def moving_window_smooth(data, width=1):
    data_copy = np.copy(data)
    #smoothed = sp.ndimage.filters.gaussian_filter(data_copy, [2,2], mode='constant')
    smoothed = np.zeros(data_copy.shape)
    for ix in range(data_copy.shape[0]):
        for iy in range(data_copy.shape[1]):
            smoothed[ix,iy] = np.nanmean(data_copy[ix-width:ix+width,iy-width:iy+width])
    return smoothed   

colors = ['cyan','blue']  # 2016, 2017
WASFILE = r'/home/wsmith/STRATOCLIM/STRATOCLIM_WAS_Data.xlsx'
PLOTDIR = '/home/wsmith/STRATOCLIM/plots/Raw_Data/'
GFSDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_RDF/GFS_rolling_archive/'
AIRFILES = sorted(glob.glob('/UTLS/Field_Data/STRATOCLIM_Data/merge_data/stratoclim*.nc'))
GFSrng16 = ['0901_00','0907_18']
GFSrng17 = ['0727_00','0810_00']
GFSvar_ind = [69,70]  # 78 is for 150 hPa
GFSvar_cont = [16770]
#GFSvar_ind = [280,281]
#GFSvar_cont = [10000.]
GFSPS_ind = 247    # From the 2017 files
GFSTROP_ind = 281  # From the 2017 files
GFSTROP_latavg = [18.,32.]  # Lat boundaries for averaging the tropopause (for XZ plot panel)
#COpts = 1  # This controls whether we will plot the locations where CO/O3 were observed, or WAS points instead
box1 = [30,130,18,40]
box2 = [75,95,18,32]

########  Read the Geophysica locations

AIRlon   = get_merge_array('AVIONIK:LON', AIRFILES)[0]
AIRlat   = get_merge_array('AVIONIK:LAT', AIRFILES)[0]
AIRprs   = get_merge_array('AVIONIK:PRESS', AIRFILES)[0]
AIRALT   = get_merge_array('AVIONIK:ALT', AIRFILES)[0]
AIRTHETA = get_merge_array('AVIONIK:THETA', AIRFILES)[0]

AIRlon = np.asarray(AIRlon)
AIRlat = np.asarray(AIRlat)
AIRprs = np.asarray(AIRprs)
AIRALT = np.asarray(AIRALT)*1.e-3
AIRTHETA = np.asarray(AIRTHETA)

# Old version using the txt file:
#datafile = '/home/wsmith/STRATOCLIM/STRATOCLIM_2017_Kathmandu.txt'
#rawdata = np.loadtxt(datafile, delimiter=',', skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
#latCO = rawdata[:,6]
#lonCO = rawdata[:,7]
#prsCO = rawdata[:,9]

########  Read the WAS data

all_data, var_headers = load_WAS(WASFILE)
mission = np.array(all_data[1][1:])

lat16 = np.array(all_data[4][1:])[mission == 1]
lon16 = np.array(all_data[5][1:])[mission == 1]
alt16 = np.array(all_data[6][1:])[mission == 1]/1.0e3
prs16 = np.array(all_data[8][1:])[mission == 1]

lat17 = np.array(all_data[4][1:])[mission == 2]
lon17 = np.array(all_data[5][1:])[mission == 2]
alt17 = np.array(all_data[6][1:])[mission == 2]/1.0e3
prs17 = np.array(all_data[8][1:])[mission == 2]
theta17 = np.array(all_data[9][1:])[mission == 2]

########  Read the GFS data

GFSFILES16 = sorted(glob.glob(GFSDIR + 'gfs_3_2016*.grb2'))
GFS16ind0 = GFSFILES16.index(GFSDIR + 'gfs_3_2016' + GFSrng16[0] + '00_000.grb2')
GFS16ind1 = GFSFILES16.index(GFSDIR + 'gfs_3_2016' + GFSrng16[1] + '00_000.grb2')
GFSFILES16 = GFSFILES16[GFS16ind0:GFS16ind1+1]

GFSFILES17 = sorted(glob.glob(GFSDIR + 'gfs_3_2017*.grb2'))
GFS17ind0 = GFSFILES17.index(GFSDIR + 'gfs_3_2017' + GFSrng17[0] + '00_000.grb2')
GFS17ind1 = GFSFILES17.index(GFSDIR + 'gfs_3_2017' + GFSrng17[1] + '00_000.grb2')
GFSFILES17 = GFSFILES17[GFS17ind0:GFS17ind1+1]

rawdata16 = pygrib.open(GFSFILES16[0])  # Read the first file to determine shape
GFSvar16 = rawdata16.select()[GFSvar_ind[0]]
latlon = GFSvar16.latlons()
GFSlat = latlon[0][:,0]
GFSlon = latlon[1][0,:]
GFSvar16 = GFSvar16.values[None,:]

rawdata17 = pygrib.open(GFSFILES17[0])  # Read the first file to determine shape
GFSvar17  = rawdata17.select()[GFSvar_ind[1]].values[None,:]
GFSPS = rawdata17.select()[GFSPS_ind].values[None,:]
#GFSTROP = rawdata17.select()[GFSTROP_ind].values[None,:]

#for f in range(len(GFSFILES16)-1):
#    rawdata16 = pygrib.open(GFSFILES16[f+1])
#    GFSvar16 = np.append(GFSvar16, rawdata16.select()[GFSvar_ind[0]].values[None,:], axis=0)
    
for f in range(len(GFSFILES17)-1):
    rawdata17 = pygrib.open(GFSFILES17[f+1])
    GFSvar17 = np.append(GFSvar17, rawdata17.select()[GFSvar_ind[1]].values[None,:], axis=0)
    GFSPS    = np.append(GFSPS,    rawdata17.select()[GFSPS_ind].values[None,:], axis=0)
#    GFSTROP  = np.append(GFSTROP,  rawdata17.select()[GFSTROP_ind].values[None,:], axis=0)
    
#GFSvar16 = np.mean(GFSvar16, axis=0)
GFSvar17 = np.mean(GFSvar17, axis=0)
GFSPS = np.mean(GFSPS, axis=0)/100.
#GFSTROP = np.mean(GFSTROP, axis=0)/100.
#print(GFSTROP)

# Restrict the GFS tropopause to the lat bounds defined above
#GFSTROP_latind0 = list(GFSlat).index(GFSTROP_latavg[1])
#GFSTROP_latind1 = list(GFSlat).index(GFSTROP_latavg[0])

#GFSTROP_avg = []
#for ln in range(len(GFSlon)):
#    curr_GFSTROP = GFSTROP[GFSTROP_latind0:GFSTROP_latind1+1,ln]
#    #print(curr_GFSTROP)
#    GFSTROP_avg.append(np.nanmean(curr_GFSTROP))

########  Read the ERA5 tropopause information

with open('/home/wsmith/STRATOCLIM/vars/trop_alt_file.pkl', 'rb') as f:
    ERA5_trop_lon, ERA5_trop_lat, ERA5_trop_jd, ERA5_trop_alt = pkl.load(f) 
with open('/home/wsmith/STRATOCLIM/vars/trop_theta_file.pkl', 'rb') as f:
    ERA5_trop_lon, ERA5_trop_lat, ERA5_trop_jd, ERA5_trop_th = pkl.load(f) 
trop_lat_ind0 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-box2[2]))
trop_lat_ind1 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-box2[3]))

#ERA5_trop_prs_mapplot = np.mean(ERA5_trop_prs, axis=0)
ERA5_trop_alt_lonplot = np.mean(ERA5_trop_alt[:,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))
ERA5_trop_th_lonplot  = np.mean( ERA5_trop_th[:,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))


########  Make the plot

fig = plt.figure(figsize = (18,10))
ax1 = fig.add_subplot(1,4,(1,2),projection = ccrs.PlateCarree())
#ax2 = fig.add_subplot(221)  # These two separate axes are for including 2016 measurements in a separate panel
#ax3 = fig.add_subplot(222)
ax3 = fig.add_subplot(143)
ax4 = fig.add_subplot(144)
#ax3.fill_between([0,360],[50,50],[200,200], color=[0.8,0.8,0.8])

ax3b = ax3.twinx()

#ax1.scatter(lon16, lat16, color = colors[0])
#if COpts == 0: ax1.scatter(lon17, lat17, color = 'blue')
#else: ax1.scatter(lonCO, latCO, color = 'blue', s=2)
#ax1.scatter(lonCO, latCO, color = 'black', s=1)
#ax1.scatter(lon17, lat17, color = 'blue', s=4)
ax1.scatter(AIRlon[::100], AIRlat[::100], color='black', s=1)

#ax1.contour(gphlon, gphlat, gph16, levels = gph_cont, colors=colors[0])
ax1.contour(GFSlon, GFSlat, GFSvar17, levels = GFSvar_cont, colors='purple', linewidths = 3)
#ERA5_trop_prs_mapplot = moving_window_smooth(ERA5_trop_prs_mapplot, width=2)
#ax1.contour(ERA5_trop_lon, ERA5_trop_lat, ERA5_trop_prs_mapplot, levels = [100.], colors = 'purple', linewidths = 2) 
ax1.coastlines()

def plot_box(ax, box):
    ax.plot([box[0],box[0]],[box[2],box[3]], color = 'gray', linewidth=3)
    ax.plot([box[1],box[1]],[box[2],box[3]], color = 'gray', linewidth=3)
    ax.plot([box[0],box[1]],[box[2],box[2]], color = 'gray', linewidth=3)
    ax.plot([box[0],box[1]],[box[3],box[3]], color = 'gray', linewidth=3)

#plot_box(ax1, box1)
plot_box(ax1, box2)

GFSPS = moving_window_smooth(GFSPS, width=1)
ax1.contour(GFSlon, GFSlat, GFSPS, [609.], colors=['gray'], linewidths = 0.5)

#ax2.scatter(lon16, prs16, color = colors[0])
#ax2.plot([0,180],[100,100], color = 'gray', linestyle = '--')
#ax3.scatter(lonCO, prsCO, color='black', s=1)

ax3b.plot(ERA5_trop_lon, ERA5_trop_alt_lonplot, color='gray', linewidth=2, linestyle='--')
ax3.scatter(AIRlon[::20], AIRprs[::20], color='black', s=1, zorder=999)
ax3.scatter(lon17, prs17, color=colors[1], zorder=9999)

#ax3.plot([0,180],[100,100], color = 'gray', linestyle = '--')  # Horizontal line at 100 hPa
#ax3.plot(GFSlon, GFSTROP_avg, color='gray',linestyle='--')
#ax3b.plot([0,180],[17.25,17.25], color='gray', linestyle='--')  # Adding the model tropopause (see Fig 3) because it's close enough

ax4.plot(ERA5_trop_lon, ERA5_trop_th_lonplot, color='gray', linewidth=2, linestyle='--')
ax4.scatter(AIRlon[::20], AIRTHETA[::20], color='black', s=1, zorder=999)
ax4.scatter(lon17, theta17, color=colors[1], zorder=9999)
#ax4.plot([0,180],[389,389], color='gray', linestyle='--')  # Adding the model tropopause (see Fig 3) because it's close enough


#ax2.set_yscale('log')
ax3.set_yscale('log')
ax1.set_xlim(15,150)
ax1.set_ylim(0,60)
#ax2.set_xlim(22,32)
#ax2.set_ylim(320,45)
ax3.set_xlim(78,92)
ax4.set_xlim(78,92)
ax3.set_ylim(1000,40)
ax3.set_yticks([])
#ax2.set_xlabel('Longitude (E)', fontsize = 12)
#ax3.set_xlabel('Longitude (E)', fontsize = 12)
ax3.set_ylabel('Pressure (hPa)', fontsize = 12)
ax3b.set_ylabel('Altitude (km)', fontsize = 12)
ax4.set_ylabel('Theta (K)', fontsize=12)
#ax2.set_yticks([60,100,200,300])
ax3.set_yticks([40,100,400,1000])
ax3.set_xticks([80,85,90])
ax4.set_xticks([80,85,90])

palt_lim = sh.prs2alt(np.asarray([1000,40]), 7.6, 1000.0)
ax3b.set_ylim(palt_lim)  # Convert to kft

ax3.tick_params(labelsize=12)
#ax3.set_yticklabels([])
ax3.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax1.set_position([0.02,0.35,0.4,0.29])

ax3pos = [0.45,0.35,0.23,0.29]
ax3.set_position(ax3pos)
ax3b.set_position(ax3pos)

ax4.set_position([0.76,0.35,0.23,0.29])

deg = u"\N{DEGREE SIGN}"
gl1 = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, color='gray', alpha=0.5, linestyle='--', \
xlocs = [30,60,90,120,150], ylocs = [0,20,40,60])
ax1.set_xticks([30,60,90,120,150])
ax1.set_xticklabels(['30'+deg+'E','60'+deg+'E','90'+deg+'E','120'+deg+'E','150'+deg+'E'], fontsize = 12)
ax1.set_yticks([0,20,40,60])
ax3.set_xticklabels(['80'+deg+'E', '85'+deg+'E', '90'+deg+'E'])
ax4.set_xticklabels(['80'+deg+'E', '85'+deg+'E', '90'+deg+'E'])

ax1.set_yticklabels(['0'+deg,'20'+deg+'N','40'+deg+'N','60'+deg+'N'], fontsize = 12) 

#plt.figtext(0.16, 0.16, str(GFSrng16), color=colors[0])
#plt.figtext(0.16, 0.15, 'Date Range = ' + str(GFSrng17), color=colors[1])
#plt.figtext(0.16, 0.17, 'GPH = ' + str(gph_cont) + 'm')

plt.figtext(0.06, 0.37, '(a)', fontsize=28, fontweight='bold')
plt.figtext(0.48, 0.37, '(b)', fontsize=28, fontweight='bold')
plt.figtext(0.78, 0.37, '(c)', fontsize=28, fontweight='bold')

plt.savefig(PLOTDIR + 'Figure1.png', dpi=300)
plt.savefig(PLOTDIR + 'Figure1.pdf', dpi=300)
plt.close('all')




