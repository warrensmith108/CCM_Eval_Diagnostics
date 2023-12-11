
#  The purpose of this script is to calculate tracer-tracer relationships for CESM variables
#  and overplot the STRATOCLIM points on top of them

import numpy as np
import matplotlib.pyplot as plt
import glob
import netCDF4 as nc4
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from read_STRATOCLIM_WAS import load_WAS
import pickle as pkl
import os
from scipy.io import readsav
import heatmap_config as hc
from read_STRATOCLIM_merge import get_merge_array
        
def heatmap_score(pixx, pixy, pixdata, var1_obs, var2_obs):
    # Currently, the approach for this is calculating the fraction of obs points that overlap the heatmap
    # At some point we may want to screen the heatmap for "small" values which might cause misleading results for this score

    var_flag = []
    for pt in range(len(var1_obs)):
        curr_var1 = var1_obs[pt]
        curr_var2 = var2_obs[pt]   
        
        # Only count a point for the score if it falls on the histogram boundaries in the first place
        if (curr_var1 >= pixx[0]) and (curr_var1 < pixx[-1]) and (curr_var2 >= pixy[0]) and (curr_var2 < pixy[-1]):
            
            # Determine the value of the pixel which the obs point lies within
            curr_pixx_bin = np.where(curr_var1 >= pixx)[0][-1]
            curr_pixy_bin = np.where(curr_var2 >= pixy)[0][-1]
            curr_pixdata = pixdata[curr_pixy_bin, curr_pixx_bin]
            
            if curr_pixdata > 0.: var_flag.append(1)
            else: var_flag.append(0)
    
    score = np.mean(var_flag)*100.
            
    return score

def get_SC_arrays(var1_inst, var1_str, var2_inst, var2_str, prs_bnds, WAS_SAVEDIR, MERGEFILES):
    
    # If one or both of the instruments is from WAS, use the WAS merge files
    if var1_inst == 'WAS' or var2_inst == 'WAS':      
        var1_SC_file = WAS_SAVEDIR + var1_inst + ':' + var1_str + '_merge.pkl'
        var2_SC_file = WAS_SAVEDIR + var2_inst + ':' + var2_str + '_merge.pkl'
        
        
        # var_data16, var_err_data16, var_data17, var_err_data17
        with open(var1_SC_file, 'rb') as f:
            var1_SC_16, var1_SC_err_16, var1_SC_17, var1_SC_err_17 = pkl.load(f)
        with open(var2_SC_file, 'rb') as f:
            var2_SC_16, var2_SC_err_16, var2_SC_17, var2_SC_err_17 = pkl.load(f)
    
    
    
    else:
        var1_SC_raw, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA = get_merge_array(var1_inst + ':' + var1_str, MERGEFILES)
        var2_SC_raw, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA = get_merge_array(var2_inst + ':' + var2_str, MERGEFILES)

        # Screen this larger set of points for the chosen pressure bounds
        var1_SC_17 = []
        var2_SC_17 = []
        var1_SC_err_17 = []
        var2_SC_err_17 = []        
        for pt in range(len(var1_SC_raw)):
            if AVIONIK_PRS[pt] >= prs_bnds[0] and AVIONIK_PRS[pt] <= prs_bnds[1]:  
                var1_SC_17.append(var1_SC_raw[pt])
                var2_SC_17.append(var2_SC_raw[pt])              
                var1_SC_err_17.append(0)
                var2_SC_err_17.append(0)
        var1_SC_16 = []
        var2_SC_16 = []      
        var1_SC_err_16 = []
        var2_SC_err_16 = []             
        
    return var1_SC_16, var2_SC_16, var1_SC_17, var2_SC_17, var1_SC_err_16, var2_SC_err_16, var1_SC_err_17, var2_SC_err_17
                
if __name__ == '__main__':

    sims = ['110L','M32L']
    fnt = 24

    #vars =    ['CO','N2O','CFC11','CFC113',  'CCL4','CFC115',  'CFC12','CH2BR2','CH3BR','CH3CCL3',  'CH3CL', \
    #            'CF2CLBR', 'HCFC141B','HCFC142B',  'HCFC22','CFC114','CF3BR','O3']
    #y_flips = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0]  # Controls whether the variables will be flipped on the yaxis or not
    
    vars = ['CO','N2O','O3']
    y_flips = [1,1,0]

    #lon_bnds = [0,180]
    #lat_bnds = [0,60]
    lon_bnds = [30,140]
    lat_bnds = [20,40]
    #lon_bnds = [257,277]   # For SEAC4RS
    #lat_bnds = [23,37]    
    #lon_bnds = [75,95]
    #lat_bnds = [18,32]
    #lon_bnds = [140,360]
    #lat_bnds = [0,60]
    #lon_bnds = [160,180]
    #lat_bnds = [0,20]
    prs_bnds = [50,200] # hPa
    
    for sim in sims: 
    
        CESM_PATH  = '/home/wsmith/forecasts/2017_' + sim + '/f*'
        ACE_PATH = '/home/wsmith/ACE_BackTraj/vars/'
        PLOTDIR = '/home/wsmith/STRATOCLIM/plots/' + sim + '/heatmaps/'
        MERGEFILES = sorted(glob.glob('/UTLS/Field_Data/STRATOCLIM_Data/merge_data/stratoclim*.nc'))
        WAS_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/WAS_merge/'
        DIST_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/heatmap_vars/'

        for v2 in range(len(vars)):
        #for v2 in range(2):                #####  Temporary!!! 
         
            var2_str = vars[v2]
            #var1_list = vars[v2+1:]    
            y_flip = y_flips[v2]
            var2_inst = hc.heatmap_info[var2_str]['instrument']
            var2_mod  = hc.heatmap_info[var2_str]['mod']
            var2_lims = hc.heatmap_info[var2_str]['lims']
            var2_bins = hc.heatmap_info[var2_str]['bins']
            var2_unit = hc.heatmap_info[var2_str]['unit']
            var2_SEAC4RS_index = hc.heatmap_info[var2_str]['SEAC4RS_index']
            
            for v1 in range(len(vars)):
                #var1_str = vars[v1+v2+1]    # Temporary to allow all instruments to plot against all others
                var1_str = vars[v1]
                
                if var2_str != var1_str: 
                
                    var1_inst = hc.heatmap_info[var1_str]['instrument']
                    var1_mod  = hc.heatmap_info[var1_str]['mod']
                    var1_lims = hc.heatmap_info[var1_str]['lims']
                    var1_bins = hc.heatmap_info[var1_str]['bins']
                    var1_unit = hc.heatmap_info[var1_str]['unit']
                    var1_SEAC4RS_index = hc.heatmap_info[var1_str]['SEAC4RS_index']
                    
                    DIST1_SAVEFILE = DIST_SAVEDIR + 'heatmap_' + sim + '_' + var1_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '-' + str(prs_bnds) + '.pkl'
                    DIST2_SAVEFILE = DIST_SAVEDIR + 'heatmap_' + sim + '_' + var2_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '-' + str(prs_bnds) + '.pkl'
                    
                    print(var1_str, var2_str)
                    
                    if os.path.exists(DIST1_SAVEFILE) and os.path.exists(DIST2_SAVEFILE):
                    
                        print('Heatmap data found... loading file...')
                        
                        with open(DIST1_SAVEFILE, 'rb') as f: var1 = pkl.load(f)
                        with open(DIST2_SAVEFILE, 'rb') as f: var2 = pkl.load(f)
                        
                    else:
                        print('Heatmap data not found... proceeding with calculation...')

                        #################  Read the information from the model
                        CESM_FILES = sorted(glob.glob(CESM_PATH))

                        var1 = []
                        var2 = []

                        for f in range(len(CESM_FILES)):
                            print('File ', f+1, ' of ', len(CESM_FILES))

                            CESM_file = nc4.Dataset(CESM_FILES[f])  

                            CESM_lon = CESM_file.variables['lon'][:]
                            CESM_lat = CESM_file.variables['lat'][:]
                            CESM_lev = CESM_file.variables['lev'][:]
                            var1_raw = CESM_file.variables[var1_str][:]*var1_mod
                            var2_raw = CESM_file.variables[var2_str][:]*var2_mod  
                            PS_raw = CESM_file.variables['PS'][:]
                            hyam     = CESM_file.variables['hyam'][:]
                            hybm     = CESM_file.variables['hybm'][:]

                            #prs_index0 = np.where(CESM_lev > prs_bnds[0])[0][0]
                            #prs_index1 = np.where(CESM_lev < prs_bnds[1])[0][-1]

                            if len(CESM_lon) > 5000.:   # If it is an unstructured grid, we need to handle it a certain way

                                for pt in range(len(CESM_lon)):
                                    
                                    curr_lon = CESM_lon[pt]
                                    curr_lat = CESM_lat[pt]
                                                                
                                    if curr_lon > lon_bnds[0] and curr_lon < lon_bnds[1] and curr_lat > lat_bnds[0] and curr_lat < lat_bnds[1]:                                
                                        curr_PS = np.mean(PS_raw[:,pt])
                                        curr_plev = (hyam*1.0e5 + hybm*curr_PS)/100.
                                        
                                        prs_index0 = np.where(curr_plev > prs_bnds[0])[0][0]
                                        prs_index1 = np.where(curr_plev < prs_bnds[1])[0][-1]
                                        
                                        var1.append(np.squeeze(var1_raw[:,prs_index0:prs_index1+1,pt]).flatten())
                                        var2.append(np.squeeze(var2_raw[:,prs_index0:prs_index1+1,pt]).flatten())
                                        
                                        #var1.append(np.squeeze(var1_raw[:,:,pt]).flatten())
                                        #var2.append(np.squeeze(var2_raw[:,:,pt]).flatten())

                            else:  # If this is a regular grid, we also handle it a certain way, but a different way, because reasons

                                lon_index0 = np.argmin(np.abs(CESM_lon - lon_bnds[0]))
                                lon_index1 = np.argmin(np.abs(CESM_lon - lon_bnds[1]))

                                lat_index0 = np.argmin(np.abs(CESM_lat - lat_bnds[0]))
                                lat_index1 = np.argmin(np.abs(CESM_lat - lat_bnds[1]))

                                for ln in range(len(CESM_lon)):
                                    curr_lon = CESM_lon[ln]
                                    for lt in range(len(CESM_lat)):
                                        curr_lat = CESM_lat[lt]
                                        if curr_lon > lon_bnds[0] and curr_lon < lon_bnds[1] and curr_lat > lat_bnds[0] and curr_lat < lat_bnds[1]:  
                                            curr_PS = np.mean(PS_raw[:,lt,ln])
                                            curr_plev = (hyam*1.0e5 + hybm*curr_PS)/100.

                                            prs_index0 = np.where(curr_plev > prs_bnds[0])[0][0]
                                            prs_index1 = np.where(curr_plev < prs_bnds[1])[0][-1]
                                            
                                            var1.append(np.squeeze(var1_raw[:,prs_index0:prs_index1+1,lt,ln]).flatten())
                                            var2.append(np.squeeze(var2_raw[:,prs_index0:prs_index1+1,lt,ln]).flatten())
                                        
                                #var1.append(CESM_file.variables[var1_str][:,prs_index0:prs_index1+1,lat_index0:lat_index1+1,lon_index0:lon_index1+1].flatten()*var1_mod)
                                #var2.append(CESM_file.variables[var2_str][:,prs_index0:prs_index1+1,lat_index0:lat_index1+1,lon_index0:lon_index1+1].flatten()*var2_mod)

                        var1 = np.asarray(var1).flatten()
                        var2 = np.asarray(var2).flatten()
                        
                        # Save the variables to have the ability to plot later
                        
                        with open(DIST1_SAVEFILE, 'wb') as f: pkl.dump((var1), f)
                        with open(DIST2_SAVEFILE, 'wb') as f: pkl.dump((var2), f)
                    
                    ################  Read the information from the STRATOCLIM data

                    var1_SC = []
                    var2_SC = []
                    for curr_var1_inst in var1_inst:
                        for curr_var2_inst in var2_inst:            
                            curr_var1_SC, curr_var2_SC = get_SC_arrays(curr_var1_inst, var1_str, curr_var2_inst, var2_str, prs_bnds, WAS_SAVEDIR, MERGEFILES)
                            var1_SC.extend(curr_var1_SC)
                            var2_SC.extend(curr_var2_SC)
                            
                    var1_SC = np.asarray(var1_SC)
                    var2_SC = np.asarray(var2_SC)
                    
                    ################  Read in SEAC4RS data, if necessary
                    
                    DATAFILE = '/home/wsmith/STRATOCLIM/SEAC4RS_Data/SEAC4RS-mrg60-ER2-merge.txt'
                    if var1_SEAC4RS_index != -999 and var2_SEAC4RS_index != -999:
                        SEAC4RSdata  = np.loadtxt(DATAFILE, delimiter=',', skiprows=218)
                        var1_SEAC4RS = SEAC4RSdata[:,var1_SEAC4RS_index]
                        var2_SEAC4RS = SEAC4RSdata[:,var2_SEAC4RS_index]
                    else:
                        var1_SEAC4RS = [-999]
                        var2_SEAC4RS = [-999]
                        
                    ################  Read in the ACE data as another layer, if it exists
                    
                    #ACEFILE1 = ACE_PATH + var1_str + '_[2017072000, 2017081000]_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(prs_bnds) + '.pkl'
                    #ACEFILE2 = ACE_PATH + var2_str + '_[2017072000, 2017081000]_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(prs_bnds) + '.pkl'       
                    
                    #if os.path.exists(ACEFILE1):
                    #    with open(ACEFILE1, 'rb') as f:
                    #        var1_ACE, date, lon, lat = pkl.load(f)
                    #else:
                    #    var1_ACE = [-999]
                        
                    #if os.path.exists(ACEFILE2):
                    #    with open(ACEFILE2, 'rb') as f:
                    #        var2_ACE, date, lon, lat = pkl.load(f)
                    #else:
                    #    var2_ACE = [-999]            
                        
                    ################  Calculate 2d histogram (heatmap)
                    pixels = np.histogram2d(var1, var2, range = [var1_lims,var2_lims], bins = [var1_bins,var2_bins]) 
                    pixx = pixels[1]
                    pixy = pixels[2]
                    pixdata = np.transpose(pixels[0])
                    pixdata = pixdata*100./np.max(pixdata)  # Normalize
                    pixdata[pixdata == 0.0] = pixdata[pixdata == 0.0]-0.01
                    
                    ################  Score the obs against the heatmap
                    #SC_score = heatmap_score(pixx, pixy, pixdata, var1_SC, var2_SC)
                    #if len(var1_ACE) > 3 and len(var2_ACE) > 3: ACE_score = heatmap_score(pixx, pixy, pixdata, var1_ACE, var2_ACE)
                    #print('StratoClim Score = ', SC_score)

                    ################  Plot the heatmap
                    cmap = plt.get_cmap('YlOrRd')
                    levels = MaxNLocator(nbins=10).tick_values(0.0, 100.0)
                    norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
                    cmap.set_under([0.8,0.8,0.8])

                    fig = plt.figure(figsize = (10,8))
                    ax = fig.add_subplot(1,1,1)

                    pixplt = ax.pcolormesh(pixx, pixy, pixdata, cmap = cmap, norm = norm)
                    cbar = plt.colorbar(pixplt, cmap = cmap, norm = norm)
                    cbar.set_label('Relative Frequency (%)', fontsize = fnt)
                    for tick in cbar.ax.yaxis.get_ticklabels():
                        tick.set_weight('medium')
                        tick.set_fontsize(fnt)

                    # Overlay the WAS observations on it
                    #plt.scatter(var1_WAS16, var2_WAS16, color = 'cyan')
                    if np.min([len(var1_SC[var1_SC >= 0]), len(var2_SC[var2_SC >= 0])]) > 10000.: plt.scatter(var1_SC[::100], var2_SC[::100], s=20, color='black')
                    else: plt.scatter(var1_SC, var2_SC, s=20, color='black')
                    
                    # Overlay the SEAC4RS observations
                    ax.scatter(var1_SEAC4RS, var2_SEAC4RS, color='blue', facecolor='none', s=0.5)            
                    
                    # Overlay the ACE observations as well, if they exist
                    #if len(var1_ACE) > 3 and len(var2_ACE) > 3:
                        #plt.scatter(var1_ACE, var2_ACE, color = 'black', s = 5, zorder=999)
                        #plt.scatter(var1_ACE, var2_ACE, color = 'black', s = 5, zorder=999)
                        
                    ax.set_xlim(var1_lims)
                    if y_flip == 1: ax.set_ylim(np.flip(var2_lims))
                    else: ax.set_ylim(var2_lims)

                    ax.set_xlabel(var1_str + ' (' + var1_unit + ')', fontsize = fnt)
                    ax.set_ylabel(var2_str + ' (' + var2_unit + ')', fontsize = fnt)
                    
                    plt.figtext(0.0,0.05, 'Model: ' + sim, fontsize = 10, color = 'orange', backgroundcolor='white')
                    plt.figtext(0.0,0.03, var1_str + ': ' + str(var1_inst), fontsize = 10, color = 'black', backgroundcolor='white')
                    plt.figtext(0.0,0.01, var2_str + ': ' + str(var2_inst), fontsize = 10, color = 'black', backgroundcolor='white')
                    #plt.figtext(0.2,0.27, 'StratoClim Score: ' + str(round(SC_score, 1)), fontsize = 10, color = 'blue', backgroundcolor='white')
                    #if len(var1_ACE) > 3 and len(var2_ACE) > 3: plt.figtext(0.6,0.24, 'ACE Score: ' + str(round(ACE_score, 1)), fontsize = 10, color = 'black', backgroundcolor='white')
                    
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label1.set_fontweight('medium')
                        tick.label1.set_fontsize(fnt)
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label1.set_fontweight('medium')
                        tick.label1.set_fontsize(fnt)

                    plt.tight_layout()
                    plt.savefig(PLOTDIR + 'heatmap_' + var2_str + 'vs' + var1_str + str(lon_bnds) + '-' + str(lat_bnds) + '.png')
                    plt.close('all')

                    
                    









