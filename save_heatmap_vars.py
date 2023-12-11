
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
from scipy.interpolate import interp1d
from joblib import Parallel, delayed

def save_heatmap_var(CESM_PATH, DIST_SAVEDIR, var_str, sim, lon_bnds, lat_bnds, prs_bnds):

    #################  Read the information from the model
    CESM_FILES = sorted(glob.glob(CESM_PATH))

    var_mod = hc.heatmap_info[var_str]['mod']
    
    var_asm   = []
    var_trop  = []
    var_strat = []

    for f in range(len(CESM_FILES)):
        print('File ', f+1, ' of ', len(CESM_FILES))

        CESM_file = nc4.Dataset(CESM_FILES[f])  

        CESM_lon = CESM_file.variables['lon'][:]
        CESM_lat = CESM_file.variables['lat'][:]
        CESM_lev = CESM_file.variables['lev'][:]
        var_raw = CESM_file.variables[var_str][:]*var_mod
        TROP_raw = CESM_file.variables['TROP_P'][:]/100.
        GPH_raw  = CESM_file.variables['Z3'][:]
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

                    for t in range(np.shape(var_raw)[0]):
                    
                        curr_GPH_prof = GPH_raw[t,:,pt]
                        curr_TROP = TROP_raw[t,pt]
                                                                                                                                            
                        curr_PS = PS_raw[t,pt]
                        curr_plev = (hyam*1.0e5 + hybm*curr_PS)/100.

                        prs_index0 = np.where(curr_plev > prs_bnds[0])[0][0]
                        prs_index1 = np.where(curr_plev < prs_bnds[1])[0][-1]
                        
                        f_GPH = interp1d(curr_plev, curr_GPH_prof)
                        curr_GPH_150 = f_GPH(150.)
                        
                        if curr_GPH_150 >= 14350.: var_asm.extend(np.squeeze(var_raw[t,prs_index0:prs_index1+1,pt]).flatten())
                        elif curr_TROP <= 150.: var_trop.extend(np.squeeze(var_raw[t,prs_index0:prs_index1+1,pt]).flatten())
                        else: var_strat.extend(np.squeeze(var_raw[t,prs_index0:prs_index1+1,pt]).flatten())                                       


        else:  # If this is a regular grid, we also handle it a certain way, but a different way, because reasons

            for ln in range(len(CESM_lon)):
                curr_lon = CESM_lon[ln]
                for lt in range(len(CESM_lat)):
                    curr_lat = CESM_lat[lt]
                    if curr_lon > lon_bnds[0] and curr_lon < lon_bnds[1] and curr_lat > lat_bnds[0] and curr_lat < lat_bnds[1]:  
                    
                        for t in range(np.shape(var_raw)[0]):
                        
                            curr_GPH_prof = GPH_raw[t,:,lt,ln]
                            curr_TROP = TROP_raw[t,lt,ln]
                            
                            curr_PS = PS_raw[t,lt,ln]
                            curr_plev = (hyam*1.0e5 + hybm*curr_PS)/100.

                            prs_index0 = np.where(curr_plev > prs_bnds[0])[0][0]
                            prs_index1 = np.where(curr_plev < prs_bnds[1])[0][-1]
                            
                            f_GPH = interp1d(curr_plev, curr_GPH_prof)
                            curr_GPH_150 = f_GPH(150.)
                            
                            if curr_GPH_150 >= 14350.: var_asm.extend(np.squeeze(var_raw[t,prs_index0:prs_index1+1,lt,ln]).flatten())
                            elif curr_TROP <= 150.: var_trop.extend(np.squeeze(var_raw[t,prs_index0:prs_index1+1,lt,ln]).flatten())
                            else: var_strat.extend(np.squeeze(var_raw[t,prs_index0:prs_index1+1,lt,ln]).flatten())                                                 
                            

    var_asm   = np.asarray(var_asm)
    var_trop  = np.asarray(var_trop)
    var_strat = np.asarray(var_strat)

    # Save the variables to have the ability to plot later

    DIST_SAVEFILE = DIST_SAVEDIR + 'heatmap_' + sim + '_' + var_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '-' + str(prs_bnds) + '.pkl'
    with open(DIST_SAVEFILE, 'wb') as f: pkl.dump((var_asm, var_trop, var_strat), f)
    
    
    
    
    
    
    
#####################  MAIN PROGRAM
    
sim = '110L'
vars =    ['CO','N2O','CFC11','CFC113', 'CH3CCL3', 'CCL4','CFC115',  'CFC12','CH2BR2','CH3BR',  'CH3CL', \
            'CF2CLBR', 'HCFC141B','HCFC142B',  'HCFC22','CFC114','CF3BR','O3']   
    
CESM_PATH  = '/data/wsmith/forecasts/2017_' + sim + '/f*'
DIST_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/heatmap_vars/'

#lon_bnds = [0,180]
#lat_bnds = [0,60]
#lon_bnds = [30,130]
#lat_bnds = [18,40]
lon_bnds = [75,95]
lat_bnds = [18,32]
prs_bnds = [40,200] # hPa

Parallel(n_jobs=10)(delayed(save_heatmap_var)(CESM_PATH, DIST_SAVEDIR, var_str, sim, lon_bnds, lat_bnds, prs_bnds) for var_str in vars)







