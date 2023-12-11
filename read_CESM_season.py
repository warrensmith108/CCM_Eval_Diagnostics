

import numpy as np
import glob
import os
import netCDF4 as nc4
import pickle as pkl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.tri as tri
import math
from scipy.interpolate import interp1d

def get_unstruc_var(FILES, var, prs_lvls, troprelz_lvls, triang, Xi, Yi):
    
    file0 = nc4.Dataset(FILES[0])
    var_raw = file0.variables[var][:]
    hyam    = file0.variables['hyam'][:]
    hybm    = file0.variables['hybm'][:]
    GPH_raw = file0.variables['Z3'][:]
    PS_raw  = file0.variables['PS'][:]
    TROP_P_raw = file0.variables['TROP_P'][:]
    
    for f in range(len(FILES)-1):    
        curr_file = nc4.Dataset(FILES[f+1])
        var_raw = np.append(var_raw, curr_file.variables[var][:], axis=0)
        GPH_raw = np.append(GPH_raw, curr_file.variables['Z3'][:], axis=0)
        PS_raw  = np.append(PS_raw,  curr_file.variables['PS'][:], axis=0)
        TROP_P_raw  = np.append(TROP_P_raw,  curr_file.variables['TROP_P'][:], axis=0)
    
    PS_raw = np.squeeze(np.mean(PS_raw, axis=0))
    NT = np.shape(var_raw)[0]
    NP = np.shape(var_raw)[2]    
    var_array_plev   = np.zeros((NT, len(prs_lvls), len(yi), len(xi)))-999  # Blank array for the variable information
    var_array_trzlev = np.zeros((NT, len(troprelz_lvls), len(yi), len(xi)))-999  # Blank array for the variable information
  
    for t in range(NT): 
        print(t+1, ' of ', NT)
        
        newvar_raw_plev   = np.zeros((NP,len(prs_lvls)))-999
        newvar_raw_trzlev = np.zeros((NP,len(troprelz_lvls)))-999
        
        for pt in range(NP):          
            
            #print('*****************')
            
            curr_PS = PS_raw[pt]
            curr_TROP_P = TROP_P_raw[t,pt]/100.
            curr_var = np.squeeze(var_raw[t,:,pt])
            curr_GPH = np.squeeze(GPH_raw[t,:,pt])
            curr_plev = (hyam*1.e5 + hybm*curr_PS)/100.
            
            # Pressure level stuff
            pvar_f = interp1d(curr_plev, curr_var)
            newvar_raw_plev[pt,:] = pvar_f(prs_lvls)
            
            # Trop-rel stuff   
            plev_f = interp1d(curr_plev, np.arange(len(curr_plev)), bounds_error=False, fill_value=-999)            
            trop_index_frac = plev_f(curr_TROP_P)
            gph_f = interp1d(np.arange(len(curr_GPH)), curr_GPH, bounds_error=False, fill_value=-999)
            curr_TROP_gph = gph_f(trop_index_frac)
            
            
            curr_deltazlev = (curr_GPH - curr_TROP_gph)/1.0e3
            f_deltaz = interp1d(curr_deltazlev, curr_var, bounds_error=False, fill_value=-999)   
            newvar_raw_trzlev[pt,:] = f_deltaz(troprelz_lvls)   
            if np.isnan(np.min(f_deltaz(troprelz_lvls))): print(f_deltaz(troprelz_lvls))

        for p in range(len(prs_lvls)):
            var_interpolator = tri.LinearTriInterpolator(triang, np.squeeze(newvar_raw_plev[:,p]))
            curr_var = var_interpolator(Xi, Yi)[None,:]
            var_array_plev[t,p,:,:] = curr_var            

        for z in range(len(troprelz_lvls)):
            var_interpolator = tri.LinearTriInterpolator(triang, np.squeeze(newvar_raw_trzlev[:,z]))
            curr_var = var_interpolator(Xi, Yi)[None,:]
            var_array_trzlev[t,z,:,:] = curr_var    
            
            
    return var_array_plev, var_array_trzlev

def get_struc_var(FILES, var, prs_lvls, troprelz_lvls, lon, lat):

    file0 = nc4.Dataset(FILES[0])
    var_array = file0.variables[var][:]    
    hyam = file0.variables['hyam'][:] 
    hybm = file0.variables['hybm'][:]    
    PS = file0.variables['PS'][:]
    TROP_P = file0.variables['TROP_P'][:]
    GPH = file0.variables['Z3'][:]
    #var_dim = len(np.shape(var_array))
    
    # Read in the full variable arrays
    for f in range(len(FILES)-1):
    
        curr_file = nc4.Dataset(FILES[f+1])
        curr_var = curr_file.variables[var][:]
        curr_PS = curr_file.variables['PS'][:]
        curr_TROP_P = curr_file.variables['TROP_P'][:]
        curr_GPH = curr_file.variables['Z3'][:]

        var_array = np.append(var_array, curr_var, axis=0)
        TROP_P = np.append(TROP_P, curr_TROP_P, axis=0)
        GPH = np.append(GPH, curr_GPH, axis=0)
        PS = np.append(PS, curr_PS, axis=0)
        
    PS = np.squeeze(np.mean(PS, axis=0))

    #if var_dim == 4:

    var_array_plev   = np.zeros((np.shape(var_array)[0], len(prs_lvls), len(lat), len(lon)))-999  # Blank array for the variable information 
    var_array_trzlev = np.zeros((np.shape(var_array)[0], len(troprelz_lvls), len(lat), len(lon)))-999  # Blank array for the variable information 
    
    for lt in range(len(lat)):
        print('lat ', lt)
        for ln in range(len(lon)):
        
            curr_PS = PS[lt,ln]           
            curr_plev = (hyam*1.e5 + hybm*curr_PS)/100.
            plev_f = interp1d(curr_plev, np.arange(len(curr_plev)), bounds_error=False, fill_value=-999)
                
            for t in range(np.shape(var_array)[0]):

                curr_varprof = np.squeeze(var_array[t,:,lt,ln])
                curr_TROP_P = TROP_P[t,lt,ln]/100.
                curr_gphprof = np.squeeze(GPH[t,:,lt,ln])
            
                # Fill the pressure level array
                pvar_f = interp1d(curr_plev, np.squeeze(curr_varprof))
                var_array_plev[t,:,lt,ln] = pvar_f(prs_lvls)

                # Fill the tropopause-relative array
                
                trop_index_frac = plev_f(curr_TROP_P)
                gph_f = interp1d(np.arange(len(curr_gphprof)), curr_gphprof, bounds_error=False, fill_value=-999)
                curr_TROP_gph = gph_f(trop_index_frac)
                #TROP_gph_array.append(curr_TROP_gph)

                curr_deltazlev = (curr_gphprof - curr_TROP_gph)/1.0e3
                f_deltaz = interp1d(curr_deltazlev, curr_varprof, bounds_error=False, fill_value=-999)   
                                
                var_array_trzlev[t,:,lt,ln] = f_deltaz(troprelz_lvls)
                
            
    #else: var_array = var_array_full
    
    return var_array_plev, var_array_trzlev
      
      
# User settings
model_name = '110L_longspin'
vars     = ['U','V','Z3','CO','CFC11','CFC12','N2O','C3H8','C2H6','HCFC22','CH3CL','CFC113','CCL4','CFC115','CH2BR2','CH3BR','CH3CCL3','CF2CLBR','HCFC141B','HCFC142B','CFC114','CF3BR','O3']
var_mods = [  1,  1,   1, 1e9,   1e12,   1e12,  1e9,  1e12,  1e12,   1e12,    1e12,    1e12,   1e12,   1e12,    1e12,   1e12,     1e12,     1e12,      1e12,      1e12,    1e12,   1e12, 1e9]
#vars = ['E90']
#var_mods = [1e9]
prs_lvls = [70.,100.,150.,200.,250.]
troprelz_lvls = [-4.0,-3.0,-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0]
SAVEDIR = '/home/wsmith/STRATOCLIM/vars/Map_Data/'

DATAROOT = '/data/wsmith/forecasts/2017_' + model_name + '/'
FILES = glob.glob(DATAROOT + 'f*')

for v in range(len(vars)):

    var = vars[v]
    var_mod = var_mods[v]
        
    print('******************** Reading ' + var + ' data')
    
    file0 = nc4.Dataset(FILES[0])
    lon = file0.variables['lon'][:]
    lat = file0.variables['lat'][:]

    if len(lon) > 5000.: # Unstructured grid
        print('Unstructured grid detected...proceeding accordingly...')
        
        xi = np.linspace(0.0, 180., 721) 
        yi = np.linspace(0, 60, 241)
        
        triang = tri.Triangulation(lon, lat)
        Xi, Yi = np.meshgrid(xi, yi)        

        var_array_plev, var_array_trzlev = get_unstruc_var(FILES, var, prs_lvls, troprelz_lvls, triang, Xi, Yi)
        
        lon = xi
        lat = yi

    else: # Structured grid
                                                     
        var_array_plev, var_array_trzlev = get_struc_var(FILES, var, prs_lvls, troprelz_lvls, lon, lat)

    var_array_plev   = var_array_plev*var_mod
    var_array_trzlev = var_array_trzlev*var_mod
        
    with open(SAVEDIR + model_name + '_' + var + '_MapData_plev.pkl', 'wb') as f:  
        pkl.dump((var_array_plev, prs_lvls, lon, lat), f)  

    with open(SAVEDIR + model_name + '_' + var + '_MapData_trzlev.pkl', 'wb') as f:  
        pkl.dump((var_array_trzlev, troprelz_lvls, lon, lat), f)
        
    #print(np.shape(var_array))
    #print(np.nanmin(var_array), np.nanmax(var_array))
            
