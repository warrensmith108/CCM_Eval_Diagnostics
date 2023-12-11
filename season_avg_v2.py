
################################################################
#
#  This code is intended to read 2017 simulation output and 
#  compute season-averaged distributions and maps of species
#
################################################################

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob
import os
import netCDF4 as nc4
from scipy.interpolate import interp1d
import math
import utils
import matplotlib.tri as tri
#from plot_season import make_hist
import pickle as pkl
from utils import moving_window_smooth as mws
import scale_hgt as sh

do_prof = 1 # Determines if the code will calculate mean profiles of model data to compare against aircraft

sims = ['NOVSL_ens']
#lon_bnds = [81,87]
#lat_bnds = [26,30]
lon_bnds = [75,95]
lat_bnds = [18,32]
#lon_bnds = [257,277]   # For SEAC4RS
#lat_bnds = [23,37]
#lon_bnds = [30,140]
#lat_bnds = [20,40]
#lon_bnds = [0,180]
#lat_bnds = [0,60]

for sim in sims:

    PLOTDIR = '/home/wsmith/STRATOCLIM/plots/' + sim + '/maps/'
    SAVEDIR = '/home/wsmith/STRATOCLIM/vars/'

    var_names = [   'O3'  ]#, 'CO'  ,   'N2O', 'CH2BR2', 'CH3CL', 'CH3BR','CH3CCL3','HCFC141B','HCFC22','HCFC142B','CF2CLBR','CCL4', \
                  #'CFC11','CF3BR','CFC113',  'CFC12','CFC114','CFC115'] 
    var_mods  =  [    1e9,    1e9,     1e9,     1e12,    1e12,    1e12,     1e12,      1e12,    1e12,      1e12,     1e12,  1e12, \
                     1e12,   1e12,    1e12,     1e12,    1e12,    1e12]
                     
    prs_lvls = [100.]
    gph_asmas = [[16730,16810]]
    #gph_asmas = [[18930,18970], [16730,16810], [14310,14390]]

    z_prof = np.arange(0,25,0.5)
    th_prof = np.arange(300,600,5)
    delta_th_prof = np.arange(-100,101,5)
    delta_z_prof = np.arange(-20,10.1,0.5)

    DATAROOT = '/data/wsmith/forecasts/2017_' + sim
    FILES = glob.glob(DATAROOT + '/f*')

    Rdcp = 287.06/1004.

    for v in range(len(var_names)):
        
        var_name = var_names[v]
        var_mod = var_mods[v]
        #rng = rngs[v]

        file0 = nc4.Dataset(FILES[0])

        var    = file0.variables[var_name][:]
        temp   = file0.variables['T'][:]
        PS     = file0.variables['PS'][:]
        TROP_P = file0.variables['TROP_P'][:]
        gph    = file0.variables['Z3'][:]
        U      = file0.variables['U'][:]
        V      = file0.variables['V'][:]

        for f in range(len(FILES)-1):
            
            file = nc4.Dataset(FILES[f+1])
            
            var    = np.append(var, file.variables[var_name][:], axis=0)
            temp   = np.append(temp, file.variables['T'][:], axis=0)
            PS     = np.append(PS, file.variables['PS'][:], axis=0)
            U      = np.append(U, file.variables['U'][:], axis=0)
            V      = np.append(V, file.variables['V'][:], axis=0)
            TROP_P = np.append(TROP_P, file.variables['TROP_P'][:], axis=0)
            gph    = np.append(gph, file.variables['Z3'][:], axis=0)
            
        var = var*var_mod
        NT = np.shape(var)[0]
        
        PS_mean = np.mean(PS, axis=0)

        for p in range(len(prs_lvls)):
        
            # Set up arrays so we can save the non-interpolated values
            var_1d = []
            troprelz_1d = []
            z_1d = []
            theta_1d = []
            tropreltheta_1d = []
        
            prs_lvl = prs_lvls[p]
            gph_asma = gph_asmas[p]
            TROP_gph_array = []
            TROP_theta_array = []

            lon  = file0.variables['lon'][:]
            lat  = file0.variables['lat'][:]
            lev  = file0.variables['lev'][:]
            hyam = file0.variables['hyam'][:]
            hybm = file0.variables['hybm'][:]
            NZ   = len(lev)

            print('**** ', var_name, ' *** ', prs_lvl)

            var_lvls_f = interp1d(lev, np.arange(len(lev)), bounds_error=None, fill_value=-999)
            var_lev_index_frac = var_lvls_f(prs_lvl)
            var_lev_index_dec  = var_lev_index_frac - math.floor(var_lev_index_frac)
            var_lev_index_flr  = math.floor(var_lev_index_frac)

            var_zprofs_in  = np.zeros((1,len(z_prof)))-999
            var_zprofs_out = np.zeros((1,len(z_prof)))-999

            var_thprofs_in  = np.zeros((1,len(th_prof)))-999
            var_thprofs_out = np.zeros((1,len(th_prof)))-999
            
            var_troprelzprofs_in  = np.zeros((1,len(delta_z_prof)))-999
            var_troprelzprofs_out = np.zeros((1,len(delta_z_prof)))-999
            
            var_troprelthprofs_in  = np.zeros((1,len(delta_th_prof)))-999
            var_troprelthprofs_out = np.zeros((1,len(delta_th_prof)))-999
            
            if len(lon) < 5000.:  # Structured grid

                # Trim the plot to just the lon/lat boundaries
                lon_ind0 = np.min(np.where(lon >= lon_bnds[0])[0])
                lon_ind1 = np.max(np.where(lon <= lon_bnds[1])[0])+1
                lat_ind0 = np.min(np.where(lat >= lat_bnds[0])[0])
                lat_ind1 = np.max(np.where(lat <= lat_bnds[1])[0])+1
                
                var_sub = var[:,:,lat_ind0:lat_ind1,lon_ind0:lon_ind1]
                TROP_P_sub = TROP_P[:,lat_ind0:lat_ind1,lon_ind0:lon_ind1]
                PS_sub = PS[:,lat_ind0:lat_ind1,lon_ind0:lon_ind1]
                temp_sub = temp[:,:,lat_ind0:lat_ind1,lon_ind0:lon_ind1]
                gph_sub = gph[:,:,lat_ind0:lat_ind1,lon_ind0:lon_ind1]
                var_lev = (1-var_lev_index_dec)*var_sub[:,var_lev_index_flr,:,:] + var_lev_index_dec*var_sub[:,var_lev_index_flr+1,:,:]   
                gph_lev = (1-var_lev_index_dec)*gph_sub[:,var_lev_index_flr,:,:] + var_lev_index_dec*gph_sub[:,var_lev_index_flr+1,:,:] 
                U_lev = (1-var_lev_index_dec)*U[:,var_lev_index_flr,lat_ind0:lat_ind1,lon_ind0:lon_ind1] + var_lev_index_dec*U[:,var_lev_index_flr+1,lat_ind0:lat_ind1,lon_ind0:lon_ind1] 
                V_lev = (1-var_lev_index_dec)*V[:,var_lev_index_flr,lat_ind0:lat_ind1,lon_ind0:lon_ind1] + var_lev_index_dec*V[:,var_lev_index_flr+1,lat_ind0:lat_ind1,lon_ind0:lon_ind1] 
                
                lon_sub = lon[lon_ind0:lon_ind1]
                lat_sub = lat[lat_ind0:lat_ind1]
                
                if do_prof:
                
                    ################ Compute average profiles within and outside the anticyclone
                    for t in range(NT):
                        print('****************  t = ', t)

                        # At each time step, we will generate a new array of profiles so we can add that to the master array all at one time
                        curr_var_zprofs_in  = np.zeros((1,len(z_prof)))-999
                        curr_var_zprofs_out = np.zeros((1,len(z_prof)))-999

                        curr_var_thprofs_in  = np.zeros((1,len(th_prof)))-999
                        curr_var_thprofs_out = np.zeros((1,len(th_prof)))-999
                        
                        curr_var_troprelzprofs_in  = np.zeros((1,len(delta_z_prof)))-999
                        curr_var_troprelzprofs_out = np.zeros((1,len(delta_z_prof)))-999

                        curr_var_troprelthprofs_in  = np.zeros((1,len(delta_th_prof)))-999
                        curr_var_troprelthprofs_out = np.zeros((1,len(delta_th_prof)))-999
                    
                        lt_in = []
                        ln_in = []
                        lt_out = []
                        ln_out = []
                        
                        for lt in range(len(lat_sub)):
                            for ln in range(len(lon_sub)):
                        
                                # Unfortunately these simulations are on "raw" model levels which means we need to account for topography.  
                                curr_TROP_P = TROP_P_sub[t,lt,ln]/100.
                                curr_PS = PS_sub[t,lt,ln] 
                                curr_plev = (hyam*1.0e5 + hybm*curr_PS)/100.
                                curr_gph_lev = gph_lev[t,lt,ln]                            
                                curr_gphprof = np.squeeze(gph_sub[t,:,lt,ln])                            
                                curr_varprof = np.squeeze(var_sub[t,:,lt,ln])
                                
                                plev_f = interp1d(curr_plev, np.arange(len(curr_plev)), bounds_error=False, fill_value=-999)
                                trop_index_frac = plev_f(curr_TROP_P)
                                gph_f = interp1d(np.arange(len(curr_gphprof)), curr_gphprof, bounds_error=False, fill_value=-999)
                                curr_TROP_gph = gph_f(trop_index_frac)
                                TROP_gph_array.append(curr_TROP_gph)

                                curr_deltazlev = (curr_gphprof - curr_TROP_gph)/1.0e3
                                f_z = interp1d(curr_gphprof*1.e-3, curr_varprof, bounds_error=False, fill_value=-999)
                                f_deltaz = interp1d(curr_deltazlev, curr_varprof, bounds_error=False, fill_value=-999)                                                                                                                                  
                                
                                curr_tempprof = np.squeeze(temp_sub[t,:,lt,ln])
                                curr_thetaprof = curr_tempprof*((1000./curr_plev)**Rdcp)
                                f_th = interp1d(curr_thetaprof, curr_varprof, bounds_error=False, fill_value=-999)

                                th_f = interp1d(np.arange(len(curr_thetaprof)), curr_thetaprof, bounds_error=False, fill_value=-999)
                                curr_TROP_theta = th_f(trop_index_frac)  
                                TROP_theta_array.append(curr_TROP_theta)
                                #curr_TROP_theta = curr_thetaprof[TROP_nrst_index]  # Current tropopause in theta coordinates
                                curr_deltathetaprof = curr_thetaprof - curr_TROP_theta
                                f_deltatheta = interp1d(curr_deltathetaprof, curr_varprof, bounds_error=False, fill_value=-999)  # The 9s are intentional to flag missing data as negative (missing)                                          

                                var_1d.extend(curr_varprof)
                                z_1d.extend(curr_gphprof*1.e-3)
                                troprelz_1d.extend(curr_deltazlev)
                                theta_1d.extend(curr_thetaprof)
                                tropreltheta_1d.extend(curr_deltathetaprof)

                                if curr_gph_lev >= np.mean(gph_asma): 
                                    #ln_in.append(ln)
                                    #lt_in.append(lt)
                                    #var_profs_in_theta  = np.append(var_profs_in_theta, f_deltatheta(delta_th_prof)[None,:], axis=0)
                                    #curr_var_zprofs_in         = np.append(, axis=0)
                                    curr_var_zprofs_in = np.append(curr_var_zprofs_in, f_z(z_prof)[None,:], axis=0)
                                    curr_var_thprofs_in = np.append(curr_var_thprofs_in, f_th(th_prof)[None,:], axis=0)
                                    curr_var_troprelthprofs_in = np.append(curr_var_troprelthprofs_in, f_deltatheta(delta_th_prof)[None,:], axis=0)
                                    curr_var_troprelzprofs_in  = np.append(curr_var_troprelzprofs_in, f_deltaz(delta_z_prof)[None,:], axis=0)
                                else:
                                    #ln_out.append(ln)
                                    #lt_out.append(lt)
                                    #var_profs_out_theta = np.append(var_profs_out_theta, f_deltatheta(delta_th_prof)[None,:], axis=0)
                                    curr_var_zprofs_out  = np.append(curr_var_zprofs_out, f_z(z_prof)[None,:], axis=0)    
                                    curr_var_thprofs_out = np.append(curr_var_thprofs_out, f_th(th_prof)[None,:], axis=0)                                    
                                    curr_var_troprelthprofs_out = np.append(curr_var_troprelthprofs_out, f_deltatheta(delta_th_prof)[None,:], axis=0)
                                    curr_var_troprelzprofs_out  = np.append(curr_var_troprelzprofs_out, f_deltaz(delta_z_prof)[None,:], axis=0)                                                                
                        
                        var_zprofs_in = np.append(var_zprofs_in, curr_var_zprofs_in[1:,:], axis=0)
                        var_zprofs_out = np.append(var_zprofs_out, curr_var_zprofs_out[1:,:], axis=0)                        

                        var_thprofs_in = np.append(var_thprofs_in, curr_var_thprofs_in[1:,:], axis=0)
                        var_thprofs_out = np.append(var_thprofs_out, curr_var_thprofs_out[1:,:], axis=0) 
                        
                        var_troprelzprofs_in = np.append(var_troprelzprofs_in, curr_var_troprelzprofs_in[1:,:], axis=0)
                        var_troprelzprofs_out = np.append(var_troprelzprofs_out, curr_var_troprelzprofs_out[1:,:], axis=0)

                        var_troprelthprofs_in  = np.append(var_troprelthprofs_in, curr_var_troprelthprofs_in[1:,:], axis=0)
                        var_troprelthprofs_out = np.append(var_troprelthprofs_out, curr_var_troprelthprofs_out[1:,:], axis=0)
                        
                        #var_profs_in  = np.append(var_profs_in, var_sub[t,:,lt_in,ln_in], axis=0)
                        #var_profs_out = np.append(var_profs_out, var_sub[t,:,lt_out,ln_out], axis=0)
                
                    #print(np.shape(var_zprofs_in))
                    #print(np.shape(var_zprofs_out))
                    
                PS_plot = np.mean(PS_sub, axis=0)
            
            else:  # Regional refinement

                print('This is a regionally refined grid - converting to rectangular grid...')
                xi = np.linspace(0.0, 359.75, 1440)  # For now, we define a quarter degree rectangular grid to interpolate to
                yi = np.linspace(-90, 90, 721)
                triang = tri.Triangulation(lon, lat)
                Xi, Yi = np.meshgrid(xi, yi)    
                
                #var_lev = np.zeros((1, len(yi), len(xi)))
                #gph_lev = np.zeros((1, len(yi), len(xi)))
                #U_lev = np.zeros((1, len(yi), len(xi)))
                #V_lev = np.zeros((1, len(yi), len(xi)))
                
                PS_interpolator = tri.LinearTriInterpolator(triang, PS_mean)   
                PS_plot = PS_interpolator(Xi, Yi)
                    
                for t in range(NT):
                    print('t = ', t)
                    
                    #var_arr = np.squeeze((1-var_lev_index_dec)*var[t,var_lev_index_flr,:] + var_lev_index_dec*var[t,var_lev_index_flr+1,:])
                    #var_interpolator = tri.LinearTriInterpolator(triang, var_arr)   
                    #var_lev = np.append(var_lev, var_interpolator(Xi, Yi)[None,:], axis=0)
                    
                    gph_arr = np.squeeze((1-var_lev_index_dec)*gph[t,var_lev_index_flr,:] + var_lev_index_dec*gph[t,var_lev_index_flr+1,:])
                    #gph_interpolator = tri.LinearTriInterpolator(triang, gph_arr)
                    #gph_lev = np.append(gph_lev, gph_interpolator(Xi, Yi)[None,:], axis=0)

                    #U_arr = np.squeeze((1-var_lev_index_dec)*U[t,var_lev_index_flr,:] + var_lev_index_dec*U[t,var_lev_index_flr+1,:])
                    #U_interpolator = tri.LinearTriInterpolator(triang, U_arr)
                    #U_lev = np.append(U_lev, U_interpolator(Xi, Yi)[None,:], axis=0)

                    #V_arr = np.squeeze((1-var_lev_index_dec)*V[t,var_lev_index_flr,:] + var_lev_index_dec*V[t,var_lev_index_flr+1,:])
                    #V_interpolator = tri.LinearTriInterpolator(triang, V_arr)   
                    #V_lev = np.append(V_lev, V_interpolator(Xi, Yi)[None,:], axis=0)
                    
                    pts_in = []  # For each time step, we will save the points inside and outside the monsoon
                    pts_out = []

                    # At each time step, we will generate a new array of profiles so we can add that to the master array all at one time
                    curr_var_zprofs_in  = np.zeros((1,len(z_prof)))-999
                    curr_var_zprofs_out = np.zeros((1,len(z_prof)))-999

                    curr_var_thprofs_in  = np.zeros((1,len(th_prof)))-999
                    curr_var_thprofs_out = np.zeros((1,len(th_prof)))-999
                    
                    curr_var_troprelzprofs_in  = np.zeros((1,len(delta_z_prof)))-999
                    curr_var_troprelzprofs_out = np.zeros((1,len(delta_z_prof)))-999

                    curr_var_troprelthprofs_in  = np.zeros((1,len(delta_th_prof)))-999
                    curr_var_troprelthprofs_out = np.zeros((1,len(delta_th_prof)))-999
                    
                    if do_prof:
                    
                        for pt in range(len(lon)):
                        
                            if pt % 10000 == 0: print(pt, ' of ', len(lon))
                           
                            curr_lon = lon[pt]
                            curr_lat = lat[pt]
                            
                            if curr_lon >= lon_bnds[0] and curr_lon <= lon_bnds[1] and curr_lat >= lat_bnds[0] and curr_lat <= lat_bnds[1]:

                                # Unfortunately these simulations are on "raw" model levels which means we need to account for topography.  
                                curr_TROP_P = TROP_P[t,pt]/100.
                                curr_PS = PS[t,pt] 
                                curr_plev = (hyam*1.0e5 + hybm*curr_PS)/100.
                                curr_gph_lev = gph_arr[pt]       
                                curr_gphprof = np.squeeze(gph[t,:,pt])                            
                                curr_varprof = np.squeeze(var[t,:,pt])
                                
                                plev_f = interp1d(curr_plev, np.arange(len(curr_plev)), bounds_error=False, fill_value=-999)
                                trop_index_frac = plev_f(curr_TROP_P)
                                gph_f = interp1d(np.arange(len(curr_gphprof)), curr_gphprof, bounds_error=False, fill_value=-999)
                                curr_TROP_gph = gph_f(trop_index_frac)
                                TROP_gph_array.append(curr_TROP_gph)
                                
                                curr_deltazlev = (curr_gphprof - curr_TROP_gph)/1.0e3
                                f_z = interp1d(curr_gphprof*1.e-3, curr_varprof, bounds_error=False, fill_value=-999)
                                f_deltaz = interp1d(curr_deltazlev, curr_varprof, bounds_error=False, fill_value=-999)  
                             
                                curr_tempprof = np.squeeze(temp[t,:,pt])                        
                                curr_thetaprof = curr_tempprof*((1000./curr_plev)**Rdcp) 
                                f_th = interp1d(curr_thetaprof, curr_varprof, bounds_error=False, fill_value=-999)
                                
                                th_f = interp1d(np.arange(len(curr_thetaprof)), curr_thetaprof, bounds_error=False, fill_value=-999)
                                curr_TROP_theta = th_f(trop_index_frac)
                                TROP_theta_array.append(curr_TROP_theta)                            
                                #curr_TROP_theta = curr_thetaprof[TROP_nrst_index]  # Current tropopause in theta coordinates
                                curr_deltathetaprof = curr_thetaprof - curr_TROP_theta
                                f_deltatheta = interp1d(curr_deltathetaprof, curr_varprof, bounds_error=False, fill_value=-999)  # The 9s are intentional to flag data below the surface as negative

                                curr_gph = gph_arr[pt]
                                
                                var_1d.extend(curr_varprof)
                                z_1d.extend(curr_gphprof*1.e-3)
                                troprelz_1d.extend(curr_deltazlev)
                                theta_1d.extend(curr_thetaprof)
                                tropreltheta_1d.extend(curr_deltathetaprof)
                                
                                if curr_gph >= np.mean(gph_asma): 
                                    #pts_in.append(pt)
                                    #var_profs_in_theta  = np.append(var_profs_in_theta, f_deltatheta(delta_theta_prof)[None,:], axis=0) 
                                    curr_var_zprofs_in = np.append(curr_var_zprofs_in, f_z(z_prof)[None,:], axis=0)
                                    curr_var_thprofs_in = np.append(curr_var_thprofs_in, f_th(th_prof)[None,:], axis=0)
                                    curr_var_troprelthprofs_in = np.append(curr_var_troprelthprofs_in, f_deltatheta(delta_th_prof)[None,:], axis=0)
                                    curr_var_troprelzprofs_in  = np.append(curr_var_troprelzprofs_in, f_deltaz(delta_z_prof)[None,:], axis=0)                                
                                else: 
                                    #pts_out.append(pt)
                                    #var_profs_out_theta = np.append(var_profs_out_theta, f_deltatheta(delta_theta_prof)[None,:], axis=0)
                                    curr_var_zprofs_out = np.append(curr_var_zprofs_out, f_z(z_prof)[None,:], axis=0)
                                    curr_var_thprofs_out = np.append(curr_var_thprofs_out, f_th(th_prof)[None,:], axis=0)
                                    curr_var_troprelthprofs_out = np.append(curr_var_troprelthprofs_out, f_deltatheta(delta_th_prof)[None,:], axis=0)
                                    curr_var_troprelzprofs_out  = np.append(curr_var_troprelzprofs_out, f_deltaz(delta_z_prof)[None,:], axis=0)                    

                        var_zprofs_in = np.append(var_zprofs_in, curr_var_zprofs_in[1:,:], axis=0)
                        var_zprofs_out = np.append(var_zprofs_out, curr_var_zprofs_out[1:,:], axis=0)  

                        var_thprofs_in = np.append(var_thprofs_in, curr_var_thprofs_in[1:,:], axis=0)
                        var_thprofs_out = np.append(var_thprofs_out, curr_var_thprofs_out[1:,:], axis=0)
                        
                        var_troprelzprofs_in = np.append(var_troprelzprofs_in, curr_var_troprelzprofs_in[1:,:], axis=0)
                        var_troprelzprofs_out = np.append(var_troprelzprofs_out, curr_var_troprelzprofs_out[1:,:], axis=0)

                        var_troprelthprofs_in  = np.append(var_troprelthprofs_in, curr_var_troprelthprofs_in[1:,:], axis=0)
                        var_troprelthprofs_out = np.append(var_troprelthprofs_out, curr_var_troprelthprofs_out[1:,:], axis=0)
                                                                                                                          
                #var_lev = var_lev[1:,:,:]
                #gph_lev = gph_lev[1:,:,:]
                #U_lev = U_lev[1:,:,:]
                #V_lev = V_lev[1:,:,:]
                
                lon = xi
                lat = yi
                
                lon_sub = lon
                lat_sub = lat
                #lon_ind0 = np.min(np.where(lon >= lon_bnds[0])[0])
                #lon_ind1 = np.max(np.where(lon <= lon_bnds[1])[0])+1
                #lat_ind0 = np.min(np.where(lat >= lat_bnds[0])[0])
                #lat_ind1 = np.max(np.where(lat <= lat_bnds[1])[0])+1

            #var_lev_mean = np.mean(var_lev, axis=0)
            #gph_lev_mean = np.mean(gph_lev, axis=0)
            #U_lev_mean = np.mean(U_lev, axis=0)
            #V_lev_mean = np.mean(V_lev, axis=0)       
            
            
            ################# Save the profile information
            
            #with open(SAVEDIR + var_name + 'profiles_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl', 'wb') as f:
            #    pkl.dump((var_profs_in[1:,:], var_profs_out[1:,:], lev, gph_asma), f)

            if do_prof:

                with open(SAVEDIR + var_name + 'zprofiles_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl', 'wb') as f:
                    pkl.dump((var_zprofs_in[1:,:], var_zprofs_out[1:,:], TROP_gph_array, z_prof, gph_asma), f)

                with open(SAVEDIR + var_name + 'thetaprofiles_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl', 'wb') as f:
                    pkl.dump((var_thprofs_in[1:,:], var_thprofs_out[1:,:], TROP_theta_array, th_prof, gph_asma), f)
                    
                with open(SAVEDIR + var_name + 'troprelthetaprofiles_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl', 'wb') as f:
                    pkl.dump((var_troprelthprofs_in[1:,:], var_troprelthprofs_out[1:,:], TROP_theta_array, delta_th_prof, gph_asma), f)

                with open(SAVEDIR + var_name + 'troprelzprofiles_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl', 'wb') as f:
                    pkl.dump((var_troprelzprofs_in[1:,:], var_troprelzprofs_out[1:,:], TROP_gph_array, delta_z_prof, gph_asma), f)
                    
                # Also save the non-interpolated values
                with open(SAVEDIR + var_name + '_1d_nointerpvals_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl', 'wb') as f:
                    pkl.dump((var_1d, z_1d, troprelz_1d, theta_1d, tropreltheta_1d), f)
                
            ################ Make a map plot
            
            #fig = plt.figure()
            #ax  = plt.subplot(1,1,1, projection=ccrs.PlateCarree())

            #levels, colormap, norm = utils.get_colorbar(var_name, prs_lvl, 1)  # Get colorbar information

            #ax.set_ylim(lat_bnds)
            #ax.set_xlim(lon_bnds)
             
            #con = ax.contourf(lon_sub, lat_sub, var_lev_mean, cmap = colormap, levels = levels, norm=norm, extend = 'both')
            #cgph = ax.contour(lon_sub, lat_sub, gph_lev_mean, levels = gph_asma, colors = 'white')
            #plt.clabel(cgph, fontsize = 6, fmt='%1.0f')
            
            #stream = ax.streamplot(lon_sub, lat_sub, U_lev_mean, V_lev_mean, density=0.5, color='grey')
                    
            #ax.coastlines()

            #cax = fig.add_axes([0.13,0.27,0.77,0.02])
            #cbar = plt.colorbar(con, cax=cax, orientation='horizontal', spacing='uniform', ticks=levels[::4], extend = 'both')
            #if var_mod == 1e9: cbar.set_label(var_name + ' (ppbv)')
            #if var_mod == 1e12: cbar.set_label(var_name + ' (pptv)')
            
            #PS_plot = mws(PS_plot, width=1)
            
            # Plot the Tibetan Plateau
            #ax.contour(lon_sub, lat_sub, PS_plot/100., [609.], colors=['blue'], linewidths = 0.5)

            #ax.set_title(sim + ' STRATOCLIM Mean ' + var_name + ' ' + str(int(prs_lvl)) + 'hPa')
            #plt.figtext(0.7,0.35,'GPH = ' + str(gph_asma) + 'm', color = 'white')

            #plt.savefig(PLOTDIR + 'Mean_JA2021_' + var_name + '_' + str(int(prs_lvl)) + '.png', dpi=300)
            #plt.close('all')

            ################  Make a distribution

            # Note that for regional refinement, we are using the post-processed regular grid for this because we are lazy right now.
            # However, because of math, it is unlikely that this would be much different than using the "unstructured" grid data

            #var_flat = var_lev.flatten()
            #gph_flat = gph_lev.flatten()

            #make_hist(var_flat, gph_flat, np.mean(gph_asma), rng, var_name, sim, prs_lvl, lat_bnds, lon_bnds, PLOTDIR)



