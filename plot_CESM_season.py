
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os
import glob
import utils
import cartopy as cp
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point


def make_hist(var_flat, gph_flat, gph_asma, rng, var, model_name, prs_lvl, lat_bnds, lon_bnds, PLOTDIR):

    var_in = []  # Inside the ASM
    var_out = [] # Outside the ASM

    # Screen for inside and outside
    for pt in range(len(var_flat)):
        
        curr_gph = gph_flat[pt]
        
        if curr_gph >= gph_asma:
            var_in.append(var_flat[pt])
        else:
            var_out.append(var_flat[pt])

    # Calculate interesting numbers
    mean_in = np.mean(var_in)
    mean_out = np.mean(var_out)    
    pct_in = len(var_in)*100./len(var_flat)

    bins = 30
            
    # Make the plot
    hist_in = plt.hist(var_in, density = True, range = rng, bins = bins, color = 'orange')
    hist_out = plt.hist(var_out, density = True, range = rng, bins = bins, histtype = 'step', color = 'blue')

    # Dress the plot up
    plt.title(model_name + ' July 15 - Aug 31 2017 ' + var + ' ' + str(int(prs_lvl)) + 'hPa')
    plt.ylabel('Relative Frequency')
    plt.xlabel(var + ' (ppbv)')

    plt.figtext(0.7,0.83,'In Mean = ' + str(round(mean_in, 1)) , color = 'orange')
    plt.figtext(0.7,0.80,'Out Mean = ' + str(round(mean_out, 1)), color = 'blue')

    plt.figtext(0.7,0.75,'Lat = ' + str(lat_bnds) + 'N')
    plt.figtext(0.7,0.72,'Lon = ' + str(lon_bnds) + 'E')
    plt.figtext(0.7,0.69,'GPH = ' + str(int(gph_asma)) + 'm')
    plt.figtext(0.7,0.66,'% in ASM = ' + str(round(pct_in, 1)))

    plt.savefig(PLOTDIR + 'Distribution_JA2021_' + var + '_' + str(int(prs_lvl)) + '.png')
    plt.close('all')


def plot_box(ax, box):
    ax.plot([box[0],box[0]],[box[2],box[3]], color = 'gray', linewidth=3)
    ax.plot([box[1],box[1]],[box[2],box[3]], color = 'gray', linewidth=3)
    ax.plot([box[0],box[1]],[box[2],box[2]], color = 'gray', linewidth=3)
    ax.plot([box[0],box[1]],[box[3],box[3]], color = 'gray', linewidth=3)
    

if __name__ == "__main__":

    model_names = ['110L']
    #vars  = ['CO','N2O','CFC11','CFC113','CCL4','CFC115','CFC12','CH2BR2','CH3BR','CH3CCL3','C2H6','C3H8',\
    #         'CH3CL','CF2CLBR','HCFC141B','HCFC142B','HCFC22','CFC114','CF3BR','O3']
    #units = [ 'b',  'b',  't',   't',     't',     't',    't',       't',    't',    't',    't',    't', \
    #             't',    't',      't',      't',       't',     't',    't',   'b']
    vars  = ['N2O','CFC12','CO','C2H6']
    var_labs = ['(d) N$_2$O','(c) CFC-12','(b) CO','(a) C$_2$H$_6$']
    units = [  'b',    't',   'b',    't', 't',    't',    't', 't',     't']
    lts   = ['15600','11600','0.3','0.2']  # Tropospheric lifetimes (y)
    #vars = ['CO']
    #lat_bnds = [0,60]
    #lon_bnds = [0,180]
    
    SAVEDIR = '/home/wsmith/STRATOCLIM/vars/Map_Data/'
    
    lev_choice = 'trzlev'

    for model_name in model_names:
    
        PLOTDIR = '/home/wsmith/STRATOCLIM/plots/' + model_name + '/maps/'
        
        try:
            os.mkdir(PLOTDIR)
        except:
            zzz = -999
        
        for v in range(len(vars)):
        
            var = vars[v]
            var_lab = var_labs[v]
            unit = units[v]
            lt = lts[v]
        
            with open(SAVEDIR + model_name + '_' + var + '_MapData_' + lev_choice + '.pkl', 'rb') as f:  
                var_array, lvls, lon, lat = pkl.load(f)

            with open(SAVEDIR + model_name + '_Z3_MapData_' + lev_choice + '.pkl', 'rb') as f:  
                gph_array, lvls, lon, lat = pkl.load(f)

            with open(SAVEDIR + model_name + '_U_MapData_' + lev_choice + '.pkl', 'rb') as f:  
                U_array, lvls, lon, lat = pkl.load(f)

            with open(SAVEDIR + model_name + '_V_MapData_' + lev_choice + '.pkl', 'rb') as f:  
                V_array, lvls, lon, lat = pkl.load(f)
            
            for p in range(len(lvls)):                               
                
                prs_lvl = lvls[p]
                print('***', model_name, '***' , var, '***', prs_lvl, ' level')
                
                ##############################################################################
                #   MAKE A HISTOGRAM OF VALUES INSIDE AND OUTSIDE THE ANTICYCLONE
                ##############################################################################
                
                # We choose manual values of GPH to distinguish between inside and outside the anticyclone
                #if prs_lvl == 150.: gph_asma = 14350.
                #if prs_lvl == 100.: gph_asma = 16770.
                #if prs_lvl == 70.:  gph_asma = 18920.
                
                #if prs_lvl > 90.: rng = (0,150)               
                #else: rng = (0,60)
        
                #var_flat = var_array.flatten()
                #gph_flat = gph_array.flatten()          
                
                #make_hist(var_flat, gph_flat, gph_asma, rng, var, model_name, prs_lvl, lat_bnds, lon_bnds, PLOTDIR)

                ##############################################################################
                #   MAKE AN AVERAGED MAP OF VALUES FOR THE SEASON
                ##############################################################################

                # Compute the time mean of the variable and GPH_arrays
                var_mean = np.nanmean(np.squeeze(var_array[:,p,:,:]), axis = 0)
                gph_mean = np.nanmean(np.squeeze(gph_array[:,p,:,:]), axis = 0)
                U_mean_raw   = np.nanmean(np.squeeze(U_array[:,p,:,:]), axis = 0)
                V_mean_raw   = np.nanmean(np.squeeze(V_array[:,p,:,:]), axis = 0)
                
                gph_mean = utils.moving_window_smooth(gph_mean, width=2)

                nlon = np.shape(var_mean)[1]
                
                # Change the wind fields to start at -180 because they won't plot right for some reason
                U_mean = np.append(U_mean_raw[:,int(nlon/2):], U_mean_raw[:,0:int(nlon/2)], axis=1)
                V_mean = np.append(V_mean_raw[:,int(nlon/2):], V_mean_raw[:,0:int(nlon/2)], axis=1)

                #lon = np.append(lon, [360.])

                var_mean, map_lon = add_cyclic_point(var_mean, coord=lon)

                if lev_choice == 'plev': levels, colormap, norm = utils.get_colorbar(var, prs_lvl, 0)  # Get colorbar information
                else: levels, colormap, norm = utils.get_colorbar(var, 100., 1)

                # Make the plot
                fig = plt.figure()
                ax  = plt.subplot(1,1,1, projection=ccrs.PlateCarree())
                #ax.set_ylim(0,60)
                #ax.set_xlim(0,180)                
                ax.set_ylim([-90,90])
                ax.set_xlim([-180,180])
                
                con = ax.contourf(map_lon, lat, var_mean, cmap = colormap, levels = levels, norm=norm, extend = 'both', transform=ccrs.PlateCarree())
                ax.contour(lon, lat, gph_mean, levels = [17270.], colors = 'white', transform=ccrs.PlateCarree())
                ax.streamplot(lon-180, lat, U_mean, V_mean, density=0.8, linewidth = 1, color='grey', transform=ccrs.PlateCarree())
                ax.coastlines()
                
                deg = u"\N{DEGREE SIGN}"
                ax.set_xticks([-180,-120,-60,0,60,120,179.9])
                ax.set_xticklabels(['180'+deg+'W','120'+deg+'W','60'+deg+'W','0'+deg,'60'+deg+'E','120'+deg+'E','180'+deg+'E'], fontsize=16, fontweight='bold')
                ax.set_yticks([-90,-60,-30,0,30,60,90])
                ax.set_yticklabels(['90'+deg+'S','60'+deg+'S','30'+deg+'S','0'+deg,'30'+deg+'N','60'+deg+'N','90'+deg+'N'], fontsize=16, fontweight='bold')                
                
                #plot_box(ax, [75,95,18,32])

                cbaraxes = [0.13,0.21,0.77,0.02]
                cax = fig.add_axes(cbaraxes)
                cbar = plt.colorbar(con, cax=cax, orientation='horizontal', spacing='uniform', ticks=levels[::4], extend = 'both')
                cbar.set_label(var_lab + ' (pp' + unit + 'v)')

                #if lev_choice == 'plev': ax.set_title(model_name + ' STRATOCLIM Mean ' + var + ' ' + str(int(prs_lvl)) + 'hPa')
                #else: ax.set_title(model_name + ' STRATOCLIM Mean ' + var + ' ' + str(prs_lvl) + 'km above local tropopause')
                #plt.figtext(0.7,0.35,'GPH = ' + str(int(gph_asma)) + 'm', color = 'white')

                plt.figtext(0.16,0.28,var_lab + ' (' + lt + ')', fontweight='bold', fontsize=16, bbox ={'edgecolor':'black','facecolor':'white'})
                
                plt.savefig(PLOTDIR + var + '_' + str(prs_lvl) + '.png')
                plt.close('all')


