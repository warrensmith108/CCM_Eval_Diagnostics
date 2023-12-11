
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
    
def make_map_plot(ax, SAVEDIR, model_name, var, var_lab, unit, prs_lvl, lt, lev_choice, cbaraxes = [0.13,0.21,0.77,0.02]):

    with open(SAVEDIR + model_name + '_' + var + '_MapData_' + lev_choice + '.pkl', 'rb') as f:  
        var_array, lvls, lon, lat = pkl.load(f)

    with open(SAVEDIR + model_name + '_Z3_MapData_' + lev_choice + '.pkl', 'rb') as f:  
        gph_array, lvls, lon, lat = pkl.load(f)

    with open(SAVEDIR + model_name + '_U_MapData_' + lev_choice + '.pkl', 'rb') as f:  
        U_array, lvls, lon, lat = pkl.load(f)

    with open(SAVEDIR + model_name + '_V_MapData_' + lev_choice + '.pkl', 'rb') as f:  
        V_array, lvls, lon, lat = pkl.load(f)
    
    print('***', model_name, '***' , var, '***', prs_lvl, ' level')
    p = list(lvls).index(prs_lvl)
    print(p)
    
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

    ax.set_ylim([-90,90])
    ax.set_xlim([-180,180])
    
    con = ax.contourf(map_lon, lat, var_mean, cmap = colormap, levels = levels, norm=norm, extend = 'both', transform=ccrs.PlateCarree())
    ax.contour(lon, lat, gph_mean, levels = [17270.], colors = 'white', transform=ccrs.PlateCarree())
    ax.streamplot(lon-180, lat, U_mean, V_mean, density=0.8, linewidth = 1, color='grey', transform=ccrs.PlateCarree())
    ax.coastlines()
    
    deg = u"\N{DEGREE SIGN}"
    ax.set_xticks([-180,-90,0,90,179.9])
    ax.set_xticklabels(['180'+deg+'W','90'+deg+'W','0'+deg,'90'+deg+'E','180'+deg+'E'], fontsize=16)
    ax.set_yticks([-90,-45,0,45,90])
    ax.set_yticklabels(['90'+deg+'S','45'+deg+'S','0'+deg,'45'+deg+'N','90'+deg+'N'], fontsize=16)                
    
    #plot_box(ax, [75,95,18,32])

    cax = fig.add_axes(cbaraxes)
    cbar = plt.colorbar(con, cax=cax, orientation='vertical', spacing='uniform', ticks=levels[::4], extend = 'both', pad=0.01, fraction=0.8) 
    cbar.set_label(var_lab + ' (pp' + unit + 'v)', fontsize=16)

    for l in cbar.ax.yaxis.get_ticklabels(): l.set_fontsize(16)
        
    #if lev_choice == 'plev': ax.set_title(model_name + ' STRATOCLIM Mean ' + var + ' ' + str(int(prs_lvl)) + 'hPa')
    #else: ax.set_title(model_name + ' STRATOCLIM Mean ' + var + ' ' + str(prs_lvl) + 'km above local tropopause')
    #plt.figtext(0.7,0.35,'GPH = ' + str(int(gph_asma)) + 'm', color = 'white')

    
    
    #plt.savefig(PLOTDIR + var + '_' + str(prs_lvl) + '.png')
    #plt.close('all')


if __name__ == "__main__":

    model_name = '110L'    
    SAVEDIR = '/home/wsmith/STRATOCLIM/vars/Map_Data/'
    PLOTDIR = '/home/wsmith/STRATOCLIM/plots/' + model_name + '/maps/'
    lev_choice = 'trzlev'

    try:
        os.mkdir(PLOTDIR)
    except:
        zzz = -999
            
    # Make the plot
    fig = plt.figure(figsize=(14,10))
    ax1  = plt.subplot(2,2,1, projection=ccrs.PlateCarree())
    ax2  = plt.subplot(2,2,2, projection=ccrs.PlateCarree())
    ax3  = plt.subplot(2,2,3, projection=ccrs.PlateCarree())
    ax4  = plt.subplot(2,2,4, projection=ccrs.PlateCarree())
    #ax.set_ylim(0,60)
    #ax.set_xlim(0,180)                

    make_map_plot(ax1, SAVEDIR, model_name, 'C2H6', 'C$_2$H$_6$', 't', 0.5, '0.2', lev_choice, cbaraxes=[0.41,0.51,0.01,0.23])
    make_map_plot(ax2, SAVEDIR, model_name, 'CO', 'CO', 'b', 0.5, '0.3', lev_choice, cbaraxes=[0.91,0.51,0.01,0.23])
    make_map_plot(ax3, SAVEDIR, model_name, 'CFC12', 'CFC-12', 't', 0.5, '11600', lev_choice, cbaraxes=[0.41,0.21,0.01,0.23])
    make_map_plot(ax4, SAVEDIR, model_name, 'N2O', 'N$_2$O', 'b', 0.5, '15600', lev_choice, cbaraxes=[0.91,0.21,0.01,0.23])
    
    ax1.set_position([0.05,0.50,0.35,0.25])
    ax2.set_position([0.55,0.50,0.35,0.25])
    ax3.set_position([0.05,0.20,0.35,0.25])
    ax4.set_position([0.55,0.20,0.35,0.25])

    plt.figtext(0.2,0.77,'WACCM 500m above the local tropopause', fontsize=24, fontweight='bold')
    plt.figtext(0.06,0.52,'(a) C$_2$H$_6$ (0.2)', fontweight='bold', fontsize=16, bbox ={'edgecolor':'black','facecolor':'white'})
    plt.figtext(0.56,0.52,'(b) CO (0.3)', fontweight='bold', fontsize=16, bbox ={'edgecolor':'black','facecolor':'white'})
    plt.figtext(0.06,0.22,'(c) CFC-12 (11600)', fontweight='bold', fontsize=16, bbox ={'edgecolor':'black','facecolor':'white'})
    plt.figtext(0.56,0.22,'(d) N$_2$O (15600)', fontweight='bold', fontsize=16, bbox ={'edgecolor':'black','facecolor':'white'})
    
    plt.savefig(PLOTDIR + 'Figure2.png', dpi=300)
    plt.savefig(PLOTDIR + 'Figure2.pdf', dpi=300)
    plt.close('all')
    
        
        

        

        
    
