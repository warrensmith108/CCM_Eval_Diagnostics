
import numpy as np
import matplotlib.pyplot as plt
import cartopy as cp
import cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

# Define a function which will turn a tracer-tracer relationship into a distribution 
def calc_TT_pix(xvals, xbins, xrange, yvals, ybins, yrange):
    pixels = np.histogram2d(xvals, yvals, bins=[xbins,ybins], range=[xrange,yrange])    
    return pixels[1], pixels[2], np.transpose(pixels[0])

# The below function plots tracer-tracer distributions as pixel distributions
def plot_TT_pix(ax, xvals, xbins, xrange, yvals, ybins, yrange, xlab, ylab, savename, lay_avg=0, boundary = -999, \
      do_yscale = 0, trop_data=[-999], n_xticks=5, yticks = [], quants=[0,100], do_save=1):
    # boundary marks the location where a curve will be added      
    fnt = 16
    
    trop_mean = np.nanmean(trop_data)
    trop_std  = np.nanstd(trop_data)
    
    if do_yscale == 1: trop_scale = trop_mean
    else: trop_scale = 0
    
    #pixx, pixy, pixdata_raw = calc_TT_pix(xvals, xbins, xrange, yvals, ybins, yrange)
    
    pixels = np.histogram2d(xvals, yvals, bins=[xbins,ybins], range=[xrange,yrange])  
    pixx = pixels[1]
    pixy = pixels[2]
    pixdata     = np.transpose(pixels[0])
    pixdata_raw = np.transpose(pixels[0])
    
    #pixdata[pixdata < 5] = 0  ############  BE CAREFUL WITH THIS!!!  THIS IS A SCREEN FOR LOW POINTS

    if lay_avg:
    
        bound_data = []
        for lev in range(np.shape(pixdata)[0]):
            curr_pixdata = pixdata[lev,:]
                   
            ####  This is for an attempt at filtering using quantiles, this is proving difficult so it 
            ####  is being abandoned for now.  Bottom line is in order to do this correctly, it will require 
            ####  a manual bypass to the calc_TT_pix function for the filtering to take place there
            
            #curr_pixdata_sum = np.sum(curr_pixdata_raw)
            #curr_lowpct_sum = curr_pixdata_sum*quants[0]/100.
            #curr_upppct_sum = curr_pixdata_sum*quants[1]/100.
            
            #curr_pixdata = []
            #for pt in range(len(curr_pixdata_raw)):

            #    curr_lowsum = np.sum(curr_pixdata_raw[:pt])
                
            #    if curr_lowsum >= curr_lowpct_sum and curr_lowsum <= curr_upppct_sum:
            #        curr_pixdata.append(curr_pixdata_raw[pt])
            #    else: curr_pixdata.append(0)            
            #curr_pixdata = np.asarray(curr_pixdata)
            
            pixdata[lev,:] = curr_pixdata*100./np.max(curr_pixdata)
            if np.isnan(pixdata[lev,0]): bound_data.append(-999)  # If there is no data at a level, all elements are nans
            else: bound_data.append(pixx[np.where(pixdata[lev,:] > boundary)[0][-1]+1])  # Add 1 to get the edge of that bin
            
    #else:
    #    pixdata = pixdata*100./np.max(pixdata)
    
    pixdata[pixdata == 0] = pixdata[pixdata == 0]-0.01
    
    rgb      = np.loadtxt('/home/wsmith/STRATOCLIM/Gray_to_Orange_cbar.txt', skiprows=5, usecols=(0,1,2))/255.
    #rgb_under = rgb[0,:]   
    #rgb_over  = rgb[-1,:]   
    rgb = rgb[1:151,:]       # Temporary to go to 150, but we only want gray to black

    levels = np.arange(1,100.1, 100./len(rgb[:,0]))
    
    cmap = ListedColormap(rgb, N=len(levels)-1)    
    cmap.set_under('white')
    
    #cmap = plt.get_cmap('YlOrRd')
    #cmap.set_under([0.85,0.85,0.85])
    #cmap = plt.get_cmap('binary')
    #cmap.set_under([1,1,1])
    #levels = MaxNLocator(nbins=10).tick_values(0, 100.)
    #levels = np.arange(0,501,50)
    #levels[0] = 5
    
    norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    

    
    #pixplt = ax.pcolormesh(pixx, pixy+trop_scale, pixdata, cmap = cmap, norm = norm)
    #cbar = plt.colorbar(pixplt, extend='both')
    #if lay_avg: cbar.ax.set_ylabel('Layer-Normalized Relative Frequency')
    ##else: cbar.ax.set_ylabel('Relative Frequency')
    #else: cbar.ax.set_ylabel('Number of measurements')
    #cbar.set_ticks([1,10,20,30,40,50,60,70,80,90,100])
    
    #print(np.min(xvals), np.max(xvals))
    #print(np.min(yvals), np.max(yvals))
        
    scatplt = ax.scatter(xvals, yvals+trop_scale, color='black', s=0.1)  # Alternative plotting method
    
    # Add tropopause error bars
    if trop_mean > -900:
        ax.fill_between(xrange, trop_mean-trop_std, trop_mean+trop_std, color='gray', alpha=0.3)
        trop_line = ax.plot(xrange,[trop_mean,trop_mean], color='k')
        #plt.plot([xrange[0], xrange[1]], [trop_mean+trop_std, trop_mean+trop_std], color='gray', linestyle='--')
        #plt.plot([xrange[0], xrange[1]], [trop_mean-trop_std, trop_mean-trop_std], color='gray', linestyle='--')
    
    # Add a convective influence boundary
    if boundary > -900:
        for lev in range(np.shape(pixdata)[0]):
            plt.plot([bound_data[lev],bound_data[lev]], [pixy[lev],pixy[lev+1]]+trop_scale, color = 'blue', linewidth=3)
    
    ax.set_xlabel(xlab, fontsize=fnt)
    ax.set_ylabel(ylab, fontsize=fnt)
    
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xticks(np.arange(xrange[0], xrange[1]+0.01, (xrange[1]-xrange[0])/n_xticks))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if len(yticks) > 2: ax.set_yticks(yticks)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontweight('medium')
        tick.label1.set_fontsize(fnt)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontweight('medium')
        tick.label1.set_fontsize(fnt)
    
    ax.set_xlim(xrange)
    ax.set_ylim([yrange[0]+trop_scale, yrange[1]+trop_scale])  # Custom y axis label
    #plt.ylim(4,21)
    
    #plt.tight_layout()
    
    if do_save:
        plt.savefig(savename + '.png')
        plt.close('all')
    
    return pixdata_raw

def plot_layerhist(xvals, xbins, xrange, yvals, ybins, yrange, xlab, ylab, savename, histbnds = [[-5,-1],[-17,-15]]):

    fnt = 14
    
    pixx, pixy, pixdata = calc_TT_pix(xvals, xbins, xrange, yvals, ybins, yrange)
    pixx_mid = (pixx[1:]+pixx[0:-1])/2.0
    pixy_mid = (pixy[1:]+pixy[0:-1])/2.0
    
    for lev in range(np.shape(pixdata)[0]):
        curr_pixdata = np.squeeze(pixdata[lev,:])
        curr_col = lev/ybins
        if pixy[lev] < 0: curr_rgb = [curr_col,0.35,0.35]
        else: curr_rgb = [curr_col, curr_col, curr_col]
        plt.hist(pixx[:-1], pixx, weights=curr_pixdata/np.sum(curr_pixdata), density=False, histtype='step', color=curr_rgb, label=str(pixy[lev]) + ' to ' + str(pixy[lev+1]) + 'km')
        
    plt.legend()    
    plt.xlabel(xlab, fontsize=fnt)
    plt.ylabel(ylab, fontsize=fnt)
    
    plt.savefig(savename + '.png')
    plt.close('all')    
 
def plot_boxwhisker(xvals, xbins, xrange, yvals, ybins, yrange, xlab, ylab, savename):

    fnt = 14
    pixx, pixy, pixdata = calc_TT_pix(xvals, xbins, xrange, yvals, ybins, yrange)
    pixx_mid = (pixx[1:]+pixx[0:-1])/2.
    pixy_mid = (pixy[1:]+pixy[0:-1])/2.
    
    low_quant = []
    up_quant = []
    low_quart = []
    up_quart = []
    medians = []
    means = []
    pixy_plot = []
    
    # Calculate the quantile values that we desire
    for lev in range(np.shape(pixdata)[0]):
        curr_var_lev = []
        curr_pixdata = np.squeeze(pixdata[lev,:])
        for val in range(len(curr_pixdata)):
            curr_val = curr_pixdata[val]
            for zzz in range(int(curr_val)):
                curr_var_lev.extend([pixx_mid[val]])
        if len(curr_var_lev) > 0:
            
            low_quant.extend([np.quantile(curr_var_lev, 0.05)])
            up_quant.extend([np.quantile(curr_var_lev, 0.95)])
            low_quart.extend([np.quantile(curr_var_lev, 0.25)])
            up_quart.extend([np.quantile(curr_var_lev, 0.75)])            
            medians.extend([np.quantile(curr_var_lev, 0.50)])
            means.extend([np.mean(curr_var_lev)])  
            #pixy_plot.extend([pixy_mid[lev]])              
            pixy_plot.extend([pixy[lev]])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    #medians_sm = []
    for lev in range(len(medians)-1):
        plt.plot([low_quant[lev],low_quant[lev]],[pixy_plot[lev],pixy_plot[lev+1]], color='gray')
        plt.plot([low_quant[lev],low_quant[lev+1]],[pixy_plot[lev+1],pixy_plot[lev+1]], color='gray')
        plt.plot([up_quant[lev],up_quant[lev]],[pixy_plot[lev],pixy_plot[lev+1]], color='gray')
        plt.plot([up_quant[lev],up_quant[lev+1]],[pixy_plot[lev+1],pixy_plot[lev+1]], color='gray')    
        
        plt.plot([low_quart[lev],low_quart[lev]],[pixy_plot[lev],pixy_plot[lev+1]], color='black')
        plt.plot([low_quart[lev],low_quart[lev+1]],[pixy_plot[lev+1],pixy_plot[lev+1]], color='black')
        plt.plot([up_quart[lev],up_quart[lev]],[pixy_plot[lev],pixy_plot[lev+1]], color='black')
        plt.plot([up_quart[lev],up_quart[lev+1]],[pixy_plot[lev+1],pixy_plot[lev+1]], color='black')  
        
        plt.plot([medians[lev],medians[lev]],[pixy_plot[lev],pixy_plot[lev+1]], color='blue')
        plt.plot([medians[lev],medians[lev+1]],[pixy_plot[lev+1],pixy_plot[lev+1]], color='blue') 
        plt.plot([means[lev],means[lev]],[pixy_plot[lev],pixy_plot[lev+1]], color='red')
        plt.plot([means[lev],means[lev+1]],[pixy_plot[lev+1],pixy_plot[lev+1]], color='red')
        #medians_sm.append(np.mean(medians[lay-1:lay+1]))
    
    #plt.plot(medians_sm, pixy_plot, color='red')

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xticks(np.arange(xrange[0], xrange[1]+0.01, (xrange[1]-xrange[0])/5))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_yticks([0,4,8,12,16,20])
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.xlabel(xlab, fontsize=fnt)
    plt.ylabel(ylab, fontsize=fnt)
    
    plt.xlim(xrange)
    plt.ylim(yrange)    

    plt.savefig(savename + '.png')
    plt.close('all') 
    

# The below function plots YvsX with a 1-to-1 line overlaid
def plot_1to1(xvals, yvals, xlab, ylab, xrange, yrange, savename):
    plt.scatter(xvals, yvals)
    plt.plot([-1000,1000], [-1000,1000], color = 'k')    
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xlim(xrange)
    plt.ylim(yrange)
    plt.savefig(savename + '.png')
    plt.close('all')    
    
# The below function produces XY, XZ and YZ plots of 1D data colored by a measured specie / variable
def xyz_survey(lon, lat, alt, var, mission, var_str, unit, savename, miss=0):

    # Filter for missing data (None type)
    lon_plt = lon[var != None]
    lat_plt = lat[var != None]
    alt_plt = alt[var != None]
    var_plt = var[var != None]
    miss_plt = mission[var != None]
    
    # Filter the data according to the mission desired    
    if miss > 0:
        lon_plt = lon_plt[miss_plt == miss]
        lat_plt = lat_plt[miss_plt == miss]
        alt_plt = alt_plt[miss_plt == miss]
        var_plt = var_plt[miss_plt == miss]    
    
    cm = plt.cm.get_cmap('viridis')
    fig = plt.figure(figsize = (12,10))
    ax1 = fig.add_subplot(projection = ccrs.PlateCarree())
    ax2 = fig.add_subplot(221)
    ax3 = fig.add_subplot(224)

    if miss == 0:
        ax1.set_position([0.125,0.08,0.35,0.41])
        ax1.set_aspect(3)         
        tstr = 'All'
    if miss == 1: 
        ax1.set_position([0.125,0.08,0.35,0.41])
        ax1.set_aspect(1)            
        tstr = '2016'
    if miss == 2:
        ax1.set_position([0.125,0.08,0.35,0.41])
        ax1.set_aspect(1.2)     
        tstr = '2017'
        
    ax3.set_position([0.53,0.105,0.2,0.36])
    
    s1 = ax1.scatter(lon_plt, lat_plt, c=var_plt, vmin = np.min(var_plt), vmax = np.max(var_plt), cmap = cm, s = 50)
    s2 = ax2.scatter(lon_plt, alt_plt, c=var_plt, vmin = np.min(var_plt), vmax = np.max(var_plt), cmap = cm)
    s3 = ax3.scatter(alt_plt, lat_plt, c=var_plt, vmin = np.min(var_plt), vmax = np.max(var_plt), cmap = cm)    

    ax2.set_ylabel('Altitude (km)')
    ax2.set_xlabel('Longitude (E)')
    ax3.set_ylabel('Latitude (N)')
    ax3.set_xlabel('Altitude (km)')
    ax1.coastlines('50m', linewidth=0.5, color='gray')
    ax1.add_feature(cp.feature.BORDERS, linestyle = '-', linewidth=0.3, facecolor = 'none', edgecolor='gray')
    
    plt.figtext(0.65, 0.70, tstr + ' Data', fontsize = 28)   
    plt.figtext(0.65, 0.65, var_str, fontsize = 28)
        
    cbax = fig.add_axes([0.51, 0.53, 0.02, 0.35])        
    cbar = plt.colorbar(s2, cax = cbax)
    cbar.ax.set_ylabel(var_str + ' ' + unit, fontsize = 16)
    cbar.ax.tick_params(labelsize = 16)
    
    plt.savefig(savename + '_' + tstr + '.png')
    plt.close('all')
    
# The below function produces TT distributions for the measured variables 
def TT(var1, var2, lab1, lab2, unit1, unit2, savename, miss=0):
    
    if miss == 0: tstr = 'All'
    if miss == 1: tstr = '2016'
    if miss == 2: tstr = '2017'
    
    #var1_plt = [(var1 != None) & (var2 != None)]
    
    plt.scatter(var1, var2)
    plt.xlabel(lab1 + ' ' + unit1)
    plt.ylabel(lab2 + ' ' + unit2)
    plt.title(tstr)
    
    plt.savefig(savename + '_' + tstr + '.png')
    plt.close('all')

    