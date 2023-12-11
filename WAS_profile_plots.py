
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import glob
import datetime as dt
import pickle as pkl
#from STRATOCLIM_plotter import plot_TT_pix, plot_layerhist, plot_boxwhisker
from read_STRATOCLIM_merge import get_model_prof_data
from read_STRATOCLIM_WAS import load_WAS
from matplotlib.ticker import ScalarFormatter
import os
import heatmap_config as hc
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from read_STRATOCLIM_merge import get_merge_array
from sklearn.utils import shuffle
import numpy.ma as ma

def get_trz_information(var_name, DATAFILES, WAS_SAVEDIR, zchoice):

    print(var_name)

    xlims = hc.heatmap_info[var_name]['lims']
    inst = str(hc.heatmap_info[var_name]['instrument'][0])
    
    if var_name == 'N2O':
    
        # Read the merge data for N2O
    
        WAS_var_raw, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_merge_array('HAGAR:N2O', DATAFILES)
        if zchoice == 'troprelz':
            
            WAS_z_raw = (np.asarray(AVIONIK_ALT)*1.e-3) - ERA5_LRT
            y_scale = np.nanmean(ERA5_LRT)
            
            #WAS_z_raw = (np.asarray(AVIONIK_ALT) - np.asarray(CPT_GPH))*1.e-3    
            #y_scale = np.nanmean(CPT_GPH)*1.e-3     
        
        # Screen the masked values in the dataset
        WAS_var = []
        WAS_z = []
        for pt in range(len(WAS_var_raw)):
            curr_WAS_var = WAS_var_raw[pt]
            if not ma.is_masked(curr_WAS_var): 
                WAS_var.append(curr_WAS_var)
                WAS_z.append(WAS_z_raw[pt])
 
        WAS_var = np.asarray(WAS_var)
        WAS_z = np.asarray(WAS_z)
        WAS_var_err = np.zeros(len(WAS_var))  # We don't have error for N2O
    
    else:
    
        # Read the tropopause information on WAS points
        
        with open(WAS_SAVEDIR + 'AVIONIK:ALT_merge.pkl', 'rb') as f:
            z, zz, WAS_ALT, zzz = pkl.load(f)
        #with open(WAS_SAVEDIR + 'CLAMS_MET:GPH_CPT_merge.pkl', 'rb') as f:
        #    z, zz, WAS_GPH_CPT, zzz = pkl.load(f)
        with open(WAS_SAVEDIR + 'ERA5_wmo_1st_z_merge.pkl', 'rb') as f:
            z, zz, WAS_ERA5_LRT, zzz = pkl.load(f)         
         
        if zchoice == 'troprelz':
            WAS_z = (np.asarray(WAS_ALT)*1.e-3) - np.asarray(WAS_ERA5_LRT)    
            y_scale = np.nanmean(WAS_ERA5_LRT)
        #elif zchoice == 'z':
        #    WAS_z = np.asarray(WAS_ALT)*1.e-3
        #    mean_trop = np.nanmean(WAS_GPH_CPT)*1.e-4
        #    y_scale = 0

        # Read the WAS data
        
        with open(WAS_SAVEDIR + inst + ':' + str(var_name) + '_merge.pkl', 'rb') as f:
            z, zz, WAS_var, WAS_var_err = pkl.load(f)

        WAS_var = np.asarray(WAS_var)
        
    return WAS_var, WAS_var_err, WAS_z, y_scale
    
    
    
WAS_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/WAS_merge/'
CESMSAVEDIR = '/home/wsmith/STRATOCLIM/vars/'
PLOTDIR = '/home/wsmith/STRATOCLIM/plots/WAS_profiles/'

DATAROOT = '/UTLS/Field_Data/STRATOCLIM_Data/merge_data/'
DATANAME = 'stratoclim*.nc'
DATAFILES = sorted(glob.glob(DATAROOT + DATANAME))

lon_bnds = [75,95]
lat_bnds = [18,32]
prs_lvl = 100.
quants = [5,95]
zchoice = 'troprelz'
z_rng = [-1,1]
#z_rng = [2,3]
#z_rng = [19.5,20.5]
var_names = ['CH2BR2','CH3CL', 'CH3BR','CH3CCL3','HCFC141B','HCFC22','HCFC142B','CF2CLBR','CCL4', \
              'CFC11','CF3BR','CFC113',  'CFC12',     'N2O','CFC114',  'CFC115']  
var_insts = ['(b) WAS',  'WAS',   'WAS',    'WAS',     'WAS','(c) WAS',     'WAS',    'WAS', 'WAS',\
                'WAS',  'WAS',   'WAS', '(a) WAS','(d) HAGAR',   'WAS',     'WAS']
var_labs = ['CH$_2$Br$_2$','CH$_3$Cl','CH$_3$Br','CH$_3$CCl$_3$','HCFC-141b','HCFC-22','HCFC-142b','CF$_2$ClBr','CCl$_4$', \
              'CFC-11','CF$_3$Br','CFC-113','CFC-12','N$_2$O','CFC-114','CFC-115'] 
var_tlife = [  '0.3',  '1.3',   '1.6',    '5.2' , '8.9',  '11.3',    '15.3',   '17.3','1230', \
               '1870', '4490',  '7620',  '11600',   '15600', '19600',  '126000']     # Tropospheric lifetime                

var1_caption = []
for vv in range(len(var_names)):
    var1_caption.append(var_labs[vv] + ' (' + var_tlife[vv] + ')')
    
score_table = np.zeros((2,len(var_names)))-999

for v in range(len(var_names)):

    var_name = var_names[v]
    WAS_var, WAS_var_err, WAS_z, y_scale = get_trz_information(var_name, DATAFILES, WAS_SAVEDIR, zchoice)

    # Read the model data
    
    lev_z, var_mean_110L, var_lowq_110L, var_upq_110L = get_model_prof_data(CESMSAVEDIR, var_name, zchoice, '110L', lon_bnds, lat_bnds, prs_lvl, quants)
    #lev_z, var_mean_110L_long, var_lowq_110L_long, var_upq_110L_long = get_model_prof_data(CESMSAVEDIR, var_name, zchoice, '110L_longspin', lon_bnds, lat_bnds, prs_lvl, quants)
    lev_z, var_mean_M32L, var_lowq_M32L, var_upq_M32L = get_model_prof_data(CESMSAVEDIR, var_name, zchoice, 'M32L', lon_bnds, lat_bnds, prs_lvl, quants)

    ######### Score the model's "stratospheric entry value"
    
    trop_WAS = []
    #trop_WAS_rand = []    
    for pt in range(len(WAS_z)):
        curr_z   = WAS_z[pt]
        curr_WAS = WAS_var[pt]
        curr_err = WAS_var_err[pt]
        if curr_z >= z_rng[0] and curr_z <= z_rng[1] and curr_WAS > -900:
            trop_WAS.append(curr_WAS)
            #trop_WAS_rand.extend(list(np.random.normal(loc=curr_WAS, scale=curr_err, size=100)))
    
    # Calculate the WAS "delta" value to use as the score denominator
    delta_WAS_var = np.nanmax(WAS_var) - np.nanmin(WAS_var[WAS_var >= 0])
    
    trop_WAS_mean = np.mean(np.asarray(trop_WAS))
    #trop_WAS_mean = np.mean(np.asarray(trop_WAS_rand))  # This is done as a crude way to account for uncertainty
    
    trop_model_index = list(lev_z).index(np.mean(z_rng))
    trop_110L_mean = var_mean_110L[trop_model_index]
    #trop_110L_long_mean = var_mean_110L_long[trop_model_index]
    trop_M32L_mean = var_mean_M32L[trop_model_index]
    
    #print(trop_WAS_mean)
    #print(np.mean(np.asarray(trop_WAS)))
    #print(trop_110L_mean)
    #print(trop_M32L_mean)
    
    score_table[0,v] = (trop_110L_mean-trop_WAS_mean)*100./delta_WAS_var
    #score_table[1,v] = np.abs(trop_110L_long_mean-trop_WAS_mean)*100./trop_WAS_mean 
    score_table[1,v] = (trop_M32L_mean-trop_WAS_mean)*100./delta_WAS_var 
    
print(score_table)

fnt = 18

# Make a final score table (scatter)

fig = plt.figure(figsize=(14,6))
ax = plt.subplot(111)

ax.hist(np.arange(0.4,len(var_names),1), len(var_names), weights=score_table[0,:], color='orange', histtype='bar', label='WACCM', rwidth=0.2)
ax.hist(np.arange(0.6,len(var_names),1), len(var_names), weights=score_table[1,:], color='red', histtype='bar', label='MUSICA', rwidth=0.2)

#plt.scatter(np.arange(0.4,len(var_names)), score_table[0,:], color='orange', s=150, marker='o', label='WACCM')
#plt.scatter(np.arange(0.6,len(var_names)), score_table[1,:], color='red', s=200, marker='x', label='MUSICA')
plt.plot([0,len(var_names)],[0,0],linestyle='--',color='k')
plt.plot([0,len(var_names)],[-40,-40],linestyle='--',color=[0.8,0.8,0.8])
plt.plot([0,len(var_names)],[-20,-20],linestyle='--',color=[0.8,0.8,0.8])
plt.plot([0,len(var_names)],[20,20],linestyle='--',color=[0.8,0.8,0.8])
plt.plot([0,len(var_names)],[40,40],linestyle='--',color=[0.8,0.8,0.8])

ax.set_title('Stratospheric Entry Mixing Ratio Errors', fontsize = fnt)

ax.set_xticks(np.arange(0.91,len(var_names),0.94))
ax.set_xticklabels(var1_caption, fontsize = fnt)
plt.xticks(rotation=90)

plt.ylabel('Percent error', fontsize=fnt)
plt.yticks(fontsize=fnt)
plt.ylim(-50,50)
plt.xlim(0,len(var_names))

plt.legend(fontsize=fnt)
plt.tight_layout()
plt.savefig(PLOTDIR + 'Figure6.png')  
plt.savefig(PLOTDIR + 'Figure6.pdf')  
plt.close('all')  
    

def make_plot(ax, var_name, var_lab, var_inst, WAS_var, WAS_z, y_scale, z_rng, zchoice):

    fnt = 16

    # Plot StratoClim data
    if var_name == 'N2O': ax.scatter(WAS_var, WAS_z+y_scale, s=1, color='k')
    else: ax.scatter(WAS_var, WAS_z+y_scale, s=10, color='k')
    #ax.plot([xlims[0],xlims[1]],[trz_rng[0],trz_rng[0]], color='gray', linestyle='--')
    #ax.plot([xlims[0],xlims[1]],[trz_rng[1],trz_rng[1]], color='gray', linestyle='--')
    
    xlims = hc.heatmap_info[var_name]['lims']
    
    lev_z, var_mean_110L, var_lowq_110L, var_upq_110L = get_model_prof_data(CESMSAVEDIR, var_name, zchoice, '110L', lon_bnds, lat_bnds, prs_lvl, quants)
    lev_z, var_mean_M32L, var_lowq_M32L, var_upq_M32L = get_model_prof_data(CESMSAVEDIR, var_name, zchoice, 'M32L', lon_bnds, lat_bnds, prs_lvl, quants)    

    # Plot StratoClim error bars
    for pt in range(len(WAS_var)):
        curr_var     = WAS_var[pt]
        curr_z       = WAS_z[pt]+y_scale
        curr_var_err = WAS_var_err[pt]
        if np.abs(curr_var_err) <= 100: ax.plot([curr_var-curr_var_err,curr_var+curr_var_err], [curr_z, curr_z], color='k', linewidth=0.5)

    ax.fill_between(np.arange(0,10000), np.zeros(10000)+z_rng[0]+y_scale, np.zeros(10000)+z_rng[1]+y_scale, color='gray', alpha = 0.5)
    ax.plot(np.arange(10000), np.zeros(10000)+y_scale, color='k')    

    ax.plot(var_mean_110L, lev_z+y_scale, color = 'orange')
    #ax.plot(var_lowq_110L, lev_z, color = 'orange', linestyle = '--')
    #ax.plot(var_upq_110L,  lev_z, color = 'orange', linestyle = '--')
    ax.fill_betweenx(lev_z+y_scale, var_lowq_110L, var_upq_110L, color='orange', alpha=0.2)

    ax.plot(var_mean_M32L, lev_z+y_scale, color = 'red')
    #ax.plot(var_lowq_M32L, lev_z, color = 'red', linestyle = '--')
    #ax.plot(var_upq_M32L,  lev_z, color = 'red', linestyle = '--')
    ax.fill_betweenx(lev_z+y_scale, var_lowq_M32L, var_upq_M32L, color='red', alpha=0.2)

    #ax.plot(var_mean_110L_long, lev_z+y_scale, color = 'blue')
    #ax.fill_betweenx(lev_z+y_scale, var_lowq_110L_long, var_upq_110L_long, color='blue', alpha=0.2)    

    if zchoice == 'z': ax.plot([0,10000], [mean_trop,mean_trop], linestyle='--', color='gray')

    ax.set_xlim(xlims)
    ax.set_ylim(12,21)

    if var_name == 'N2O': ax.set_xlabel(var_lab + ' (ppbv)', fontsize = fnt)
    else: ax.set_xlabel(var_lab + ' (pptv)', fontsize = fnt)

    if zchoice == 'troprelz': ax.set_ylabel('Adj Trop Rel Alt (km)', fontsize = fnt)
    elif zchoice == 'z':      ax.set_ylabel('Altitude (km)', fontsize = fnt)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontweight('medium')
        tick.label1.set_fontsize(fnt)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontweight('medium')
        tick.label1.set_fontsize(fnt)    

    
# Make Figure 5
fig = plt.figure(figsize=(12,10))
ax1  = plt.subplot(2,2,1)
ax2  = plt.subplot(2,2,2)
ax3  = plt.subplot(2,2,3)
ax4  = plt.subplot(2,2,4)

WAS_var, WAS_var_err, WAS_z, y_scale = get_trz_information('CFC12', DATAFILES, WAS_SAVEDIR, zchoice)
make_plot(ax1, 'CFC12', 'CFC-12', 'WAS', WAS_var, WAS_z, y_scale, z_rng, zchoice)
    
WAS_var, WAS_var_err, WAS_z, y_scale = get_trz_information('CH2BR2', DATAFILES, WAS_SAVEDIR, zchoice)
make_plot(ax2, 'CH2BR2', 'CH$_2$Br$_2$', 'WAS', WAS_var, WAS_z, y_scale, z_rng, zchoice)    

WAS_var, WAS_var_err, WAS_z, y_scale = get_trz_information('HCFC22', DATAFILES, WAS_SAVEDIR, zchoice)
make_plot(ax3, 'HCFC22', 'HCFC-22', 'WAS', WAS_var, WAS_z, y_scale, z_rng, zchoice)    

WAS_var, WAS_var_err, WAS_z, y_scale = get_trz_information('N2O', DATAFILES, WAS_SAVEDIR, zchoice)
make_plot(ax4, 'N2O', 'N$_2$O', 'HAGAR', WAS_var, WAS_z, y_scale, z_rng, zchoice)   

plt.figtext(0.08, 0.59, '(a) WAS CFC-12', fontsize=16, fontweight='bold')
plt.figtext(0.57, 0.59, '(b) WAS CH$_2$Br$_2$', fontsize=16, fontweight='bold')
plt.figtext(0.08, 0.09, '(c) WAS HCFC-22', fontsize=16, fontweight='bold')
plt.figtext(0.57, 0.09, '(d) HAGAR N$_2$O', fontsize=16, fontweight='bold')

plt.tight_layout()
plt.savefig(PLOTDIR + 'Figure5.png')
plt.savefig(PLOTDIR + 'Figure5.pdf')
plt.close('all')



# Make a final score table (colormap)
    
#cmap = plt.get_cmap('YlOrRd')
#cmap.set_under('white')
##levels = MaxNLocator(nbins=7).tick_values(5, 41.)
#levels = [5,10,15,20,25,30,35,40,45]
#norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
#
#fig = plt.figure(figsize=(14,5))
#ax = plt.subplot(111)
#    
#scorepix = ax.pcolormesh(np.abs(score_table), norm = norm, cmap = cmap)
#plt.xticks(rotation=90)
#
#cax = fig.add_axes([0.90,0.15,0.03,0.7])
#cbar = plt.colorbar(scorepix, extend = 'both', orientation = 'vertical', ticks=[5,10,20,30,40], cax = cax)
#cbar.set_label('Percent error', fontsize = fnt)
#cbar.ax.tick_params(labelsize=fnt)
#
#ax.set_xticks(np.arange(0.5,len(var_names)))
#ax.set_xticklabels(var1_caption, fontsize = fnt)
#ax.set_yticks([0.5,1.5])
#ax.set_yticklabels(['WACCM','MUSICA'], fontsize = fnt)
#
##ax.set_title('Lower Stratospheric Mixing Ratio Errors', fontsize = fnt)
#ax.set_title('Stratospheric Entry Mixing Ratio Errors', fontsize = fnt)
#ax.set_position([0.15,0.55,0.72,0.3])
#
#plt.savefig(PLOTDIR + 'error_scores_colormap.png')  
#plt.close('all')    


