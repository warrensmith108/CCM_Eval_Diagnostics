
#############################################
# This script is intended to overplot histograms
# of species in tropopause-relative bins
#############################################

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import glob
from read_STRATOCLIM_merge import get_merge_array
from STRATOCLIM_plotter import calc_TT_pix
import heatmap_config as hc

def get_model_hist_data(var, sim, ztype, alt_rng, lon_bnds, lat_bnds, prs_lvl):
    #simfile_z = CESMSAVEDIR + var + ztype + 'profiles_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl'
    #with open(simfile_z, 'rb') as f:
    #    var_zprofs_in, var_zprofs_out, TROP_array_z, lev_z, gph_asma_z = pkl.load(f) 
    #var_profs_z = np.append(var_zprofs_in, var_zprofs_out, axis=0)
    
    simfile = CESMSAVEDIR + sim + '_' + var + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '_' + str(alt_rng[0]) + '-' + str(alt_rng[1]) + '.pkl'
    with open(simfile, 'rb') as f:
        var_data = pkl.load(f)     

    #lev_index0 = list(lev_z).index(alt_rng[0])
    #lev_index1 = list(lev_z).index(alt_rng[1])

    #var_data_in  = var_zprofs_in[:,lev_index0:lev_index1+1]
    #var_data_out = var_zprofs_out[:,lev_index0:lev_index1+1]
    
    return var_data
    
def get_STRATOCLIM_hist_data(var_str):
 
    DATAFILES = sorted(glob.glob('/UTLS/Field_Data/STRATOCLIM_Data/merge_data/stratoclim*.nc'))
    
    if isinstance(var_str, str): 
        if var_str[0:3] == 'WAS': 
            with open('/home/wsmith/STRATOCLIM/vars/WAS_merge/CLAMS_MET:GPH_CPT_merge.pkl','rb') as f: CPT_GPH = np.asarray(pkl.load(f))*0.1  # 0.1 Necessary for some reason
            with open('/home/wsmith/STRATOCLIM/vars/WAS_merge/AVIONIK:ALT_merge.pkl','rb') as f: AVIONIK_ALT = np.asarray(pkl.load(f))
            with open('/home/wsmith/STRATOCLIM/vars/WAS_merge/' + var_str + '_merge.pkl','rb') as f: var_data = pkl.load(f) 
        else: var_data, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_merge_array(var_str, DATAFILES)
    if isinstance(var_str, list): 
        var_data, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_merge_array(var_str[0], DATAFILES)
        for sub_v in range(len(var_str)-1):
            var_data_sub, var_unit_sub, dateint_sub, jd, AVIONIK_LON_sub, AVIONIK_LAT_sub, AVIONIK_ALT_sub, AVIONIK_PRS_sub, CPT_GPH_sub, CLAMS_THETA_sub, CPT_THETA_sub, LRT_THETA_sub, LRT_PRS_sub, ERA5_LRT_sub = get_merge_array(var_str[sub_v+1], DATAFILES)
            var_data.extend(var_data_sub)
            AVIONIK_ALT.extend(AVIONIK_ALT_sub)
            AVIONIK_PRS.extend(AVIONIK_PRS_sub)
            CPT_GPH.extend(CPT_GPH_sub)
            CLAMS_THETA.extend(CLAMS_THETA_sub)
            CPT_THETA.extend(CPT_THETA_sub)
            LRT_THETA.extend(LRT_THETA_sub)
            LRT_PRS.extend(LRT_PRS_sub)
            ERA5_LRT.extend(ERA5_LRT_sub)
            dateint.extend(dateint_sub)
    
    CPT_REL_ALT = (np.asarray(AVIONIK_ALT)*1.e-3) - np.asarray(CPT_GPH)
    CPT_REL_THETA = []
    
    #LRT_ALT = np.log(1000/np.asarray(LRT_PRS))*(287.06*255)/9.8  # Hypsometric equation, assuming a layer averaged temperature
    LRT_REL_ALT   = (np.asarray(AVIONIK_ALT)*1.e-3) - np.asarray(ERA5_LRT)
    
    #CPT_REL_THETA = (np.asarray(CLAMS_THETA) - np.asarray(CPT_THETA))
    
    return var_data, CPT_REL_ALT, LRT_REL_ALT, CPT_REL_THETA, CLAMS_THETA
            
    
# Define the relevant information
var_infos = [[['COLD:CO','AMICA:CO'],'CO',[0,140],40,0]]#,['FOZAN:O3','O3',[0,1000],40,1]]
#var_info = ['WAS:CH3BR','CH3BR',[4,9],50]
#var_info = ['WAS:HCFC22','HCFC22',[200,260],60]
#var_info = ['HAGAR:N2O','N2O',[260,360],20]
#lon_bnds = [75,95]
#lat_bnds = [18,32]
#lon_bnds = [0,180]
#lat_bnds = [0,60]
#lon_bnds = [257,277]
#lat_bnds = [23,37]
prs_lvl = 100
#alt_info = [-4,-3,-2,-1,0,1,2]    # Note these must lie on exact lev_z points for this script to work
#alt_info = [-10,0,10,20,30,40,50]
#alt_info = list(np.arange(-1.,2.1,2))

# Define the relevant file paths
CESMSAVEDIR = '/home/wsmith/STRATOCLIM/vars/z_bins/'
PLOTDIR = '/home/wsmith/STRATOCLIM/plots/alt_hists/'

#for a in alt_info:
#    alt_rng = [a,a+2.0]

#if ztype == 'troprelz': 
#    alt_rngs = [[-2.0,0.0],[0.0,2.0]]
#    alt_unit = ' km'
#if ztype == 'theta': 
#    alt_rngs = [[350.0,370.0],[390.0,410.0]]
#    alt_unit = ' K'

alt_rngs = [[0.0,2.0],[-2.0,0.0],[350.0,370.0]]
ztypes = ['troprelz','troprelz','theta']
alt_units = [' km',' km',' K']
bin_ys = [[0.7,0.8,0.9],[0.35,0.45,0.55],[0.0,0.1,0.2]]
markers = ['^','s','*']

for var_info in var_infos:

    fig, ax = plt.subplots()
        
    for aaa in range(len(alt_rngs)):

        alt_rng = alt_rngs[aaa]
        bin_y = bin_ys[aaa]
        ztype = ztypes[aaa]
        alt_unit = alt_units[aaa]
        marker = markers[aaa]

        # Get the StratoClim dataset and compute the histogram data
        var_data_obs_raw, CPT_REL_ALT_raw, LRT_REL_ALT_raw, CPT_REL_THETA_raw, THETA_raw = get_STRATOCLIM_hist_data(var_info[0])
        var_data_STRATOCLIM = []
        #CPT_REL_DATA = []
        for pt in range(len(var_data_obs_raw)):
            curr_var = var_data_obs_raw[pt]
            #if ztype == 'troprelz': curr_CPT = CPT_REL_ALT_raw[pt]
            if ztype == 'troprelz': curr_CPT = LRT_REL_ALT_raw[pt]
            elif ztype == 'tropreltheta': curr_CPT = CPT_REL_THETA_raw[pt]
            elif ztype == 'theta': curr_CPT = THETA_raw[pt]
            if curr_CPT >= alt_rng[0] and curr_CPT <= alt_rng[1]: 
                var_data_STRATOCLIM.append(curr_var)
                #CPT_REL_DATA.append(curr_CPT)

        # Get the SEAC4RS data
        #DATAFILE = '/home/wsmith/STRATOCLIM/SEAC4RS_Data/SEAC4RS-mrg60-ER2-merge.txt'
        #SEAC4RSdata  = np.loadtxt(DATAFILE, delimiter=',', skiprows=218)
        #var_SEAC4RS_index = hc.heatmap_info[var_info[1]]['SEAC4RS_index']        
        #var_SEAC4RS_raw = SEAC4RSdata[:,var_SEAC4RS_index]
        #alt = SEAC4RSdata[:,8]
        #TROP_z_raw = SEAC4RSdata[:,171]

        #var_data_SEAC4RS = []
        #for pt in range(len(TROP_z_raw)):
        #    curr_alt = alt[pt]
        #    curr_TROP_z = TROP_z_raw[pt]            
        #    if curr_TROP_z > 0 and var_SEAC4RS_raw[pt] > 0 and curr_alt-curr_TROP_z >= alt_rng[0] and curr_alt-curr_TROP_z <= alt_rng[1]:            
        #        #print(pt, curr_alt, curr_TROP_z, var_SEAC4RS_raw[pt])
        #        var_data_SEAC4RS.append(var_SEAC4RS_raw[pt])

        # Get the model datasets
        var_data_110L = get_model_hist_data(var_info[1], '110L', ztype, alt_rng, [75,95], [18,32], prs_lvl)
        #var_data_110L_nama_in, var_data_110L_nama_out = get_model_hist_data(var_info[1], '110L', ztype, alt_rng, [257,277], [23,37], prs_lvl)
        var_data_M32L = get_model_hist_data(var_info[1], 'M32L', ztype, alt_rng, [75,95], [18,32], prs_lvl)
        #var_data_M32L_nama_in, var_data_M32L_nama_out = get_model_hist_data(var_info[1], 'M32L', ztype, alt_rng, [257,277], [23,37], prs_lvl)
        #var_data_M58L = get_model_hist_data(var_info[1], 'M58L', ztype, alt_rng, [75,95], [18,32], prs_lvl)

        # Calculate metrics we wish to include
        var_mean_STRATOCLIM  = np.nanmean(var_data_STRATOCLIM)
        #var_mean_SEAC4RS  = np.nanmean(var_data_SEAC4RS)
        var_mean_M32L  = np.nanmean(var_data_M32L)
        #var_mean_M58L  = np.nanmean(var_data_M58L)
        var_mean_110L = np.nanmean(var_data_110L)
        #var_mean_M32L_out  = np.nanmean(var_data_M32L_nama_out)
        #var_mean_110L_out = np.nanmean(var_data_110L_nama_out)
        
        print(var_mean_STRATOCLIM, var_mean_M32L, var_mean_110L)

        var_mean_STRATOCLIM_str = str(round(var_mean_STRATOCLIM, 1))
        #var_mean_SEAC4RS_str    = str(round(var_mean_SEAC4RS, 1))
        var_mean_M32L_str  = str(round(var_mean_M32L, 1))
        #var_mean_M58L_str  = str(round(var_mean_M58L, 1))
        var_mean_110L_str  = str(round(var_mean_110L, 1))
        #var_mean_M32L_out_str = str(round(var_mean_M32L_out, 1))
        #var_mean_110L_out_str = str(round(var_mean_110L_out, 1))

        ########### Make a box and whisker plot
        def get_percentiles(data, pcts):
            data_pcts = []
            for pct in pcts:
                data_pcts.append(np.nanpercentile(data, pct))
                
            return data_pcts

        def box_plot(ax, pcts, mean, yval, color, marker):
            d_y = 0.05
                
            plt.plot([pcts[1],pcts[1]],[yval-d_y,yval+d_y], color = color)
            plt.plot([pcts[2],pcts[2]],[yval-d_y,yval+d_y], color = color)
            plt.plot([pcts[3],pcts[3]],[yval-d_y,yval+d_y], color = color)
            plt.plot([pcts[1],pcts[3]],[yval+d_y,yval+d_y], color = color)
            plt.plot([pcts[1],pcts[3]],[yval-d_y,yval-d_y], color = color)
            plt.plot([pcts[0],pcts[1]],[yval,yval], color = color)
            plt.plot([pcts[3],pcts[4]],[yval,yval], color = color)
            plt.scatter([mean],[yval],marker=marker,color=color)
            
            return pcts

        STRATOCLIM_pcts = get_percentiles(var_data_STRATOCLIM, [5,25,50,75,95])
        W110L_pcts = get_percentiles(var_data_110L, [5,25,50,75,95])
        M32L_pcts = get_percentiles(var_data_M32L, [5,25,50,75,95])
             
        obs_pcts = box_plot(ax, STRATOCLIM_pcts, var_mean_STRATOCLIM, bin_y[2], 'black', marker=marker)
        box_plot(ax, W110L_pcts, var_mean_110L, bin_y[1], 'orange', marker=marker)    
        box_plot(ax, M32L_pcts, var_mean_M32L, bin_y[0], 'red', marker=marker)
        
        #plt.text(obs_pcts[1], 0.22+bin_y[2], str(alt_rng) + alt_unit, color='gray', fontsize=12)

    plt.yticks([0.1,0.45,0.8])
    ax.set_yticklabels(['350-370K','0-2km below LRT','0-2km above LRT'], fontsize=12)
    ax.tick_params(axis='x', which='major', labelsize=12)

    if var_info[4] == 1: plt.xscale('log')
    plt.xlim(var_info[2])
    plt.ylim(-0.1,1)
    plt.title(var_info[1] + ' distributions over South Asia', fontsize = 12)
    plt.xlabel(var_info[1] + ' (ppbv)', fontsize=12)

    plt.tight_layout()
    #plt.savefig(PLOTDIR + var_info[1] + '_' + str(alt_rngs[0]) + '-' + str(alt_rngs[1]) + '_boxwhisker.png')
    #plt.savefig(PLOTDIR + var_info[1] + '_' + str(alt_rngs[0]) + '-' + str(alt_rngs[1]) + '_boxwhisker.pdf')
    
    plt.savefig(PLOTDIR + 'Figure4.png')
    plt.savefig(PLOTDIR + 'Figure4.pdf')
    
    plt.close('all')

    ########### Make a histogram plot
    fnt = 10

    plt.hist(var_data_STRATOCLIM, range = var_info[2], bins = var_info[3], density = True, color='black', histtype='step', linestyle='solid',\
             label = 'STRATOCLIM Obs (' + var_mean_STRATOCLIM_str + ')', linewidth=2)
    plt.hist(var_data_M32L, range = var_info[2], bins = var_info[3], density = True, color='red', histtype='step', linestyle='solid',\
             label = 'MUSICA 32L ASM (' + var_mean_M32L_str + ')', linewidth=2)
    #plt.hist(var_data_M58L, range = var_info[2], bins = var_info[3], density = True, color='blue', histtype='step', linestyle='solid',\
    #         label = 'MUSICA 58L ASM (' + var_mean_M58L_str + ')', linewidth=2)
    plt.hist(var_data_110L, range = var_info[2], bins = var_info[3], density = True, color='orange', histtype='step', linestyle='solid',\
             label = 'WACCM ASM (' + var_mean_110L_str + ')', linewidth=2)

    #plt.hist(var_data_SEAC4RS, range = var_info[2], bins = var_info[3], density = True, color='black', histtype='step', linestyle='dashed',\
    #         label = 'SEAC4RS Obs (' + var_mean_SEAC4RS_str + ')', linewidth=2)      
    #plt.hist(var_data_M32L_nama_out.flatten(), range = var_info[2], bins = var_info[3], density = True, color='red', histtype='step', linestyle='dashed',\
    #         label = 'MUSICA Bkgrnd (' + var_mean_M32L_out_str + ')', linewidth=2)
    #plt.hist(var_data_110L_nama_out.flatten(), range = var_info[2], bins = var_info[3], density = True, color='orange', histtype='step', linestyle='dashed',\
    #         label = 'WACCM Bkgrnd (' + var_mean_110L_out_str + ')', linewidth=2)
    plt.legend(fontsize=fnt)

    # Mean lines
    #plt.plot([var_mean_obs,var_mean_obs],[0,1], color='black', linestyle='--', linewidth=2)
    #plt.plot([var_mean_M32L,var_mean_M32L],[0,1], color='red', linestyle='--', linewidth=2)
    #plt.plot([var_mean_110L,var_mean_110L],[0,1], color='orange', linestyle='--', linewidth=2)

    #plt.xlim(var_info[2])
    #plt.xlim(10,100)
    #if alt_rng[0] >= 1.5: plt.ylim(0,0.1) 
    #else: plt.ylim(0,0.13)

    plt.title(var_info[1] + ' between ' + str(alt_rng[0]) + ' and ' + str(alt_rng[1]) + ' km relative to the local tropopause', fontsize = fnt)
    plt.xlabel(var_info[1] + ' (ppbv)', fontsize = fnt)
    plt.ylabel('Relative Frequency', fontsize = fnt)
    plt.xticks(fontsize=fnt)
    plt.yticks(fontsize=fnt)

    #plt.ylim(0,0.042)  # For the -5 to 3.5 plot

    plt.savefig(PLOTDIR + var_info[1] + '_' + str(alt_rng[0]) + '-' + str(alt_rng[1]) + '_hist.png', dpi=300)
    plt.close('all')

