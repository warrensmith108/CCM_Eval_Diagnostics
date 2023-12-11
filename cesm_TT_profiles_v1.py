
#############################################################
#   This script is intended to compute profiles of a tracer, 
#   plotted in another "tracer altitude" space
#############################################################

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import heatmap_config as hc
from cesm_TT_heatmap_v2 import get_SC_arrays
import glob
import pickle as pkl
import os
from sklearn.metrics import r2_score, mean_squared_error

def get_chem_prof(MODEL_SAVEDIR, var1_str, var2_str, var2_binbnds, var2_binmeans, sim, lon_bnds, lat_bnds, prs_bnds, var_choice='all'):

    MODEL_SAVEFILE = MODEL_SAVEDIR + var2_str + 'vs' + var1_str + '_' + sim + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '_' + var_choice + '_chemprofinfo.pkl'

    DIST_SAVEFILE = DIST_SAVEDIR + 'heatmap_' + sim + '_' + var1_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '-' + str(prs_bnds) + '.pkl'
    with open(DIST_SAVEFILE, 'rb') as f: var1_asm, var1_trop, var1_strat = pkl.load(f)

    DIST_SAVEFILE = DIST_SAVEDIR + 'heatmap_' + sim + '_' + var2_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '-' + str(prs_bnds) + '.pkl'
    with open(DIST_SAVEFILE, 'rb') as f: var2_asm, var2_trop, var2_strat = pkl.load(f)

    if var_choice == 'asm':
        var1 = var1_asm
        var2 = var2_asm
    elif var_choice == 'trop':
        var1 = var1_trop
        var2 = var2_trop
    elif var_choice == 'strat':
        var1 = var1_strat
        var2 = var2_strat
    elif var_choice == 'all':
        var1 = []
        var2 = []         
        var1.extend(var1_asm)
        var1.extend(var1_trop)
        var1.extend(var1_strat)
        var2.extend(var2_asm)
        var2.extend(var2_trop)
        var2.extend(var2_strat)    

    var1 = np.asarray(var1)
    var2 = np.asarray(var2)
    
    if os.path.exists(MODEL_SAVEFILE):
        
        print('File exists.....reading file....')
        
        with open(MODEL_SAVEFILE, 'rb') as f:               
            var2_binbnds, var2_binmeans, model_mean, model_std, model_10q, model_90q = pkl.load(f)
        
    else: 
    
        print('File not found.....proceeding with calculating the profile...')

        #############  Calculate profiles for each bin

        model_mean = []
        model_std = []
        model_10q = []
        model_90q = []
        #model_counts = []

        for z in range(len(var2_binmeans)):    
            
            curr_binmin = var2_binbnds[z]
            curr_binmax = var2_binbnds[z+1]
            
            print(curr_binmin, curr_binmax)
            
            MODELBIN_SAVEFILE = MODEL_SAVEDIR + var2_str + '_' + sim + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '_bininfo_' + str(curr_binmin) + '-' + str(curr_binmax) + '_' + var_choice + '.pkl'
            
            if os.path.exists(MODELBIN_SAVEFILE):
                with open(MODELBIN_SAVEFILE, 'rb') as f: curr_model_inds = pkl.load(f)
            else: 
                curr_model_inds = []            
                for pt in range(len(var1)):
                    curr_var2 = var2[pt]
                    if curr_var2 >= curr_binmin and curr_var2 <= curr_binmax:
                        curr_model_inds.append(pt)
                #model_counts.append(len(curr_model_inds))           
                curr_model_inds = np.asarray(curr_model_inds)        
                
                with open(MODELBIN_SAVEFILE, 'wb') as f: pkl.dump((curr_model_inds), f)
                
            if len(curr_model_inds) > 0: 
                model_mean.append(np.nanmean(var1[curr_model_inds]))
                model_std.append(np.nanstd(var1[curr_model_inds]))
                model_10q.append(np.nanquantile(var1[curr_model_inds], 0.1))
                model_90q.append(np.nanquantile(var1[curr_model_inds], 0.9))
            else: 
                model_mean.append(np.nan)
                model_std.append(np.nan)
                model_10q.append(np.nan)               
                model_90q.append(np.nan)  

        model_mean = np.asarray(model_mean)
        model_std = np.asarray(model_std)
        model_10q = np.asarray(model_10q)
        model_90q = np.asarray(model_90q)

        with open(MODEL_SAVEFILE, 'wb') as f: 
            pkl.dump((var2_binbnds, var2_binmeans, model_mean, model_std, model_10q, model_90q), f)
            
    return model_mean, model_std, model_10q, model_90q, var1, var2


#lon_bnds = [0,180]
#lat_bnds = [0,60]
#lon_bnds = [30,140]
#lat_bnds = [20,40]
lon_bnds = [75,95]
lat_bnds = [18,32]
prs_bnds = [50,200] # hPa
do_plot = 1   # This should be disabled when running a nohup command

MERGEFILES = sorted(glob.glob('/UTLS/Field_Data/STRATOCLIM_Data/merge_data/stratoclim*.nc'))
PLOTDIR = '/home/wsmith/STRATOCLIM/plots/chemical_profiles/'
MODEL_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/chemprofinfo/'
WAS_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/WAS_merge/'
DIST_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/heatmap_vars/'
  
#var2_list = ['CFC11','CFC12','N2O'] # Var 2 will be on the y axis, as the "chemical vertical coordinate"
var2_list = ['CFC12','CO']
var2_labs = ['CFC-12','CO']

var1_list = ['CH2BR2','CH3CL','CH3BR','CH3CCL3','HCFC141B','HCFC22','HCFC142B','CF2CLBR', 'CCL4', \
              'CFC11','CF3BR','CFC113','N2O','CFC114','CFC115']  
var1_tlife = [    '?',  '1.3', '1.6',    '5.2',      '8.9',  '11.3',    '15.3',    '17.3','1230', \
               '1870', '4490',  '7620','15600','19600','126000']     # Tropospheric lifetime
var1_slife = [    '?', '30.4', '26.3',  '37.7',     '72.3',   '161',     '212',   '33.5', '50.6', \
                 '57', '73.5',  '88.4','116',   '191',   '997']     # Stratospheric lifetime
var1_labs = ['CH$_2$Br$_2$','CH$_3$Cl','CH$_3$Br','CH$_3$CCl$_3$','HCFC-141b','HCFC-22','HCFC-142b','CF$_2$ClBr','CCl$_4$', \
              'CFC-11','CF$_3$Br','CFC-113','N$_2$O','CFC-114','CFC-115']
              
#var1_list = ['CH2BR2','CH3BR','CH3CL','CH3CCL3','HCFC141B','HCFC22','CF2CLBR','HCFC142B','CCL4', \
#              'CFC11','CF3BR','CFC113','N2O','CFC114','CFC115']  
#var1_tlife = [    '?',  '1.8', '1.57',    '6.1',    '10.7',  '13.0',      '-',    '19.3',   '-', \
#                  '-',    '-',     '-',    '-',  '-',     '-',     '-']     # Tropospheric lifetime
#var1_slife = [    '?', '26.3', '30.4',     '38',    '72.3',   '161',     '41',     '212',  '44', \
#                 '55', '73.5',  '94.5','123',   '191',   '664']     # Stratospheric lifetime
#var1_labs = ['CH$_2$Br$_2$','CH$_3$Br','CH$_3$Cl','CH$_3$CCl$_3$','HCFC-141b','HCFC-22','CF$_2$ClBr','HCFC-142b','CCl$_4$', \
#              'CFC-11','CF$_3$Br','CFC-113','N$_2$O','CFC-114','CFC-115']
              
                         
#fit_degs  = [       2,      1,      1,        1,         1,       1,        1,         1,     1, \
#                    1,      1,      1,      1,        1,   1,       1,      1]

#  These are the variable total lifetimes from WMO (2018), printed on 2022 AGU poster
#var1_life = [   '0.3',  '0.8',  '0.9',    '5.0',     '9.4',  '11.9',     '16',      '18',  '32',
#                 '52',   '72',   '93',   '102','123',   '189',   '540']

#error_score_table = np.zeros((3,len(var1_list)))
#
#var1_caption = []
#for vv in range(len(var1_list)):
#    var1_caption.append(var1_labs[vv] + ' (' + var1_tlife[vv] + ')')# (' + var1_slife[vv] + ')')

def make_TT_plot(ax, var2_str, var2_lab, var1_str, var1_lab):

    print('**************** ' + var2_str + ' vs ' + var1_str)

    var1_inst = hc.heatmap_info[var1_str]['instrument']
    var1_lims = hc.heatmap_info[var1_str]['lims']
    var1_unit = hc.heatmap_info[var1_str]['unit']
    var1_SEAC4RS_index = hc.heatmap_info[var1_str]['SEAC4RS_index']
    var1_SEAC4RS_scale = hc.heatmap_info[var1_str]['SEAC4RS_scale']
    var2_inst = hc.heatmap_info[var2_str]['instrument']
    var2_lims = hc.heatmap_info[var2_str]['lims']
    var2_bins = hc.heatmap_info[var2_str]['bins']
    var2_unit = hc.heatmap_info[var2_str]['unit']
    var2_SEAC4RS_index = hc.heatmap_info[var2_str]['SEAC4RS_index']
    var2_SEAC4RS_scale = hc.heatmap_info[var2_str]['SEAC4RS_scale']

    #fit_deg = fit_degs[v1]

    var2_binbnds  = np.arange(var2_lims[0],var2_lims[1]+0.01,(var2_lims[1]-var2_lims[0])/var2_bins)
    var2_binmeans = (var2_binbnds[1:]+var2_binbnds[0:-1])/2.

    ############# Get the model variable files

    #ne30_mean, ne30_std_all, ne30_10q_all, ne30_90q_all, ne30_var1_all, ne30_var2_all = get_chem_prof(MODEL_SAVEDIR, var1_str, var2_str, var2_binbnds, var2_binmeans, 'ne30',  lon_bnds, lat_bnds, prs_bnds)
    #ne120_mean, ne120_std, ne120_var1, ne120_var2 = get_chem_prof(MODEL_SAVEDIR, var1_str, var2_str, var2_binbnds, var2_binmeans, 'ne120', lon_bnds, lat_bnds, prs_bnds)
    WACCM_mean_all, WACCM_std_all, WACCM_10q_all, WACCM_90q_all, WACCM_var1_all, WACCM_var2_all = get_chem_prof(MODEL_SAVEDIR, var1_str, var2_str, var2_binbnds, var2_binmeans, '110L',  lon_bnds, lat_bnds, prs_bnds, var_choice='all')                      
    M32L_mean_all, M32L_std_all, M32L_10q_all, M32L_90q_all, M32L_var1_all, M32L_var2_all = get_chem_prof(MODEL_SAVEDIR, var1_str, var2_str, var2_binbnds, var2_binmeans, 'M32L',  lon_bnds, lat_bnds, prs_bnds, var_choice='all')
    #M58L_mean_all, M58L_std_all, M58L_10q_all, M58L_90q_all, M58L_var1_all, M58L_var2_all = get_chem_prof(MODEL_SAVEDIR, var1_str, var2_str, var2_binbnds, var2_binmeans, 'M58L',  lon_bnds, lat_bnds, prs_bnds, var_choice='all')

    #############  Get the StratoClim data

    var1_SC_16 = []
    var2_SC_16 = []
    var1_SC_17 = []
    var2_SC_17 = []
    var1_SC_err_16 = []
    var2_SC_err_16 = []
    var1_SC_err_17 = []
    var2_SC_err_17 = []            
    for curr_var1_inst in var1_inst:
        for curr_var2_inst in var2_inst:            
            curr_var1_SC_16, curr_var2_SC_16, curr_var1_SC_17, curr_var2_SC_17, curr_var1_SC_err_16, curr_var2_SC_err_16, curr_var1_SC_err_17, curr_var2_SC_err_17 = \
                                                      get_SC_arrays(curr_var1_inst, var1_str, curr_var2_inst, var2_str, prs_bnds, WAS_SAVEDIR, MERGEFILES)
            var1_SC_16.extend(curr_var1_SC_16)
            var2_SC_16.extend(curr_var2_SC_16)
            var1_SC_17.extend(curr_var1_SC_17)
            var2_SC_17.extend(curr_var2_SC_17)                    
            var1_SC_err_16.extend(curr_var1_SC_err_16)
            var2_SC_err_16.extend(curr_var2_SC_err_16)
            var1_SC_err_17.extend(curr_var1_SC_err_17)
            var2_SC_err_17.extend(curr_var2_SC_err_17)                    
            
    var1_SC_16 = np.asarray(var1_SC_16)
    var2_SC_16 = np.asarray(var2_SC_16)
    var1_SC_17 = np.asarray(var1_SC_17)
    var2_SC_17 = np.asarray(var2_SC_17)
    var1_SC_err_16 = np.asarray(var1_SC_err_16)
    var2_SC_err_16 = np.asarray(var2_SC_err_16)
    var1_SC_err_17 = np.asarray(var1_SC_err_17)
    var2_SC_err_17 = np.asarray(var2_SC_err_17)

    print(np.nanmin(var1_SC_err_17), np.nanmax(var1_SC_err_17))

    #####  Read SEAC4RS Data    
    DATAFILE = '/home/wsmith/STRATOCLIM/SEAC4RS_Data/SEAC4RS-mrg60-ER2-merge.txt'
    if var1_SEAC4RS_index != -999 and var2_SEAC4RS_index != -999:
        SEAC4RSdata  = np.loadtxt(DATAFILE, delimiter=',', skiprows=218)
        var1_SEAC4RS_raw = SEAC4RSdata[:,var1_SEAC4RS_index]
        var2_SEAC4RS_raw = SEAC4RSdata[:,var2_SEAC4RS_index]
        SEAC4RS_prs      = SEAC4RSdata[:,9]
        var1_SEAC4RS = []
        var2_SEAC4RS = []
        for SV in range(len(var1_SEAC4RS_raw)):
            curr_prs = SEAC4RS_prs[SV]
            if curr_prs >= prs_bnds[0] and curr_prs <= prs_bnds[1]:
                var1_SEAC4RS.append(var1_SEAC4RS_raw[SV]*var1_SEAC4RS_scale)
                var2_SEAC4RS.append(var2_SEAC4RS_raw[SV]*var2_SEAC4RS_scale)
    else:
        var1_SEAC4RS = [-999]
        var2_SEAC4RS = [-999]
        
        
    #if os.path.exists(SC_SAVEFILE):
    #
    #    with open(SC_SAVEFILE, 'rb') as f:               
    #        var2_binbnds, var2_binmeans, SC_mean, SC_std = pkl.load(f)
    #
    #else:

    SC_mean = []
    SC_std = []
    #ne30_var2binmeans = []
    #ne120_var2binmeans = []
    WACCM_var2binmeans = []
    M32L_var2binmeans = []
    #M58L_var2binmeans = []
    #ne30_mean_nonan = []
    #ne120_mean_nonan = []
    WACCM_mean_nonan = []
    M32L_mean_nonan = []
    #M58L_mean_nonan = []

    for z in range(len(var2_binmeans)): 

        curr_binmean = var2_binmeans[z]
        curr_binmin = var2_binbnds[z]
        curr_binmax = var2_binbnds[z+1]
        
        #curr_SC_inds = []
        #for msr in range(len(var1_SC)):
        #    curr_var1_SC = var1_SC[msr]
        #    curr_var2_SC = var2_SC[msr]
        #    if curr_var2_SC >= curr_binmin and curr_var2_SC <= curr_binmax and curr_var1_SC >= 0:
        #        curr_SC_inds.append(msr)            

        #curr_SC_inds = np.asarray(curr_SC_inds)            

        #if len(curr_SC_inds) > 0: 
        #    SC_mean.append(np.nanmean(var1_SC[curr_SC_inds]))
        #    SC_std.append(np.nanstd(var1_SC[curr_SC_inds]))
        #else: 
        #    SC_mean.append(np.nan)
        #    SC_std.append(np.nan)            

        # We only want to score the profiles if they are constrained by measurements.  Remove the below statement if that doesn't matter.
        if curr_binmax >= np.nanmin(var2_SC_17) and curr_binmin <= np.nanmax(var2_SC_17):
        
            #curr_ne30_mean = ne30_mean[z]
            #curr_ne120_mean = ne120_mean[z]
            curr_WACCM_mean = WACCM_mean_all[z]
            curr_M32L_mean = M32L_mean_all[z]
            #curr_M58L_mean = M58L_mean_all[z]
            #if ~np.isnan(curr_ne30_mean): 
            #    ne30_mean_nonan.append(curr_ne30_mean)
            #    ne30_var2binmeans.append(curr_binmean)
            #if ~np.isnan(curr_ne120_mean):
            #    ne120_mean_nonan.append(curr_ne120_mean)
            #    ne120_var2binmeans.append(curr_binmean)
            if ~np.isnan(curr_M32L_mean): 
                M32L_mean_nonan.append(curr_M32L_mean)
                M32L_var2binmeans.append(curr_binmean)
            #if ~np.isnan(curr_M58L_mean): 
            #    M58L_mean_nonan.append(curr_M58L_mean)
            #    M58L_var2binmeans.append(curr_binmean)
            if ~np.isnan(curr_WACCM_mean): 
                WACCM_mean_nonan.append(curr_WACCM_mean)
                WACCM_var2binmeans.append(curr_binmean)

    #print(np.nanmin(var2_SC), np.nanmax(var2_SC))
    #print(len(ne120_mean_nonan))                        

    SC_mean = np.asarray(SC_mean)
    SC_std = np.asarray(SC_std)
    #ne30_var2binmeans = np.asarray(ne30_var2binmeans)
    #ne120_var2binmeans = np.asarray(ne120_var2binmeans)
    WACCM_var2binmeans = np.asarray(WACCM_var2binmeans)
    M32L_var2binmeans  = np.asarray(M32L_var2binmeans)
    #M58L_var2binmeans  = np.asarray(M58L_var2binmeans)
    #ne30_mean_nonan = np.asarray(ne30_mean_nonan)
    #ne120_mean_nonan = np.asarray(ne120_mean_nonan)
    WACCM_mean_nonan = np.asarray(WACCM_mean_nonan)
    M32L_mean_nonan = np.asarray(M32L_mean_nonan)
    #M58L_mean_nonan = np.asarray(M58L_mean_nonan)

    #with open(SC_SAVEFILE, 'wb') as f: 
    #    pkl.dump((var2_binbnds, var2_binmeans, SC_mean, SC_std), f)

    var1_SC_nonan = []
    var2_SC_nonan = []            
    for xx in range(len(var1_SC_17)):
        curr_var1 = var1_SC_17[xx]
        curr_var2 = var2_SC_17[xx]
        if (~np.isnan(curr_var1)) and (~np.isnan(curr_var2)) and curr_var1 > -900 and curr_var2 > -900:
            var1_SC_nonan.append(curr_var1)
            var2_SC_nonan.append(curr_var2)
            
    ##  For the exponential fits, we will use the 2016 observations to help
    #if fit_deg == 2 and var2_lab != 'CO': 
    #    for xx in range(len(var1_SC_16)):
    #        curr_var1 = var1_SC_16[xx]
    #        curr_var2 = var2_SC_16[xx]
    #        if (~np.isnan(curr_var1)) and (~np.isnan(curr_var2)) and curr_var1 > -900 and curr_var2 > -900:
    #            var1_SC_nonan.append(curr_var1)
    #            var2_SC_nonan.append(curr_var2)           

    #####  Fit the observations using a curve
           
    #if fit_deg == 1:   ### Linear
    #    print('Fitting with a linear function...')
    #    SC_polyfit = np.polyfit(np.asarray(var2_SC_nonan), np.asarray(var1_SC_nonan), fit_deg)             
    #    #ne30_var1_predict = ne30_var2binmeans*SC_polyfit[0]+SC_polyfit[1]
    #    WACCM_var1_predict = WACCM_var2binmeans*SC_polyfit[0]+SC_polyfit[1]
    #    M32L_var1_predict  = M32L_var2binmeans*SC_polyfit[0]+SC_polyfit[1]
    #    #M58L_var1_predict  = M58L_var2binmeans*SC_polyfit[0]+SC_polyfit[1]
    #    best_r2 = r2_score(np.asarray(var1_SC_nonan), np.asarray(var2_SC_nonan)*SC_polyfit[0]+SC_polyfit[1])               

    #elif fit_deg == 2:    
    #    ne30_var1_predict   = (ne30_var2binmeans**2)*SC_polyfit[0]+ne30_var2binmeans*SC_polyfit[1]+SC_polyfit[2]
    #    ne120_var1_predict  = (ne120_var2binmeans**2)*SC_polyfit[0]+ne120_var2binmeans*SC_polyfit[1]+SC_polyfit[2]
    #    WACCM_var1_predict  = (WACCM_var2binmeans**2)*SC_polyfit[0]+WACCM_var2binmeans*SC_polyfit[1]+SC_polyfit[2]
    #    r2_measfit = str(round(r2_score(np.asarray(var1_SC_nonan), (np.asarray(var2_SC_nonan)**2)*SC_polyfit[0]+np.asarray(var2_SC_nonan)*SC_polyfit[1]+SC_polyfit[2]), 2))            
    #elif fit_deg == 2:  ### Exponential
    #    print('Fitting the curve with an exponential function...')
    #    a_guess = np.concatenate((np.arange(0.0000001,0.0001,0.0000001), np.arange(0.0001,0.1,0.0001), np.arange(0.1,10,0.1)))         
    #    b_guess = np.arange(0,10,0.1)
    #    c_value = np.quantile(var1_SC_nonan, 0.10)            
    #    best_a = -99999
    #    best_b = -99999
    #    best_r2 = -99999
    #    for ag in a_guess:                
    #        for bg in b_guess:
    #            var1_guess = ag * ((np.asarray(var2_SC_nonan)-np.min(var2_SC_nonan)) ** bg) + c_value                                               
    #            obs_r2 = r2_score(np.asarray(var1_guess), np.asarray(var1_SC_nonan))
    #            if obs_r2 > best_r2:
    #                best_r2 = obs_r2
    #                best_a = ag
    #                best_b = bg
    #                print(best_r2, best_a, best_b)
    #    var1_guess = best_a * ((np.asarray(var2_SC_nonan)-np.min(var2_SC_nonan)) ** best_b) + c_value 
    #    #ne30_var1_predict = best_a * ((np.asarray(ne30_var2binmeans)-np.nanmin(var2_SC_nonan)) ** best_b) + c_value    
    #    WACCM_var1_predict = best_a * ((np.asarray(WACCM_var2binmeans)-np.nanmin(var2_SC_nonan)) ** best_b) + c_value    
    #    M32L_var1_predict  = best_a * ((np.asarray(M32L_var2binmeans) -np.nanmin(var2_SC_nonan)) ** best_b) + c_value 
    #    #M58L_var1_predict  = best_a * ((np.asarray(M58L_var2binmeans) -np.nanmin(var2_SC_nonan)) ** best_b) + c_value 

    #score_ne30 = str(round(np.sum(np.abs(ne30_var1_predict-ne30_mean_nonan))*100./(len(ne30_mean_nonan)*np.nanmean(var1_SC_nonan)),1))
    #score_M32L = str(round(np.sum(np.abs(M32L_var1_predict-M32L_mean_nonan))*100./(len(M32L_mean_nonan)*np.nanmean(var1_SC_nonan)),1))
    #score_110L = str(round(np.sum(np.abs(WACCM_var1_predict-WACCM_mean_nonan))*100./(len(WACCM_mean_nonan)*np.nanmean(var1_SC_nonan)),1))   
    #score_M58L = str(round(np.sum(np.abs(M58L_var1_predict-M58L_mean_nonan))*100./(len(M58L_mean_nonan)*np.nanmean(var1_SC_nonan)),1))   

    #print(M32L_var2binmeans)
    #print(M32L_var1_predict)
    #print(var2_SC_nonan)
    #print(np.asarray(M32L_var2binmeans) -np.nanmin(var2_SC_nonan))
    #print(best_a)
    #print(best_b)
    #print(c_value)

    #error_score_table[0,v1] = score_110L
    #error_score_table[1,v1] = score_M32L
    #error_score_table[2,v1] = score_M58L

    #score_ne30   = str(round(np.sum(np.abs(ne30_var1_predict-ne30_mean_nonan))*100./(len(ne30_mean_nonan)*np.nanmax(var1_SC)),1))
    #score_ne120  = str(round(np.sum(np.abs(ne120_var1_predict-ne120_mean_nonan))*100./(len(ne120_mean_nonan)*np.nanmax(var1_SC)),1))

    # Test plot of the predicted and expected scores
    #plt.scatter(ne120_var1, ne120_var1_predict, s=0.5, color='blue')
    #plt.plot([-9999,9999],[-9999,9999],color='k')
    #plt.xlabel('MUSICA ' + var1_str + ' ACTUAL')
    #plt.ylabel('MUSICA ' + var1_str + ' PREDICTED based on observed relationship')
    #plt.xlim(var1_lims)
    #plt.ylim(var1_lims)
    #plt.figtext(0.6,0.2,'Mean Abs Error = ' + err_ne120)
    #plt.savefig('/home/wsmith/STRATOCLIM/plots/predict_scatter_' + var1_str + '.png')
    #plt.close('all')

    ##############   Make a profile plot

    #fig, ax = plt.subplots(figsize=(6,4))            
    #ax.set_position([0.15,0.2,0.5,0.6])

    #ax.plot(ne30_mean, var2_binmeans, color='blue')
    #plt.plot(ne30_10q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'blue')
    #plt.plot(ne30_90q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'blue')
    #ax.plot(ne120_mean, var2_binmeans, color='red')
    #plt.plot(ne120_mean-ne120_std, var2_binmeans, linestyle='--', color='red')
    #plt.plot(ne120_mean+ne120_std, var2_binmeans, linestyle='--', color='red')
    #ax.plot(WACCM_mean_all, var2_binmeans, color='orange')
    ax.fill_betweenx(var2_binmeans, WACCM_10q_all, WACCM_90q_all, color='orange', alpha=0.5)
    #plt.plot(WACCM_10q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'orange')
    #plt.plot(WACCM_90q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'orange')
    #plt.plot(WACCM_mean-WACCM_std, var2_binmeans, linewidth = 0.5, linestyle='--', color='orange')
    #plt.plot(WACCM_mean+WACCM_std, var2_binmeans, linewidth = 0.5, linestyle='--', color='orange')
    #ax.plot(M32L_mean_all, var2_binmeans, color='red')        
    ax.fill_betweenx(var2_binmeans, M32L_10q_all, M32L_90q_all, color='red', alpha=0.5)                
    #plt.plot(M32L_10q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'red')
    #plt.plot(M32L_90q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'red')
    #plt.plot(M32L_mean-M32L_std, var2_binmeans, linewidth = 0.5, linestyle='--', color='red')
    #plt.plot(M32L_mean+M32L_std, var2_binmeans, linewidth = 0.5, linestyle='--', color='red')
    #ax.plot(M58L_mean_all, var2_binmeans, color='blue')                
    #plt.plot(M58L_10q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'blue')
    #plt.plot(M58L_90q_all, var2_binmeans, linewidth = 0.5, linestyle = '--', color = 'blue')

    # Plot StratoClim
    #if fit_deg == 1: ax.plot(np.arange(np.nanmin(var2_SC_17),np.nanmax(var2_SC_17))*SC_polyfit[0]+SC_polyfit[1],np.arange(np.nanmin(var2_SC_17),np.nanmax(var2_SC_17)), color='black')
    #if fit_deg == 2: ax.plot([x for _,x in sorted(zip(var2_SC_nonan,var1_guess))], sorted(var2_SC_nonan), color='black')
    #if var1_str != 'N2O': ax.scatter(var1_SC_16, var2_SC_16, color='gray', s=5, zorder = 999)
    ax.scatter(var1_SC_17, var2_SC_17, color='black', s=5, zorder = 999)

    # Plot StratoClim error bars
    for pt in range(len(var1_SC_17)):
        curr_var1     = var1_SC_17[pt]
        curr_var1_err = var1_SC_err_17[pt]
        curr_var2     = var2_SC_17[pt]
        curr_var2_err = var2_SC_err_17[pt]
                            
        if np.abs(curr_var1_err) <= 100: ax.plot([curr_var1-curr_var1_err,curr_var1+curr_var1_err], [curr_var2, curr_var2], color='k', linewidth=0.3)
        if np.abs(curr_var2_err) <= 100: ax.plot([curr_var1, curr_var1], [curr_var2-curr_var2_err,curr_var2+curr_var2_err], color='k', linewidth=0.3)                    

    # Plot SEAC4RS
    ax.scatter(var1_SEAC4RS, var2_SEAC4RS, color='gray', facecolor='none', s=0.5, zorder=999) 

    ax.set_xlim(var1_lims)                    
    if var2_str == 'N2O': ax.set_ylim([340,260])
    elif var2_str == 'CFC11': ax.set_ylim([240,130])
    elif var2_str == 'CO': ax.set_ylim([160,0])
    else: ax.set_ylim(np.flip(var2_lims))

    ax.set_xlabel(var1_lab + ' (' + var1_unit + ')', fontsize = 16)
    ax.set_ylabel(var2_lab + ' (' + var2_unit + ')', fontsize = 16)

    ax.tick_params(labelsize=16)

    #plt.figtext(0.17,0.30,'\u03C4$_{trop}$: ' + var1_tlt + ' y', fontsize=10)
    #plt.figtext(0.17,0.26,'\u03C4$_{strat}$: ' + var1_slt + ' y', fontsize=10)
    #plt.figtext(0.17,0.22,'Obs r$^2$ = ' + str(round(best_r2, 2)), fontsize = 10, color='black')
    #plt.figtext(0.7,0.48,'CAM-Chem (Score=' + score_ne30 + ')', fontsize = 8, color = 'blue')
    #plt.figtext(0.7,0.45,'MUSICA (Score=' + score_ne120 + ')', fontsize = 8, color = 'red')
    #plt.figtext(0.7,0.45,'MUSICA 32L (Score=' + score_M32L + ')', fontsize = 8, color = 'red')
    #plt.figtext(0.7,0.48,'MUSICA 58L (Score=' + score_M58L + ')', fontsize = 8, color = 'blue')
    #plt.figtext(0.7,0.42,'WACCM (Score=' + score_110L + ')', fontsize = 8, color = 'orange')
    #plt.figtext(0.7,0.39,var2_str + ' from ' + str(var2_inst), fontsize = 8, color = 'black')
    #plt.figtext(0.7,0.36,var1_str + ' from ' + str(var1_inst), fontsize = 8, color = 'black')

    #for tick in ax.xaxis.get_major_ticks():
    #    tick.label1.set_fontweight('medium')
    #    tick.label1.set_fontsize(fnt)
    #for tick in ax.yaxis.get_major_ticks():
    #    tick.label1.set_fontweight('medium')
    #    tick.label1.set_fontsize(fnt)
        
    #plt.savefig(PLOTDIR + str(var2_inst) + var2_str + 'vs' + str(var1_inst) + var1_str + '_profiles_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(prs_bnds) + '.png', dpi=300)
    #plt.close('all')


# Make Figure 7
fig = plt.figure(figsize=(12,10))
ax1  = plt.subplot(2,2,1)
ax2  = plt.subplot(2,2,2)
ax3  = plt.subplot(2,2,3)
ax4  = plt.subplot(2,2,4)


make_TT_plot(ax1, 'CFC12', 'CFC-12', 'CH3CL', 'CH$_3$Cl')
make_TT_plot(ax2, 'CO', 'CO', 'CH3BR', 'CH$_3$Br')
make_TT_plot(ax3, 'CFC12', 'CFC-12', 'CH3CL', 'CH$_3$Cl')
make_TT_plot(ax4, 'CO', 'CO', 'CH3BR', 'CH$_3$Br')

plt.figtext(0.10,0.60,'(a)', fontweight='bold', fontsize=16)
plt.figtext(0.60,0.60,'(b)', fontweight='bold', fontsize=16)
plt.figtext(0.10,0.10,'(c)', fontweight='bold', fontsize=16)
plt.figtext(0.60,0.10,'(d)', fontweight='bold', fontsize=16)

plt.tight_layout()
plt.savefig(PLOTDIR + 'Figure7.png', dpi=300)
plt.savefig(PLOTDIR + 'Figure7.pdf', dpi=300)
plt.close('all')


    
    ##############  Make a pixel plot with all the scores combined
    
    #fnt = 18
    #
    #cmap = plt.get_cmap('YlOrRd')
    #levels = MaxNLocator(nbins=150).tick_values(0, 15)
    #norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    #
    #fig = plt.figure(figsize=(14,5))
    #ax = plt.subplot(111)
    #    
    #scorepix = ax.pcolormesh(error_score_table, norm = norm, cmap = cmap)
    #plt.xticks(rotation=90)
    #
    #cax = fig.add_axes([0.90,0.15,0.03,0.7])
    #cbar = plt.colorbar(scorepix, extend = 'both', orientation = 'vertical', ticks=[0,5,10,15,20], cax = cax)
    #cbar.set_label('Mean absolute percent error', fontsize = fnt)
    #cbar.ax.tick_params(labelsize=fnt)
    #
    #ax.set_xticks(np.arange(0.5,len(var1_labs)))
    #ax.set_xticklabels(var1_caption, fontsize = fnt)    
    #ax.set_yticks([0.5,1.5,2.5])
    #ax.set_yticklabels(['WACCM','MUSICA 32L','MUSICA 58L'], fontsize = fnt)
    #
    #ax.set_title('Errors for ' + str(var2_str) + ' relationships over '  + str(lon_bnds) + 'E, ' + str(lat_bnds) + 'N, ' + str(prs_bnds) + ' hPa', fontsize = fnt)
    #
    #ax.set_position([0.15,0.55,0.72,0.3])
    #
    #plt.savefig(PLOTDIR + 'error_scores_' + var2_str + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(prs_bnds) + '.png')







            #############   Make a 1-to-1 plot

            #print(SC_mean)
            #print(model_mean)
                
            #SC_mean_cc = []
            #model_mean_cc = []

            #plt.scatter(SC_mean, model_mean)
            #plt.plot([-9999,9999],[-9999,9999],color='k')

            #plt.xlim(var1_lims)
            #plt.ylim(var1_lims)

            #plt.xlabel(var1_str + ' (' + var1_unit + ') Observations Mean')#, fontsize = fnt)
            #plt.ylabel(var1_str + ' (' + var1_unit + ') Model Mean')#, fontsize = fnt)

            #curr_corr = str(round(np.ma.corrcoef(np.ma.masked_invalid(SC_mean), np.ma.masked_invalid(model_mean))[0,1], 3))
            #curr_corr = str(round(np.corrcoef(SC_mean_cc, model_mean_cc)[0,1], 3))
            #print(curr_corr)

            #plt.figtext(0.3,0.23,'Model: ' + sim, fontsize = 10, color = 'orange')
            #plt.figtext(0.3,0.20,'Coordinate = ' + var2_str + ' from ' + str(var2_inst), fontsize = 10, color = 'blue')
            #plt.figtext(0.3,0.17,'vs ' + var1_str + ' from ' + str(var1_inst), fontsize = 10, color = 'blue')
            #plt.figtext(0.3,0.14,'Correlation = ' + curr_corr)

            #plt.savefig(PLOTDIR + var2_str + 'vs' + var1_str + '_MeanScatter.png')
            #plt.close('all')


