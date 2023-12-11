
#########################################################
#
#  The following script is designed to calculate linear fits for the tracer relationships instead of binning the points in var2 space (cesm_TT_profiles)
#
#########################################################

# Import packages

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
from sklearn.utils import shuffle

# Define settings
sims = ['M32L'] #110L
cols = [['lightsalmon','pink','red']]#['orange','wheat','brown']]

var2_strs = ['CFC12','N2O']#,'CFC12','N2O']
var2_labs = ['CFC-12','N$_2$O','CFC-12','N$_2$O']
model_boxes = ['Sm','Sm','Lg','Lg']
var2_val_rngs = [[394,442],[265,292],[394,442],[265,292]]  # For now, just screen for the tropopause at a particular value
do_wgts = [1,1,1,1]  # For var2, indicate whether weighted regression should be used

var1_strs = [ 'CH3CCL3','HCFC141B','HCFC22','HCFC142B','CF2CLBR','CCL4', \
              'CFC11','CF3BR','CFC113',  'CFC12',     'N2O','CFC114',  'CFC115']
#var1_tlife = [  '1.3',   '1.6',    '5.2',     '8.9',  '11.3',    '15.3',   '17.3','1230', \
#               '1870', '4490',  '7620',  '11600',   '15600',   '19600','126000']     # Tropospheric lifetime      # 'CH3CL','CH3BR',
var1_slife = [  '37.7',     '72.3',  '161',   '212',   '33.5','50.6', \
                 '57', '73.5',  '88.4',   '95.5',     '116',   '191',     '997']     # Stratospheric lifetime      # '30.4', '26.3',  
var1_labs = ['CH$_3$CCl$_3$','HCFC-141b','HCFC-22','HCFC-142b','CF$_2$ClBr','CCl$_4$', \
              'CFC-11','CF$_3$Br','CFC-113','CFC-12','N$_2$O','CFC-114','CFC-115']                                 # 'CH$_3$Cl','CH$_3$Br',   
var1_us   = [     't',       't',     't',      't',       't',   't',\
                  't',    't',     't',      't',          'b',     't',     't']
fit_degs  = [       1,      1,        1,        1,        1,        1,         1,     1, \
                    1,      1,       1,        1,            1,       1,       1]
use_alls  = [       0,      0,        0,        0,        0,        0,         0,     0, \
                    0,      0,            0,        0,       0,       0,        0]
stars     = [    '* ',   '* ',       '',       '',     '* ',    '* ',\
                 '* ',   '* ',     '* ',     '* ',     ''     , '',     '' ]      # This controls whether an * shall be included to denote photolysis dominant loss                    

var1_caption = []
for vv in range(len(var1_strs)):
    var1_caption.append(stars[vv] + var1_labs[vv] + ' (' + var1_slife[vv] + ')')
    
DIST_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/heatmap_vars/'
PLOTDIR = '/home/wsmith/STRATOCLIM/plots/tracer_fits/'
WAS_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/WAS_merge/'
TXT_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/model_txt/'
MERGEFILES = sorted(glob.glob('/UTLS/Field_Data/STRATOCLIM_Data/merge_data/stratoclim*.nc'))


#['salmon','mistyrose','red']]# 

#sims = ['WITHVSL','NOVSL']
#cols = [['cyan','lightsteelblue','blue'],['pink','mistyrose','magenta']]

slope_score_table  = np.zeros((len(var2_strs), len(var1_strs)))  # Or replace with len(sims) if more than one is desired
offset_score_table = np.zeros((len(var2_strs), len(var1_strs)))
  
for s in range(len(sims)):

    sim = sims[s]
    col = cols[s]
    
    for v2 in range(len(var2_strs)):

        var2_str = var2_strs[v2]
        var2_lab = var2_labs[v2]
        var2_val_rng = var2_val_rngs[v2]
        do_wgt = do_wgts[v2]   
        model_box = model_boxes[v2]

        if model_box == 'Lg':
            lon_bnds = [30,130]
            lat_bnds = [18,40]
            model_sub = 500      
            #model_sub = 1         # For the VSL sensitivity runs
        elif model_box == 'Sm':
            lon_bnds = [75,95]
            lat_bnds = [18,32]
            model_sub = 500     #50
            #model_sub = 1         # For the VSL sensitivity runs
        prs_bnds = [50,200]
         
        for v1 in range(len(var1_strs)):
        
            var1_str = var1_strs[v1]
            var1_lab = var1_labs[v1]
            var1_u = var1_us[v1]
            fit_deg = fit_degs[v1]
            use_all = use_alls[v1]
        
            print('****', sim, '****', var2_str, 'vs', var1_str)

            ##################################################
            #  Read in the model values
            ##################################################        

            # Read in info for the variable choices
            var1_lims = hc.heatmap_info[var1_str]['lims']
            var2_lims = hc.heatmap_info[var2_str]['lims']
            var1_inst = hc.heatmap_info[var1_str]['instrument']
            var2_inst = hc.heatmap_info[var2_str]['instrument']
            var1_SEAC4RS_index = hc.heatmap_info[var1_str]['SEAC4RS_index']
            var1_SEAC4RS_scale = hc.heatmap_info[var1_str]['SEAC4RS_scale']
            var2_SEAC4RS_index = hc.heatmap_info[var2_str]['SEAC4RS_index']
            var2_SEAC4RS_scale = hc.heatmap_info[var2_str]['SEAC4RS_scale']            

            # Read in these files, from the save_heatmap_vars.py script
            DIST_SAVEFILE = DIST_SAVEDIR + 'heatmap_' + sim + '_' + var1_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '-' + str(prs_bnds) + '.pkl'
            with open(DIST_SAVEFILE, 'rb') as f: var1_asm, var1_trop, var1_strat = pkl.load(f)

            DIST_SAVEFILE = DIST_SAVEDIR + 'heatmap_' + sim + '_' + var2_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '-' + str(prs_bnds) + '.pkl'
            with open(DIST_SAVEFILE, 'rb') as f: var2_asm, var2_trop, var2_strat = pkl.load(f)

            # We don't need to separate for different regimes, so combine them all
            var1_raw = []
            var2_raw = []         
            var1_raw.extend(var1_asm)
            var1_raw.extend(var1_trop)
            var1_raw.extend(var1_strat)  # Note that these are not just the stratosphere points, they are supposed to be the points north of the jet
            var2_raw.extend(var2_asm)
            var2_raw.extend(var2_trop)
            var2_raw.extend(var2_strat)  # Theoretically, we could rewrite the reader to grab points above the tropopause if we wanted to
            
            var1_sim_strat = []
            var2_sim_strat = []
            var1_sim_trop  = []
            var2_sim_trop  = []
            for v in range(len(var1_raw)):
                curr_var2 = var2_raw[v]
                if (curr_var2 >= var2_val_rng[0] and curr_var2 <= var2_val_rng[1]) or use_all == 1:
                    var1_sim_strat.append(var1_raw[v])
                    var2_sim_strat.append(curr_var2)
                else:
                    var1_sim_trop.append(var1_raw[v])
                    var2_sim_trop.append(curr_var2)
                 
            # Save the model variables in txt files for Elliot           
            #var2_txt = open(TXT_SAVEDIR + var2_str + '.txt', 'w')
            #print(var2, file=var2_txt)
            #var1_txt = open(TXT_SAVEDIR + var1_str + '.txt', 'w')
            #print(var1, file=var1_txt)            
           
            # We need to subscript these arrays because their size prevents adequate fitting from happening
            # However, the issue is that subscripting seems to omit certain chunks of the data.  So we shuffle them 
            # while keeping the same structure
            var1_sim_strat, var2_sim_strat = shuffle(var1_sim_strat, var2_sim_strat, random_state=0)            

            # Fit the model relationship
            if fit_deg == 1:          m_sim, b_sim = np.polyfit(var2_sim_strat[::model_sub], var1_sim_strat[::model_sub], 1)
            elif fit_deg == 2: a_sim, b_sim, c_sim = np.polyfit(var2_sim_strat[::model_sub], var1_sim_strat[::model_sub], 2)
            
            #xplot_vals_model = np.arange(np.min(var2), np.max(var2), 0.01)
            #if var2_str == 'CFC12': xplot_vals_model = np.arange(var2_val_rng)
            #if var2_str == 'N2O': xplot_vals_model = np.arange(265,var2_tropp_val)
            xplot_vals_model = np.arange(var2_val_rng[0],var2_val_rng[1])

            ##################################################
            #  Get the StratoClim data
            ##################################################
            
            var1_SC_16_r = []
            var2_SC_16_r = []
            var1_SC_17_r = []
            var2_SC_17_r = []
            var1_SC_err_16_r = []
            var2_SC_err_16_r = []
            var1_SC_err_17_r = []
            var2_SC_err_17_r = []
            for curr_var1_inst in var1_inst:
                for curr_var2_inst in var2_inst:            
                    curr_var1_SC_16_r, curr_var2_SC_16_r, curr_var1_SC_17_r, curr_var2_SC_17_r, curr_var1_SC_err_16_r, curr_var2_SC_err_16_r, curr_var1_SC_err_17_r, curr_var2_SC_err_17_r = \
                                                              get_SC_arrays(curr_var1_inst, var1_str, curr_var2_inst, var2_str, prs_bnds, WAS_SAVEDIR, MERGEFILES)
                    var1_SC_16_r.extend(curr_var1_SC_16_r)
                    var2_SC_16_r.extend(curr_var2_SC_16_r)
                    var1_SC_17_r.extend(curr_var1_SC_17_r)
                    var2_SC_17_r.extend(curr_var2_SC_17_r)
                    var1_SC_err_16_r.extend(curr_var1_SC_err_16_r)
                    var2_SC_err_16_r.extend(curr_var2_SC_err_16_r)
                    var1_SC_err_17_r.extend(curr_var1_SC_err_17_r)
                    var2_SC_err_17_r.extend(curr_var2_SC_err_17_r)
            
            # Also read the altitude and tropopause arrays so we can screen for trop or strat
            #with open(WAS_SAVEDIR + 'AVIONIK:ALT_merge.pkl', 'rb') as f:
            #    z, zz, WAS_ALT, zzz = pkl.load(f)
            #with open(WAS_SAVEDIR + 'CLAMS_MET:GPH_CPT_merge.pkl', 'rb') as f:
            #    z, zz, WAS_GPH_CPT, zzz = pkl.load(f)
                
            # Process the observations according to desires

            var1_SC_16_all = []
            var2_SC_16_all = []
            var1_SC_16_trop = []
            var2_SC_16_trop = []
            var1_SC_16_strat = []
            var2_SC_16_strat = []
            var1_SC_17_all = []
            var2_SC_17_all = []
            var1_SC_17_trop = []
            var2_SC_17_trop = []
            var1_SC_17_strat = []
            var2_SC_17_strat = []
            var1_SC_17_strat_err = []
            var2_SC_17_strat_err = []
            var_SC_17_strat_err_wgt = []
            
            for msr in range(len(var1_SC_17_r)):     # 2017       
                curr_var1 = var1_SC_17_r[msr]
                curr_var2 = var2_SC_17_r[msr]
                curr_var1_err = var1_SC_err_17_r[msr]
                curr_var2_err = var2_SC_err_17_r[msr]
                #curr_trz  = WAS_ALT[msr] - (WAS_GPH_CPT[msr]*0.1)  # 0.1 necessary for some goofy reason
                if curr_var1 > -900 and curr_var2 > -900:  # Make sure there are no negative values
                    var1_SC_17_all.append(curr_var1)
                    var2_SC_17_all.append(curr_var2)
                    
                    if (curr_var2 >= var2_val_rng[0] and curr_var2 <= var2_val_rng[1]) or use_all == 1: 
                        var1_SC_17_strat.append(curr_var1)
                        var2_SC_17_strat.append(curr_var2)
                        var1_SC_17_strat_err.append(curr_var1_err)
                        var2_SC_17_strat_err.append(curr_var2_err)                        
                                                
                        if (curr_var1_err == 0 and curr_var2_err == 0) or do_wgt == 0: var_SC_17_strat_err_wgt.append(0.000001) # This will only happen when no errors are present for the entire specie, so the number shouldn't matter
                        else: var_SC_17_strat_err_wgt.append(1./((curr_var1_err**2)+(curr_var2_err**2)))
                        
                    else:
                        var1_SC_17_trop.append(curr_var1)
                        var2_SC_17_trop.append(curr_var2)                    


            #for msr in range(len(var1_SC_16_r)):     # 2016           
            #    curr_var1 = var1_SC_16_r[msr]
            #    curr_var2 = var2_SC_16_r[msr]
            #    #curr_trz  = WAS_ALT[msr] - (WAS_GPH_CPT[msr]*0.1)  # 0.1 necessary for some goofy reason
            #    if curr_var1 > -900 and curr_var2 > -900:
            #        var1_SC_16_all.append(curr_var1)
            #        var2_SC_16_all.append(curr_var2)
            #        if curr_var2 >= var2_tropp_val: 
            #            var1_SC_16_trop.append(curr_var1)
            #            var2_SC_16_trop.append(curr_var2)
            #        else:
            #            var1_SC_16_strat.append(curr_var1)
            #            var2_SC_16_strat.append(curr_var2) 
            
            # Make lists into numpy arrays so they play nice with other functions
            #var1_SC_16_all = np.asarray(var1_SC_16_all)
            #var2_SC_16_all = np.asarray(var2_SC_16_all)
            #var1_SC_16_trop = np.asarray(var1_SC_16_trop)
            #var2_SC_16_trop = np.asarray(var2_SC_16_trop)
            #var1_SC_16_strat = np.asarray(var1_SC_16_strat)
            #var2_SC_16_strat = np.asarray(var2_SC_16_strat)
            var1_SC_17_all = np.asarray(var1_SC_17_all)
            var2_SC_17_all = np.asarray(var2_SC_17_all)
            var1_SC_17_trop = np.asarray(var1_SC_17_trop)
            var2_SC_17_trop = np.asarray(var2_SC_17_trop)
            var1_SC_17_strat = np.asarray(var1_SC_17_strat)
            var2_SC_17_strat = np.asarray(var2_SC_17_strat)
            var_SC_17_strat_err_wgt = np.asarray(var_SC_17_strat_err_wgt)
            
            #var1_SC_1617_strat = np.append(var1_SC_16_strat, var1_SC_17_strat)
            #var2_SC_1617_strat = np.append(var2_SC_16_strat, var2_SC_17_strat)

            #xplot_vals_obs   = np.arange(np.min(var2_SC_1617_strat), var2_tropp_val,0.01)
            xplot_vals_obs   = np.arange(np.min(var2_SC_17_strat), np.max(var2_SC_17_strat),0.01)           

            ##################################################
            #  Perform a fit for observations
            ##################################################

            #a_obs, b_obs, c_obs = np.polyfit(var1_SC_17, var2_SC_17, 2)            
            if fit_deg == 1:          m_obs, b_obs = np.polyfit(var2_SC_17_strat, var1_SC_17_strat, 1, w=var_SC_17_strat_err_wgt)
            elif fit_deg == 2: a_obs, b_obs, c_obs = np.polyfit(var2_SC_17_strat, var1_SC_17_strat, 2, w=var_SC_17_strat_err_wgt)
            
            ##################################################
            #  Get SEAC4RS data
            ##################################################
   
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
                
            ##################################################
            #  Calculate the slope score
            ##################################################     

            #if fit_deg == 1: slope_error = m_sim/m_obs                               
            #elif fit_deg == 2: slope_error = (a_sim/b_sim)/(a_obs/b_obs)
            
            if fit_deg == 1: slope_error = (m_sim-m_obs)*100./m_obs
            #elif fit_deg == 2: slope_error = np.abs((a_sim/b_sim - a_obs/b_obs)*100./(a_obs/b_obs))

            ### New attempt - assign a negative number if the observed slope is bigger
            #if m_sim >= m_obs: slope_error = ((m_sim/m_obs)-1)*100.
            #if m_sim <= m_obs: slope_error = -(m_sim/m_obs)*100.

            slope_score_table[v2,v1] = slope_error                        

            ##################################################
            #  Calculate the stratosphere offset score
            ##################################################  
            
            if fit_deg == 1: 
                sim_midpt = m_sim*np.mean(var2_val_rng)+b_sim
                obs_midpt = m_obs*np.mean(var2_val_rng)+b_obs
                #print(sim_midpt)
                #print(obs_midpt)
                offset_error = (sim_midpt-obs_midpt)*100./(np.max(var1_SC_17_all)-np.min(var1_SC_17_all))
            
            #print(np.mean(var1_sim_strat))
            #print(np.mean(var1_SC_17_strat))
            #print(np.max(var2_SC_17_strat))
            #print(np.min(var2_SC_17_strat))
            #print(offset_error)
            
            offset_score_table[v2,v1] = offset_error
            
            ##################################################
            #  Calculate the boundary condition (tropp value) score
            ##################################################                 
                
            #if fit_deg == 1:                
            #    tropp_obs = m_obs*var2_tropp_val+b_obs
            #    tropp_sim = m_sim*var2_tropp_val+b_sim
            #elif fit_deg == 2: 
            #    tropp_obs = a_obs*(var2_tropp_val**2)+b_obs*var2_tropp_val+c_obs
            #    tropp_sim = a_sim*(var2_tropp_val**2)+b_sim*var2_tropp_val+c_sim            
            #    
            #score_table[1,v1] = np.abs(tropp_sim-tropp_obs)*100./tropp_obs                
                
            ##################################################
            #  Make a plot
            ##################################################        

            fnt = 14
            fig, ax = plt.subplots(figsize=(5,5))            
            #ax.set_position([0.15,0.2,0.5,0.6])                    

            # Plot model points
            ax.scatter(var1_sim_strat[::model_sub], var2_sim_strat[::model_sub], color = col[0], s=1)
            ax.scatter(var1_sim_trop[::model_sub], var2_sim_trop[::model_sub], color = col[1], s=1)
            if fit_deg == 1:   ax.plot(m_sim*xplot_vals_model+b_sim, xplot_vals_model, color=col[2], zorder=999)
            elif fit_deg == 2: ax.plot(a_sim*(xplot_vals_model**2)+b_sim*xplot_vals_model+c_sim, xplot_vals_model, color=col[2])     # Blue line for model fit
            
            # Plot StratoClim Obs
            #ax.scatter(var1_SC_17_all, var2_SC_17_all, color='k', s=2, zorder = 999)
            #ax.scatter(var1_SC_16_trop,   var2_SC_16_trop, color='gray', marker='s', s=3, zorder = 999)
            #ax.scatter(var1_SC_16_strat, var2_SC_16_strat, color='gray', marker='s', s=3, zorder = 999)            
            ax.scatter(var1_SC_17_trop,   var2_SC_17_trop, color='gray', s=5, zorder = 999)
            ax.scatter(var1_SC_17_strat, var2_SC_17_strat, color='black', s=5, zorder = 999)

            # Plot StratoClim error bars
            for pt in range(len(var1_SC_17_strat)):
                curr_var1     = var1_SC_17_strat[pt]
                curr_var1_err = var1_SC_17_strat_err[pt]
                curr_var2     = var2_SC_17_strat[pt]
                curr_var2_err = var2_SC_17_strat_err[pt]                                   
                if np.abs(curr_var1_err) <= 100: ax.plot([curr_var1-curr_var1_err,curr_var1+curr_var1_err], [curr_var2, curr_var2], color='k', linewidth=0.5)
                if np.abs(curr_var2_err) <= 100: ax.plot([curr_var1, curr_var1], [curr_var2-curr_var2_err,curr_var2+curr_var2_err], color='k', linewidth=0.5)
                    
            # Plot SEAC4RS Obs
            #ax.scatter(var1_SEAC4RS, var2_SEAC4RS, color='cyan', s=1, zorder = 999)
            
            if fit_deg == 1:   ax.plot(m_obs*xplot_vals_obs+b_obs, xplot_vals_obs, color='k')
            elif fit_deg == 2: ax.plot(a_obs*(xplot_vals_obs**2)+b_obs*xplot_vals_obs+c_obs, xplot_vals_obs, color='k')

            ax.set_ylim(np.flip(var2_lims))
            ax.set_xlim(var1_lims)
            #ax.set_ylim([550,350])
            ax.set_ylabel(var2_lab + ' (pptv)', fontsize = fnt)
            ax.set_xlabel(var1_lab + ' (pp' + var1_u + 'v)', fontsize = fnt)
            
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontweight('medium')
                tick.label1.set_fontsize(fnt)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontweight('medium')
                tick.label1.set_fontsize(fnt) 
        
            plt.figtext(0.18,0.24,'Slope score = ' + str(round(slope_error, 3)), fontsize=10, color='gray') 
            plt.figtext(0.18,0.21,'Offset score = ' + str(round(offset_error, 3)), fontsize=10, color='gray') 
            plt.figtext(0.50,0.85,var1_lab, fontsize=24, fontweight='bold', color='k')
            if fit_deg == 1:            
                if b_obs > 0: oper_obs = '+'
                else: oper_obs = ''
                if b_sim > 0: oper_sim = '+'
                else: oper_sim = ''
                plt.figtext(0.18,0.15,'Obs: ' + str(round(m_obs, 6)) + '(' + var2_lab + ')' + oper_obs + str(round(b_obs, 2)), fontsize=8, color='k')
                plt.figtext(0.18,0.18,'WACCM: ' + str(round(m_sim, 6)) + '(' + var2_lab + ')' + oper_sim + str(round(b_sim, 2)), fontsize=8, color=col[2])
            #elif fit_deg == 2:
            #    plt.figtext(0.18,0.15,'Obs: ' + str(round(a_obs, 6)) + '(' + var2_str + ')$^2$+' + str(round(b_obs, 2)) + '(' + var2_str + ')+' + str(round(c_obs, 2)), fontsize=8, color='k')
            #    plt.figtext(0.18,0.18,sim + ': ' + str(round(a_sim, 6)) + '(' + var2_str + ')$^2$+' + str(round(b_sim, 2)) + '(' + var2_str + ')+' + str(round(c_sim, 2)), fontsize=8, color=col[2])
            
            plt.tight_layout()
            plt.savefig(PLOTDIR + sim + '_' + var2_str + 'vs' + var1_str + '_' + model_box + 'box.png', dpi=300)
            plt.close('all')
                                
    fnt = 18
    
    ##################################################
    #  Make a slope score table
    ##################################################          
    print(slope_score_table)

    #### SCATTER VERSION

    fig = plt.figure(figsize=(14,6))
    ax = plt.subplot(111)

    ##### These make the table for the StratoClim paper
    #plt.scatter(np.arange(0.26,len(var1_strs)), slope_score_table[0,:], color='orange', edgecolors='k', s=150, marker='o', label=var2_labs[0] + ', ' + model_boxes[0] + ' box')
    #plt.scatter(np.arange(0.42,len(var1_strs)), slope_score_table[1,:], color='orange', edgecolors='k', s=150, marker='^', label=var2_labs[1] + ', ' + model_boxes[1] + ' box')
    #plt.scatter(np.arange(0.58,len(var1_strs)), slope_score_table[2,:], color='orange', edgecolors='k', s=150, marker='*', label=var2_labs[2] + ', ' + model_boxes[2] + ' box')
    #plt.scatter(np.arange(0.74,len(var1_strs)), slope_score_table[3,:], color='orange', edgecolors='k', s=150, marker='s', label=var2_labs[3] + ', ' + model_boxes[3] + ' box')

    ##### These were used to make a table for the NASA ECIP-ES proposal
    plt.scatter(np.arange(0.4,len(var1_strs)), slope_score_table[0,:], color='orange', edgecolors='k', s=200, marker='o', label=var2_labs[0] + ', ' + model_boxes[1] + ' box')
    plt.scatter(np.arange(0.6,len(var1_strs)), slope_score_table[1,:], color='orange', edgecolors='k', s=200, marker='^', label=var2_labs[1] + ', ' + model_boxes[2] + ' box')


    plt.plot([0,len(var1_strs)],[0,0],linestyle='--',color='k')
    plt.plot([0,len(var1_strs)],[20,20],linestyle='--',color=[0.8,0.8,0.8])
    plt.plot([0,len(var1_strs)],[40,40],linestyle='--',color=[0.8,0.8,0.8])
    plt.plot([0,len(var1_strs)],[-20,-20],linestyle='--',color=[0.8,0.8,0.8])
    plt.plot([0,len(var1_strs)],[-40,-40],linestyle='--',color=[0.8,0.8,0.8])    

    ax.set_title('(b) Lower Stratospheric Slope Errors', fontsize = fnt)

    ax.set_xticks(np.arange(0.5,len(var1_strs)))
    ax.set_xticklabels(var1_caption, fontsize = fnt)
    plt.xticks(rotation=90)

    plt.ylabel('Slope Error (%)', fontsize=fnt)
    plt.yticks(fontsize=fnt)
    plt.ylim(-60,60)
    plt.xlim(0,len(var1_strs))

    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(PLOTDIR + 'Figure9b.png')  
    plt.savefig(PLOTDIR + 'Figure9b.pdf')  
    plt.close('all') 

    
    ####  COLORMAP VERSION
    
    cmap = plt.get_cmap('bwr')
    #cmap.set_under('white')
    powers = [-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0]
    levels = np.power(np.zeros(len(powers))+1.6, powers)
    norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)

    fig = plt.figure(figsize=(14,5))
    ax = plt.subplot(111)
        
    scorepix = ax.pcolormesh(slope_score_table, norm = norm, cmap = cmap)
    plt.xticks(rotation=90)

    cax = fig.add_axes([0.89,0.15,0.03,0.7])
    cbar = plt.colorbar(scorepix, extend = 'both', orientation = 'vertical', ticks=[0.625,1,1.6], cax = cax)
    cbar.set_label('Slope Ratio', fontsize = fnt)
    cbar.ax.tick_params(labelsize=fnt)

    ax.set_xticks(np.arange(0.5,len(var1_strs)))
    ax.set_xticklabels(var1_caption, fontsize = fnt)
    ax.set_yticks([0.5])
    ax.set_yticklabels([var2_str + ', ' + model_box + ' box'], fontsize = fnt)
    #ax.set_yticklabels([var2_str], fontsize = fnt)

    ax.set_title('Errors for attenuation rates against ' + str(var2_str) + ' over '  + str(lon_bnds) + 'E, ' + str(lat_bnds) + 'N', fontsize = fnt)
    ax.set_position([0.15,0.62,0.72,0.15])

    plt.savefig(PLOTDIR + 'slope_error' + var2_str + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(prs_bnds) + '.png')  
    plt.close('all')
    
    
    ##################################################
    #  Make an offset score table
    ##################################################          
    
    
    #### SCATTER VERSION

    fig = plt.figure(figsize=(14,6))
    ax = plt.subplot(111)

    plt.scatter(np.arange(0.26,len(var1_strs)), offset_score_table[0,:], color='orange', s=150, marker='o', edgecolors='k', label=var2_labs[0] + ', ' + model_boxes[0] + ' box')
    plt.scatter(np.arange(0.42,len(var1_strs)), offset_score_table[1,:], color='orange', s=150, marker='^', edgecolors='k', label=var2_labs[1] + ', ' + model_boxes[1] + ' box')
    plt.scatter(np.arange(0.58,len(var1_strs)), offset_score_table[2,:], color='orange', s=150, marker='*', edgecolors='k', label=var2_labs[2] + ', ' + model_boxes[2] + ' box')
    plt.scatter(np.arange(0.74,len(var1_strs)), offset_score_table[3,:], color='orange', s=150, marker='s', edgecolors='k', label=var2_labs[3] + ', ' + model_boxes[3] + ' box')  

    plt.plot([0,len(var1_strs)],[0,0],linestyle='--',color='k')
    plt.plot([0,len(var1_strs)],[-40,-40],linestyle='--',color=[0.8,0.8,0.8])
    plt.plot([0,len(var1_strs)],[-20,-20],linestyle='--',color=[0.8,0.8,0.8])
    plt.plot([0,len(var1_strs)],[20,20],linestyle='--',color=[0.8,0.8,0.8])
    plt.plot([0,len(var1_strs)],[40,40],linestyle='--',color=[0.8,0.8,0.8])

    ax.set_title('(a) Lower Stratospheric Mixing Ratio Errors', fontsize = fnt)

    ax.set_xticks(np.arange(0.5,len(var1_strs)))
    ax.set_xticklabels(var1_caption, fontsize = fnt)
    plt.xticks(rotation=90)

    plt.ylabel('Mixing Ratio Error (%)', fontsize=fnt)
    plt.yticks(fontsize=fnt)
    plt.ylim(-50,50)
    plt.xlim(0,len(var1_strs))

    plt.legend(loc=3, fontsize=10)
    plt.tight_layout()
    plt.savefig(PLOTDIR + 'Figure9a.png')  
    plt.savefig(PLOTDIR + 'Figure9a.pdf')  
    plt.close('all')  


    #### COLORMAP VERSION
    
    cmap = plt.get_cmap('YlOrRd')
    cmap.set_under('white')
    levels = [5,10,15,20,25,30,35,40,45]
    norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)

    fig = plt.figure(figsize=(14,5))
    ax = plt.subplot(111)
        
    scorepix = ax.pcolormesh(offset_score_table, norm = norm, cmap = cmap)
    plt.xticks(rotation=90)

    #cax = fig.add_axes([0.90,0.15,0.03,0.7])
    cax = fig.add_axes([0.90,0.3,0.03,0.7])    
    cbar = plt.colorbar(scorepix, extend = 'both', orientation = 'vertical', ticks=[1,10,20,30,40], cax = cax)
    cbar.set_label('Mixing Ratio Error (%)', fontsize = fnt)
    cbar.ax.tick_params(labelsize=fnt)

    ax.set_xticks(np.arange(0.5,len(var1_strs)))
    ax.set_xticklabels(var1_caption, fontsize = fnt)
    ax.set_yticks([0.5])
    ax.set_yticklabels([var2_str + ', ' + model_box + ' box'], fontsize = fnt)
    #ax.set_yticklabels([var2_str], fontsize = fnt)        

    ax.set_title('Offset errors against ' + str(var2_str) + ' over '  + str(lon_bnds) + 'E, ' + str(lat_bnds) + 'N', fontsize = fnt)
    ax.set_position([0.15,0.62,0.72,0.15])

    plt.savefig(PLOTDIR + 'offset_error' + var2_str + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(prs_bnds) + '.png')  
    plt.close('all')    
