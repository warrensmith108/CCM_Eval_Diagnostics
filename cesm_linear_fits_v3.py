
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

def get_model_vals(DIST_SAVEDIR, sim, var1_str, var2_str, var2_val_rng, model_box, prs_bnds):

    if model_box == 'Lg':
        lon_bnds = [30,130]
        lat_bnds = [18,40]
        model_sub = 500      
    elif model_box == 'Sm':
        lon_bnds = [75,95]
        lat_bnds = [18,32]
        if sim == '110L': model_sub = 100      
        elif sim == 'M32L': model_sub = 1500

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
        if (curr_var2 >= var2_val_rng[0] and curr_var2 <= var2_val_rng[1]):
            var1_sim_strat.append(var1_raw[v])
            var2_sim_strat.append(curr_var2)
        else:
            var1_sim_trop.append(var1_raw[v])
            var2_sim_trop.append(curr_var2)
                         
    # We need to subscript these arrays because their size prevents adequate fitting from happening
    # However, the issue is that subscripting seems to omit certain chunks of the data.  So we shuffle them 
    # while keeping the same structure
    var1_sim_strat, var2_sim_strat = shuffle(var1_sim_strat, var2_sim_strat, random_state=0)         

    # Fit the model relationship
    m_sim, b_sim = np.polyfit(var2_sim_strat[::model_sub], var1_sim_strat[::model_sub], 1)
    #elif fit_deg == 2: a_sim, b_sim, c_sim = np.polyfit(var2_sim_strat[::model_sub], var1_sim_strat[::model_sub], 2)
 
    return var1_sim_strat, var2_sim_strat, var1_sim_trop, var2_sim_trop, m_sim, b_sim, model_sub
    
    
def get_obs_vals(WAS_SAVEDIR, MERGEFILES, prs_bnds, var1_str, var2_str, var2_val_rng):

    var1_inst = hc.heatmap_info[var1_str]['instrument']
    var2_inst = hc.heatmap_info[var2_str]['instrument']
    var1_SEAC4RS_index = hc.heatmap_info[var1_str]['SEAC4RS_index']
    var1_SEAC4RS_scale = hc.heatmap_info[var1_str]['SEAC4RS_scale']
    var2_SEAC4RS_index = hc.heatmap_info[var2_str]['SEAC4RS_index']
    var2_SEAC4RS_scale = hc.heatmap_info[var2_str]['SEAC4RS_scale']
    
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
        
    # Process the observations according to desires

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
        if curr_var1 > -900 and curr_var2 > -900:  # Make sure there are no negative values
            var1_SC_17_all.append(curr_var1)
            var2_SC_17_all.append(curr_var2)
            
            if (curr_var2 >= var2_val_rng[0] and curr_var2 <= var2_val_rng[1]): 
                var1_SC_17_strat.append(curr_var1)
                var2_SC_17_strat.append(curr_var2)
                var1_SC_17_strat_err.append(curr_var1_err)
                var2_SC_17_strat_err.append(curr_var2_err)                        
                                        
                if (curr_var1_err == 0 and curr_var2_err == 0): var_SC_17_strat_err_wgt.append(0.000001) # This will only happen when no errors are present for the entire specie, so the number shouldn't matter
                else: var_SC_17_strat_err_wgt.append(1./((curr_var1_err**2)+(curr_var2_err**2)))
                
            else:
                var1_SC_17_trop.append(curr_var1)
                var2_SC_17_trop.append(curr_var2)                    
    
    # Make lists into numpy arrays so they play nice with other functions
    var1_SC_17_all = np.asarray(var1_SC_17_all)
    var2_SC_17_all = np.asarray(var2_SC_17_all)
    var1_SC_17_trop = np.asarray(var1_SC_17_trop)
    var2_SC_17_trop = np.asarray(var2_SC_17_trop)
    var1_SC_17_strat = np.asarray(var1_SC_17_strat)
    var2_SC_17_strat = np.asarray(var2_SC_17_strat)
    var_SC_17_strat_err_wgt = np.asarray(var_SC_17_strat_err_wgt)
              
    ##################################################
    #  Perform a fit for observations
    ##################################################
          
    m_obs, b_obs = np.polyfit(var2_SC_17_strat, var1_SC_17_strat, 1, w=var_SC_17_strat_err_wgt)
    #elif fit_deg == 2: a_obs, b_obs, c_obs = np.polyfit(var2_SC_17_strat, var1_SC_17_strat, 2, w=var_SC_17_strat_err_wgt)
    
    ##################################################
    #  Get SEAC4RS data - not currently used
    ##################################################

    #DATAFILE = '/home/wsmith/STRATOCLIM/SEAC4RS_Data/SEAC4RS-mrg60-ER2-merge.txt'
    #if var1_SEAC4RS_index != -999 and var2_SEAC4RS_index != -999:
    #    SEAC4RSdata  = np.loadtxt(DATAFILE, delimiter=',', skiprows=218)
    #    var1_SEAC4RS_raw = SEAC4RSdata[:,var1_SEAC4RS_index]
    #    var2_SEAC4RS_raw = SEAC4RSdata[:,var2_SEAC4RS_index]
    #    SEAC4RS_prs      = SEAC4RSdata[:,9]
    #    var1_SEAC4RS = []
    #    var2_SEAC4RS = []
    #    for SV in range(len(var1_SEAC4RS_raw)):
    #        curr_prs = SEAC4RS_prs[SV]
    #        if curr_prs >= prs_bnds[0] and curr_prs <= prs_bnds[1]:
    #            var1_SEAC4RS.append(var1_SEAC4RS_raw[SV]*var1_SEAC4RS_scale)
    #            var2_SEAC4RS.append(var2_SEAC4RS_raw[SV]*var2_SEAC4RS_scale)
    #else:
    #    var1_SEAC4RS = [-999]
    #    var2_SEAC4RS = [-999]

    return var1_SC_17_trop, var2_SC_17_trop, var1_SC_17_strat, var2_SC_17_strat, var1_SC_17_strat_err, var2_SC_17_strat_err, m_obs, b_obs 


def plot_on_axis(ax, var1_str, var2_str, var1_lab, var2_lab, var1_let, var1_SC_17_trop, var2_SC_17_trop, var1_SC_17_strat, var2_SC_17_strat, var1_SC_17_strat_err, var2_SC_17_strat_err,\
                var1_sim1_trop, var2_sim1_trop, var1_sim1_strat, var2_sim1_strat, var1_sim2_trop, var2_sim2_trop, var1_sim2_strat, var2_sim2_strat,\
                var2_val_rng, var1_u, var2_u, model_sub_sim1, model_sub_sim2, m_obs, b_obs, m_sim1, b_sim1, m_sim2, b_sim2, col1, col2):

    ##################################################
    #  Make a plot
    ##################################################        

    fnt = 14
    #fig, ax = plt.subplots(figsize=(5,5))            
    #ax.set_position([0.15,0.2,0.5,0.6])

    # Read in info for the variable choices
    var1_lims = hc.heatmap_info[var1_str]['lims']
    var2_lims = hc.heatmap_info[var2_str]['lims']  

    xplot_vals_obs   = np.arange(np.min(var2_SC_17_strat), np.max(var2_SC_17_strat),0.01)
    xplot_vals_model = np.arange(var2_val_rng[0],var2_val_rng[1])    

    # Plot model points
    ax.scatter(var1_sim1_strat[::model_sub_sim1], var2_sim1_strat[::model_sub_sim1], color = col1[0], s=0.5)
    ax.scatter(var1_sim1_trop[::model_sub_sim1], var2_sim1_trop[::model_sub_sim1], color = col1[1], s=0.5)
    ax.plot(m_sim1*xplot_vals_model+b_sim1, xplot_vals_model, color=col1[2], zorder=999)
    
    ax.scatter(var1_sim2_strat[::model_sub_sim2], var2_sim2_strat[::model_sub_sim2], color = col2[0], s=0.5)
    ax.scatter(var1_sim2_trop[::model_sub_sim2], var2_sim2_trop[::model_sub_sim2], color = col2[1], s=0.5)
    ax.plot(m_sim2*xplot_vals_model+b_sim2, xplot_vals_model, color=col2[2], zorder=999)    
    #elif fit_deg == 2: ax.plot(a_sim*(xplot_vals_model**2)+b_sim*xplot_vals_model+c_sim, xplot_vals_model, color=col[2])     # Blue line for model fit

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

    ax.plot(m_obs*xplot_vals_obs+b_obs, xplot_vals_obs, color='k')
    #elif fit_deg == 2: ax.plot(a_obs*(xplot_vals_obs**2)+b_obs*xplot_vals_obs+c_obs, xplot_vals_obs, color='k')

    ax.set_ylim(np.flip(var2_lims))
    ax.set_xlim(var1_lims)
    #ax.set_ylim([550,350])
    ax.set_ylabel(var2_lab + ' (pp' + var2_u + 'v)', fontsize = fnt)
    ax.set_xlabel(var1_lab + ' (pp' + var1_u + 'v)', fontsize = fnt)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontweight('medium')
        tick.label1.set_fontsize(fnt)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontweight('medium')
        tick.label1.set_fontsize(fnt) 
        
    ax.text(x=0.41,y=0.88,s=var1_let + ' ' + var1_lab, fontsize=fnt, fontweight='bold', color='k', transform=ax.transAxes)

    #plt.figtext(0.50,0.85,var1_lab, fontsize=14, fontweight='bold', color='k')

    #plt.figtext(0.18,0.24,'Slope score = ' + str(round(slope_error, 3)), fontsize=10, color='gray') 
    #plt.figtext(0.18,0.21,'Offset score = ' + str(round(offset_error, 3)), fontsize=10, color='gray') 
    
    #if fit_deg == 1:            
    #    if b_obs > 0: oper_obs = '+'
    #    else: oper_obs = ''
    #    if b_sim > 0: oper_sim = '+'
    #    else: oper_sim = ''
    #    plt.figtext(0.18,0.15,'Obs: ' + str(round(m_obs, 6)) + '(' + var2_lab + ')' + oper_obs + str(round(b_obs, 2)), fontsize=8, color='k')
    #    plt.figtext(0.18,0.18,'WACCM: ' + str(round(m_sim, 6)) + '(' + var2_lab + ')' + oper_sim + str(round(b_sim, 2)), fontsize=8, color=col[2])
    ##elif fit_deg == 2:
    ##    plt.figtext(0.18,0.15,'Obs: ' + str(round(a_obs, 6)) + '(' + var2_str + ')$^2$+' + str(round(b_obs, 2)) + '(' + var2_str + ')+' + str(round(c_obs, 2)), fontsize=8, color='k')
    ##    plt.figtext(0.18,0.18,sim + ': ' + str(round(a_sim, 6)) + '(' + var2_str + ')$^2$+' + str(round(b_sim, 2)) + '(' + var2_str + ')+' + str(round(c_sim, 2)), fontsize=8, color=col[2])



if __name__ == '__main__':

    DIST_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/heatmap_vars/'
    PLOTDIR = '/home/wsmith/STRATOCLIM/plots/tracer_fits/'
    WAS_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/WAS_merge/'
    TXT_SAVEDIR = '/home/wsmith/STRATOCLIM/vars/model_txt/'
    MERGEFILES = sorted(glob.glob('/UTLS/Field_Data/STRATOCLIM_Data/merge_data/stratoclim*.nc'))

    # Define settings
    #sims = ['M32L'] #110L
    #cols = [['salmon','mistyrose','red']]#['orange','wheat','brown']]

    var2_strs = ['CFC12','N2O']#,'CFC12','N2O']
    var2_labs = ['CFC-12','N$_2$O','CFC-12','N$_2$O']
    model_boxes = ['Sm','Sm','Lg','Lg']
    var2_us   = [   't',   'b',   't',   'b']
    var2_eq_inds = [9, 10]  # The index where var1 is equal to var2
    var2_val_rngs = [[394,442],[265,292],[394,442],[265,292]]  # For now, just screen for the tropopause at a particular value

    var1_strs = [ 'CH3CCL3','HCFC141B','HCFC22','HCFC142B','CF2CLBR','CCL4', \
                  'CFC11','CF3BR','CFC113',  'CFC12',     'N2O','CFC114',  'CFC115']
    var1_slife = [  '37.7',     '72.3',  '161',   '212',   '33.5','50.6', \
                     '57', '73.5',  '88.4',   '95.5',     '116',   '191',     '997']     # Stratospheric lifetime      # '30.4', '26.3',  
    var1_labs = ['CH$_3$CCl$_3$','HCFC-141b','HCFC-22','HCFC-142b','CF$_2$ClBr','CCl$_4$', \
                  'CFC-11','CF$_3$Br','CFC-113','CFC-12','N$_2$O','CFC-114','CFC-115']                                 # 'CH$_3$Cl','CH$_3$Br',   
    var1_lets = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(j)','(k)','(l)']
    var1_us   = [     't',       't',     't',      't',       't',   't',\
                      't',    't',     't',      't',          'b',     't',     't']
    stars     = [    '* ',   '* ',       '',       '',     '* ',    '* ',\
                     '* ',   '* ',     '* ',     '* ',     ''     , '',     '' ]      # This controls whether an * shall be included to denote photolysis dominant loss                    
                     

    slope_score_table  = np.zeros((len(var1_strs), len(var2_strs)*2))  # Or replace with len(sims) if more than one is desired
    #offset_score_table = np.zeros((len(var2_strs)*2, len(var1_strs)))
    
    var1_caption = []
    for vv in range(len(var1_strs)):
        var1_caption.append(stars[vv] + var1_labs[vv] + ' (' + var1_slife[vv] + ')')

    for v2 in range(len(var2_strs)):

        var2_str = var2_strs[v2]
        var2_lab = var2_labs[v2]
        var2_val_rng = var2_val_rngs[v2] 
        var2_u = var2_us[v2]
        var2_eq_ind = var2_eq_inds[v2]
        model_box = model_boxes[v2]

        fig, axs = plt.subplots(4,3, figsize=(12,16))
        
        for v1 in range(len(var1_strs)):
        
            if v1 == var2_eq_ind: slope_score_table[v1,v2]=0
            else:                    
                if v1 < var2_eq_ind: v1_ind = v1
                if v1 > var2_eq_ind: v1_ind = v1-1          

                xax = v1_ind % 3
                yax = int(v1_ind/3)                

                var1_str = var1_strs[v1]
                var1_lab = var1_labs[v1]
                var1_u = var1_us[v1]
                var1_let = var1_lets[v1]
                
                var1_110L_strat, var2_110L_strat, var1_110L_trop, var2_110L_trop, m_110L, b_110L, model_sub_110L = get_model_vals(DIST_SAVEDIR, '110L', var1_str, var2_str, var2_val_rng, model_box, [50,200])
                var1_M32L_strat, var2_M32L_strat, var1_M32L_trop, var2_M32L_trop, m_M32L, b_M32L, model_sub_M32L = get_model_vals(DIST_SAVEDIR, 'M32L', var1_str, var2_str, var2_val_rng, model_box, [50,200])
                var1_SC_17_trop, var2_SC_17_trop, var1_SC_17_strat, var2_SC_17_strat, var1_SC_17_strat_err, var2_SC_17_strat_err, m_obs, b_obs = get_obs_vals(WAS_SAVEDIR, MERGEFILES, [50,200], var1_str, var2_str, var2_val_rng)
                
                slope_error_110L = (m_110L-m_obs)*100./m_obs      
                slope_error_M32L = (m_M32L-m_obs)*100./m_obs      
                slope_score_table[v1,v2] = slope_error_110L              
                slope_score_table[v1,v2+2] = slope_error_M32L              
                                       
                plot_on_axis(axs[yax, xax], var1_str, var2_str, var1_lab, var2_lab, var1_let, var1_SC_17_trop, var2_SC_17_trop, var1_SC_17_strat, var2_SC_17_strat, var1_SC_17_strat_err, var2_SC_17_strat_err, \
                    var1_110L_trop, var2_110L_trop, var1_110L_strat, var2_110L_strat, var1_M32L_trop, var2_M32L_trop, var1_M32L_strat, var2_M32L_strat, var2_val_rng, var1_u, var2_u, \
                    model_sub_110L, model_sub_M32L, m_obs, b_obs, m_110L, b_110L, m_M32L, b_M32L, ['khaki','wheat','orange'], ['salmon','mistyrose','red'])


        plt.tight_layout()
        plt.savefig(PLOTDIR + 'Figure8_' + var2_str + '_' + model_box + 'box.png', dpi=300)
        plt.savefig(PLOTDIR + 'Figure8_' + var2_str + '_' + model_box + 'box.pdf', dpi=300)
        plt.close('all')

    
    with open('/home/wsmith/STRATOCLIM/vars/Fig_9_data.pkl', 'wb') as f:
        pkl.dump((var1_caption, var1_strs, var2_labs, slope_score_table),f)
        
    print(slope_score_table)

      
    #fnt = 18
    #
    ##################################################
    #  Make a slope score table
    ##################################################          
    #
    #### SCATTER VERSION
    #
    #fig = plt.figure(figsize=(14,6))
    #ax = plt.subplot(111)
    #
    ##### These were used to make a table for the NASA ECIP-ES proposal
    #plt.scatter(np.arange(0.35,len(var1_strs)), slope_score_table[:,0], color='orange', edgecolors='k', s=100, marker='o', label = 'WACCM ' + var2_labs[0] + ', ' + model_boxes[0] + ' box')
    #plt.scatter(np.arange(0.45,len(var1_strs)), slope_score_table[:,1], color='orange', edgecolors='k', s=100, marker='^', label = 'WACCM ' + var2_labs[1] + ', ' + model_boxes[1] + ' box')    
    #plt.scatter(np.arange(0.55,len(var1_strs)), slope_score_table[:,2], color='red', edgecolors='k', s=100, marker='o', label = 'MUSICA ' + var2_labs[0] + ', ' + model_boxes[0] + ' box')
    #plt.scatter(np.arange(0.65,len(var1_strs)), slope_score_table[:,3], color='red', edgecolors='k', s=100, marker='^', label = 'MUSICA ' + var2_labs[1] + ', ' + model_boxes[1] + ' box')
    #
    #plt.plot([0,len(var1_strs)],[0,0],linestyle='--',color='k')
    #plt.plot([0,len(var1_strs)],[20,20],linestyle='--',color=[0.8,0.8,0.8])
    #plt.plot([0,len(var1_strs)],[40,40],linestyle='--',color=[0.8,0.8,0.8])
    #plt.plot([0,len(var1_strs)],[-20,-20],linestyle='--',color=[0.8,0.8,0.8])
    #plt.plot([0,len(var1_strs)],[-40,-40],linestyle='--',color=[0.8,0.8,0.8])    
    #
    #ax.set_title('Lower Stratospheric Chemistry Errors', fontsize = fnt)
    #
    #ax.set_xticks(np.arange(0.5,len(var1_strs)))
    #ax.set_xticklabels(var1_caption, fontsize = fnt)
    #plt.xticks(rotation=90)
    #
    #plt.ylabel('Slope Error (%)', fontsize=fnt)
    #plt.yticks(fontsize=fnt)
    #plt.ylim(-60,60)
    #plt.xlim(0,len(var1_strs))
    #
    #plt.legend(fontsize=10)
    #plt.tight_layout()
    #plt.savefig(PLOTDIR + 'Figure9.png')  
    #plt.savefig(PLOTDIR + 'Figure9.pdf')  
    #plt.close('all') 

    
    



