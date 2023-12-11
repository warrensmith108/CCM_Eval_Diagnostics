
#############################################################################
#
#  This script is intended to read the merge datasets provided to us by 
#  Marc von Hobe from the StratoClim campaign
#
#############################################################################

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import glob
import datetime as dt
import pickle as pkl
from STRATOCLIM_plotter import plot_TT_pix, plot_layerhist, plot_boxwhisker
from read_STRATOCLIM_WAS import load_WAS
from matplotlib.ticker import ScalarFormatter
import os
import jdcal

def get_merge_array(var_str, DATAFILES):

    var_data = []
    AVIONIK_LON = []
    AVIONIK_LAT = []
    AVIONIK_ALT = []
    AVIONIK_PRS = []
    CLAMS_THETA = []
    CPT_GPH = []
    CPT_THETA = []
    LRT_THETA = []
    LRT_PRS = []
    dateint = []
    jd = []
    time_ref = dt.datetime(2000, 1, 1, 0)

    for f in range(len(DATAFILES)):
        
        filef = nc4.Dataset(DATAFILES[f])
        NF = len(filef.variables['AMICA:CO'][:])
                     
        AVIONIK_ALT.extend(filef.variables['AVIONIK:ALT'][:])
        AVIONIK_LAT.extend(filef.variables['AVIONIK:LAT'][:])
        AVIONIK_LON.extend(filef.variables['AVIONIK:LON'][:])
        AVIONIK_PRS.extend(filef.variables['AVIONIK:PRESS'][:])
        CLAMS_THETA.extend(filef.variables['CLAMS_MET:THETA'][:])
        #CPT_GPH.extend(filef.variables['CLAMS_MET:GPH_CPT'][:]*0.1)  # NOTE!!!!!!  It is clear that this 0.1 scaling is required to give realistic tropopause values  
        CPT_THETA.extend(filef.variables['CLAMS_MET:THETA_CPT'][:])
        LRT_THETA.extend(filef.variables['CLAMS_MET:THETA_TROP1'][:])  
        LRT_PRS.extend(filef.variables['CLAMS_MET:PRESS_TROP1'][:])  
        curr_time = filef.variables['time'][:]
        if var_str in filef.variables.keys(): 
            var_data.extend(filef.variables[var_str][:])
            var_unit = getattr(filef.variables[var_str], 'units')                        
        else:
            nanarray = np.empty(NF)
            nanarray[:] = np.NaN
            var_data.extend(nanarray)  
        for t in range(NF):
            curr_dt = time_ref + dt.timedelta(seconds=curr_time[t])
            jd.append(sum(jdcal.gcal2jd(curr_dt.year, curr_dt.month, curr_dt.day)) + curr_dt.hour/24. + curr_dt.minute/1440. + curr_dt.second/86400.)
            dateint.extend([int(str(curr_dt.year) + str(curr_dt.month).zfill(2) + str(curr_dt.day).zfill(2) + \
                           str(curr_dt.hour).zfill(2) + str(curr_dt.minute).zfill(2) + str(curr_dt.second).zfill(2))])

        # As of 7-31-23, we no longer want to use the CPT from the merge files.  We instead use the ERA5 merge (produced by ERA5_trop_merge.py)
        with open('/home/wsmith/STRATOCLIM/vars/clp_z_merge.pkl','rb') as f:
            ERA5_CPT = pkl.load(f)
            
        # This is kind of silly, but to keep everything on the same playing field we need to open the ERA5 LRT here (produced by ERA5_trop_merge.py)
        with open('/home/wsmith/STRATOCLIM/vars/wmo_1st_z_merge.pkl','rb') as f:
            ERA5_LRT = pkl.load(f)
            
    return var_data, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, ERA5_CPT, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT

def get_model_prof_data(CESMSAVEDIR, var, ztype, sim, lon_bnds, lat_bnds, prs_lvl, quants):

    simfile_z = CESMSAVEDIR + var + ztype + 'profiles_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_' + str(int(prs_lvl)) + 'GPHfilter.pkl'       
    with open(simfile_z, 'rb') as f:
        var_zprofs_in, var_zprofs_out, TROP_array_z, lev_z, gph_asma_z = pkl.load(f) 
        
    var_profs_z = np.append(var_zprofs_in, var_zprofs_out, axis=0)
            
    var_mean = []
    var_lowq = []
    var_upq  = []
    
    for z in range(len(lev_z)):
        curr_var_z = np.array(var_profs_z[:,z])
        curr_var_z = curr_var_z[curr_var_z > 0]  # Correct for the fact that some points may have been flagged as under the surface
        if len(curr_var_z) > 1:
            var_mean.append(np.mean(curr_var_z))
            var_lowq.append(np.percentile(curr_var_z, quants[0]))
            var_upq.append(np.percentile(curr_var_z, quants[1]))
        else: 
            var_mean.append(-999)
            var_lowq.append(-999)
            var_upq.append(-999)
        
    return lev_z, var_mean, var_lowq, var_upq
    
def get_model_pix_data(CESMSAVEDIR, var_str, ztype, sim, ybins, yrng, lon_bnds, lat_bnds, quants, yscale, col):

    print(sim, var_str)

    ybin_bnds = np.arange(yrng[0], yrng[1]+0.1, (yrng[1]-yrng[0])/ybins)
    
    with open(CESMSAVEDIR + var_str + '_1d_nointerpvals_' + sim + '_' + str(lon_bnds) + '_' + str(lat_bnds) + '_100GPHfilter.pkl', 'rb') as f:
        var_1d, z_1d, troprelz_1d, theta_1d, tropreltheta_1d = pkl.load(f)
        
    if ztype == 'z': vert_coord_1d = z_1d
    if ztype == 'troprelz': vert_coord_1d = troprelz_1d
    if ztype == 'theta': vert_coord_1d = theta_1d
    if ztype == 'tropreltheta': vert_coord_1d = tropreltheta_1d
                       
    var_mean = []
    var_lowq = []
    var_upq  = []
    
    for bin in range(ybins):
        curr_var_in_bin = []
        #curr_vert_in_bin = []
        curr_binmin = ybin_bnds[bin]
        curr_binmax = ybin_bnds[bin+1]

        print('****', curr_binmin, 'to', curr_binmax)        
        curr_file = '/home/wsmith/STRATOCLIM/vars/z_bins/' + sim + '_' + var_str + '_' + str(lon_bnds) + '-' + str(lat_bnds) + '_' + str(curr_binmin) + '-' + str(curr_binmax) + '.pkl'
        
        if os.path.exists(curr_file):
            print('...recovered from a file')
            with open(curr_file, 'rb') as f:
                curr_var_in_bin = pkl.load(f)
        
        else: 
            for pt in range(len(var_1d)):
                curr_vert = vert_coord_1d[pt]
                if curr_vert >= curr_binmin and curr_vert <= curr_binmax:
                    curr_var_in_bin.append(var_1d[pt])
            
            with open(curr_file, 'wb') as f:
                pkl.dump((curr_var_in_bin), f)
        

        curr_var_in_bin = np.asarray(curr_var_in_bin)

        if len(curr_var_in_bin) > 0: 
            var_mean.append(np.mean(curr_var_in_bin))
            var_lowq.append(np.percentile(curr_var_in_bin, quants[0]))
            var_upq.append(np.percentile(curr_var_in_bin, quants[1]))
        else: 
            var_mean.append(np.nan)
            var_lowq.append(np.nan)
            var_upq.append(np.nan)
            
        
        plt.plot([var_mean[-1],var_mean[-1]], [curr_binmin+yscale,curr_binmax+yscale], color=col)
        plt.plot([var_lowq[-1],var_lowq[-1]], [curr_binmin+yscale,curr_binmax+yscale], color=col, linewidth=2)
        plt.plot([var_upq[-1], var_upq[-1]],  [curr_binmin+yscale,curr_binmax+yscale], color=col, linewidth=2)



def plot_model_overlay(ax, var_data, var_bin, var_rng, var_str, var_lab, lon_bnds, lat_bnds, ALT_DATA, CPT_DATA, ztype, ybins, yrng, quants, n_xticks, ylab, yticks, var_unit, do_yscale, do_xlog, CESMSAVEDIR, PLOTDIR, do_modelpix = 0):

    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.set_position([0.1,0.1,0.8,0.8])
    
    if do_yscale == 1: yscale = np.nanmean(CPT_DATA)
    else: yscale = 0.
    plot_TT_pix(ax, np.asarray(var_data), var_bin, var_rng, ALT_DATA, ybins, yrng, str(var_lab) + ' (' + var_unit + ')', ylab, \
                PLOTDIR + str(var_str).replace(':', '_') + '_' + ztype, lay_avg = 1, do_yscale = do_yscale, \
                trop_data=CPT_DATA, n_xticks=n_xticks, yticks=yticks, do_save=0, quants=quants)

    if do_modelpix:
    
        # Call function
        get_model_pix_data(CESMSAVEDIR, var_str, ztype, '110L', ybins, yrng, lon_bnds, lat_bnds, quants, yscale, 'orange')
        get_model_pix_data(CESMSAVEDIR, var_str, ztype, 'M32L', ybins, yrng, lon_bnds, lat_bnds, quants, yscale, 'red')
        #get_model_pix_data(CESMSAVEDIR, var_str, ztype, '110L_longspin', ybins, yrng, lon_bnds, lat_bnds, quants, yscale, 'blue')
    
    else:

        lev_z, var_mean_M32L, var_lowq_M32L, var_upq_M32L = get_model_prof_data(CESMSAVEDIR, var_str, ztype, 'M32L', lon_bnds, lat_bnds, 100, quants) 
        ax.fill_betweenx(lev_z+yscale, var_lowq_M32L, var_upq_M32L, color='red', alpha=0.3)          
        ax.plot(var_mean_M32L, lev_z+yscale, color='red', linewidth = 2)
        #plt.plot(var_lowq_M32L, lev_z+yscale, color='red', linestyle='--')
        #plt.plot(var_upq_M32L,  lev_z+yscale, color='red', linestyle='--')
        
        lev_z, var_mean_110L, var_lowq_110L, var_upq_110L = get_model_prof_data(CESMSAVEDIR, var_str, ztype, '110L', lon_bnds, lat_bnds, 100, quants)        
        ax.fill_betweenx(lev_z+yscale, var_lowq_110L, var_upq_110L, color='orange', alpha=0.3)  
        ax.plot(var_mean_110L, lev_z+yscale, color='orange', linewidth = 2)
        #plt.plot(var_lowq_110L, lev_z+yscale, color='orange', linestyle='--')
        #plt.plot(var_upq_110L,  lev_z+yscale, color='orange', linestyle='--')

        # This curve goes in the supplement
        #lev_z, var_mean_lnox, var_lowq_lnox, var_upq_lnox = get_model_prof_data(CESMSAVEDIR, var_str, ztype, 'Simone_32L_noLNOx', lon_bnds, lat_bnds, 100, quants)        
        #plt.plot(var_mean_lnox, lev_z+yscale, color='cyan')
        #plt.plot(var_lowq_lnox, lev_z+yscale, color='cyan', linestyle='--')
        #plt.plot(var_upq_lnox,  lev_z+yscale, color='cyan', linestyle='--')
    
    if do_xlog == 1: ax.set_xscale('log')
    
    
    #plt.savefig(PLOTDIR + str(var_str).replace(':', '_') + '_' + ztype + '_modeloverlay.png')
    #plt.close('all')

def get_plot_data(DATAFILES, var_str):

    if isinstance(var_str, str): var_data, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_merge_array(var_str, DATAFILES)
    if isinstance(var_str, list): 
        var_data, var_unit, dateint, jd, AVIONIK_LON, AVIONIK_LAT, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_merge_array(var_str[0], DATAFILES)
        for sub_v in range(len(var_str)-1):
            var_data_sub, var_unit_sub, dateint_sub, jd_sub, AVIONIK_LON_sub, AVIONIK_LAT_sub, AVIONIK_ALT_sub, AVIONIK_PRS_sub, CPT_GPH_sub, CLAMS_THETA_sub, CPT_THETA_sub, LRT_THETA_sub, LRT_PRS_sub, ERA5_LRT_sub = get_merge_array(var_str[sub_v+1], DATAFILES)
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
            jd.extend(jd_sub)
                
    # Convert all the lists to arrays for modification and plotting
    var_data = np.asarray(var_data)
    AVIONIK_ALT = np.asarray(AVIONIK_ALT)*1.e-3
    AVIONIK_PRS = np.asarray(AVIONIK_PRS)
    CPT_GPH = np.asarray(CPT_GPH)
    CLAMS_THETA = np.asarray(CLAMS_THETA)
    CPT_THETA = np.asarray(CPT_THETA)    
    LRT_THETA = np.asarray(LRT_THETA)
    LRT_PRS = np.asarray(LRT_PRS)
    ERA5_LRT = np.asarray(ERA5_LRT)

    
    return var_data, jd, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT
        

#########################################################
#         MAIN PROGRAM
#########################################################

if __name__ == '__main__':

    # Define the relevant paths
    DATAROOT = '/UTLS/Field_Data/STRATOCLIM_Data/merge_data/'
    DATANAME = 'stratoclim*.nc'
    DATAFILES = sorted(glob.glob(DATAROOT + DATANAME))
    PLOTDIR = '/home/wsmith/STRATOCLIM/plots/merge_data_survey/'
    SAVEDIR = '/home/wsmith/STRATOCLIM/vars/WAS_merge/'
    CESMSAVEDIR = '/home/wsmith/STRATOCLIM/vars/'
    
    sub = 10

    # Define the variable information
    #var_strs = ['AMICA:CO','AMICA:H2O','AMICA:OCS','CCP-CIPGS:CONC_tot','CLAMS_AGE:E90','CLAMS_AGE:MODE','NIXECAPS:IWC','FOZAN:O3','COLD:CO','COLD:N2O','ERICA:NITRATE','ERICA:SULPHATE','FLASH:H2O',  \
    #            'HAGAR:N2O','HAGAR:CFC12','HAGAR:CFC11','HAGAR:H1211','HAGAR:SF6','HAGAR:CO2','UHSAS:N_tot','SIOUX:NO','SIOUX:NOY']
    #var_names = ['CO','H2O','OCS','CONC','
    #var_strs = ['AMICA:CO']
    #var_rngs = [[0,200],     [1,100000],  [0,1050],     [0,16],              [0,1.1],        [0,11],         [0,3000],  [0,1000],  [0,200],  [230,360],   [0,2.1],          [0,4.1],     [1,6000],    \
    #            [230,360],     [370,520],  [100,300],    [1,5.1],        [5,14],      [395,410],  [0,7000],    [0,20],    [0,25]]
    #var_bins = [    40,       100,        105,             16,                  11,             11,             30,         20,        40,      26,          21,               41,          60,  \
    #              26,             15,          20,         41,             9,            15,        70,         20,         25]  
    
    #var_strs = [,'HAGAR:N2O','FOZAN:O3','HAGAR:CFC12','HAGAR:CFC11']
    var_strs =  [['COLD:CO','AMICA:CO'],'COLD:CO','AMICA:CO','FOZAN:O3']#,'HAGAR:N2O','HAGAR:CFC12','HAGAR:CFC11','FLASH:H2O','AVIONIK:ALT','CLAMS_MET:GPH_CPT','AVIONIK:PRESS','AVIONIK:THETA','CLAMS_MET:PRESS_TROP1']
    var_names = ['CO','CO','CO','O3','N2O','CFC12','CFC11','H2O','Altitude','CPT_Altitude','p','Theta','LRT_Prs']
    var_labs  = ['CO','CO','CO','O$_3$','N2O','CFC12','CFC11','H2O','Altitude','CPT_Altitude','p','Theta','LRT_Prs']
    var_bins = [np.arange(0,201,5),np.arange(0,201,5),np.arange(0,201,5),np.logspace(1,4,num=50),np.arange(0,201,5),np.arange(0,1201,40),np.arange(260,361,5),np.arange(370,531,5),np.arange(100,301,10),np.logspace(-1,4,num=100),np.arange(0,25000),np.arange(10000,20001,1000),np.arange(10000,20001,1000)]   # ozone: ,np.arange(0,2001,20)
    var_rngs = [[0,200],[0,200],[0,200],[10,2000],[0,200],[260,360],[370,530],[100,300],[0.1,10000],[0,25000],[10000,20000],[10000,20000],[10000,20000]]
    var_logs = [  0,     0,  0,  1,   0,    0,    0,     0, 1, 0, 0,0,0,0]
    var_units = ['b','b','b','b','b','t','t','m','','','','','']

    # Make Figure 3
    fig = plt.figure(figsize=(12,10))
    ax1  = plt.subplot(2,2,1)
    ax2  = plt.subplot(2,2,2)
    ax3  = plt.subplot(2,2,3)
    ax4  = plt.subplot(2,2,4)
    
    CO_data, jd, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_plot_data(DATAFILES, ['COLD:CO','AMICA:CO'])
        
    CPT_REL_ALT   = AVIONIK_ALT - CPT_GPH  # CPT-rel geopotential height (in km)
    CPT_REL_THETA = CLAMS_THETA - CPT_THETA  # CPT-rel theta
    LRT_REL_THETA = CLAMS_THETA - LRT_THETA # LRT-rel theta
    LRT_REL_ALT   = AVIONIK_ALT - ERA5_LRT    
    
    CO_data[np.isnan(CO_data)] = -999
    CPT_REL_ALT[np.isnan(CPT_REL_ALT)] = -999
    CPT_REL_THETA[np.isnan(CPT_REL_THETA)] = -999
    
    plot_model_overlay(ax1, CO_data[::sub], np.arange(0,201,5), [0,200], 'CO', 'CO', [75,95], [18,32], LRT_REL_ALT[::sub], ERA5_LRT[::sub], 'troprelz', 17, [-13,4], [5,95], 4, 'Adj. Tropopause-rel. Alt (km)', [0,5,10,15,20], \
                       'ppbv', 1, 0, CESMSAVEDIR, PLOTDIR, do_modelpix=0)  #####  TROPOPAUSE RELATIVE ALTITUDE        
    plot_model_overlay(ax3, CO_data[::sub], np.arange(0,201,5), [0,200], 'CO', 'CO', [75,95], [18,32], CLAMS_THETA[::sub], LRT_THETA[::sub], 'theta', 32, [320,480], [5,95], 4, 'Potential Temperature (K)', [], \
                       'ppbv', 0, 0, CESMSAVEDIR, PLOTDIR, do_modelpix=0)                         #####  THETA    
    
    O3_data, jd, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_plot_data(DATAFILES, 'FOZAN:O3')

    CPT_REL_ALT   = AVIONIK_ALT - CPT_GPH  # CPT-rel geopotential height (in km)
    CPT_REL_THETA = CLAMS_THETA - CPT_THETA  # CPT-rel theta
    LRT_REL_THETA = CLAMS_THETA - LRT_THETA # LRT-rel theta
    LRT_REL_ALT   = AVIONIK_ALT - ERA5_LRT    
    
    O3_data[np.isnan(O3_data)] = -999
    CPT_REL_ALT[np.isnan(CPT_REL_ALT)] = -999
    CPT_REL_THETA[np.isnan(CPT_REL_THETA)] = -999
    
    plot_model_overlay(ax2, O3_data[::sub], np.logspace(1,4,num=50), [10,2000], 'O3', 'O$_3$', [75,95], [18,32], LRT_REL_ALT[::sub], ERA5_LRT[::sub], 'troprelz', 17, [-13,4], [5,95], 4, 'Adj. Tropopause-rel. Alt (km)', [0,5,10,15,20], \
                       'ppbv', 1, 1, CESMSAVEDIR, PLOTDIR, do_modelpix=0)  #####  TROPOPAUSE RELATIVE ALTITUDE        
    plot_model_overlay(ax4, O3_data[::sub], np.logspace(1,4,num=50), [10,2000], 'O3', 'O$_3$', [75,95], [18,32], CLAMS_THETA[::sub], LRT_THETA[::sub], 'theta', 32, [320,480], [5,95], 4, 'Potential Temperature (K)', [], \
                       'ppbv', 0, 1, CESMSAVEDIR, PLOTDIR, do_modelpix=0)                         #####  THETA  

    fnt = 16
    plt.figtext(0.26,0.94,'AMICA&COLD2 CO', fontweight='bold', fontsize=fnt)
    plt.figtext(0.82,0.75,'FOZAN-II O$_3$', fontweight='bold', fontsize=fnt)
    plt.figtext(0.84,0.72,'WACCM', fontweight='bold', fontsize=fnt, color='orange')
    plt.figtext(0.84,0.69,'MUSICA',fontweight='bold', fontsize=fnt, color='red')
    
    plt.figtext(0.44,0.60,'(a)', fontweight='bold', fontsize=fnt)
    plt.figtext(0.94,0.60,'(b)', fontweight='bold', fontsize=fnt)
    plt.figtext(0.44,0.10,'(c)', fontweight='bold', fontsize=fnt)
    plt.figtext(0.94,0.10,'(d)', fontweight='bold', fontsize=fnt)
                     
    plt.tight_layout()
    plt.savefig(PLOTDIR + 'Figure3.png', dpi=200)
    plt.savefig(PLOTDIR + 'Figure3.pdf', dpi=200)
    plt.close('all')

          
    ###########################################  MERGE THE MERGE DATA ONTO THE WAS TIME GRID
    
    all_data, var_headers = load_WAS(r'/home/wsmith/STRATOCLIM/STRATOCLIM_WAS_Times.xlsx')
    
    for v in range(len(var_strs)):
    
        var_str = var_strs[v]
        var_name = var_names[v]
    
        var_data, jd, AVIONIK_ALT, AVIONIK_PRS, CPT_GPH, CLAMS_THETA, CPT_THETA, LRT_THETA, LRT_PRS, ERA5_LRT = get_plot_data(DATAFILES, var_str)
        
        merge_WAS_arr = []  # Create a blank array to be filled with the merged array
        merge_prs_arr = []
        nWAS = (np.shape(all_data)[1]-1)/2  # The number of canisters
        #print(nWAS)
        
        for msr in range(int(nWAS)):
            
            curr_ind = msr*2+1    # 1, 3, 5, etc
            
            curr_open_date    = all_data[1][curr_ind]
            curr_open_time    = all_data[2][curr_ind]            
            curr_open_dt      = curr_open_date + dt.timedelta(hours = curr_open_time.hour, minutes = curr_open_time.minute, seconds = curr_open_time.second)
            curr_open_jd      = sum(jdcal.gcal2jd(curr_open_dt.year, curr_open_dt.month, curr_open_dt.day)) + curr_open_dt.hour/24. + curr_open_dt.minute/1440. + curr_open_dt.second/86400.
            #curr_open_dateint = int(str(curr_open_dt.year) + str(curr_open_dt.month).zfill(2) + str(curr_open_dt.day).zfill(2) + \
            #                   str(curr_open_dt.hour).zfill(2) + str(curr_open_dt.minute).zfill(2) + str(curr_open_dt.second).zfill(2))

            curr_close_date    = all_data[1][curr_ind+1]
            curr_close_time    = all_data[2][curr_ind+1]
            curr_close_dt      = curr_close_date + dt.timedelta(hours = curr_close_time.hour, minutes = curr_close_time.minute, seconds = curr_close_time.second)
            curr_close_jd      = sum(jdcal.gcal2jd(curr_close_dt.year, curr_close_dt.month, curr_close_dt.day)) + curr_close_dt.hour/24. + curr_close_dt.minute/1440. + curr_close_dt.second/86400.            
            #curr_close_dateint = int(str(curr_close_dt.year) + str(curr_close_dt.month).zfill(2) + str(curr_close_dt.day).zfill(2) + \
            #                     str(curr_close_dt.hour).zfill(2) + str(curr_close_dt.minute).zfill(2) + str(curr_close_dt.second).zfill(2))                    
            
            curr_var_arr = []  # This will be filled with all the merge data falling within each open/close, to be averaged at the end        
            curr_prs_arr = []
            for merge_t in range(len(jd)):
                curr_merge_jd = jd[merge_t]
                if curr_merge_jd >= curr_open_jd and curr_merge_jd <= curr_close_jd:
                    curr_prs_arr.append(AVIONIK_PRS[merge_t])
                    curr_var_arr.append(var_data[merge_t])
            
            #print(np.nanmean(curr_var_arr), np.nanstd(curr_var_arr))
            
            merge_prs_arr.append(np.nanmean(curr_prs_arr))
            merge_WAS_arr.append(np.nanmean(curr_var_arr))
            
        if var_name == 'CO':  merge_WAS_error = np.asarray(merge_WAS_arr)*0.03  # This is based on the references for the CO instruments
        if var_name == 'O3':  merge_WAS_error = np.zeros(len(merge_WAS_arr))
        if var_name == 'N2O': merge_WAS_error = np.asarray(merge_WAS_arr)*0.01  # This is based on the reference for HAGAR
        

        savefile = SAVEDIR + str(var_str) + '_merge.pkl'    
        #print(len(merge_WAS_arr))
        with open(savefile, 'wb') as f:
            pkl.dump(([], [], merge_WAS_arr, merge_WAS_error), f)
         
         
    ####### GENERATE SOME PLOTS
    

    
    
    
    
    
    # Make a trop-rel altitude plot
    #plot_TT_pix(np.asarray(var_data), var_bin, var_rng, CPT_REL_ALT, 40, [-15,5], str(var_str) + ' (' + var_unit + ')', 'Adj. Tropopause-rel. Alt (km)', \
    #            PLOTDIR + str(var_str).replace(':', '_') + '_CPT_GPH', lay_avg = 1, do_yscale = 1, trop_data=CPT_GPH, n_xticks=n_xticks, yticks=[0,4,8,12,16,20])        

    # FOR SAVING THE VARIABLES
    #plot_model_overlay(var_data, var_bin, var_rng, var_name, var_lab, [75,95], [18,32], LRT_REL_ALT, ERA5_LRT, 'troprelz', 2, [-2,2], [5,95], 4, 'Adj. Tropopause-rel. Alt (km)', [0,5,10,15,20], \
    #                   'pp' + unit + 'v', 1, var_log, CESMSAVEDIR, PLOTDIR, do_modelpix=1)  #####  TROPOPAUSE RELATIVE ALTITUDE
    #plot_model_overlay(var_data, var_bin, var_rng, var_name, var_lab, [75,95], [18,32], CLAMS_THETA, LRT_THETA, 'theta', 2, [350,390], [5,95], 4, 'Potential Temperature (K)', [], \
    #                   'pp' + unit + 'v', 0, var_log, CESMSAVEDIR, PLOTDIR, do_modelpix=1)                         #####  THETA
                       
            
    # Make a trop-rel altitude plot with the model data overlain on the observations 
    #plot_model_overlay(var_data, var_bin, var_rng, var_name, var_lab, [75,95], [18,32], LRT_REL_ALT, ERA5_LRT, 'troprelz', 17, [-13,4], [5,95], 4, 'Adj. Tropopause-rel. Alt (km)', [0,5,10,15,20], \
    #                   'pp' + unit + 'v', 1, var_log, CESMSAVEDIR, PLOTDIR, do_modelpix=0)  #####  TROPOPAUSE RELATIVE ALTITUDE
    #plot_model_overlay(var_data, var_bin, var_rng, var_name, [75,95], [18,32], AVIONIK_ALT, CPT_GPH, 'z', 34, [4,21], [5,95], 4, 'Altitude (km)', [0,5,10,15,20], \
    #                   'pp' + unit + 'v', 0, var_log, CESMSAVEDIR, PLOTDIR)      #####  ALTITUDE
    #plot_model_overlay(var_data, var_bin, var_rng, var_name, [75,95], [18,32], CPT_REL_THETA, CPT_THETA, 'tropreltheta', 34, [-70,100], [5,95], 4, 'Adj. Tropopause-rel. Theta (K)', [], \
    #                   'pp' + unit + 'v', 1, var_log, CESMSAVEDIR, PLOTDIR)  #####  TROPOPAUSE RELATIVE THETA
    #plot_model_overlay(var_data, var_bin, var_rng, var_name, var_lab, [75,95], [18,32], CLAMS_THETA, LRT_THETA, 'theta', 32, [320,480], [5,95], 4, 'Potential Temperature (K)', [], \
    #                   'pp' + unit + 'v', 0, var_log, CESMSAVEDIR, PLOTDIR, do_modelpix=0)                         #####  THETA
    


    
    #plot_boxwhisker(np.asarray(var_data), var_bin, var_rng, AVIONIK_ALT, 22, [0,22], str(var_str) + ' (' + var_unit + ')', 'Altitude (km)', \
    #            PLOTDIR + str(var_str).replace(':', '_') + '_BoxWhisker')
                
    #plot_TT_pix(np.asarray(var_data), var_bin, var_rng, CPT_REL_THETA, 24, [-60,60], str(var_str) + ' (' + var_unit + ')', 'Adj. Tropopause-rel. Alt (km)', \
    #            PLOTDIR + str(var_str).replace(':', '_') + '_CPT_THETA', lay_avg = 1, do_yscale=1, trop_data=CPT_THETA, yticks=[340,360,380,400,420])

    #plot_layerhist(np.asarray(var_data), var_bin, var_rng, CPT_REL_ALT, 8, [-12,4], str(var_str) + ' (' + var_unit + ')', 'Relative Frequency', \
    #            PLOTDIR + str(var_str).replace(':', '_') + '_CPT_GPH_layerhist')
    
