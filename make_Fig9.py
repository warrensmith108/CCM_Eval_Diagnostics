
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

PLOTDIR = '/home/wsmith/STRATOCLIM/plots/tracer_fits/'

with open('/home/wsmith/STRATOCLIM/vars/Fig_9_data.pkl', 'rb') as f:
    var1_caption, var1_strs, var2_labs, slope_score_table = pkl.load(f)

####  HISTOGRAM VERSION
fnt = 14

fig = plt.figure(figsize=(14,6))
ax = plt.subplot(111)

#score_inds = np.zeros(np.shape(slope_score_table))
#score_inds[:,0] = np.arange(0,len(var1_strs)+0.1)
#score_inds[:,1] = np.arange(0,len(var1_strs)+0.1)
#score_inds[:,2] = np.arange(0,len(var1_strs)+0.1)
#score_inds[:,3] = np.arange(0,len(var1_strs)+0.1)

ax.hist(np.arange(0.3,len(var1_strs),1), len(var1_strs), weights=slope_score_table[:,0], color='orange', histtype='bar', label='WACCM ' + var2_labs[0], rwidth=0.1)
ax.hist(np.arange(0.4,len(var1_strs),1), len(var1_strs), weights=slope_score_table[:,1], color='orange', hatch='///', histtype='bar', label='WACCM ' + var2_labs[1], rwidth=0.1)
ax.hist(np.arange(0.5,len(var1_strs),1), len(var1_strs), weights=slope_score_table[:,2], color='red', histtype='bar', label='MUSICA ' + var2_labs[0], rwidth=0.1)
ax.hist(np.arange(0.6,len(var1_strs),1), len(var1_strs), weights=slope_score_table[:,3], color='red', hatch='///', histtype='bar', label='MUSICA ' + var2_labs[1], rwidth=0.1)

#ax.hist(score_inds, len(var1_strs), weights=slope_score_table, color=['orange','khaki','red','salmon'], histtype='bar')

plt.plot([0,len(var1_strs)],[0,0],linestyle='--',color='k')
plt.plot([0,len(var1_strs)],[20,20],linestyle='--',color=[0.8,0.8,0.8])
plt.plot([0,len(var1_strs)],[40,40],linestyle='--',color=[0.8,0.8,0.8])
plt.plot([0,len(var1_strs)],[-20,-20],linestyle='--',color=[0.8,0.8,0.8])
plt.plot([0,len(var1_strs)],[-40,-40],linestyle='--',color=[0.8,0.8,0.8])    

ax.set_title('Lower Stratospheric Chemistry Errors', fontsize = fnt)

ax.set_xticks(np.arange(0.85,13,0.928))
ax.set_xticklabels(var1_caption, fontsize = fnt)
plt.xticks(rotation=90)

plt.ylabel('Slope Error (%)', fontsize=fnt)
plt.yticks(fontsize=fnt)
plt.ylim(-60,60)
plt.xlim(0,len(var1_strs))

plt.legend(fontsize=fnt)
plt.tight_layout()
plt.savefig(PLOTDIR + 'Figure9.png')  
plt.savefig(PLOTDIR + 'Figure9.pdf')  
plt.close('all') 
    