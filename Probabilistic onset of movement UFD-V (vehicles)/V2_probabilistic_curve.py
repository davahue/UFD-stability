# -*- coding: utf-8 -*-
"""
19.01.2022
Onset of stability for urban debris.
@author: D. Valero
"""

import time
#import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# In-house library to support auxiliary calculations
import faux as fa

start = time.time()

# For reproducibility
np.random.seed(120589)
Nr = 1000 # Number of Monte Carlo draws

folder = "V2-input"

Lx_list = fa.get_var(folder, Nr, "Lx")

# The Lx value defines the rest of parameters. Step-by-step:
# 1. Obtain Ly:
A, B, RMSE = 2.030, 0.476, 0.088 # From curve fitting (Simple Geom Rel)
Ly_list = fa.OneMinusExp(Lx_list, A, B) + np.random.normal(0,RMSE,Nr) # + np.random.uniform(-dev,+dev,Nr) 
Ly_list[Ly_list<1.6] = 1.6 # Cars are hardly ever shorter.

# 2. Obtain Lz:
A, B, RMSE = 0.086, 1.156, 0.117
Lz_list = fa.Linear(Lx_list, A, B) + np.random.normal(0,RMSE,Nr)
Lz_list[Lz_list<1.45] = 1.45 # Cars are hardly ever shorter.

# 3. Obtain rhob:
A, B, RMSE = 56.856, 20.198, 14.858
rho_bulk = fa.Linear(Ly_list, A, B) + np.random.normal(0,RMSE,Nr)
M_list = rho_bulk*Lx_list*Ly_list*Lz_list 

# 4. Obtain clearance zc:
zc_list =  np.random.uniform(0.125,0.32,Nr) # fa.Linear(Lz_list, A, B) + np.random.normal(0,RMSE,Nr)


# 5. Other independent, random quantities:
Cd_list = fa.get_var(folder, Nr, "Cd")
mu_list = fa.get_var(folder, Nr, "mu")
fA_RMSE_list = fa.get_var(folder, Nr, "fA_RMSE")
fVmax_list = fa.get_var(folder, Nr, "fV")


# # h, vx as key parameters ----------------------------------------------------
nh = 200; nv = 400

# ----------------------------------------------------------------------------
# Start calculations ---------------------------------------------------------

hcrit_flot = np.zeros(Nr)
vcrit_sliding = np.zeros([Nr,nv])
hcrit_sliding = np.zeros([Nr,nv])
vcrit_toppling = np.zeros([Nr,nv])
hcrit_toppling = np.zeros([Nr,nv])

for k in range(0, Nr):
    [hcrit_flot[k], vcrit_sliding[k,:], hcrit_sliding[k,:], vcrit_toppling[k,:], hcrit_toppling[k,:]] = \
        fa.get_stability_curves_V2(M_list[k],Lx_list[k],Ly_list[k],Lz_list[k],
                                  zc_list[k],Cd_list[k],mu_list[k],
                                  fA_RMSE_list[k], fVmax_list[k], nh, nv)
    
    if (k % 100) == 0:
        print("Stability curve number: ", k)

# We do not know if sliding or toppling will dominate. Thus, check the most unstable one...
hcrit = np.zeros_like(hcrit_sliding)
for i in range(0, nv):
    for k in range(0, Nr):
        hcrit[k, i] = np.nanmin([hcrit_sliding[k,i], hcrit_toppling[k,i]])

# PLOTTING --------------------------------------------------------------------        

# Fig. opts.
plt.figure(figsize=(4,4))

plt.plot(vcrit_sliding[k,0], hcrit[k,0], ls='-', c='k', zorder=7, alpha=0.05,
         lw=2, label='Stability curve')
# plt.plot(vcrit_toppling[k,0], hcrit_toppling[k,0], ls='--', c='b', zorder=7, alpha=0.05,
#           label='Toppling')

for k in range(1,Nr):
    # Flotation ----
    # plt.scatter(0, hcrit_flot[k], marker= 'o', s=80, zorder=7, alpha=0.2,
    #             edgecolor='k', facecolor='w', label='Flotation')
    # Sliding ----
    plt.plot(vcrit_sliding[k,:], hcrit[k,:], ls='-', c='k', zorder=7, alpha=0.05,
             lw=2)
    # Toppling ----
    # plt.plot(vcrit_toppling[k,:], hcrit_toppling[k,:], ls='--', c='b', zorder=7, alpha=0.05)

# imed = fa.find_median_curve(vcrit_sliding, hcrit)
# plt.plot(vcrit_sliding[imed,:], hcrit[imed,:], ls='--', c='r', zorder=7, 
#          lw=2, label="Median curve")

#hperc500 = np.percentile(hcrit, 50.0, axis=0)
hperc500 = fa.nan_percentiles_bias_corr(vcrit_sliding, hcrit, 50)
plt.plot(vcrit_sliding[0,:], hperc500, ls='-.', c='r', zorder=7, 
         lw=2, label="50 % percentile")
# hperc025 = np.percentile(hcrit, 2.5, axis=0)
hperc025 = fa.nan_percentiles_bias_corr(vcrit_sliding, hcrit, 2.5)
#hperc975 = np.percentile(hcrit, 97.5, axis=0)
hperc975 = fa.nan_percentiles_bias_corr2(vcrit_sliding, hcrit, 97.5)
plt.plot(vcrit_sliding[0,:], hperc025, ls=':', c='r', zorder=7, 
         lw=2)
plt.plot(vcrit_sliding[0,:], hperc975, ls=':', c='r', zorder=7, 
         lw=2, label="95 % uncertainty bounds")



# plt.plot(np.mean(vcrit_sliding,0), np.mean(hcrit,0), ls='-.', c='g', zorder=7, 
#          lw=2, label="Sliding - Mean curve")



plt.ylabel('$h$ (m)')
plt.xlabel('$v$ (m/s)')
plt.legend(loc='best')

plt.xticks(np.arange(0, 9.0+1, 1.0))
plt.xlim(0, 9)
plt.ylim([0,1.5])

# Save the plot and the data now!
folder_s = "V2-results"
np.save(folder_s +"\\"+"V2_vcrit.npy",vcrit_sliding)
np.save(folder_s +"\\"+"V2_hcrit.npy",hcrit)
np.save(folder_s +"\\"+"V2_vcrit_sliding.npy",vcrit_sliding)
np.save(folder_s +"\\"+"V2_hcrit_sliding.npy",hcrit_sliding)
np.save(folder_s +"\\"+"V2_vcrit_toppling.npy",vcrit_toppling)
np.save(folder_s +"\\"+"V2_hcrit_toppling.npy",hcrit_toppling)
plt.savefig(folder_s +"\\"+'V2_curve.pdf', dpi=600, format='pdf',  bbox_inches='tight')
plt.savefig(folder_s +"\\"+'V2_curve.png', dpi=600, format='png',  bbox_inches='tight')
plt.savefig(folder_s +"\\"+'V2_curve.svg', dpi=600, format='svg',  bbox_inches='tight')

end = time.time()
print("Elapsed time: ", end-start)

plt.show()