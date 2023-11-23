# -*- coding: utf-8 -*-

try:
    from IPython import get_ipython
    get_ipython().magic('reset -sf')
except:
    pass

# from PIL import Image
import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
# import glob
# from scipy.optimize import curve_fit
# import pandas as pd
# import os
# import scipy.stats

# In-house libraries:
import faux as fa

###################
### IMPORT DATA ###
###################

folder_DC = "DC-results - 1000"
folder_DM = "DM-results - 1000"
folder_DP = "DP-results - 1000"
folder_DW = "DW-results - 1000"
folder_DO = "DO-results - 1000"

Xth = 50 # The threshold at which we define a critical move of vehicles.

# DC --
DC_vcrit = np.load(folder_DC +"/"+"DC_vcrit.npy")
DC_hcrit = np.load(folder_DC +"/"+"DC_hcrit.npy")

"""
Note that these curves are heavily populated with np.nan. These np.nan present
and intrinsic bias: for a given velocity "v", a np.nan in depth "h" mean that
the vehicles is stable for any h. This presents a bias to be considered when
extracting the percentiles. For instance, if 50 % of the values are np.nan, 
then the 50 % stabiity curve would be the last curve, not the one in the middle
of the available values.
"""

DC_hperc= fa.nan_percentiles_bias_corr(DC_vcrit, DC_hcrit, Xth)
    # bounds
DC_h25=  fa.nan_percentiles_bias_corr(DC_vcrit, DC_hcrit, 25)
DC_h75=  fa.nan_percentiles_bias_corr(DC_vcrit, DC_hcrit, 75)
# DC_hperc= np.nanpercentile(DC_hcrit, Xth, axis=0)
#     # bounds
# DC_h25= np.nanpercentile(DC_hcrit, 25, axis=0)
# DC_h75= np.nanpercentile(DC_hcrit, 75, axis=0)

# DM --
DM_vcrit = np.load(folder_DM +"/"+"DM_vcrit.npy")
DM_hcrit = np.load(folder_DM +"/"+"DM_hcrit.npy")
DM_hperc= fa.nan_percentiles_bias_corr(DM_vcrit, DM_hcrit, Xth)
    # bounds
DM_h25=  fa.nan_percentiles_bias_corr(DM_vcrit, DM_hcrit, 25)
DM_h75=  fa.nan_percentiles_bias_corr(DM_vcrit, DM_hcrit, 75)

# DO --
DO_vcrit = np.load(folder_DO +"/"+"DO_vcrit.npy")
DO_hcrit = np.load(folder_DO +"/"+"DO_hcrit.npy")

DO_hperc= fa.nan_percentiles_bias_corr(DO_vcrit, DO_hcrit, Xth)
    # bounds
DO_h25=  fa.nan_percentiles_bias_corr(DO_vcrit, DO_hcrit, 25)
DO_h75=  fa.nan_percentiles_bias_corr(DO_vcrit, DO_hcrit, 75)
# DP --
DP_vcrit = np.load(folder_DP +"/"+"DP_vcrit.npy")
DP_hcrit = np.load(folder_DP +"/"+"DP_hcrit.npy")
DP_hperc= fa.nan_percentiles_bias_corr(DP_vcrit, DP_hcrit, Xth)
    # bounds
DP_h25=  fa.nan_percentiles_bias_corr(DP_vcrit, DP_hcrit, 25)
DP_h75=  fa.nan_percentiles_bias_corr(DP_vcrit, DP_hcrit, 75)
# DW --
DW_vcrit = np.load(folder_DW +"/"+"DW_vcrit.npy")
DW_hcrit = np.load(folder_DW +"/"+"DW_hcrit.npy")
DW_hperc= fa.nan_percentiles_bias_corr(DW_vcrit, DW_hcrit, Xth)
    # bounds
DW_h25=  fa.nan_percentiles_bias_corr(DW_vcrit, DW_hcrit, 25)
DW_h75=  fa.nan_percentiles_bias_corr(DW_vcrit, DW_hcrit, 75)



###############
### MANNING ###
###############

"""
Table 5-6 Ven Te Chow (1959)

0.010 cement neat
0.013 smooth asphalt
0.016 rough asphalt
0.023 gravel (random stone in mortar)
0.025 gravel (excavated: uniform section clean)
0.050 flood plains, light brushes in Winter
0.060 flood plains, light brushes in Summer

"""

nmin = 0.011 # smooth cement, table 5-6 Ven Te Chow 1959.
nmax = 0.025 # gravel
n = nmin

maxh = 2.0
maxV = 7
# h = np.linspace(0, maxh, 500)
h = 10**np.linspace(-4,1,500)


s0 = 1e-4
s1 = 1e-3
s2 = 1e-2
s3 = 1e-1

# Central value:
v0 = n**-1 * h**(2/3) * s0**(1/2)
v1 = n**-1 * h**(2/3) * s1**(1/2)
v2 = n**-1 * h**(2/3) * s2**(1/2)
v3 = n**-1 * h**(2/3) * s3**(1/2)

# Fill area:

v     = np.linspace(0, 10, 100)
v     = 10**np.linspace(-2, 1, 500)
h_min0 = (v * nmax * s0**(-1/2))**(3/2)
h_max0 = (v * nmin * s0**(-1/2))**(3/2)
h_min1 = (v * nmax * s1**(-1/2))**(3/2)
h_max1 = (v * nmin * s1**(-1/2))**(3/2)
h_min2 = (v * nmax * s2**(-1/2))**(3/2)
h_max2 = (v * nmin * s2**(-1/2))**(3/2)
h_min3 = (v * nmax * s3**(-1/2))**(3/2)
h_max3 = (v * nmin * s3**(-1/2))**(3/2)



############
### PLOT ###
############


font = {'family' : 'Arial', # Ubuntu / we have all in Arial. The text is in Times New Roman. 
        'weight' : 'normal',
        'size'   : 12}


#plt.rcParams["legend.facecolor"] = "k"
#plt.rcParams["legend.facealpha"] = 1

plt.rc('font', **font)

#fig, axs = plt.subplots(1, sharex=True, sharey=False, figsize=(7,5))

###############
### FILLING ###
###############

fig = plt.figure(figsize=(9, 7))
spec = fig.add_gridspec(3, 3)

ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[0, 2])
ax4 = fig.add_subplot(spec[1, 2])
ax5 = fig.add_subplot(spec[2, 2])
axC = fig.add_subplot(spec[1:3, 0:2])

###############
###    DC    ##
###############
for k in range(1,len(DC_vcrit)):
    ax1.plot(DC_vcrit[k,:], DC_hcrit[k,:], ls='-', c='green', zorder=7, alpha=0.01,
             lw=2, rasterized=True)

# DC_hperc50 = np.nanpercentile(DC_hcrit, 50.0, axis=0)
DC_hperc50 = fa.nan_percentiles_bias_corr(DC_vcrit, DC_hcrit, 50)
ax1.plot(DC_vcrit[0,:], DC_hperc50, ls=':', c='green', zorder=7, 
         lw=2, label="50 % percentile")

# ax1.loglog()
ax1.semilogy()
ax1.set_xlim(1E-2, maxV)
ax1.set_ylim(1E-3, maxh)
ax1.set_ylabel('$h$ (m)')
ax1.set_xlabel('$v$ (m/s)')

###############
###    DM    ##
###############
for k in range(1,len(DM_vcrit)):
    ax2.plot(DM_vcrit[k,:], DM_hcrit[k,:], ls='-', c='b', zorder=7, alpha=0.01,
             lw=2, rasterized=True)

DM_hperc50 = fa.nan_percentiles_bias_corr(DM_vcrit, DM_hcrit, 50)
ax2.plot(DM_vcrit[0,:], DM_hperc50, ls='-', c='b', zorder=7, 
         lw=2, label="50 % percentile")

# ax2.loglog()
ax2.semilogy()
ax2.set_xlim(1E-2, maxV)
ax2.set_ylim(1E-3, maxh)
ax2.set_ylabel('$h$ (m)')
ax2.set_xlabel('$v$ (m/s)')


###############
###    DP    ##
###############
for k in range(1,len(DP_vcrit)):
    ax3.plot(DP_vcrit[k,:], DP_hcrit[k,:], ls='-', c='orange', zorder=7, alpha=0.01,
             lw=2, rasterized=True)

DP_hperc50 = fa.nan_percentiles_bias_corr(DP_vcrit, DP_hcrit, 50)
ax3.plot(DP_vcrit[0,:], DP_hperc50, ls='--', c='orange', zorder=7, 
         lw=2, label="50 % percentile")

ax3.loglog()
ax3.set_xlim(1E-2, maxV)
ax3.set_ylim(1E-3, maxh)
ax3.set_ylabel('$h$ (m)')
ax3.set_xlabel('$v$ (m/s)')

###############
###    DW    ##
###############
for k in range(1,len(DW_vcrit)):
    ax4.plot(DW_vcrit[k,:], DW_hcrit[k,:], ls='-', c='chocolate', zorder=7, alpha=0.01,
             lw=2, rasterized=True)

DW_hperc50 = fa.nan_percentiles_bias_corr(DW_vcrit, DW_hcrit, 50)
ax4.plot(DW_vcrit[0,:], DW_hperc50, ls='-', c='chocolate', zorder=7, 
         lw=2, label="50 % percentile")

ax4.loglog()
ax4.set_xlim(1E-2, maxV)
ax4.set_ylim(1E-3, maxh)
ax4.set_ylabel('$h$ (m)')
ax4.set_xlabel('$v$ (m/s)')



###############
###    DO    ##
###############
for k in range(1,len(DO_vcrit)):
    ax5.plot(DO_vcrit[k,:], DO_hcrit[k,:], ls='-', c='darkmagenta', zorder=7, alpha=0.01,
             lw=2, rasterized=True)

DO_hperc50 = fa.nan_percentiles_bias_corr(DO_vcrit, DO_hcrit, 50)
ax5.plot(DO_vcrit[0,:], DO_hperc50, ls='-.', c='darkmagenta', zorder=7, 
         lw=2, label="50 % percentile")

ax5.loglog()
ax5.set_xlim(1E-2, maxV)
ax5.set_ylim(1E-3, maxh)
ax5.set_ylabel('$h$ (m)')
ax5.set_xlabel('$v$ (m/s)')

###############
### COMPLETE ##
###############
# Curves:
step = 1
axC.plot([],[], lw=0, label=f"UFD-H's mobility: {Xth/100:.0%}")
axC.plot(DC_vcrit[0,::step], DC_hperc[::step], ls=':', c='g', marker='', markersize=5, zorder=5, lw=2, label="DC: construction")
# axC.fill_between(DC_vcrit[0,::step], DC_h25[::step], DC_h75[::step], zorder=0, color = "k", alpha=0.05)

axC.plot(DM_vcrit[0,::step], DM_hperc[::step], ls='-', c='b', marker='', markersize=5, zorder=5, lw=2, label="DM: metal")
# axC.fill_between(DM_vcrit[0,::step], DM_h25[::step], DM_h75[::step], zorder=0, color = "b", alpha=0.05)


axC.plot(DP_vcrit[0,::step], DP_hperc[::step], ls='--', c='orange', marker='', markersize=5, zorder=5, lw=2, label="DP: plastic")
# axC.fill_between(DP_vcrit[0,::step], DP_h25[::step], DP_h75[::step], zorder=0, color = "k", alpha=0.05)

axC.plot(DW_vcrit[0,::step], DW_hperc[::step], ls='-', c='chocolate', marker='', markersize=5, zorder=5, lw=2, label="DW: wood")
# axC.fill_between(DW_vcrit[0,::step], DW_h25[::step], DW_h75[::step], zorder=0, color = "k", alpha=0.05)

axC.plot(DO_vcrit[0,::step], DO_hperc[::step], ls='-.', c='darkmagenta',  marker='', markersize=5, zorder=5, lw=2, label="DO: other")
# axC.fill_between(DO_vcrit[0,::step], DO_h25[::step], DO_h75[::step], zorder=0, color = "k", alpha=0.05)

# Manning:
# axC.plot([],[], lw=0, label="Uniform Flow, slope")
axC.plot(v0, h, ls='-', c='k', zorder=5, lw=1, label=f"")
axC.plot(v1, h, ls='-', c='k', zorder=5, lw=1, label=f"")
axC.plot(v2, h, ls='-', c='k', zorder=5, lw=1, label=f"")
axC.plot(v3, h, ls='-', c='k', zorder=5, lw=1, label=f"")

# Fill:
axC.fill_between(v, h_min0, h_max0, zorder=0, color = "k", alpha=0.15) #, label=f"Flow, $S$ = 1/{1/s0:.0f}")
axC.fill_between(v, h_min1, h_max1, zorder=0, color = "k", alpha=0.25) #, label=f"Flow, $S$ = 1/{1/s1:.0f}")
axC.fill_between(v, h_min2, h_max2, zorder=0, color = "k", alpha=0.35) #, label=f"Flow, $S$ = 1/{1/s2:.0f}")
axC.fill_between(v, h_min3, h_max3, zorder=0, color = "k", alpha=0.45) #, label=f"Flow, $S$ = 1/{1/s3:.0f}")

axC.plot([],[], ls='-',  c='black', alpha=1, lw=1, zorder=10) #,label="Flow over smooth cement")


axC.loglog()
# axC.semilogy()
axC.set_xlim(1E-2, maxV)
axC.set_ylim(1E-3, maxh)
axC.set_ylabel('$h$ (m)')
axC.set_xlabel('$v$ (m/s)')
axC.legend(loc='upper left', fancybox=True, shadow=False, ncol=1, facecolor='white', framealpha=0.8, title="" )# bbox_to_anchor=(0.5, -0.1)
#axs.legend(framealpha=1)


fig.tight_layout()

plt.gcf().text(0.02, 0.94, "A", fontsize=14, weight='bold')
plt.gcf().text(0.03+0.09, 0.93-0.15, "Construction debris", c="g", fontsize=11)
plt.gcf().text(0.02+0.33, 0.94, "B", fontsize=14, weight='bold')
plt.gcf().text(0.03+0.33+0.09, 0.93-0.15, "Metal debris", c="b", fontsize=11)
plt.gcf().text(0.02+0.66, 0.94, "C", fontsize=14, weight='bold')
plt.gcf().text(0.03+0.66+0.09, 0.93, "Plastic debris", c="orange", fontsize=11)
plt.gcf().text(0.02+0.66, 0.94-0.33, "D", fontsize=14, weight='bold')
plt.gcf().text(0.03+0.66+0.09, 0.93-0.45, "Wood drifters", c="chocolate", fontsize=11)
plt.gcf().text(0.02+0.66, 0.94-0.66, "E", fontsize=14, weight='bold')
plt.gcf().text(0.03+0.66+0.09, 0.93-0.65-0.0, "Other litter debris", c="darkmagenta", fontsize=11)
plt.gcf().text(0.02, 0.94-0.33, "F", fontsize=14, weight='bold')



fig.savefig("Complete_plot_stability_debris.svg", dpi=600)
fig.savefig("Complete_plot_stability_debris.pdf", dpi=600)
fig.savefig("Complete_plot_stability_debris.png", dpi=600)

