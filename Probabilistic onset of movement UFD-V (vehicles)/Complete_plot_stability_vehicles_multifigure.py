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

folder_V1 = "V1-results - 1000"
folder_V2 = "V2-results - 1000"
folder_V3 = "V3-results - 1000"
folder_V4 = "V4-results - 1000"
folder_V5 = "V5-results - 1000"

Xth = 50 # The threshold at which we define a critical move of vehicles.

# V1 --
V1_vcrit = np.load(folder_V1 +"/"+"V1_vcrit.npy")
V1_hcrit = np.load(folder_V1 +"/"+"V1_hcrit.npy")

"""
Note that these curves are heavily populated with np.nan. These np.nan present
and intrinsic bias: for a given velocity "v", a np.nan in depth "h" mean that
the vehicles is stable for any h. This presents a bias to be considered when
extracting the percentiles. For instance, if 50 % of the values are np.nan, 
then the 50 % stabiity curve would be the last curve, not the one in the middle
of the available values.
"""

V1_hperc= fa.nan_percentiles_bias_corr(V1_vcrit, V1_hcrit, Xth)
    # bounds
V1_h25=  fa.nan_percentiles_bias_corr(V1_vcrit, V1_hcrit, 25)
#V1_h75=  fa.nan_percentiles_bias_corr(V1_vcrit, V1_hcrit, 75)
V1_h75=  np.nanpercentile(V1_hcrit, 75, axis=0)
# V1_hperc= np.nanpercentile(V1_hcrit, Xth, axis=0)
#     # bounds
# V1_h25= np.nanpercentile(V1_hcrit, 25, axis=0)
# V1_h75= np.nanpercentile(V1_hcrit, 75, axis=0)

# V2 --
V2_vcrit = np.load(folder_V2 +"/"+"V2_vcrit.npy")
V2_hcrit = np.load(folder_V2 +"/"+"V2_hcrit.npy")
V2_hperc= np.nanpercentile(V2_hcrit, Xth, axis=0)
    # bounds
V2_h25= np.nanpercentile(V2_hcrit, 25, axis=0)
V2_h75= np.nanpercentile(V2_hcrit, 75, axis=0)
# V3 --
V3_vcrit = np.load(folder_V3 +"/"+"V3_vcrit.npy")
V3_hcrit = np.load(folder_V3 +"/"+"V3_hcrit.npy")
V3_hperc= np.nanpercentile(V3_hcrit, Xth, axis=0)
    # bounds
V3_h25= np.nanpercentile(V3_hcrit, 25, axis=0)
V3_h75= np.nanpercentile(V3_hcrit, 75, axis=0)
# V4 --
V4_vcrit = np.load(folder_V4 +"/"+"V4_vcrit.npy")
V4_hcrit = np.load(folder_V4 +"/"+"V4_hcrit.npy")
V4_hperc= np.nanpercentile(V4_hcrit, Xth, axis=0)
    # bounds
V4_h25= np.nanpercentile(V4_hcrit, 25, axis=0)
V4_h75= np.nanpercentile(V4_hcrit, 75, axis=0)
# V5 --
V5_vcrit = np.load(folder_V5 +"/"+"V5_vcrit.npy")
V5_hcrit = np.load(folder_V5 +"/"+"V5_hcrit.npy")
V5_hperc= np.nanpercentile(V5_hcrit, Xth, axis=0)
    # bounds
V5_h25= np.nanpercentile(V5_hcrit, 25, axis=0)
V5_h75= np.nanpercentile(V5_hcrit, 75, axis=0)

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
h = np.linspace(0, maxh, 500)

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

colorV, colorF, colorD  = [plt.cm.Oranges, plt.cm.Blues, plt.cm.Greens]
inner_colors = [colorV(.9), colorV(.8), colorV(.7), colorV(.6), colorV(.8),
                colorF(.6), colorF(.4),
                colorD(.3), colorD(.4), colorD(.6), colorD(.7), colorD(.8) ]

colorV1 = colorV(.3)# colorV(.9)
colorV2 = colorV(.4)# colorV(.8)
colorV3 = colorV(.5)# colorV(.7)
colorV4 = colorV(.6)# colorV(.6)
colorV5 = colorV(.8)# colorV(.8)


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
###    V1    ##
###############
for k in range(1,len(V1_vcrit)):
    ax1.plot(V1_vcrit[k,:], V1_hcrit[k,:], ls='-', c=colorV1, zorder=7, alpha=0.01,
             lw=2, rasterized=True)

# V1_hperc50 = np.nanpercentile(V1_hcrit, 50.0, axis=0)
V1_hperc50 = fa.nan_percentiles_bias_corr(V1_vcrit, V1_hcrit, 50)
ax1.plot(V1_vcrit[0,:], V1_hperc50, ls=':', c=colorV1, zorder=7, 
         lw=2, label="50 % percentile")

ax1.set_xlim(0, maxV)
ax1.set_ylim(0, maxh)
ax1.set_ylabel('$h$ (m)')
ax1.set_xlabel('$v$ (m/s)')

###############
###    V2    ##
###############
for k in range(1,len(V2_vcrit)):
    ax2.plot(V2_vcrit[k,:], V2_hcrit[k,:], ls='-', c=colorV2, zorder=7, alpha=0.01,
             lw=2, rasterized=True)

V2_hperc50 = np.nanpercentile(V2_hcrit, 50.0, axis=0)
ax2.plot(V2_vcrit[0,:], V2_hperc50, ls='-', c=colorV2, zorder=7, 
         lw=2, label="50 % percentile")

ax2.set_xlim(0, maxV)
ax2.set_ylim(0, maxh)
ax2.set_ylabel('$h$ (m)')
ax2.set_xlabel('$v$ (m/s)')

###############
###    V3    ##
###############
for k in range(1,len(V3_vcrit)):
    ax3.plot(V3_vcrit[k,:], V3_hcrit[k,:], ls='-', c=colorV3, zorder=7, alpha=0.01,
             lw=2, rasterized=True)

V3_hperc50 = np.nanpercentile(V3_hcrit, 50.0, axis=0)
ax3.plot(V3_vcrit[0,:], V3_hperc50, ls='-.', c=colorV3, zorder=7, 
         lw=2, label="50 % percentile")

ax3.set_xlim(0, maxV)
ax3.set_ylim(0, maxh)
ax3.set_ylabel('$h$ (m)')
ax3.set_xlabel('$v$ (m/s)')

###############
###    V4    ##
###############
for k in range(1,len(V4_vcrit)):
    ax4.plot(V4_vcrit[k,:], V4_hcrit[k,:], ls='-', c=colorV4, zorder=7, alpha=0.01,
             lw=2, rasterized=True)

V4_hperc50 = np.nanpercentile(V4_hcrit, 50.0, axis=0)
ax4.plot(V4_vcrit[0,:], V4_hperc50, ls='--', c=colorV4, zorder=7, 
         lw=2, label="50 % percentile")

ax4.set_xlim(0, maxV)
ax4.set_ylim(0, maxh)
ax4.set_ylabel('$h$ (m)')
ax4.set_xlabel('$v$ (m/s)')

###############
###    V5    ##
###############
for k in range(1,len(V5_vcrit)):
    ax5.plot(V5_vcrit[k,:], V5_hcrit[k,:], ls='-', c=colorV5, zorder=7, alpha=0.01,
             lw=2, rasterized=True)

V5_hperc50 = np.nanpercentile(V5_hcrit, 50.0, axis=0)
ax5.plot(V5_vcrit[0,:], V5_hperc50, ls='-', c=colorV5, zorder=7, 
         lw=2, label="50 % percentile")

ax5.set_xlim(0, maxV)
ax5.set_ylim(0, maxh)
ax5.set_ylabel('$h$ (m)')
ax5.set_xlabel('$v$ (m/s)')

###############
### COMPLETE ##
###############
# Curves:
step = 1
axC.plot([],[], lw=0, label=f"Vehicles' mobility: {Xth/100:.0%}")
axC.plot(V1_vcrit[0,::step], V1_hperc[::step], ls=':', c=colorV1, marker='', markersize=5, zorder=7, lw=2, label="V1: Two-wheelers")
axC.fill_between(V1_vcrit[0,::step], V1_h25[::step], V1_h75[::step], zorder=6, color = colorV1, alpha=0.35)

axC.plot(V2_vcrit[0,::step], V2_hperc[::step], ls='-', c=colorV2, marker='', markersize=5, zorder=7, lw=2, label="V2: Cars")
axC.fill_between(V2_vcrit[0,::step], V2_h25[::step], V2_h75[::step], zorder=6, color = colorV2, alpha=0.35)

axC.plot(V3_vcrit[0,::step], V3_hperc[::step], ls='-.', c=colorV3,  marker='', markersize=5, zorder=7, lw=2, label="V3: Vans")
axC.fill_between(V3_vcrit[0,::step], V3_h25[::step], V3_h75[::step], zorder=6, color = colorV3, alpha=0.35)

axC.plot(V4_vcrit[0,::step], V4_hperc[::step], ls='--', c=colorV4, marker='', markersize=5, zorder=7, lw=2, label="V4: Caravans & RVs")
axC.fill_between(V4_vcrit[0,::step], V4_h25[::step], V4_h75[::step], zorder=6, color = colorV4, alpha=0.35)

axC.plot(V5_vcrit[0,::step], V5_hperc[::step], ls='-', c=colorV5, marker='', markersize=5, zorder=7, lw=3.5, label="V5: Large heavy veh.")
axC.fill_between(V5_vcrit[0,::step], V5_h25[::step], V5_h75[::step], zorder=6, color = colorV5, alpha=0.35)

# Manning:
axC.plot([],[], lw=0, label="Uniform Flow, slope")
axC.plot(v0, h, ls='-', c='k', zorder=5, lw=1, label=f"")
axC.plot(v1, h, ls='-', c='k', zorder=5, lw=1, label=f"")
axC.plot(v2, h, ls='-', c='k', zorder=5, lw=1, label=f"")
axC.plot(v3, h, ls='-', c='k', zorder=5, lw=1, label=f"")

# Fill:
axC.fill_between(v, h_min0, h_max0, zorder=0, color = "k", alpha=0.15, label=f"Flow, $S$ = 1/{1/s0:.0f}")
axC.fill_between(v, h_min1, h_max1, zorder=0, color = "k", alpha=0.25, label=f"Flow, $S$ = 1/{1/s1:.0f}")
axC.fill_between(v, h_min2, h_max2, zorder=0, color = "k", alpha=0.35, label=f"Flow, $S$ = 1/{1/s2:.0f}")
axC.fill_between(v, h_min3, h_max3, zorder=0, color = "k", alpha=0.45, label=f"Flow, $S$ = 1/{1/s3:.0f}")

axC.plot([],[], ls='-',  c='black', alpha=1, lw=1, zorder=10,label="Flow over smooth cement")

axC.set_xlim(0, maxV)
axC.set_ylim(0, maxh)
axC.set_ylabel('$h$ (m)')
axC.set_xlabel('$v$ (m/s)')
axC.legend(loc='upper right', fancybox=True, shadow=False, ncol=2, facecolor='white', framealpha=0.8, title="" )# bbox_to_anchor=(0.5, -0.1)
#axs.legend(framealpha=1)

fig.tight_layout()

plt.gcf().text(0.02, 0.94, "A", fontsize=14, weight='bold')
plt.gcf().text(0.02+0.33, 0.94, "B", fontsize=14, weight='bold')
plt.gcf().text(0.02+0.66, 0.94, "C", fontsize=14, weight='bold')
plt.gcf().text(0.02+0.66, 0.94-0.33, "D", fontsize=14, weight='bold')
plt.gcf().text(0.02+0.66, 0.94-0.66, "E", fontsize=14, weight='bold')
plt.gcf().text(0.02, 0.94-0.33, "F", fontsize=14, weight='bold')


fig.savefig("Complete_plot_stability_vehicles.svg", dpi=600)
fig.savefig("Complete_plot_stability_vehicles.pdf", dpi=600)
fig.savefig("Complete_plot_stability_vehicles.png", dpi=600)

