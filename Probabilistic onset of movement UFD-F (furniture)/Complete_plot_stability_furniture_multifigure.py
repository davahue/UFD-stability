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

folder_F1 = "F1-results - 1000"
folder_F2 = "F2-results - 1000"


Xth = 50 # The threshold at which we define a critical move of vehicles.

# V1 --
F1_vcrit = np.load(folder_F1 +"/"+"F1_vcrit.npy")
F1_hcrit = np.load(folder_F1 +"/"+"F1_hcrit.npy")

"""
Note that these curves are heavily populated with np.nan. These np.nan present
and intrinsic bias: for a given velocity "v", a np.nan in depth "h" mean that
the vehicles is stable for any h. This presents a bias to be considered when
extracting the percentiles. For instance, if 50 % of the values are np.nan, 
then the 50 % stabiity curve would be the last curve, not the one in the middle
of the available values.
"""

F1_hperc= fa.nan_percentiles_bias_corr(F1_vcrit, F1_hcrit, Xth)
    # bounds
F1_h25=  fa.nan_percentiles_bias_corr(F1_vcrit, F1_hcrit, 25)
#F1_h75=  fa.nan_percentiles_bias_corr(V1_vcrit, V1_hcrit, 75)
F1_h75=  np.nanpercentile(F1_hcrit, 75, axis=0)
# V1_hperc= np.nanpercentile(V1_hcrit, Xth, axis=0)
#     # bounds
# V1_h25= np.nanpercentile(V1_hcrit, 25, axis=0)
# V1_h75= np.nanpercentile(V1_hcrit, 75, axis=0)

# V2 --
F2_vcrit = np.load(folder_F2 +"/"+"F2_vcrit.npy")
F2_hcrit = np.load(folder_F2 +"/"+"F2_hcrit.npy")
F2_hperc= np.nanpercentile(F2_hcrit, Xth, axis=0)
    # bounds
F2_h25= np.nanpercentile(F2_hcrit, 25, axis=0)
F2_h75= np.nanpercentile(F2_hcrit, 75, axis=0)


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

maxh = 0.60
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

colorF1 = colorF(.8)
colorF2 = colorF(.4)
#fig, axs = plt.subplots(1, sharex=True, sharey=False, figsize=(7,5))

###############
### FILLING ###
###############

fig = plt.figure(figsize=(9, 4.5))
spec = fig.add_gridspec(2, 3)

ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[1, 0])
#ax3 = fig.add_subplot(spec[0, 2])
#ax4 = fig.add_subplot(spec[1, 2])
#ax5 = fig.add_subplot(spec[2, 2])
axC = fig.add_subplot(spec[0:2, 1:3])

###############
###    F1    ##
###############
for k in range(1,len(F1_vcrit)):
    ax1.plot(F1_vcrit[k,:], F1_hcrit[k,:], ls='-', c=colorF1, zorder=7, alpha=0.015,
             lw=2, rasterized=True)

# V1_hperc50 = np.nanpercentile(V1_hcrit, 50.0, axis=0)
F1_hperc50 = fa.nan_percentiles_bias_corr(F1_vcrit, F1_hcrit, 50)
ax1.plot(F1_vcrit[0,:], F1_hperc50, ls='-', c=colorF1, zorder=7, 
         lw=2, label="50 % percentile")

ax1.semilogy()
ax1.set_xlim(0, maxV)
ax1.set_ylim(1E-3, maxh)
ax1.set_ylabel('$h$ (m)')
ax1.set_xlabel('$v$ (m/s)')

###############
###    F2    ##
###############
for k in range(1,len(F2_vcrit)):
    ax2.plot(F2_vcrit[k,:], F2_hcrit[k,:], ls='-', c=colorF2, zorder=7, alpha=0.015,
             lw=2, rasterized=True)

F2_hperc50 = np.nanpercentile(F2_hcrit, 50.0, axis=0)
ax2.plot(F2_vcrit[0,:], F2_hperc50, ls='-', c=colorF2, zorder=7, 
         lw=2, label="50 % percentile")

ax2.semilogy()
ax2.set_xlim(0, maxV)
ax2.set_ylim(1E-3, maxh)
ax2.set_ylabel('$h$ (m)')
ax2.set_xlabel('$v$ (m/s)')


###############
### COMPLETE ##
###############
# Curves:
step = 1
axC.plot([],[], lw=0, label=f"Furniture' mobility: {Xth/100:.0%}")
axC.plot(F1_vcrit[0,::step], F1_hperc[::step], ls=':', c=colorF1, marker='', markersize=5, zorder=7, lw=2, label="F1: urban furniture")
axC.fill_between(F1_vcrit[0,::step], F1_h25[::step], F1_h75[::step], zorder=4, color = colorF1, alpha=0.55)

axC.plot(F2_vcrit[0,::step], F2_hperc[::step], ls='--', c=colorF2, marker='', markersize=5, zorder=7, lw=2, label="F2: household fixtures")
axC.fill_between(F2_vcrit[0,::step], F2_h25[::step], F2_h75[::step], zorder=4, color = colorF2, alpha=0.55)

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
axC.legend(loc='upper right', fancybox=True, shadow=False, ncol=1,
           facecolor='white', framealpha=0.8, title="" )# bbox_to_anchor=(0.5, -0.1)
#axs.legend(framealpha=1)

fig.tight_layout()

plt.gcf().text(0.02, 0.94, "A", fontsize=14, weight='bold')
plt.gcf().text(0.02+0.10, 0.94-0.29, "F1: urban furniture", c=colorF(.90), fontsize=11)

plt.gcf().text(0.02, 0.94-0.5, "B", fontsize=14, weight='bold')
plt.gcf().text(0.02+0.10, 0.94-0.75, "F2: household fixtures", c=colorF(.90), fontsize=11)


plt.gcf().text(0.02+0.33, 0.94, "C", fontsize=14, weight='bold')


fig.savefig("Complete_plot_stability_furniture.svg", dpi=600)
fig.savefig("Complete_plot_stability_furniture.pdf", dpi=600)
fig.savefig("Complete_plot_stability_furniture.png", dpi=600)

