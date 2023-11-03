# -*- coding: utf-8 -*-
"""
Created: 19.01.2022
Onset of stability for urban flood drifters (UFD).
Exemplary code for the calculation of a given car (Audi Q7).
The code can be easily adapted to other types of cars by
modifying their physical characteristics (section: Params definition).
"Other relevant parameters" should only be considered when
there is new information.
@author: D. Valero

"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# in-house library to support auxiliary calculations
import faux as fa

# Params definition ----------------------------------------------------------
M =  2345 # kg, car mass
Lx =  5.089 # m, car length
Ly =  1.983 # m, car width
Lz =  1.737 # m, car height
zc =  0.204 # m, distance from street floor to car underfloor
# Other relevant parameters --------------------------------------------------
rhow = 998 # kg/m3, water
g = 9.81 # m/s2, gravity
CD = 0.85 # 1.38 # drag coefficient # FIG.11 Smith et al. 2019
fV1 = 0.05
fV2 = 0.50

mu = 0.3 # FIG.11 Smith et al. 2019
# h, vx as key parameters ----------------------------------------------------
nh = 200; nv = 400
h = np.linspace(0, Lz, nh)
v = np.linspace(0, 9.5, nv)
# ----------------------------------------------------------------------------
# Start calculations ---------------------------------------------------------

# GET LIMIT STATES

# Flotation limit state ------------------------------------------------------

Rz = M*g # Resistance
fV = fa.get_fV(fV1, fV2, h, zc, Lz) # Get an eta for each h?
Sz = rhow*g*fV*h*Lx*Ly # Load

hflot = Rz/(rhow*g*fV*Lx*Ly)
ZC = np.where(np.diff(np.sign(hflot-h)))[0] # zero crossings
hcrit_flot = h[ZC[0]]

# Sliding limit state ------------------------------------------------------
fA = fa.get_fA_V2(h,Lz)


H, V = np.meshgrid(h,v)
LHSs = np.zeros_like(H); RHSs = np.zeros_like(H)
Sy = np.zeros_like(H)
for i in range(0, nv): # - velocity
    for j in range(0, nh): # - depths

        # Alternative 2 - nested in y
        Sy[i,j] = 0.5 *(rhow*CD*fA[j]*Lx*H[i,j]*V[i,j]*V[i,j])
        LHSs[i,j] = Sy[i,j] # Sx
        RHSs[i,j] = mu*(Rz-Sz[j]) # Rx
        

# css = plt.contour(V, H, (LHSs-RHSs), levels=[-0.0001,0,0.0001], colors='k')
# ps = css.collections[0].get_paths()[0]
# vs = ps.vertices
# vxcrit_sliding = vs[:,0]
# hcrit_sliding = vs[:,1]

# Toppling limit state ------------------------------------------------------

LHSt = np.zeros_like(H); RHSt = np.zeros_like(H)
for i in range(0, nv):
    for j in range(0, nh):

        # Alternative 2 - nested in x and y
        
        if H[i,j] > zc:
            zp = (H[i,j]-zc)/2
        else:
            zp = (H[i,j])/2

        yp = Ly/2
        
        LHSt[i,j] = Sy[i,j]*zp + Sz[j]*yp
        RHSt[i,j] = Rz*Ly/2
        
        
"""
I need a zero-finding routine to map:
    - sliding: LHSs-RHSs
    - toppling: LHSt-RHSt
"""


arrs = LHSs-RHSs        
hcrit_sliding = fa.zero_finding(arrs, v, h)

arrt = LHSt-RHSt
hcrit_toppling = fa.zero_finding(arrt, v, h)


# PLOTTING --------------------------------------------------------------------        

# Quick plot ------------------------------------------------------------------


# Fig opts
plt.figure(figsize=(6,4))

# plt.scatter(0, hcrit_flot, marker= 'o', s=80, zorder=7,
#             edgecolor='r', facecolor='w', label='Flotation limit state  (Q7)')
plt.plot(v, hcrit_sliding, ls='--', c='r', zorder=7,
         lw=2, label='Sliding limit state (Q7)')

plt.plot(v, hcrit_toppling, ls=':', c='r', label='Toppling limit state (Q7)')

# np.savetxt("Flotation limit state (Q7).txt", [0, hcrit_flot])
# np.savetxt("Sliding limit state (Q7).txt", [v, hcrit_sliding])
# np.savetxt("Toppling limit state (Q7).txt", [v, hcrit_toppling])

plt.xlim(v.min(), v.max())

plt.ylabel('$h$ (m)')
plt.xlabel('$v$ (m/s)')

plt.legend(loc='best', ncol=1)

plt.xticks(np.arange(0, 9.0+1, 1.0))
plt.xlim(0, 9)
plt.ylim([0,1])

# plt.savefig('Stability.pdf', dpi=600, format='pdf',  bbox_inches='tight')
# plt.savefig('Stability.png', dpi=600, format='png',  bbox_inches='tight')
# plt.savefig('Stability.svg', dpi=600, format='svg',  bbox_inches='tight')

# plt.savefig('Stability_legend.pdf', dpi=600, format='pdf',  bbox_inches='tight')
# plt.savefig('Stability_legend.png', dpi=600, format='png',  bbox_inches='tight')
# plt.savefig('Stability_legend.svg', dpi=600, format='svg',  bbox_inches='tight')

plt.show()
