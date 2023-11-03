# -*- coding: utf-8 -*-
"""
Created: 19.01.2022
Onset of stability for urban flood drifters (UFD).
Verification case: vehicles from literature.
This code imports calculations done for different cars, and plots them together.
For new calculations, consider the Onset_test_Q7.py code.
@author: D. Valero
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# in-house library to support auxiliary calculations
# import faux as fa

# PLOTTING --------------------------------------------------------------------        
# Quick plot ------------------------------------------------------------------

# Import data
folder = "Protop car tests"
print("Data digitized from original studies:")
# try: # Windows
#     dataYaris1 = pd.read_excel(folder+"\\"+"Yaris_red_triangles.xlsx")
#     dataYaris2 = pd.read_excel(folder+"\\"+"Yaris_blue_triangles.xlsx")
#     dataFestiva = pd.read_excel(folder+"\\"+"Festiva_brown_markers.xlsx")
#     dataPatrol = pd.read_excel(folder+"\\"+"Patrol_green_markers.xlsx")
#     dataGolf = pd.read_excel(folder+"\\"+"2016 Kramer Golf.xlsx")
#     dataAccord = pd.read_excel(folder+"\\"+"2013 Xia Honda Accord.xlsx")
#     dataQ7 = pd.read_excel(folder+"\\"+"2013 Xia Audi Q7.xlsx")
# except: # Linux
#     dataYaris1 = pd.read_excel(folder+"/"+"Yaris_red_triangles.xlsx")
#     dataYaris2 = pd.read_excel(folder+"/"+"Yaris_blue_triangles.xlsx")
#     dataFestiva = pd.read_excel(folder+"/"+"Festiva_brown_markers.xlsx")
#     dataPatrol = pd.read_excel(folder+"/"+"Patrol_green_markers.xlsx")
#     dataGolf = pd.read_excel(folder+"/"+"2016 Kramer Golf.xlsx")
#     dataAccord = pd.read_excel(folder+"/"+"2013 Xia Honda Accord.xlsx")
#     dataQ7 = pd.read_excel(folder+"/"+"2013 Xia Audi Q7.xlsx")

# Fig opts
plt.figure(figsize=(6,4))



[_, hcrit_flot_Yaris] = np.loadtxt("Flotation limit state (Yaris).txt")
[v, hcrit_sliding_Yaris] = np.loadtxt("Sliding limit state (Yaris).txt")

[_, hcrit_flot_Festiva] = np.loadtxt("Flotation limit state (Festiva).txt")
[v, hcrit_sliding_Festiva] = np.loadtxt("Sliding limit state (Festiva).txt")

[_, hcrit_flot_Golf] = np.loadtxt("Flotation limit state (Golf).txt")
[v, hcrit_sliding_Golf] = np.loadtxt("Sliding limit state (Golf).txt")

[_, hcrit_flot_Accord] = np.loadtxt("Flotation limit state (Accord).txt")
[v, hcrit_sliding_Accord] = np.loadtxt("Sliding limit state (Accord).txt")

[_, hcrit_flot_Patrol] = np.loadtxt("Flotation limit state (Patrol).txt")
[v, hcrit_sliding_Patrol] = np.loadtxt("Sliding limit state (Patrol).txt")

[_, hcrit_flot_Q7] = np.loadtxt("Flotation limit state (Q7).txt")
[v, hcrit_sliding_Q7] = np.loadtxt("Sliding limit state (Q7).txt")

# Header
plt.scatter([], [], lw=0, marker='', label=r'$\mathbf{Bibliography~datasets}$')    


# plt.scatter(0, hcrit_flot_Yaris, marker= 'o', s=80, zorder=7,
#             edgecolor='b', facecolor='w', label='Flotation limit state  (Yaris)')
plt.plot(v, hcrit_sliding_Yaris, ls='--', c='b', zorder=7,
         lw=2, label='Sliding limit state (Yaris)')

# plt.scatter(0, hcrit_flot_Festiva, marker= 'o', s=80, zorder=7,
#             edgecolor='hotpink', facecolor='w', label='Flotation limit state  (Festiva)')
plt.plot(v, hcrit_sliding_Festiva, ls='--', c='hotpink', zorder=7,
         lw=2, label='Sliding limit state (Festiva)')
#plt.plot(vxcrit_toppling, hcrit_toppling, ls='-', c='k', label='Toppling')

# plt.scatter(0, hcrit_flot_Golf, marker= 'o', s=80, zorder=7,
#             edgecolor='orange', facecolor='w', label='Flotation limit state  (Golf)')
plt.plot(v, hcrit_sliding_Golf, ls='--', c='orange', zorder=7,
         lw=2, label='Sliding limit state (Golf)')

# plt.scatter(0, hcrit_flot_Accord, marker= 'o', s=80, zorder=7,
#             edgecolor='yellowgreen', facecolor='w', label='Flotation limit state  (Accord)')
plt.plot(v, hcrit_sliding_Accord, ls='--', c='yellowgreen', zorder=7,
         lw=2, label='Sliding limit state (Accord)')

# plt.scatter(0, hcrit_flot_Patrol, marker= 'o', s=80, zorder=7,
#             edgecolor='darkolivegreen', facecolor='w', label='Flotation limit state  (Patrol)')
plt.plot(v, hcrit_sliding_Patrol, ls='--', c='darkolivegreen', zorder=7,
         lw=2, label='Sliding limit state (Patrol)')

# plt.scatter(0, hcrit_flot_Q7, marker= 'o', s=80, zorder=7,
#             edgecolor='r', facecolor='w', label='Flotation limit state  (Q7)')
plt.plot(v, hcrit_sliding_Q7, ls='--', c='r', zorder=7,
         lw=2, label='Sliding limit state (Q7)')


# Header
plt.scatter([], [], lw=0, marker='', label=' '), plt.scatter([], [], lw=0, marker='', label=r'$\mathbf{Limit~states~of~stability}$')    

    
# Prototype data:
# plt.scatter(dataYaris1["V"], dataYaris1["h"], c='b',marker='^', label="Yaris, dataset 1 of Smith et al. (2019)")    
# plt.scatter(dataYaris2["V"], dataYaris2["h"], c='b',marker='v', label="Yaris, dataset 2 of Smith et al. (2019)")    
# plt.scatter(dataFestiva["V"], dataFestiva["h"], s=20, c='hotpink',marker='s', label="Festiva, Smith et al. (2019)")    
# plt.scatter(dataGolf["V"], dataGolf["h"], s=20, c='orange',marker='D', label="Golf, Kramer et al. (2016)")
# plt.scatter(dataAccord["V"], dataAccord["h"], s=20, c='yellowgreen',marker='>', label="Accord, Xia et al. (2014)")    
# plt.scatter(dataPatrol["V"], dataPatrol["h"], s=20, c='darkolivegreen',marker='d', label="Patrol, Smith et al. (2019)")    
# plt.scatter(dataQ7["V"], dataQ7["h"], s=20, c='r',marker='<', label="Q7, Xia et al. (2014)")    
   

# Text
plt.text(0.5, 0.1, 'Stable', fontsize=12, color='g', weight='bold')
plt.text(5.5, 0.6, 'Unstable', fontsize=12, color='r', weight='bold')

plt.text(3., 0.7, 'Buoyancy\ndominated', fontsize=12, color='k', rotation=90)
plt.arrow(4., 0.74, 0.0, 0.05, width=0.01, head_width=0.15, head_length=0.05, fc='k', ec='k')

plt.text(7., 0.35, '   Drag\ndominated', fontsize=12, color='k', rotation=0)
plt.arrow(7.4, 0.28, 0.45, 0.00, width=0.005, head_width=0.028, head_length=0.35, fc='k', ec='k')


# plt.legend()

# Axis
plt.xlim(v.min(), v.max())
plt.ylabel('$h$ (m)')
plt.xlabel('$v$ (m/s)')
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

