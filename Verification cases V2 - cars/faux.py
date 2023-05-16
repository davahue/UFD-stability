# -*- coding: utf-8 -*-
"""
19.01.2022

@author: D. Valero

"""
import csv
import numpy as np
from scipy import interpolate
from scipy.special import comb

from matplotlib import pyplot as plt

# auxiliary functions --------------------------------------------------------

def zero_finding(arr, v, h):
    """
    Inputs a 2D array with RHS-LHS and finds the points where its zero.
    
    The x-variable is velocity:
        v = np.linspace(0, 9.5, nv)
    and then we find h at which RHS-LHS == 0.
    
    The method should be robust vs nan.
    
    Tested against plt.contour function.
    """
    
    [nx, ny] = np.shape(arr)
    h_stab = np.zeros(nx)
    # print(nx) # This is v-direction
    # print(ny)
    
    for i in range(nx):
        jleft = np.where(np.diff(np.sign(arr[i,:])))[0][0]
        jright = jleft+1 # next index, arr is negative
        
        y1 = arr[i,jleft]
        y2 = arr[i,jright]
        
        X = y1/(y1-y2)
        
        h_stab[i] = h[jleft] +X*(h[jright] - h[jleft])
    
    
    return h_stab

def find_median_curve(XX,YY):
    """
    Takes a list of lines, and finds the median one based on the
    concept of depth of data.
    Input:
        XX: x-data
        YY: y-data
    This function transforms data on equally sized and then computes the depth
    The line with larger depth represents the median, as stated by:
        Sara López-Pintado & Juan Romo (2009) On the
        Concept of Depth for Functional Data, Journal of the American
        Statistical Association, 104:486, 718-734, DOI: 10.1198/jasa.2009.0108
    """
    N = len(YY)
    M = [] # number of elements of each hydrograph
    for i in range(0,N):
        M.append(len(YY[i]))
    
    Data_grid = np.zeros([N,np.max(M)])
    xgrid = np.linspace(0,1,np.max(M)) # poitns to interpolate in.
    for i in range(0,N):
        xdata = XX[i]
        ydata = YY[i]
        
        finterp = interpolate.interp1d(xdata, ydata, kind='linear')
        for j in range(0,np.max(M)):
            Data_grid[i,j] = finterp(xgrid[j])
    """
    # Verification plot -- verified 03/08/2021
    plt.figure()
    for i in range(0,N):
        
        plt.plot(xgrid, NHydrograph_grid[i,:],
                 c='k', alpha=0.075)
    """
    """
    Data is now ready (same length), we can obtain the depth. Depth defined
    as per:    
    statsmodels.graphics.functional.banddepth(data, method='MBD')
    """
    # The lines below must be credited to statsmodels.
    # See band depth: https://www.statsmodels.org/dev/_modules/statsmodels/graphics/functional.html#banddepth
    n, p = Data_grid.shape
    rv = np.argsort(Data_grid, axis=0)
    rmat = np.argsort(rv, axis=0) + 1
    # method: BD2
    down = np.min(rmat, axis=1) - 1
    up = n - np.max(rmat, axis=1)
    depth = (up * down + n - 1) / comb(n, 2)
    
    # The median curve is that holding maximum depth, see: Neto The Concept of
    # Depth in Statistics
    # This is also stated in: Sara López-Pintado & Juan Romo (2009) On the
    # Concept of Depth for Functional Data, Journal of the American
    # Statistical Association, 104:486, 718-734, DOI: 10.1198/jasa.2009.0108
    i_med = np.argmax(depth)
    
    return i_med



def get_var(folder, Nr, var_str):
    """
    Access the $var_str$.txt file in the specified folder and generates Nr samples
    following the predefined probability distribution.
    """
    
    file = folder + "//" + var_str + ".txt"
    with open(file) as f:
        csv.reader(file, delimiter='\t')
        lines = f.readlines()
    
    [Min, Mode, Max, Std, distr] = lines[1].split("\t")
    
    try:
        Min = float(Min); Max = float(Max)
    except Exception:
        Min = []; Max = []
        
    try:
        Std = float(Std)
        Mode = float(Mode)
    except Exception:
        Std = []
    
    print("Parameter: ", var_str)
    print("Min \t Mode \t Max \t STD\n", [Min, Mode, Max, Std])

    if distr == "Uniform":
        var_rand = np.random.uniform(Min,Max,Nr)
    elif distr == "Triangular":
        var_rand = np.random.triangular(Min,Mode,Max,Nr)
    elif distr == "Normal":
        var_rand = np.random.normal(Mode,Std,Nr)
    else:
        print("Distribution not found for variable: ", var_str)
        var_rand = []
    
    return var_rand



# SOLVER ----------------------------------------------------------------------

# def get_stability_curves(M,Lx,Ly,Lz,zuf,Cd,mu,PSI_std,nh,nv):
    
#     g = 9.8
#     rhow = 998
#     # h, vx as key parameters
#     # nh = 200; nv = 400
#     h = np.linspace(0, Lz, nh)
#     v = np.linspace(0, 9.5, nv)
#     # 
#     # Start calculations 

#     # Flotation limit state 

#     Ry = M*g # Resistance
#     eta = get_eta(h, zuf, Lz) # Get an eta for each h?
#     Sy = rhow*g*eta*h*Lx*Ly # Load

#     hflot = Ry/(rhow*g*eta*Lx*Ly)
#     ZC = np.where(np.diff(np.sign(hflot-h)))[0] # zero crossings
#     hcrit_flot = h[ZC[0]]

#     # Sliding limit state
#     fA = get_PSI_car(h,Lz) + PSI_std # MEAN + DEV.


#     H, V = np.meshgrid(h,v)
#     LHSs = np.zeros_like(H); RHSs = np.zeros_like(H)
#     Sx = np.zeros_like(H)
#     for i in range(0, nv): # - velocity
#         for j in range(0, nh): # - depths

#             Sx[i,j] = 0.5 *(rhow*Cd*fA[j]*Lx*H[i,j]*V[i,j]*V[i,j])
#             LHSs[i,j] = Sx[i,j] # Sx
#             RHSs[i,j] = mu*(Ry-Sy[j]) # Rx
            
            
#     # Toppling limit state

#     LHSt = np.zeros_like(H); RHSt = np.zeros_like(H)
#     for i in range(0, nv):
#         for j in range(0, nh):
           
#             if H[i,j] > zuf:
#                 yp = (H[i,j]-zuf)/2
#             else:
#                 yp = (H[i,j])/2
#             #xp1 = Ly*(1-eta[j]/3)
#             xp2 = Ly/2
#             xp = xp2 # (xp1+xp2)/2
            
#             LHSt[i,j] = Sx[i,j]*yp + Sy[j]*xp
#             RHSt[i,j] = Ry*Ly/2
            
#     """
#     I need a zero-finding routine to map:
#         - sliding: LHSs-RHSs
#         - toppling: LHSt-RHSt
#     """
#     # Slidding:
#     arr_sliding = LHSs-RHSs        
#     hcrit_sliding = zero_finding(arr_sliding, v, h)
#     # Toppling:
#     arr_toppling = LHSt-RHSt        
#     hcrit_toppling = zero_finding(arr_toppling, v, h)

#     return hcrit_flot, v, hcrit_sliding, v, hcrit_toppling

# Limit state: flotation -----------------------------------------------------

def get_fV(fV1, fV2, h, zc, Lz):
    """
    gets the eta paramerter, related to the percentage of underfloor area wetted.

    h: water level (vector).
    z_uf: elevation of the underfloor relative to the street level.
    Lz: height of the boundaing box.
    """
    # fV1 = 0.05 # MIN percentage of bounding box area submerged, below clearance.
    # 0.075 full wheels, if they were a solid squared block
    
    
    fV = np.zeros_like(h)
    for i in range(0, len(h)):
        if h[i]<zc: # below the car underfloor
            fV[i] = fV1 # wheels.
        elif h[i]<Lz: # over the car underfloor, and up to the car height
            # fV[i] = fVmax
            fV[i] = (fV1*zc + fV2*(h[i]-zc)) / h[i]
            
        else: # over the car height
            fV[i] = np.nan
    
    return fV


# Limit state: sliding 

def get_fA_V2(h,Lz):
    """
    This function obtains the projected wetted area in the frontal direction
    for the case of a car. This is based on analysis of side images of cars.
    Total area:
        Axz = fA*Lx*Lz
        
        fA = a0 + a1*x + a2*x**2 + a3*x**3 
    """
        
    a0 = -0.010491 # zero-order param
    a1 = 2.761454 # first-order param
    a2 = -2.904696 # second-order param
    a3 = 0.845684 # third-order param
    # a4 = 0 # fourth-order param
    
    fA = a0 + a1*(h/Lz) + a2*(h/Lz)**2 + a3*(h/Lz)**3 #+ a4*(h/Lz)**4 
    
    return fA

# def get_PSI_car(h,Lz):
#     """
#     This function obtains the projected wetted area in the frontal direction
#     for the case of a car. This is based on analysis of side images of cars.
#     Total area:
#         Axz = PSI*Lx*Lz
#     """
#     # PSI = np.ones_like(h)   
    
#     a0 = -0.022209790806162293 # zero-order param
#     a1 = 2.7013864027501824 # first-order param
#     a2 = -2.738099023102069 # second-order param
#     a3 = 0.7524523329486932 # third-order param
#     a4 = 0 # fourth-order param
    
#     PSI = a0 + a1*(h/Lz) + a2*(h/Lz)**2 + a3*(h/Lz)**3 + a4*(h/Lz)**4 
    
#     #PSI = 0.85*np.ones_like(h)  
#     return PSI




# Limit state: toppling ------------------------------------------------------