# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 18:41:47 2015
@author: Asif Arain

Description: 
Sensor placement problem (SPP) for robot assisted gas tomography (RAGT):

G = ...
U = beta*sum(V,2)
max. C'(GC.*U)
s.t. |C| <= number and C in {0,1}

"""

# from spp import solveSPP_hotspot
from spp import *

# from spp import solveSPP_fusion
import time  # time
import math  # math
from gurobipy import *  # gurobi
import numpy as np  # np
from numpy import linalg as LA  # norm
import scipy.io  # matlab files

# import statistics as stat

# initialize computation time
#time_start = time.clock()
time_start = time.time()

# n = 3 # number of allowed sensing configutions

# alpha = 0.75
# beta  = 0.25
# gamma = 0.50

#gdo = 0.1

# problem = 'spp hotspots'
# problem = 'spp fusion'


titlefile = scipy.io.loadmat('SPP_HOTSPOT_DATA.mat')

EXP_TITLE = titlefile['EXP_TITLE']
EXP_TITLE = str(''.join(EXP_TITLE))

strategy = titlefile['strategy']
strategy = str(''.join(strategy))

hotspot_num = titlefile['hotspot_num']
hotspot_num = int(hotspot_num)

# number of allowed conf
n = titlefile['NumOfAllowedConf']
n = int(n)

# alpha
alpha = titlefile['alpha']
alpha = float(alpha)

# beta
beta = titlefile['beta']
beta = float(beta)

# print EXP_TITLE
# print hotspot_num

# gamma
gamma = titlefile['gamma']
gamma = float(gamma)

print '*************************************************************'
print ' Optimization for the local solution of hotspot# %g' % (hotspot_num)
print '*************************************************************'

if strategy in 'two-step-exploration-jfr':

    print '****** plan type: JFR17 ******'

    # gamma
    #gamma = titlefile['gamma']
    #gamma = float(gamma)
    
    print 'alpha: %g' % (alpha)
    print 'beta:  %g' % (beta)
    print 'gamma: %g' % (gamma)


    # for hotspot_num in range(1,2,1):

    # --- path to read/write solutions    
    #FilePath = 'results/{0}/tomography_jfr17/local-solutions'.format(EXP_TITLE)
    #FilePath = ("results/{0}/tomography_jfr17/"
    #            "alpha{1:03d}_beta{2:03d}_gamma{3:03d}/local-solutions").\
    #             format(EXP_TITLE,int(alpha*100),int(beta*100),int(gamma*100))

    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/local-solutions").\
                 format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))
                 
    print FilePath

    # Data from matlab files
    matfile = scipy.io.loadmat('{0}/preprocess_XVD_hotspot{1}.mat'.\
                               format(FilePath, hotspot_num))    
    # V   = matfile['V']
    V = matfile['wV']
    G = matfile['conf_crossAngles_G']
    D = matfile['D']
    Vdo = matfile['Conf_DistOrnGainV']

    G = G.tolist()
    V = V.tolist()
    # Vdo = Vdo.tolist()

    # number of conf to be selected (float)
    nc = float(n)    
    print 'number of conf to select: %g' % (nc)
    
    

    # alpha times G
    for i in range(np.shape(G)[0]):
        for j in range(np.shape(G)[1]):
            G[i][j] = (alpha * beta) * (1 / (nc * (nc - 1))) * G[i][j]
    # alpha times G
    # for i in range(np.shape(G)[0]):
    #    for j in range(np.shape(G)[1]):
    #        G[i][j] = (alpha*beta)*G[i][j]



    # Uc is V gain    
    # Uc = [  sum(x) for x in zip(*V) ]
    # Uc = [float(i) for i in Uc]
    Uc_w = [sum(i) for i in zip(*V)]
    Uc_w = [float(i) for i in Uc_w]

    Vdo = [float(j) for j in Vdo]

    Uc = [(gamma * i) + ((1 - gamma) * j) for i, j in zip(Uc_w, Vdo)]

    # gamma times normalized Uc        
    # Uc[:] = [x/max(Uc)*gamma for x in Uc]
    # Uc[:] = [x/max(Uc) for x in Uc] # temp
    # print Uc
    # Uc[:] = [alpha*(1-beta)*x for x in Uc]
    Uc[:] = [alpha * (1 - beta) * x / nc for x in Uc]

    # Ud is D gain    
    Ud = [float(i) for i in D]
    # normalized (1-gamma) times Ud
    # Ud[:] = [(1-(x/max(Ud)))*(1-gamma) for x in Ud]
    Ud[:] = [(1 - (x / max(Ud))) for x in Ud]
    # Ud[:] = [(1-alpha)*x for x in Ud]
    Ud[:] = [(1 - alpha) * x / nc for x in Ud]
    # print Ud


    # U is beta times (Uc + Ud)
    # U = [beta*(i+j) for i, j in zip(Uc,Ud)]
    U = [(i + j) for i, j in zip(Uc, Ud)]

    # solve the optimization problem
    [C, m] = solveSPP_hotspot(n, G, U, V)

    # computation time
    #time_elapsed = (time.clock() - time_start)
    time_elapsed = (time.time() - time_start)
    

    # --------------------------------------------------------------------------
    # PRINT RESULTS
    # --------------------------------------------------------------------------

    num = 0
    akhriC = {}
    for v in m.getVars():
        # print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1

    print 'The optimal vale: %g' % (m.ObjVal)

    print 'Computation time: %g sec' % (time_elapsed)

    print C_x_GC(akhriC.values(), GC_x_U(G_x_C(G, akhriC.values()), U))
    # ERQ = C_x_GC(akhriC,GC_x_U(G_x_C(G,akhriC),U))
    # print 'ERQ calculated: %g sec' %( ERQ )



    # --------------------------------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------------------------------

    with open('{0}/selectedConf_hotspot{1}.dat'.\
              format(FilePath, hotspot_num), 'w') as f:
        for v in m.getVars():
            f.write(str(akhriC[v]))
            f.write(str("\n"))

    with open('{0}/computationTime_hotspot{1}.dat'.\
              format(FilePath,hotspot_num), 'w') as f:
        f.write(str(time_elapsed))

if strategy in 'two-step-exploration-icra16' or \
    strategy in 'two-step-exploration-icra16-newHc':

    print '**** plan type: ICRA16 ****'

    # --- path to read/write solutions    
    #FilePath = 'results/{0}/tomography_{1}/local-solutions'.\
    #            format(EXP_TITLE,plan_type)
    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/local-solutions").\
                 format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))
    
    print FilePath

    # Data from matlab files
    matfile = scipy.io.loadmat('{0}/preprocess_XVD_hotspot{1}.mat'.\
                               format(FilePath,hotspot_num))
    V = matfile['V']
    G = matfile['conf_crossAngles_G']
    D = matfile['D']

    G = G.tolist()
    V = V.tolist()

    # alpha times G
    for i in range(np.shape(G)[0]):
        for j in range(np.shape(G)[1]):
            G[i][j] = alpha * G[i][j]

    # Uc is V gain
    Uc = [sum(x) for x in zip(*V)]
    Uc = [float(i) for i in Uc]
    # gamma times normalized Uc        
    Uc[:] = [x / max(Uc) * beta for x in Uc]

    # Ud is D gain    
    Ud = [float(i) for i in D]
    # normalized (1-gamma) times Ud
    Ud[:] = [(1 - (x / max(Ud))) * (1 - beta) for x in Ud]

    # U is beta times (Uc + Ud)
    U = [(1 - alpha) * (i + j) for i, j in zip(Uc, Ud)]

    # solve the optimization problem
    [C,m] = solveSPP_hotspot(n,G,U,V)

    # computation time
    time_elapsed = (time.clock() - time_start)

    # --------------------------------------------------------------------------
    # PRINT RESULTS
    # --------------------------------------------------------------------------

    num = 0
    akhriC = {}
    for v in m.getVars():
        # print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1

    print 'The optimal vale: %g' % (m.ObjVal)

    print 'Computation time: %g sec' % (time_elapsed)

    # ERQ = C_x_GC(akhriC,GC_x_U(G_x_C(G,akhriC),U))
    # print 'ERQ calculated: %g sec' %( ERQ )



    # --------------------------------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------------------------------

    with open('{0}/selectedConf_hotspot{1}.dat'.\
              format(FilePath, hotspot_num), 'w') as f:
        for v in m.getVars():
            f.write(str(akhriC[v]))
            f.write(str("\n"))

            # with open('{0}/spp_computationtime_hotspot{1}.dat'.format(FilePath,hotspot_num), 'w') as f:
            #    f.write(str(time_elapsed))
