# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 18:41:47 2015
@author: Asif Arain

Description: 
Sensor placement problem (SPP) for robot emission monitoring (REM):

G = ...
U = beta*sum(V,2)
max. C'(GC.*U)
s.t. |C| <= number and C in {0,1}

"""

from spp import solveSPP_replanning
#from spp import solveSPP_fusion
import time # time
import math # math
#from gurobipy import * # gurobi
import numpy as np # np
from numpy import linalg as LA # norm
import scipy.io # matlab files
#import statistics as stat

# initialize computation time
#time_start = time.clock()
time_start = time.time()


titlefile = scipy.io.loadmat('DATA_SPP_HOTSPOT_GEOMETRY_1t.mat')

EXP_TITLE = titlefile['EXP_TITLE']
EXP_TITLE = str(''.join(EXP_TITLE))

strategy = titlefile['strategy']
strategy = str(''.join(strategy))

FilePath = titlefile['FilePath']
FilePath = str(''.join(FilePath))

sensing_system = titlefile['sensing_system']
sensing_system = str(''.join(sensing_system))

hotspot_num = titlefile['hotspot_num']
hotspot_num = int(hotspot_num)

#conf_num = titlefile['conf_num']
#conf_num = int(conf_num)

print(' **********************************************')
print('                OPTIMIZATION:                  ')
print('   for optimal hotspot sensing geometry (hc#%g)' % (hotspot_num))
print(' **********************************************')

# alpha
alpha = titlefile['alpha']
alpha = float(alpha)

# beta
beta = titlefile['beta']
beta = float(beta)

# beta
gamma = titlefile['gamma']
gamma = float(gamma)

# number of allowed conf
n = titlefile['NumOfAllowedConf']
n = int(n)

# number of local conf
#nc = titlefile['NumOfLocalConf']
#nc = float(nc)
nc = float(n)

print('--- Number of confs to select: %g' %(n))

#if strategy == '1t-armex':
#    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/solutions/{5}").\
#             format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100),sensing_system)
#else:
#    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/replanning").\
#                 format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))
             
#print FilePath

matfile = scipy.io.loadmat('{0}/optimal_geometry_hc{1:02d}_xvd.mat'.\
                           format(FilePath,hotspot_num))

#V   = matfile['V']
V   = matfile['wV']
G   = matfile['conf_crossAngles_G']
D   = matfile['D']
executedConf = matfile['executedConfsHot_num']
Vdo = matfile['Conf_DistOrnGainV']

G = G.tolist()
V = V.tolist()
executedConf = executedConf.tolist()
#executedConf = executedConf.todense()
#executedConf = np.array(executedConf)
#executedConf = int(executedConf)

# alpha times G
for i in range(np.shape(G)[0]):
    for j in range(np.shape(G)[1]):
        G[i][j] = (alpha*beta)*(1/(nc*(nc-1)))*G[i][j]
        #G[i][j] = alpha*G[i][j]


Uc_w = [  sum(i) for i in zip(*V)]
Uc_w = [float(i) for i in Uc_w]
Vdo  = [float(j) for j in Vdo]
Uc   = [(gamma*i)+((1-gamma)*j) for i, j in zip(Uc_w,Vdo)]
          
# Uc is V gain
#Uc = [  sum(x) for x in zip(*V) ]
#Uc = [float(i) for i in Uc]
# gamma times normalized Uc        
#Uc[:] = [x/max(Uc)*gamma for x in Uc]
Uc[:] = [alpha*(1-beta)*x/nc for x in Uc]

# Ud is D gain    
Ud = [float(i) for i in D]
# normalized (1-gamma) times Ud
#Ud[:] = [(1-(x/max(Ud)))*(1-gamma) for x in Ud]
Ud[:] = [(1-(x/max(Ud))) for x in Ud]
Ud[:] = [(1-alpha)*x/nc for x in Ud]

# U is beta times (Uc + Ud)
#U = [beta*(i+j) for i, j in zip(Uc,Ud)]
U = [(i+j) for i, j in zip(Uc,Ud)]

#--------------------------------------------------------------------------
#                       SOLVE OPTIMIZATION PROBLEM
#--------------------------------------------------------------------------
     
# solve the optimization problem
[C,m] = solveSPP_replanning(n,G,U,V,executedConf)

# computation time
#time_elapsed = (time.clock() - time_start)
time_elapsed = (time.time() - time_start)

# --------------------------------------------------------------------------
#                           PRINT RESULTS
# --------------------------------------------------------------------------

num = 0
akhriC = {}
for v in m.getVars():
    #print('%s %g' % (v.varName, v.x))
    akhriC[num] = v.x
    num += 1

print('The optimal vale: %g' %(m.ObjVal))
print('Computation time: %g sec' %(time_elapsed))
    
# --------------------------------------------------------------------------
# SAVE RESULTS
# --------------------------------------------------------------------------

   
with open('{0}/optimal_geometry_hc{1:02d}_spp.dat'.\
          format(FilePath,hotspot_num), 'w') as f:
    for v in m.getVars():
        f.write(str(akhriC[v]))
        f.write(str("\n"))

#with open('{0}/spp_computationtime_hotspot{1}.dat'.format(FilePath,hotspot_num), 'w') as f:
#    f.write(str(time_elapsed)) 


