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

#n = 3     # number of allowed sensing configutions
#hotspot_num  = 5

#alpha = 0.75
#beta  = 0.25
#gamma = 0.50


titlefile = scipy.io.loadmat('SPP_REPLANNING_DATA.mat')

EXP_TITLE = titlefile['EXP_TITLE']
EXP_TITLE = str(''.join(EXP_TITLE))


strategy = titlefile['strategy']
strategy = str(''.join(strategy))


hotspot_num = titlefile['hotspot_num']
hotspot_num = int(hotspot_num)


conf_num = titlefile['conf_num']
conf_num = int(conf_num)




print '###########################################################################'
print ' Optimization for the replanning of hotspot# %g and conf# %g' %(hotspot_num,conf_num)
print '###########################################################################'



# alpha
alpha = titlefile['alpha']
alpha = float(alpha)

# beta
beta = titlefile['beta']
beta = float(beta)

# beta
gamma = titlefile['gamma']
gamma = float(gamma)

# number of local conf
nc = titlefile['NumOfLocalConf']
nc = float(nc)



# number of allowed conf
n = titlefile['NumOfAllowedConf']
n = int(n)

print '--- Number of allowed conf %g' %(n)



# --- path to read/write solutions
#FilePath = ("results/{0}/tomography_{1}/"
#            "alpha{2:03d}_beta{3:03d}_gamma{4:03d}/replanning").\
#             format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))


if strategy == 'one-step-exploration-jfr':
    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/solutions").\
             format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))
             
else:
    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/replanning").\
                 format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))
             
print FilePath

# Data from matlab files
#matfile = scipy.io.loadmat('{0}/VXD_ReEstimatedHotspot_Conf{1}_Hotspot{2}.mat'.\
#                           format(FilePath,conf_num,hotspot_num))


matfile = scipy.io.loadmat('{0}/tomo_adaptive_hc{1:02d}_conf{2:02d}_xvd.mat'.\
                           format(FilePath,hotspot_num,conf_num))

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

# solve the optimization problem
[C,m] = solveSPP_replanning(n,G,U,V,executedConf)

# computation time
#time_elapsed = (time.clock() - time_start)
time_elapsed = (time.time() - time_start)

# --------------------------------------------------------------------------
# PRINT RESULTS
# --------------------------------------------------------------------------

num = 0
akhriC = {}
for v in m.getVars():
    #print('%s %g' % (v.varName, v.x))
    akhriC[num] = v.x
    num += 1

print 'The optimal vale: %g' %(m.ObjVal)

print 'Computation time: %g sec' %(time_elapsed)
    
# --------------------------------------------------------------------------
# SAVE RESULTS
# --------------------------------------------------------------------------

   
#with open('{0}/selectedConf_Replanning_Conf{1}_Hotspot{2}.dat'.\
#          format(FilePath,conf_num,hotspot_num), 'w') as f:
#    for v in m.getVars():
#        f.write(str(akhriC[v]))
#        f.write(str("\n"))
with open('{0}/tomo_adaptive_hc{1:02d}_conf{2:02d}_spp.dat'.\
          format(FilePath,hotspot_num,conf_num), 'w') as f:
    for v in m.getVars():
        f.write(str(akhriC[v]))
        f.write(str("\n"))

#with open('{0}/spp_computationtime_hotspot{1}.dat'.format(FilePath,hotspot_num), 'w') as f:
#    f.write(str(time_elapsed)) 


