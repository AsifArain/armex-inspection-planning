# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 19:22:34 2015
@author: Asif Arain

Description: 

"""

from spp import *
import time # time
import math # math
from gurobipy import * # gurobi
import numpy as np # np
from numpy import linalg as LA # norm
import scipy.io # matlab files
#import statistics as stat

# initialize computation time
#time_start = time.clock()
time_start = time.time()


print(' **********************************************')
print('                OPTIMIZATION                   ')
print('     for the fusion of hotspot geometries      ')
print(' **********************************************')

titlefile = scipy.io.loadmat('DATA_SPP_FUSION_1t.mat')

EXP_TITLE = titlefile['EXP_TITLE']
EXP_TITLE = str(''.join(EXP_TITLE))

strategy = titlefile['strategy']
strategy = str(''.join(strategy))

sensing_system = titlefile['sensing_system']
sensing_system = str(''.join(sensing_system))

FilePath = titlefile['FilePath']
FilePath = str(''.join(FilePath))


if strategy in '1t-armex':
    
    print('----- plan type: 1t-armex')
    
    # -- number of allowed conf
    n = titlefile['NumOfAllowedConfFusion']
    n = int(n)
    
    
    hc_block_num = titlefile['hc_block_num']
    hc_block_num = int(hc_block_num)
    
    
    # -- number of local conf
    #nc = titlefile['NumOfAllowedConfLocal']
    #nc = nc.astype(float)
    
    # -- alpha
    alpha = titlefile['alpha']
    alpha = float(alpha)
    
    # -- beta
    beta = titlefile['beta']
    beta = float(beta)
    
    # -- gamma
    gamma = titlefile['gamma']
    gamma = float(gamma)
    
    G = titlefile['G']
    U = titlefile['U']

    # --- path to read/write solutions
    #FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/solutions/{5}").\
    #             format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100),sensing_system)
                     
    # Data from matlab files
    matfilename = ("{0}/pairwise_fusion_block{1:02d}_preprocess.mat").format(FilePath,hc_block_num)
    
    matfile = scipy.io.loadmat(matfilename)
    
    
    
    fixedConfs = matfile['fixedConfsAll_num']
    

    G = G.tolist()
    fixedConfs = fixedConfs.tolist()
    U = U.tolist()
    
    #print G
    #print V
    #print fixedConfs
    #print num_of_hc
    #print Uc
    #print Ud
    #print U
    #print n
    
    c = np.shape(G[0][0])[1]
    #c = np.shape(V[0][0])[1]
    #c = np.shape(V[0][1])[1]
    #c = np.shape(V[0][n])[1]
    
    #--------------------------------------------------------------------------
    # OPTIMIZATION
    #--------------------------------------------------------------------------
    # solve the optimization problem
    [C,m] = solveSPP_fusionTMP(n,c,G,U,fixedConfs)
        
    # computation time
    #time_elapsed = (time.clock() - time_start)
    time_elapsed = (time.time() - time_start)
    
    #--------------------------------------------------------------------------
    # PRINT RESULTS
    #--------------------------------------------------------------------------
    
    num = 0
    akhriC = {}
    for v in m.getVars():
        #print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1
    
    print('The optimal vale: %g' %(m.ObjVal))
    print('Computation time: %g sec' %(time_elapsed))
    
    #print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G[0][0],akhriC.values()),U[0][0]))
    #print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G[0][1],akhriC.values()),U[0][1]))
    
    # --------------------------------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------------------------------
    
    sppfilename = ("{0}/pairwise_fusion_block{1:02d}_spp.dat").format(FilePath,hc_block_num)
    # with open(sppfilename,'w') as f:
    #     for v in m.getVars():
    #         f.write(str(akhriC[v]))
    #         f.write(str("\n"))
    with open(sppfilename,'w') as f:
        #for v in len(m.getVars()):
        #for i in xrange(len(m.getVars())):
        for v in range(0, len(m.getVars())):
            #print(v)
            #print(akhriC[v])
            #f.write("some %g"%akhriC[v])
            f.write(str(akhriC[v]))
            f.write(("\n"))


        
if strategy in 'two-step-exploration-icra16' \
    or strategy in 'two-step-exploration-icra16-newHc':
    
    print('**** plan type: ICRA16 ****')
    
    PairNumToHandle = titlefile['PairNumToHandle']
    PairNumToHandle = int(PairNumToHandle)
    
    # -- alpha
    alpha = titlefile['alpha']
    alpha = float(alpha)
    
    # -- beta
    beta = titlefile['beta']
    beta = float(beta)
    
    # -- gamma
    gamma = titlefile['gamma']
    gamma = float(gamma)
    
    # compromise on the gain
    deltaGain = titlefile['deltaGain']
    deltaGain = float(deltaGain) #0.75
    
    # Pair number from matlab files
    #pairfile = scipy.io.loadmat('PairNumToHandle.mat')
    #PairNumToHandle = pairfile['PairNumToHandle']
        
    #print plan_type
    # --- path to read/write solutions
    #FilePath = 'results/{0}/tomography_{1}/fusion'.format(EXP_TITLE,plan_type)
    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/fusion").\
                 format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))
    
    #print FilePath
    
    
    # Data from matlab files
    matfile = scipy.io.loadmat('{0}/preprocess_pairwise_fusion{1}.mat'.\
                               format(FilePath,PairNumToHandle))
    Z = matfile['Z']
    Z = Z.tolist()
    
    infoGainSubH_sppH = matfile['infoGainSubH_sppH']
    infoGainSubH_sppH = infoGainSubH_sppH.tolist()
    
    constV = matfile['constV']
    constV = constV.tolist()
    
        
    # solve the optimization problem
    [C,m] = solveSPP_fusion_old(Z,infoGainSubH_sppH,deltaGain,constV)
    
    # computation time
    time_elapsed = (time.clock() - time_start)
    
    # --------------------------------------------------------------------------
    # PRINT RESULTS
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
             
    with open('{0}/selectedConf_PairwiseFusion{1}.dat'.\
              format(FilePath,PairNumToHandle), 'w') as f:
        for v in m.getVars():
            f.write(str(akhriC[v]))
            f.write(str("\n"))
        
    
