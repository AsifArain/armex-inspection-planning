# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 19:22:34 2015
@author: Asif Arain

Description: 

"""

#from spp import solveSPP_hotspot
#from spp import solveSPP_fusion
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


print '*************************************************************'
print '             Optimization for fusion'
print '*************************************************************'

#PairNumToHandle = 5;
# Pair number from matlab files
#pairfile = scipy.io.loadmat('PairNumToHandle.mat')
#PairNumToHandle = pairfile['PairNumToHandle']
#PairNumToHandle = int(PairNumToHandle)

titlefile = scipy.io.loadmat('SPP_FUSION_DATA.mat')

EXP_TITLE = titlefile['EXP_TITLE']
EXP_TITLE = str(''.join(EXP_TITLE))


strategy = titlefile['strategy']
strategy = str(''.join(strategy))


#print PairNumToHandle

#problem = 'spp hotspots'
#problem = 'spp fusion'



if strategy in 'two-step-exploration-jfr':
    
    print '**** plan type: JFR17 ****'
    
    # -- number of allowed conf
    n = titlefile['NumOfAllowedConf']
    n = int(n)
    
    # -- number of local conf
    nc = titlefile['NumOfLocalConf']
    nc = float(nc)
    
    # -- alpha
    alpha = titlefile['alpha']
    alpha = float(alpha)
    
    # -- beta
    beta = titlefile['beta']
    beta = float(beta)
    
    # -- gamma
    gamma = titlefile['gamma']
    gamma = float(gamma)
    
    # -- DisiredERQ1
    #ERQ_Desired1 = titlefile['ERQ_Desired1']
    #ERQ_Desired1 = float(ERQ_Desired1)
    
    # -- DisiredERQ2
    #ERQ_Desired2 = titlefile['ERQ_Desired2']
    #ERQ_Desired2 = float(ERQ_Desired2)
    
    #print alpha
    #print beta
        
    #for hotspot_num in range(1,2,1):
        
    # --- path to read/write solutions    
    #FilePath = 'results/{0}/tomography_jfr17/fusion'.format(EXP_TITLE)
    #FilePath = ("results/{0}/tomography_jfr17/"
    #            "alpha{1:03d}_beta{2:03d}_gamma{3:03d}/fusion").\
    #             format(EXP_TITLE,int(alpha*100),int(beta*100),int(gamma*100))
    FilePath = ("results/{0}/{1}/alpha{2:03d}_beta{3:03d}_gamma{4:03d}/fusion").\
                 format(EXP_TITLE,strategy,int(alpha*100),int(beta*100),int(gamma*100))
                     
    # Data from matlab files
    #matfile = scipy.io.loadmat('{0}/pairwise_fusion{1}.mat'.format(FilePath,PairNumToHandle))
    matfile = scipy.io.loadmat('{0}/preprocess_pairwise_fusion.mat'.format(FilePath))
    
    #V   = matfile['V']
    V1 = matfile['wV1']
    V2 = matfile['wV2']
    G1 = matfile['conf_crossAngles1_G']
    G2 = matfile['conf_crossAngles2_G']
    D1 = matfile['D1']
    D2 = matfile['D2']
    fixedConf1 = matfile['fixedConf1_num']
    fixedConf2 = matfile['fixedConf2_num']
    Vdo1 = matfile['Conf_DistOrnGainV1']
    Vdo2 = matfile['Conf_DistOrnGainV2']
    
    G1 = G1.tolist()
    G2 = G2.tolist()
    V1 = V1.tolist()
    V2 = V2.tolist()
    fixedConf1 = fixedConf1.tolist()
    fixedConf2 = fixedConf2.tolist()
    
    
    # (alpha*beta) times G
    for i in range(np.shape(G1)[0]):
        for j in range(np.shape(G1)[1]):
            G1[i][j] = (alpha*beta)*(1/(nc*(nc-1)))*G1[i][j]
    
    
    # (alpha*beta) times G
    for i in range(np.shape(G2)[0]):
        for j in range(np.shape(G2)[1]):
            G2[i][j] = (alpha*beta)*(1/(nc*(nc-1)))*G2[i][j]

    
    Uc_w1 = [  sum(i) for i in zip(*V1) ]
    Uc_w1 = [float(i) for i in Uc_w1]
             
    Uc_w2 = [  sum(i) for i in zip(*V2) ]
    Uc_w2 = [float(i) for i in Uc_w2]
            
    Vdo1 = [float(j) for j in Vdo1]
    Vdo2 = [float(j) for j in Vdo2]
    
    Uc1 = [(gamma*i)+((1-gamma)*j) for i, j in zip(Uc_w1,Vdo1)]
    Uc2 = [(gamma*i)+((1-gamma)*j) for i, j in zip(Uc_w2,Vdo2)] 
          
    # Uc1 is V1 gain
    #Uc1 = [  sum(x) for x in zip(*V1) ]
    #Uc1 = [float(i) for i in Uc1]
             
          
    # Uc2 is V2 gain
    #Uc2 = [  sum(x) for x in zip(*V2) ]
    #Uc2 = [float(i) for i in Uc2]
    
    # gamma times normalized Uc        
    #Uc1[:] = [x/max(Uc1) for x in Uc1]
    Uc1[:] = [alpha*(1-beta)*x/nc for x in Uc1]

    # gamma times normalized Uc        
    #Uc2[:] = [x/max(Uc2) for x in Uc2]
    Uc2[:] = [alpha*(1-beta)*x/nc for x in Uc2]
    
    # Ud is D gain    
    Ud1 = [float(i) for i in D1]

    # Ud is D gain    
    Ud2 = [float(i) for i in D2]
          
    # normalized (1-alpha) times Ud
    Ud1[:] = [(1-(x/max(Ud1))) for x in Ud1]
    Ud1[:] = [(1-alpha)*x/nc for x in Ud1]


    # normalized (1-alpha) times Ud
    Ud2[:] = [(1-(x/max(Ud2))) for x in Ud2]
    Ud2[:] = [(1-alpha)*x/nc for x in Ud2]
    
    # U is (Uc + Ud)
    U1 = [(i+j) for i, j in zip(Uc1,Ud1)]
         
    # U is (Uc + Ud)
    U2 = [(i+j) for i, j in zip(Uc2,Ud2)]
          
    # solve the optimization problem
    [C,m] = solveSPP_fusion(n,G1,G2,U1,U2,V1,V2,fixedConf1,fixedConf2)
    
    #[C,m] = solveSPP_fusion(n,XX1,XX2,VV1,VV2,TT1,TT2,fixedConf1,fixedConf2)

    
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
    
    print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G1,akhriC.values()),U1))
    print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G2,akhriC.values()),U2))
    
    # --------------------------------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------------------------------
    with open('{0}/selectedConf_PairwiseFusion.dat'.format(FilePath), 'w') as f:
        for v in m.getVars():
            f.write(str(akhriC[v]))
            f.write(str("\n"))
    


        
if strategy in 'two-step-exploration-icra16' \
    or strategy in 'two-step-exploration-icra16-newHc':
    
    print '**** plan type: ICRA16 ****'
    
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
    
    print 'The optimal vale: %g' %(m.ObjVal)
    print 'Computation time: %g sec' %(time_elapsed)
        
    # --------------------------------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------------------------------
             
    with open('{0}/selectedConf_PairwiseFusion{1}.dat'.\
              format(FilePath,PairNumToHandle), 'w') as f:
        for v in m.getVars():
            f.write(str(akhriC[v]))
            f.write(str("\n"))
        
    
