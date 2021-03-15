# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 03:19:28 2015
@author: Asif Arain

Description: 

Sensor placement problem
G = ...
U = beta*sum(V,2)
max. C'(GC.*U)
s.t. |C| <= number and C in {0,1}


"""


import math # math
from gurobipy import * # gurobi
import numpy as np # np
from numpy import linalg as LA # norm
# --------------------------------------------------------------------------
# LOCAL FUNCTIONS -- SPP HOTSPOT
# --------------------------------------------------------------------------

# G times C
def G_x_C(G,C):
    GC = {}
    #print C
    for i in range(np.shape(G)[0]):
        #GC[i] = quicksum( [G[i][j]*C[j] for j in range(np.shape(G)[1])] )
        GC[i] = quicksum( [G[i][j]*C[j] for j in range(len(C))] )
        #GC[i] =     sum( [G[i][j]*C[j] for j in range(len(C))] )
    return GC

## G times C
#def G_x_C(G,C):
#    GC = {}
#    for i in range(np.shape(G)[0]):
#        tmp_vec = {}
#        for j in range(np.shape(G)[1]):
#            tmp_vec[j] = G[i][j]*C[j]
#        #GC[i] = sum(tmp_vec)
#        GC[i] = quicksum(tmp_vec.values())
#    return GC

    
#def G_x_C(G,C):
#    GC = {}
#    print np.shape(C)
#    print np.shape(G)
#    for i in range(np.shape(G)[0]):
#        tmp_vec = {};
#        for j in range(np.shape(G)[1]):
#            tmp_vec[j] = G[i][j]*C[j]
#        nonzeros = [i for i, e in enumerate(tmp_vec) if e != 0]
#        print nonzeros
#        new_vec = {}
#        for k in range(len(nonzeros)):
#            new_vec[k] = tmp_vec[k]
#        GC[i] = np.min(new_vec.values())
#    print GC.values()
#    return GC.values()

# GC times U
def GC_x_U(GC,U):
    GCU = {};
    for i in range(len(GC)):
        #GCU[i] = GC[i]*U[i]
        GCU[i] = GC[i]+U[i]
    return GCU

    
# C times GC
def C_x_GC(C,GC):
    CGC = quicksum( [ C[i]*GC[i] for i in range(len(GC)) ] )
    #CGC =     sum( [ C[i]*GC[i] for i in range(len(GC)) ] )
    return CGC
   
    
# C times GCU
#-----------------------
def C_x_GCU(C,GCU):
    CGC = quicksum( [ C[i]*GCU[i] for i in range(len(GCU)) ] )
    return CGCU
    
#     

# FusedGain
#-----------------
def FusedCGC(C,GC):
    CGC = quicksum( [ C[i]*GC[i] for i in range(len(GC)) ] )
    return CGC
    
    
# Fused ERG (expected reconstruction gain)
#---------------------------------------------
def fFusedERG(C,G,U):
    fusedERG = quicksum( [(C_x_GC(C,GC_x_U(G_x_C(G[0][i],C),U[0][i]))) for i in range(np.shape(G[0])[0])])
    return fusedERG
    
    #-- Multiplication of ERGs does not work.
    # possible reason: "Gurobi only solves convex optimization problems and multiplying two variables is not necessarily convex (so that Gurobi does not allow it)!" at https://stackoverflow.com/questions/45579555/python-gurobi-pi-multiplication-for-binary-variables-invalid-argument-to-quade
    #----------------------------------------------------------------------------
    #fusedERG = float(0.00)
    #for i in range(np.shape(G[0])[0]):
    #    fusedERG = fusedERG + ( (C_x_GC(C,GC_x_U(G_x_C(G[0][i],C),U[0][i]))) )
    #    print fusedERG
    #return fusedERG


# Fused ERG (expected reconstruction gain)
#---------------------------------------------
def fFusedERG2(C,G,U):
    fusedERG = quicksum( \
                        [(C_x_GC(C,GC_x_U(G_x_C(G[0][i],C),U[0][i]))) \
                         for i in range(np.shape(G[0])[0])] \
                        )
    #-- Multiplication of ERGs does not work.
    # possible reason: "Gurobi only solves convex optimization problems and multiplying two variables is not necessarily convex (so that Gurobi does not allow it)!" at https://stackoverflow.com/questions/45579555/python-gurobi-pi-multiplication-for-binary-variables-invalid-argument-to-quade
    #----------------------------------------------------------------------------    
    #fusedERG = 1
    #for i in range(np.shape(G[0])[0]):
    #    fusedERG = fusedERG * (C_x_GC(C,GC_x_U(G_x_C(G[0][i],C),U[0][i])))
        
    return fusedERG




# variance of coverage    
def IIc_VisVar(V,C):
    VC = G_x_C(V,C)
    #print "--- hon VC"
    #print VC.values()
    negVC = np.dot(-1,VC.values())
    #print "------- manfi VC"
    #print negVC
    #print negVC.tolist()
    threeNegVC = np.power(3.1,float(negVC))
    #print threeNegVC
    OneMinusThreeNegVC = 1-threeNegVC
    #print negVC
    VarOneMinusThreeNegVC = np.var(OneMinusThreeNegVC)
    return 1-VarOneMinusThreeNegVC
    
# mean coverage
def IIc_VisMean(V,C):
    VC = G_x_C(V,C)
    #print "--- mean idher"
    #print VC.values()
    meanVC = np.mean(VC)
    meanVC3 = meanVC/3
    return meanVC3


# V times C
def V_x_C(V,C):
    VC = {}
    for i in range(np.shape(V)[0]):
        ij_sum = 0;
        for j in range(len(C)):
            ij_sum +=  V[i][j]*C[j]
        VC[i] = ij_sum
    return VC

# ObjMean
def objMean(V,C):
    VC = V_x_C(V,C)
    sumVC = 0
    for i in range(len(VC)):
        sumVC += VC[i]
    meanVC = sumVC/len(VC)
    #meanVC = float(sumVC)/len(VC)
    return meanVC
    
# Mean Coverage
def meanCov(V,C):
    VC = V_x_C(V,C)
    sumVC = 0
    for i in range(len(VC)):
        sumVC += VC[i]
    mean_coverage = sumVC/len(VC)
    #meanVC = float(sumVC)/len(VC)
    return mean_coverage

# Even Coverage
def evenCov(V,C):
    VC = V_x_C(V,C)
    mean_coverage = meanCov(V,C)    
    #even_coverage = quicksum(abs(VC-mean_coverage))    
    sumVC = 0
    for i in range(len(VC)):
        #sumVC += abs(VC[i]-mean_coverage)
        #sumVC += VC[i]-mean_coverage
        #sumVC += IloCplexModeler.abs(VC[i]-mean_coverage)
        sumVC += (VC[i]-mean_coverage) #*(VC[i]-mean_coverage)
    even_coverage = sumVC/len(VC)
    #meanVC = float(sumVC)/len(VC)
    return even_coverage


# --------------------------------------------------------------------------
# OPTIMIZATION - SPP REPLANNING
# --------------------------------------------------------------------------

def solveSPP_replanning(n,G,U,V,executedConf):
    
    c = np.shape(V)[1]
    
    # Create a new model
    m = Model("SPP-REM")
    
    #-- progress display
    m.Params.OutputFlag = 0

    
    # Create variables    
    conf_vector = range(c)
    C = {}    
    for conf_num in conf_vector:
        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)    
    
    
    # The objective is to maximize the gain
    m.modelSense = GRB.MAXIMIZE
        
    # Integrate new variables
    m.update()
    
    # --- OBJECTIVE FUNCTION
    m.setObjective( C_x_GC(C,GC_x_U(G_x_C(G,C),U)) )
    
    
    
    # --- CONSTRAINTS
    # -- visibility
    #for cell_num in range(np.shape(V)[0]):
    #    m.addConstr(
    #        quicksum(V[cell_num][conf_num] * C[conf_num] for conf_num in range(np.shape(V)[1])) >= 2,
    #        "visibility%d" % conf_num)
    
    # -- number of conf
    #m.addConstr( quicksum (C) <= 2, "conf_num")
    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= 2, "conf_num")
    m.addConstr( quicksum ([C[i] for i in range(c)]) <= n)
    #m.addConstr( C[0]+C[1]+C[2]+C[3]+C[4]+C[5] <= 2, "conf_num")
    #m.addConstr( quicksum(C) <= 1, "conf_num" )
    #for i in range(c):
    #    m.addConstr( C[i] >= 0 )
    #m.addConstr( C[i] >=0 for i in range(c))
    
    # -- fixed conf
    for x in range(np.shape(executedConf)[0]):
        num = executedConf[x][0]-1
        print('--- fixed conf #%g'%(num))
        m.addConstr( C[num] >= 1 )
        
    
    # optimize the model
    m.optimize()
    
    num = 0
    akhriC = {}
    for v in m.getVars():
        #print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1
        
       
    #print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G,akhriC.values()),U))
    
    return [C,m]


# --------------------------------------------------------------------------
# OPTIMIZATION - SPP HOTSPOT
# --------------------------------------------------------------------------

def solveSPP_hotspot(n,G,U,V):
    
    c = np.shape(V)[1]
    
    # Create a new model
    m = Model("SPP-RAGT")
    
    #-- progress display
    m.Params.OutputFlag = 0
    
    
    # Create variables    
    conf_vector = range(c)
    C = {}
    for conf_num in conf_vector:
        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)
    
    
    # The objective is to maximize the gain
    m.modelSense = GRB.MAXIMIZE
        
    # Integrate new variables
    m.update()
    
    # --- OBJECTIVE FUNCTION
    m.setObjective( C_x_GC(C,GC_x_U(G_x_C(G,C),U)) )
    
    
    # --- CONSTRAINTS    
    m.addConstr( quicksum ([C[i] for i in range(c)]) <= n)
        
    # optimize the model
    m.optimize()
    
    num = 0
    akhriC = {}
    for v in m.getVars():
        #print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1
        
        
    
    #print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G,akhriC.values()),U))
    
    return [C,m]


## --------------------------------------------------------------------------
## OPTIMIZATION - SPP HOTSPOT
## --------------------------------------------------------------------------
#
#def solveSPP_hotspot(alpha,beta,n,G,U,V):
#    
#    c = np.shape(V)[1]
#    
#    # Create a new model
#    m = Model("SPP-RAGT")
#    
#    # Create variables
#    #C = m.addVar(vtype=GRB.BINARY, name="Conf")
#    #conf_vector = range(c)
#    #conf = {}
#    #for conf_num in conf_vector:
#    #    conf[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
#    
#    conf_vector = range(c)
#    C = {}
#    #for conf_num in conf_vector:
#    #    C[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
#    for conf_num in conf_vector:
#        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)    
#    
#    #C = {}
#    #for conf_num in conf_vector:
#    #    C[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
#    #y = m.addVar(ub=1.0, name="y")
#    #z = m.addVar(ub=1.0, name="z")
#    
#    # The objective is to minimize the costs
#    m.modelSense = GRB.MAXIMIZE
#        
#    # Integrate new variables
#    m.update()
#    
#    # --- OBJECTIVE FUNCTION
#    
#    #m.setObjective( alpha*(C_x_GC(C,G_x_C(G,C))/quicksum(C)) )
#    
#    m.setObjective( C_x_GC(C,GC_x_U(G_x_C(G,C),U)) )
#    #m.setObjective( C_x_GC(C,G_x_C(G,C)) )
#    
#    
#    # --- CONSTRAINTS
#    #for cell_num in range(np.shape(V)[0]):
#    #    m.addConstr(
#    #        quicksum(V[cell_num][conf_num] * C[conf_num] for conf_num in range(np.shape(V)[1])) >= 2,
#    #        "visibility%d" % conf_num)
#    
#    
#    #m.addConstr( quicksum (C) <= 2, "conf_num")
#    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= 2, "conf_num")
#    m.addConstr( quicksum ([C[i] for i in range(c)]) <= n)
#    #m.addConstr( C[0]+C[1]+C[2]+C[3]+C[4]+C[5] <= 2, "conf_num")
#    #m.addConstr( quicksum(C) <= 1, "conf_num" )
#    #for i in range(c):
#    #    m.addConstr( C[i] >= 0 )
#    #m.addConstr( C[i] >=0 for i in range(c))
#    
#    # optimize the model
#    m.optimize()
#    
#    num = 0
#    akhriC = {}
#    for v in m.getVars():
#        #print('%s %g' % (v.varName, v.x))
#        akhriC[num] = v.x
#        num += 1
#        
#        
#    
#    print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G,akhriC.values()),U))
#    
#    return [C,m]


## --------------------------------------------------------------------------
## LOCAL FUNCTIONS - SPP FUSION
## --------------------------------------------------------------------------
#
## Z times C
#def Z_x_C(Z,C):
#    ZC = {}
#    #print C
#    for i in range(np.shape(Z)[0]):
#        #GC[i] = quicksum( [G[i][j]*C[j] for j in range(np.shape(G)[1])] )
#        ZC[i] = quicksum( [Z[i][j]*C[j] for j in range(len(C))] )
#    return ZC
#
## sum of ZC
#def sumZC(ZC):
#    #sZC = quicksum( [ZC[i] for i in range(ZC)] )    
#    sZC = 0;
#    for i in range(len(ZC)):
#        sZC +=  ZC[i]    
#    return sZC
#
#
#
## --------------------------------------------------------------------------
## OPTIMIZATION - SPP FUSION
## --------------------------------------------------------------------------
#
#def solveSPP_fusion(Z):
#    
#    c = np.shape(Z)[1]
#    
#    # Create a new model
#    m = Model("SPP-RAGT")
#    
#    conf_vector = range(c)
#    C = {}
#    #for conf_num in conf_vector:
#    #    C[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
#    for conf_num in conf_vector:
#        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)    
#    
#   
#    # The objective is to minimize the costs
#    m.modelSense = GRB.MAXIMIZE
#        
#    # Integrate new variables
#    m.update()
#    
#    # --- OBJECTIVE FUNCTION ---
#    
#    #m.setObjective( C_x_GC(C,GC_x_U(G_x_C(G,C),U)) )
#    m.setObjective( sumZC(Z_x_C(Z,C)) )
#    
#    #m.setObjective( C_x_GC(C,G_x_C(G,C)) )
#    
#    
#    # --- CONSTRAINTS ---
#    #for cell_num in range(np.shape(V)[0]):
#    #    m.addConstr(
#    #        quicksum(V[cell_num][conf_num] * C[conf_num] for conf_num in range(np.shape(V)[1])) >= 2,
#    #        "visibility%d" % conf_num)
#    
#    
#    #m.addConstr( quicksum (C) <= 2, "conf_num")
#    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= 2, "conf_num")
#    m.addConstr( quicksum ([C[i] for i in range(c)]) <= 1)
#    #m.addConstr( C[0]+C[1]+C[2]+C[3]+C[4]+C[5] <= 2, "conf_num")
#    #m.addConstr( quicksum(C) <= 1, "conf_num" )
#    #for i in range(c):
#    #    m.addConstr( C[i] >= 0 )
#    #m.addConstr( C[i] >=0 for i in range(c))
#    
#    # optimize the model
#    m.optimize()
#    
#    return [C,m]


# --------------------------------------------------------------------------
# LOCAL FUNCTIONS - SPP FUSION
# --------------------------------------------------------------------------

# Z times C
def Z_x_C(Z,C):
    ZC = {}
    #print C
    for i in range(np.shape(Z)[0]):
        #GC[i] = quicksum( [G[i][j]*C[j] for j in range(np.shape(G)[1])] )
        ZC[i] = quicksum( [Z[i][j]*C[j] for j in range(len(C))] )
    return ZC

# sum of ZC
def sumZC(ZC):
    #sZC = quicksum( [ZC[i] for i in range(ZC)] )    
    sZC = 0;
    for i in range(len(ZC)):
        sZC +=  ZC[i]    
    return sZC

# V times C
def V_x_C(V,C):
    VC = {}
    for i in range(np.shape(V)[0]):
        ij_sum = 0;
        for j in range(len(C)):
            ij_sum +=  V[i][j]*C[j]
        VC[i] = ij_sum
    return VC

    
# --------------------------------------------------------------------------
# OPTIMIZATION - SPP FUSION (NEW TMP)
# --------------------------------------------------------------------------

def solveSPP_fusionTMP(n,c,G,U,fixedConfs):
    
    #c = np.shape(V1)[1]
    
    # Create a new model
    m = Model("SPP-REM")
    
    #-- progress display
    m.Params.OutputFlag = 0

    conf_vector = range(c)
    C = {}
    for conf_num in conf_vector:
        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)
    
    
    # The objective is to maximize the gain
    m.modelSense = GRB.MAXIMIZE
        
    # Integrate new variables
    m.update()
    
    
    
    # --- OBJECTIVE FUNCTION
    #m.setObjective( (C_x_GC(C,GC_x_U(G_x_C(G1,C),U1))) + (C_x_GC(C,GC_x_U(G_x_C(G2,C),U2))) )
    m.setObjective( fFusedERG(C,G,U) )
    #m.setObjective(  (C_x_GC(C,GC_x_U(G_x_C(G2,C),U2))) )
    #m.setObjective( (C_x_GC(C,GC_x_U(G_x_C(G1,C),U1))) )
    
    
    # --- CONSTRAINTS
    m.addConstr( quicksum ([C[i] for i in range(c)]) <= n)
    
    # -- Fixed Conf all hotspot
    #for x in range(np.shape(fixedConfs)[0]):
    for x in range(np.shape(fixedConfs)[1]):
        #num = fixedConfs[x][0]-1
        num = fixedConfs[0][x]-1
        #print num
        m.addConstr( C[num] >= 1 )
        
        
        
    # -- 
    #m.addConstr( (C_x_GC(C,GC_x_U(G_x_C(G1,C),U1))) >= ERQ_Desired1)
    #m.addConstr( (C_x_GC(C,GC_x_U(G_x_C(G2,C),U2))) >= ERQ_Desired2)
        
 
    #print '..... optimization begins'
    # optimize the model
    m.optimize()
    #print '..... optimization finished'
    
    num = 0
    akhriC = {}
    for v in m.getVars():
        #print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1
 
        
    
    #print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G[0][0],akhriC.values()),U[0][0]))
    #print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G[0][1],akhriC.values()),U[0][1]))
    
    
    
    return [C,m]
    
    
# --------------------------------------------------------------------------
# OPTIMIZATION - SPP FUSION
# --------------------------------------------------------------------------

def solveSPP_fusion(n,G1,G2,U1,U2,V1,V2,fixedConf1,fixedConf2):
    
    c = np.shape(V1)[1]
    
    # Create a new model
    m = Model("SPP-REM")
    

    conf_vector = range(c)
    C = {}
    for conf_num in conf_vector:
        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)
    
    
    # The objective is to minimize the costs
    m.modelSense = GRB.MAXIMIZE
        
    # Integrate new variables
    m.update()
    
    # --- OBJECTIVE FUNCTION
    m.setObjective( (C_x_GC(C,GC_x_U(G_x_C(G1,C),U1))) + (C_x_GC(C,GC_x_U(G_x_C(G2,C),U2))) )
    #m.setObjective(  (C_x_GC(C,GC_x_U(G_x_C(G2,C),U2))) )
    #m.setObjective( (C_x_GC(C,GC_x_U(G_x_C(G1,C),U1))) )
    
    
    # --- CONSTRAINTS
    m.addConstr( quicksum ([C[i] for i in range(c)]) <= n)
    
    # -- Fixed Conf hotspot 1
    for x in range(np.shape(fixedConf1)[0]):
        num1 = fixedConf1[x][0]-1
        #print num1
        m.addConstr( C[num1] >= 1 )
        
    # -- Fixed Conf hotspot 2
    for x in range(np.shape(fixedConf2)[0]):
        num2 = fixedConf2[x][0]-1
        #print num2
        m.addConstr( C[num2] >= 1 )
        
    # -- 
    #m.addConstr( (C_x_GC(C,GC_x_U(G_x_C(G1,C),U1))) >= ERQ_Desired1)
    #m.addConstr( (C_x_GC(C,GC_x_U(G_x_C(G2,C),U2))) >= ERQ_Desired2)
        
 
    # optimize the model
    m.optimize()
    
    num = 0
    akhriC = {}
    for v in m.getVars():
        #print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1
 
        
    print(C_x_GC(akhriC.values(),GC_x_U(G_x_C(G1,akhriC.values()),U1)))
    print(C_x_GC(akhriC.values(),GC_x_U(G_x_C(G2,akhriC.values()),U2)))
    
    
    # --- test
    #C1 = akhriC
    #for x in range(np.shape(fixedConf2)[0]):
    #    num1 = fixedConf2[x][0]-1
    #    print num1
    #    C1[num1] = 0
    
    # --- test    
    #for x in range(np.shape(fixedConf2)[0]):
    #    num1 = fixedConf2[x][0]-1
    #    print num1
    #    akhriC[num1] = 0    

    # --- test
    #totalC = 0;
    #for v in m.getVars():
    #    totalC = totalC + akhriC[v]
    #print totalC
        
    #print C1
    #print akhriC

    # --- test
    #C2 = akhriC
    #for x in range(np.shape(fixedConf1)[0]):
    #    num2 = fixedConf1[x][0]-1
    #    print num2
    #    C2[num2] = 0

    #print C_x_GC(C1.values(),GC_x_U(G_x_C(G1,C1.values()),U1)) 
    #print C_x_GC(C2.values(),GC_x_U(G_x_C(G2,C2.values()),U2)) 
    
    
    print(C_x_GC(akhriC.values(),GC_x_U(G_x_C(G1,akhriC.values()),U1)))
    print(C_x_GC(akhriC.values(),GC_x_U(G_x_C(G2,akhriC.values()),U2)))
    
    return [C,m]


# --------------------------------------------------------------------------
# OPTIMIZATION - SPP FUSION
# --------------------------------------------------------------------------

def solveSPP_fusion_old(Z,infoGainSubH_sppH,deltaGain,constV):
    
    c = np.shape(Z)[1]
    
    # Create a new model
    m = Model("SPP-RAGT")
    
    conf_vector = range(c)
    C = {}
    #for conf_num in conf_vector:
    #    C[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
    for conf_num in conf_vector:
        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)    
    
   
    # The objective is to minimize the costs
    m.modelSense = GRB.MINIMIZE
        
    # Integrate new variables
    m.update()
    
    # --- OBJECTIVE FUNCTION ---
    
    #m.setObjective( C_x_GC(C,GC_x_U(G_x_C(G,C),U)) )
    #m.setObjective( sumZC(Z_x_C(Z,C)) )
    m.setObjective( quicksum ([C[i] for i in range(c)]) )
    
    #m.setObjective( C_x_GC(C,G_x_C(G,C)) )
    
    
    # --- CONSTRAINTS ---
    for cell_num in range(len(Z_x_C(Z,C))):
        m.addConstr( Z_x_C(Z,C)[cell_num] >= infoGainSubH_sppH[cell_num][0]-deltaGain )
        
            #quicksum(V[cell_num][conf_num] * C[conf_num] for conf_num in range(np.shape(V)[1])) >= 2,
            #"visibility%d" % conf_num)
    
    
    for cell_num in range(np.shape(constV)[0]):
            m.addConstr(
                quicksum(constV[cell_num][conf_num] * C[conf_num] for conf_num in range(np.shape(constV)[1])) >= 1,"visibility%d" % conf_num)    
    
    # distance constraint
    #m.addConstr( C_x_GC(C,G_x_C(mDg,C)) >= 0.2, "dist_const")
    
    
    #m.addConstr( quicksum (C) <= 2, "conf_num")
    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= 2, "conf_num")
    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= 1)
    
    
    
    #m.addConstr( C[0]+C[1]+C[2]+C[3]+C[4]+C[5] <= 2, "conf_num")
    #m.addConstr( quicksum(C) <= 1, "conf_num" )
    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= 1)
    #for i in range(c):
    #    m.addConstr( C[i] >= 0 )
    #m.addConstr( C[i] >=0 for i in range(c))
    
    # optimize the model
    m.optimize()
    
    num = 0
    akhriC = {}
    for v in m.getVars():
        #print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1
    
    #print C_x_GC(akhriC.values(),GC_x_U(G_x_C(G,akhriC.values()),U))
    print( Z_x_C(Z,akhriC.values()))
    print( infoGainSubH_sppH)
    
    return [C,m]


# --------------------------------------------------------------------------
# OPTIMIZATION - SPP REPLANNING
# --------------------------------------------------------------------------

def solveSPP_replanning_tmp(G,U,V,executedConf):
    
    c = np.shape(V)[1]
    
    # Create a new model
    m = Model("SPP-REM")
    
    # Create variables
    #C = m.addVar(vtype=GRB.BINARY, name="Conf")
    #conf_vector = range(c)
    #conf = {}
    #for conf_num in conf_vector:
    #    conf[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
    
    conf_vector = range(c)
    C = {}
    #for conf_num in conf_vector:
    #    C[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
    for conf_num in conf_vector:
        C[conf_num] = m.addVar(vtype=GRB.BINARY,lb=0.0,ub=1.0)    
    
    #C = {}
    #for conf_num in conf_vector:
    #    C[conf_num] = m.addVar(vtype=GRB.BINARY, name="conf_num%d" % conf_num)
    #y = m.addVar(ub=1.0, name="y")
    #z = m.addVar(ub=1.0, name="z")
    
    # The objective is to minimize the costs
    #m.modelSense = GRB.MAXIMIZE
    m.modelSense = GRB.MINIMIZE
        
    # Integrate new variables
    m.update()
    
    # --- OBJECTIVE FUNCTION
       
    #m.setObjective( alpha*(C_x_GC(C,G_x_C(G,C))/quicksum(C)) )
    
    #m.setObjective( C_x_GC(C,GC_x_U(G_x_C(G,C),U)) )
    #m.setObjective( C_x_GC(C,G_x_C(G,C)) )
    m.setObjective( quicksum ([C[i] for i in range(c)]) )
    
    # --- CONSTRAINTS
    #for cell_num in range(np.shape(V)[0]):
    #    m.addConstr(
    #        quicksum(V[cell_num][conf_num] * C[conf_num] for conf_num in range(np.shape(V)[1])) >= 2,
    #        "visibility%d" % conf_num)
    
    
    #m.addConstr( quicksum (C) <= 2, "conf_num")
    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= 2, "conf_num")
    #m.addConstr( quicksum ([C[i] for i in range(c)]) <= n)
        
    for x in range(np.shape(executedConf)[0]):
        num = executedConf[x][0]-1
        print('fixed conf #%g'%(num))
        m.addConstr( C[num] >= 1 )
    
    
    m.addConstr( (C_x_GC(C,GC_x_U(G_x_C(G,C),U))) >= 0.9 )
       

    # optimize the model
    m.optimize()
    
    num = 0
    akhriC = {}
    for v in m.getVars():
        #print('%s %g' % (v.varName, v.x))
        akhriC[num] = v.x
        num += 1
        
        
    
    print( C_x_GC(akhriC.values(),GC_x_U(G_x_C(G,akhriC.values()),U)))
    
    return [C,m]



