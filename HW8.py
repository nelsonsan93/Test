# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 15:48:34 2019

@author: toshiba
"""

import numpy as np
import csv

def mesh_selector(n):
    if n == 0:
        return  r"C:\Users\toshiba\Documents\Sim_Sciences_1_Semester\1_Sem_Materias\NMPDE\HW8\str10_test.csv"
    if n == 1:
        return r"C:\Users\toshiba\Documents\Sim_Sciences_1_Semester\1_Sem_Materias\NMPDE\HW8\unstr8_test.csv"

def f(x,y):
    return  np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)

def fsol(x,y):
    return  np.sin(np.pi*x)*np.sin(np.pi*y)

def matmat(m,n,mn,u,v):
    result = np.zeros((m,n))
    for i in range(0,m):
       for j in range(0,n):
           for k in range(0,mn):
               result[i][j] += u[i][k] * v[k][j]
    return result

def transpose(m,n,A):
    B = np.zeros((n,m))
    for j in range(m):
        for k in range(n):
            B[k][j] = A[j][k]
    return B

def chunkdet(a,b,c,d):
    return np.abs((a*d) - (c*b))
    

##Select the desired mesh
n = 0
filepath = mesh_selector(n)
#Place holders for the needed quantities
x,y = [], []  #Coordinates
v1,v2,v3 = [], [], [] #Connectivity
Nt, N = [] , [] #Total nodes and inner nodes
index = [] #Number of triangles
##Keep track in the file
count = 0
##After some modification to the initial data, extracting the info is easier now
for d in csv.DictReader(open(filepath), delimiter=','):
    if count == 0:
        Nt.append(int(d['x']))
        N.append(int(d['y']))
    if count > 0 and count < 1 + Nt[0]:
        x.append(float(d['x']))
        y.append(float(d['y']))
    if count == 1 + Nt[0]:
        index.append(int(d['x']))
    if count > 1 + Nt[0]:
        v1.append(int(d['x']))
        v2.append(int(d['y']))
        v3.append(int(d['z']))
    count+=1
##Convert into a integer for extra convenience
index = np.int(index[0])
##Number of triangles = 200 so to match the 0 indexing we need to substract 1
for j in range(0,index):
    v1[j] = v1[j] - 1
    v2[j] = v2[j] - 1
    v3[j] = v3[j] - 1
 
##Now with the values extracted, nodes and coordinates, we can continue
##Create the stifness matrix and rhs vector, p arrays and everything to assembly the system of equations    
STM = np.zeros([N[0],N[0]])
rhs = np.zeros([N[0]])
alpha_A = np.zeros((2,2))
alpha_B = np.zeros((2,2))
p = np.array([[-1,-1],[1,0],[0,1]])
p = transpose(3,2,p)
#Assembly the matrix
for k in range(0,index):
    npt = [v1[k], v2[k], v3[k]]
    BT = np.array([[x[npt[1]] - x[npt[0]], x[npt[2]] - x[npt[0]]],[y[npt[1]] - y[npt[0]], y[npt[2]] - y[npt[0]]]])
    ##Let's get the determinant of this bad boy
    ##Try the new detchunk function, selecting the chunks from the BT matrix
    detBT = chunkdet(BT[0,0],BT[0,1],BT[1,0],BT[1,1])
    ##Following the parallelogram law it's time to extract the areas of the elements.
    Trarea = detBT / 2.0
    ##Create the B matrix
    B = np.array([[y[npt[2]]- y[npt[0]],y[npt[0]]- y[npt[1]]],[x[npt[0]]- x[npt[2]],x[npt[1]]- x[npt[0]]]])
    ##Due to the definition of the problem itself, this matrices will always be the same size, alpha and B and BT.
    mm = 2
    nn = 3
    mn = 2
    alpha = matmat(mm,nn,mn,B,p)
    for l in range(0,3):
        i = npt[l]
        if i >= N[0]: ##Purge the boundary nodes
            continue
        for m in range(0,3):
            j = npt[m]
            if j >= N[0]: ##Purge the boundary nodes
                continue
            ##Create the alpha coefficients according to option no.3
            ##Select matrices
            alphaij = 0
            for t in range(0,2):
                alpha_A = alpha[t][l]
                alpha_B = alpha[t][m]
                alphaij += alpha_A* alpha_B
            STM[i][j] = STM[i][j] + alphaij / (4.0*Trarea)
        ##Once that the boundary nodes are out, there's a chance to calculate the right hand side integrals
        rhs[i] = rhs[i] + 2.0*(f(x[i],y[i])*(Trarea)*(1 / 3.0))
        
##Solve the system
Sol = np.linalg.solve(STM,rhs)
##Calculate the error
freal = np.zeros([N[0]])
for j in range(0,N[0]):
        freal[j] =  fsol(x[j],y[j])
norm = 0
error = np.abs(Sol-freal)
##Integrate the norm to obtain the l2 estimate
for k in range(0,index):
    npt = [v1[k], v2[k], v3[k]]
    BT = np.array([[x[npt[1]] - x[npt[0]], x[npt[2]] - x[npt[0]]],[y[npt[1]] - y[npt[0]], y[npt[2]] - y[npt[0]]]])
    detBT = chunkdet(BT[0,0],BT[0,1],BT[1,0],BT[1,1])
    Trarea = detBT / 2.0
    for l in range(0,3):
        i = npt[l]
        if i >= N[0]:
            continue
        norm = norm + Trarea * error[i]**2* (1/3.0)
##Collect the error in l2 norm
print(np.sqrt(norm))



