# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 20:41:24 2019

@author: toshiba
"""

import numpy as np
import csv

def mesh_selector(n):
    if n == 0:
        return  r"C:\Users\toshiba\Documents\Sim_Sciences_1_Semester\1_Sem_Materias\NMPDE\HW8\str10_test.csv"
    if n == 1:
        return r"C:\Users\toshiba\Documents\Sim_Sciences_1_Semester\1_Sem_Materias\NMPDE\HW8\MESH.unstr.8"
    
x,y = [], [] 
v1,v2,v3 = [], [], []
Nt, N = [] , []
index = []
count = 0
filepath = mesh_selector(0)
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
        ##index.append(float(d['x']))
        v1.append(int(d['x']))
        v2.append(int(d['y']))
        v3.append(int(d['z']))
    count+=1


##Assemblying the matrix working stuff
for l in range(0,3):
        i = npt[l]
        if i >= N[0]: ##Purge the boundary nodes
            continue
        
        for m in range(l,3):
            j = npt[m]
            if j >= N[0]: ##Purge the boundary nodes
                continue
            
            alpha_lm = alpha[:,l].dot(alpha[:,m])
            STM[i][j] = STM[i][j] + alpha_lm / (4.0*Trarea)
            if i is not j:
                STM[j][i] = STM[j][i] + alpha_lm / (4.0*Trarea)
        

        rhs[i] = rhs[i] + 2.0*(f(x[i],y[i])*(Trarea)*(1 / 3.0))
        
    




for t in range(0,2):
                alpha_A[][] = alpha[j][l]
                alpha_B[][] =alpha[j][m]
                

#Plotting the error                
x = np.linspace(0,1,N[0])
plt.plot(x,freal)
plt.plot(x,Sol)
plt.grid()
plt.show()



## Determinant calculated withouth the chunks
 # print("This is the ddet %f",ddetBT)
    #detBT = abs(np.linalg.det(BT))
   # print("This is the det %f",detBT)