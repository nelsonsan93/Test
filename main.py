#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt
#import scipy.sparse as sps
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve
import matplotlib.tri as tri
import re
import pdb
import time as ti

class Mesh:
    def __init__(self,filename):
        self.filename = filename
        mesh = open(meshes[mesh_choosen],"r")
        mesh_properties = mesh.readlines(1)
        match = re.match(r"\s+(\d+)\s+(\d+)",mesh_properties[0])
        num_points = int(match.group(1))
        num_inner_points = int(match.group(2))
        #print('Number of total points: {}'.format(num_points))
        #print('Number of triangles: {}'.format(num_triangles))
        self.num_points = num_points
        self.num_inner_points = num_inner_points
        with open(meshes[mesh_choosen]) as f:
          content = f.readlines()
        properties = re.compile(r"\s+(\d+)\s+(\d+)\s*$")
        pattern_triangles = re.compile(r"\s+(\d+)\s+(\d+)\s+(\d+)")
        pattern_points = re.compile(r"\s+(-?\d\.\d+(E\-\d{3})?)\s+(-?\d\.\d+(E\-\d{3})?)\s+$")
        list_properties = list(filter(properties.match, content))
        list_points = list(filter(pattern_points.match, content))
            #print(list_points[44])
        mama = pattern_points.match(list_points[44])
        #print(mama.group(1))
        #print(mama.group(2))
        #print(mama.group(3))
        #print(mama.group(4))
        list_triangles = list(filter(pattern_triangles.match, content))
        triangles = np.zeros((len(list_triangles),3))
        self.num_triangles = len(triangles)
        i = 0
        for i_triangle in list_triangles:
            match = pattern_triangles.match(i_triangle)
            triangles[i,0]=int(match.group(1))
            triangles[i,1]=int(match.group(2))
            triangles[i,2]=int(match.group(3))
            i = i + 1
        points = np.zeros((len(list_points),2))
        i = 0
        for iCoordinate in list_points:
            match = pattern_points.match(iCoordinate)
            points[i,0]=float(match.group(1))
            points[i,1]=float(match.group(3))
            i = i + 1
        self.triangles = triangles
        self.points = points

    def print_properties(self):
        print('Number of grid points: {}'.format(self.num_points))
        print('Number of inner points: {}'.format(self.num_inner_points))
        print('Number of triangles: {}'.format(self.num_triangles))

    def plotTriangle(self,triangle_id):
        vertice1 = self.points[int(self.triangles[triangle_id,0]-1)]
        vertice2 = self.points[int(self.triangles[triangle_id,1]-1)]
        vertice3 = self.points[int(self.triangles[triangle_id,2]-1)]
        tri_x = np.array([float(vertice1[0]),float(vertice2[0]),float(vertice3[0]),float(vertice1[0])])
        tri_y = np.array([float(vertice1[1]), float(vertice2[1]), float(vertice3[1]), float(vertice1[1]) ])
        plt.plot(tri_x,tri_y)

    def plot_inner_points(self):
        for i_inner_point in range(self.num_inner_points):
            vertice = self.points[i_inner_point]
            plt.plot(vertice[0],vertice[1],'bo')
        for i_boundary_points in range(self.num_inner_points,self.num_points):
            vertice = self.points[i_boundary_points]
            plt.plot(vertice[0],vertice[1],'ro')
        plt.show()

    def plot(self):
        length = len(self.triangles)
        print('Plotting mesh...')
        for i in range(len(self.triangles)):
            print('Plotting triangle number {} of a total of {}...'.format(i+1,length))
            self.plotTriangle(i)
            #plt.pause(0.05)
        plt.show()
        print('Finished plotting.')

class Problem(Mesh):

    def __init__(self,filename):
        Mesh.__init__(self,filename)
        self.basis = np.array([[-1,-1],[1,0],[0,1]])
        print('Matrix A size: {} x {}'.format(self.num_inner_points,self.num_inner_points))
        self.A = sps.lil_matrix((self.num_inner_points,self.num_inner_points),dtype=np.float64)
        #self.A = np.zeros((self.num_inner_points,self.num_inner_points))
        self.f = np.zeros(self.num_inner_points)

    def get_vertices(self,triangle_id):
        self.vertice0 = self.points[int(self.triangles[triangle_id,0]-1)]
        self.vertice1 = self.points[int(self.triangles[triangle_id,1]-1)]
        self.vertice2 = self.points[int(self.triangles[triangle_id,2]-1)]
        self.vertices = np.stack((self.vertice0,self.vertice1,self.vertice2))

    def cal_B(self,triangle_id):
        self.get_vertices(triangle_id)
        self.B = np.array([[ self.vertice1[0] - self.vertice0[0], self.vertice2[0] - self.vertice0[0] ],\
            [ self.vertice1[1] - self.vertice0[1], self.vertice2[1] - self.vertice0[1]]])
        self.B_tilde = np.linalg.det(self.B) * np.linalg.inv(np.transpose(self.B))
        self.B_tilde_manual = (1/np.linalg.det(self.B)) *np.array([[self.vertice2[1]-self.vertice0[1],    -self.vertice1[1]+self.vertice0[1]],\
                [-self.vertice2[0]+self.vertice2[0],    self.vertice1[0]-self.vertice1[0]]])
        self.triangle_area_manual = 0.5 *( self.vertice0[0]* ( self.vertice1[1] -self.vertice2[1] ) + self.vertice1[0] * ( self.vertice2[1]  - self.vertice0[1] ) + self.vertice2[0] * ( self.vertice0[1] - self.vertice1[1] ) )
        self.triangle_area = abs(np.linalg.det(self.B))/2
        #print(self.B_tilde)
        #print(self.B_tilde_manual)
        #print('Triangle area: {}'.format(self.triangle_area))
        #print(self.triangle_area_manual)
        #diff = ( self.triangle_area - self.triangle_area_manual ) 
        #if ( diff > 1e-15 ):
        #    print(diff)

    def cal_numerical_solution(self):
        for k in range(self.num_triangles):
            self.cal_B(k)
            #if ( k == 0 ):
               # print(self.B_tilde)
            for l in range(3):
                i = int( self.triangles[k,l] - 1 )
                #print('Value of i={}, k={}, l={}'.format(i,k,l))
                #print('k={},i={},vertice={}'.format(k,i,self.points[i]))
                for m in range(3):
                    j = int( self.triangles[k,m] - 1 )
                    #print('Value of j: {}'.format(j))
                    if ( i < self.num_inner_points ) and ( j < self.num_inner_points ):
                        self.alpha1 = ( self.B_tilde[0,0]*self.basis[l,0] + self.B_tilde[0,1]*self.basis[l,1] ) * ( self.B_tilde[0,0]*self.basis[m,0] + self.B_tilde[0,1]*self.basis[m,1] ) 
                        self.alpha2 = ( self.B_tilde[1,0]*self.basis[l,0] + self.B_tilde[1,1]*self.basis[l,1] ) * ( self.B_tilde[1,0]*self.basis[m,0] + self.B_tilde[1,1]*self.basis[m,1] )
                        self.alpha = self.alpha1 + self.alpha2
                        self.A[i,j] = self.A[i,j] + self.alpha/( 4 * abs(self.triangle_area) )
                if ( i < self.num_inner_points ):
                    #print(i)
                    self.f[i] = self.f[i] + (self.triangle_area/3)*(2*np.pi**2*np.sin(np.pi*self.vertices[l,0])*np.sin(np.pi*self.vertices[l,1]))
        #print(self.A)
        #print(self.f)
        self.A = self.A.tocsr()
        self.numerical_solution = linsolve.spsolve(self.A,self.f)
        #self.numerical_solution = np.linalg.solve(self.A,self.f)
        print('Max: {}, min: {}'.format(np.max(self.f),np.min(self.f)))
        print('Max: {}, min: {}'.format(np.max(self.A),np.min(self.A)))
        #print('Max: {}, min: {}'.format(np.max(self.numerical_solution),np.min(self.numerical_solution)))
                                                                                
    def cal_analytical_solution(self):
        analytical_solution = np.zeros((self.num_inner_points,1))
        for i in range(self.num_inner_points):
            vertice = self.points[i]
            #print('x = {}, y = {}'.format(vertice[0],vertice[1]))
            analytical_solution[i] = np.sin(np.pi*vertice[0])*np.sin(np.pi*vertice[1])
        self.analytical_solution = analytical_solution

    def cal_error(self):
        error = 0
        for i_triangle in range(self.num_triangles):
            self.cal_B(i_triangle)
            for l in range(3):
                index = int( self.triangles[i_triangle,l] - 1 )
                if index < self.num_inner_points:
                    error = error + (self.triangle_area/3) * ( self.numerical_solution[index] - self.analytical_solution[index] ) ** 2
        error = np.sqrt(error)
        self.error = error

    def solve(self):
        self.cal_numerical_solution()
        self.cal_analytical_solution()
        self.cal_error()
        plt.plot(self.numerical_solution,label='Numerical')
        plt.plot(self.analytical_solution,label='Analytical')
        plt.legend()
        plt.show()
        print('Error: {}'.format(self.error))
    
mesh_choosen = 0
meshes = ['MESH.str.10','MESH.str.100','MESH.unstr.256','MESH.unstr.8']
pde = Problem(meshes[mesh_choosen])
pde.solve()
