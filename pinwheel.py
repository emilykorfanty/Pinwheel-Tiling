# --------------------------------
# Pinwheel Substitution
#
# Created by: Emily Rose Korfanty
# Last updated: 2022-03-05
#---------------------------------

import os
os.environ['path'] += r';C:\msys64\mingw64\bin' # add cairo library to path
import numpy as np
import drawSvg as draw
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

def rot(alpha):
    """tile orientation"""
    return np.exp(complex(0,alpha))
        
def to_plane(V):
    """complex number to plane"""
    V = np.array(V,dtype=np.complex64)      # complex array
    W = V.view(np.float32)                  # real array
    W = W.reshape(V.shape + (2,))           # points on the plane
    return(W)

def get_tile(p):
    """control point to vertices"""

    x = p[0]        # distinguished point
    theta = p[1]    # orientation
    chi = p[2]      # chirality
        
    u = x + (theta + theta*rot(-chi*np.pi/2))               # right-angle vertex
    v = x + (theta*rot(np.pi) + theta*rot(-chi*np.pi/2))    # larger-angle vertex
    w = x + (theta*rot(np.pi) + 3*theta*rot(chi*np.pi/2))   # smaller-angle vertex
    
    return to_plane([u,v,w])
    
def get_vertices(Lambda):
    return [get_tile(p) for p in Lambda]
    
def fit_canvas(points,scale=10,ratio=1.5):

    # find the convex hull of the patch
    pts = np.concatenate(points)
    hull = np.array([[pts[vertex,0], pts[vertex,1]] for vertex in ConvexHull(pts).vertices])

    # find the smallest rectangle enclosing the patch
    xlen = hull[:,0].max() - hull[:,0].min()
    ylen = hull[:,1].max() - hull[:,1].min()
    canvas = ratio*scale*np.array([xlen, ylen])
    hull = scale*hull
    return([canvas,hull])
    
def draw_polygons(polygons,canvas,scale=10):
    
    # initialize the svg
    d = draw.Drawing(*canvas.tolist(), origin='center', displayInline=False)
    
    # inflate the polygons and canvas
    polygons = scale*polygons
    
    # draw each polygon on the canvas
    for polygon in polygons:
        d.append(draw.Lines(*polygon.flatten().tolist(),
                             close=True,
                             fill='none',
                             stroke='black'))
    return(d)
    
def F(p,j):
    """Pinwheel substitution for control points"""
    omega = -np.arctan(0.5)                 # rotation angle        
    x = p[0]                                # current distinguished point
    theta = p[1]                            # current orientation
    chi = p[2]                              # current chirality
        
    if(j==0):      
        x = np.sqrt(5)*rot(omega)*x
        theta = theta*rot(omega - chi*omega)
        chi = chi
              
    elif(j==1):            
        x = 2*theta*rot(omega - chi*omega + chi*np.pi/2) + np.sqrt(5)*rot(omega)*x
        theta = theta*rot(omega - chi*omega + np.pi)
        chi = chi
        
    elif(j==2):
        x = 4*theta*rot(omega - chi*omega + chi*np.pi/2) + np.sqrt(5)*rot(omega)*x
        theta = theta*rot(omega - chi*omega + np.pi)
        chi = -chi
        
    elif(j==3):    
        x = 2*theta*rot(omega - chi*omega + chi*np.pi) + np.sqrt(5)*rot(omega)*x
        theta = theta*rot(omega - chi*omega + np.pi)
        chi = -chi

    elif(j==4):
        x = 2*theta*rot(omega - chi*omega - chi*np.pi/2) + np.sqrt(5)*rot(omega)*x
        theta = theta*rot(omega - chi*omega - chi*np.pi/2)
        chi = -chi
        
    return [x,theta,chi]    # new control point
    
def get_patch(vertices, scale=10):
    """draw tiles on svg"""
    # create the drawing
    fit = fit_canvas(vertices)
    canvas = fit[0]
    hull = fit[1]
    d = draw.Drawing(*canvas.tolist(), origin='center', displayInline=False)

    # scale and center the patch
    ctr = np.array(Polygon(hull).centroid.coords).flatten()
    tiles = [scale*tile - ctr for tile in vertices]
        
    # draw each tile
    for tile in tiles:
        d.append(draw.Lines(*tile.flatten().tolist(),
                            close=True,
                            fill='none',
                            stroke='black'))
    return(d)

class pinwheel:
    """Pinwheel tiling control points"""  
    
    def __init__(self,x=0,theta=1,chi=1):
        Lambda = [[x,theta,chi]]
        vertices = get_vertices(Lambda)
        self.Lambda = Lambda
        self.vertices = vertices
        self.patch = get_patch(vertices)
        self.level = 0
        
    def sub(self, N=1):     
        
        n = 0        
        while(n < N):
        
            prev_Lambda = self.Lambda
            new_Lambda = []
            
            for p in prev_Lambda:
                for j in range(5):
                    new_Lambda.append(F(p,j))      
                    
            self.Lambda = new_Lambda
            self.vertices = get_vertices(self.Lambda)
            self.patch = get_patch(self.vertices)
            self.level = self.level+1
            
            n = n+1
            
    def draw(self):
        return(self.patch)
    
    def save(self, path):
        d = self.patch
        d.saveSvg(path)
        
        