#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 14:35:25 2020

@author: dimitrisgiovanis
"""
import numpy as np


def include_atoms_gr(x, y, s, t, arclen, l):
    """
    A function used to find which atoms to include in the graphene sheet.
    Input:
        x: coordinates of atoms in the x-axis
        y: corrdinates of atoms in the y-axis
        s: distance along the C-direction.
        t: distance along the T-direction.
        arclen: Length of the Ch vector.
        l: length of the perpendicular edge T.
    Output:
        pos: A numpy array containing the coordinates of the atoms
               included in the graphene sheet, (pos_nt.shape = (N, 2)).
        
    """
    x = np.array(x).T
    y = np.array(y).T
    s = np.array(s).T
    t = np.array(t).T
    
    tol=0.1;
    cond1 = s+tol>0
    cond2 = s-tol<arclen
    cond3 = t+tol>0
    cond4 = t-tol<l

    include = np.ones(shape=(s.shape[0], s.shape[1]), dtype=bool)
    for i in range(s.shape[0]):
        for j in range(s.shape[1]):
            if (cond1[i, j] == False) or (cond2[i, j] == False) or (cond3[i, j] == False) or (cond4[i, j] == False):
                include[i, j] = False
    
    return np.hstack([x[include].reshape(-1, 1), y[include].reshape(-1, 1)])
    
 
def include_atoms_nt(pos, c_hat, arclen, tubrad):
    """
    A function used to find which atoms to include in the nanotube
    Input:
        pos: A numpy array containing the coordinates of the atoms, calculated
             with the Graphene(n, m, l) function (pos.shape = (N, 2)).
        c_hat: Numpy array  holding the coordinates of normalized Ch vector.
        arclen: Norm of vector Ch (type: float).
        tubrad: Radius of the nanotube.
    Output:
          pos_nt: A numpy array containing the coordinates of the atoms
          included in the nanotube, (pos_nt.shape = (N, 3)).
    """
    
    tol=0.1;
    s = c_hat[0]*pos[:,0] + c_hat[1]*pos[:,1]
    t = -c_hat[1]*pos[:,0] + c_hat[0]*pos[:,1]
    
    tol=0.1;
    cond1 = s+tol>0
    cond2 = s+tol<arclen
    
    include = np.full((s.shape[0]), True, dtype=bool)
    
    for i in range(s.shape[0]):
        if (cond1[i] == False) or (cond2[i] == False):
            include[i] = False
    
    pos_ = [tubrad*np.cos(s[include]/tubrad), tubrad*np.sin(s[include]/tubrad),
            t[include]]
    pos_nt = np.vstack((pos_[0],pos_[1],pos_[2])).T
    
    return pos_nt



def grid_pq(p, q):
    """
    A function used to create the rectiliear grid of points given the 
    maximum and minimum values of p and q.
    Input:
        p: A numpy array holding the minimum and maximum value for p, i.e.
        p = np.array([p_min, p_max]), p_min, p_max are integers.
        q: A numpy array holding the minimum and maximum value for q, i.e.
        q = np.array([q_min, q_max]), q_min, q_max are integers.
    Output:
        pgrid: Numpy array holding the integers used identify each atom 
        along the  x-direction, i.e. pgrid.shape(2, N)
        qgrid: Numpy array holding the integers used identify each atom 
        along the y-direction, i.e. qgrid.shape(2, N)
    """
    P = np.arange(p[0], p[1] + 1, 1)
    Q = np.arange(q[0], q[1] + 1, 1)
    pgrid, qgrid = np.meshgrid(P,Q)
    return pgrid, qgrid