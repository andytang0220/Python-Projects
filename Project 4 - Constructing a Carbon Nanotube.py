import numpy as np
import math
from atomplot import plot
from misc import include_atoms_gr, include_atoms_nt, grid_pq


def vector(n, m):
    """
    A function used to find the coordinates of vector Ch.
    Input:
        n: number of n hops (type: integer)
        m: number of m hops (type: integer)
    Output:
        Ch: A numpy array holding the coordinates of vector Ch
        i.e. Ch.shape = (2,)
    """
    vec = np.array([round(math.sqrt(3), 3) * n + round(math.sqrt(3) / 2, 3) * m, 0 - 1.5 * m])
    return vec


def normVector(vec):
    """
    A function used to find the coordinates of the normalized Ch vector.
    Input:
        Ch: A numpy array holding the coordinates of vector Ch.
    Output:
        norm_vec: A numpy array holding the coordinates of normalized vector.
        i.e. norm_vec.shape = (2,)
        norm: Length of vector Ch (type: float)
    """
    norm = round(math.sqrt(math.pow(vec[0],2) + math.pow(vec[1],2)), 3)
    norm_vec = np.array([round(vec[0] / norm, 3), round(vec[1] / norm, 3)])
    return norm_vec, norm


def normTvector(c_hat):
    """
    A function used to find the coordinates of vector t_hat, perpendicular
    to c_hat. Dot product of c_hat and t_hat = 0.
    Input:
        c_hat: A numpy array holding the coordinates of normalized Ch vector.
    Output:
        t_hat: A numpy array holding the coordinates of normalized
        perpendicular vector t_hat
        i.e. t_hat.shape = (2,)
    """
    t_hat = np.array([-1 * c_hat[1], c_hat[0]])
    return t_hat


def TVector(Ch, l):
    """
    A function used to find the coordinates of vector T, perpendicular
    to Ch.
    Input:
        Ch: A numpy array holding the coordinates of vector Ch.
        l: Length of the perpendicular edge in units of carbon-carbon bond
        lengths (type: float).
    Output:
        T: A numpy array holding the coordinates of perpendicular vector T.
        i.e. T.shape = (2,)
    """
    temp1 = normTvector(Ch)
    a = normVector(temp1)[1]
    temp2 = np.array([temp1[0] / a, temp1[1] / a])
    T = np.array([round(l*temp2[0], 3), round(l * temp2[1], 3)])
    return T


def pq(Ch, T):
    """
    A function used to find the maximum and minimum values of p and q.
    Input:
        Ch: A numpy array holding the coordinates of vector Ch.
        T: A numpy array holding the coordinates of perpendicular vector T.
    Output:
        p: A numpy array holding the minimum and maximum value for p, i.e.
        p = np.array([p_min, p_max]), p_min, p_max are integers.
        q: A numpy array holding the minimum and maximum value for q, i.e.
        q = np.array([q_min, q_max]), q_min, q_max are integers.
    """
    p_max = math.ceil((Ch[0] + T[0]) / (math.sqrt(3)/2))
    p_min = int(0)
    p = np.array([p_min, p_max])
    q = np.array([math.floor(Ch[1] / 1.5), math.ceil(T[1] / 1.5)])
    return p, q


def coordinates(pg, qg):
    """
    A function used to return the coordinates of each atom.
    Input:
        p: A 2-d numpy array holding the integers used identify each atom 
        along the  x-direction.
        q: A 2-d numpy array holding the integers used identify each atom 
        along the y-direction.
    Output:
        x: A 2-d numpy array holding the coordinates of each atom 
        along the  x-direction.
        y: A 2-d numpy array holding the coordinates of each atom  
        along the y-direction.
    """
    x_raw = [[] for list in pg]
    y_raw = [[] for list in qg]
  
    for a in range(len(pg)):
        for b in range(len(pg[a])):
            i = round(pg[a][b] * (math.sqrt(3) / 2), 3)
            x_raw[a].append(i)
        
    for j in range(len(qg)):
        for k in range(len(qg[j])):
            u = round(1.5*qg[j][k],3)
            if (pg[j][k] + qg[j][k]) % 2 != 0:
                u = u - (0.5*((qg[j][k] + pg[j][k]) % 2))
            u = round(u,3)
            y_raw[j].append(u)
                
    x = np.array(x_raw)
    y = np.array(y_raw)
    return x, y


def distance(x, y, c_hat):
    """
    A function used to calculate the distance along  Ch- and T- direction
    of each atom.
    Input:
        x: A 2-d numpy array holding the coordinates of each atom
        along the  x-direction.
        y: A 2-d numpy array holding the coordinates of each atom
        along the y-direction.
        c_hat: A numpy array holding the coordinates of normalized Ch vector.
    Output:
        s: A 2-d numpy array holding the distance of each atom
        along the  Ch-direction.
        t: A 2-d numpy array holding the distance of each atom
        along the T-direction.
    """
    s_raw = [[] for list in x]
    t_raw = [[] for list in y]
    
    for a in range(len(x)):
        for b in range(len(x[a])):
            t_dot = round((-1 * c_hat[1] * x[a][b]) + (c_hat[0] * y[a][b]), 3)
            t_raw[a].append(t_dot)
            s_dot = round(c_hat[0] * x[a][b] + c_hat[1] * y[a][b], 3)
            s_raw[a].append(s_dot)
    
    t = np.array(t_raw)
    s = np.array(s_raw)
    return s, t


def Graphene(n, m, l):
    # Graphene determines the x, y and z coordinates
    # for a graphene sheet.  This set of coordinates
    # represents a rectangular slice of a carbon sheet where
    # the lower left edge is represented by the vector
    # n*a1+m*a2, where a1 and a2 are vectors given in the
    # project description and n>=m are positive integers.
    # The integer len gives the length of the perpendicular edge
    # in units of carbon-carbon bond lengths.
    vec = vector(n, m)
    T = TVector(vec, l)
    p, q = pq(vec, T)
    pg, qg = grid_pq(p, q)
    x, y = coordinates(pg, qg)
    c_hat, arclen = normVector(vec)
    s, t = distance(x, y, c_hat)
    pos_gr = include_atoms_gr(x, y, s, t, arclen, l)
    tuberad = float(arclen / (2*(math.pi)))
    pos_nt = include_atoms_nt(pos_gr, c_hat, arclen, tuberad)
    return pos_gr, pos_nt


if __name__ == "__main__": 
    pos_gr, pos_nt = Graphene(6,6,6)
    plot(pos_gr)
    plot(pos_nt)
    

