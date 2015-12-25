'''
This set of functions fills in holes in a data array. There must be data surrounding the 
hole for it to be filled in; if not, the holes are left in their present state.

LMK, UW, 10/10/2015
'''

import numpy as np

# -------------------------------------
def interior_missing_points(q, nodatavalue=float('NaN')):
    '''
    Return a mask of points that have nodatavalue data, but whose surroundings
    have data.
    
    Inputs:
    q: data array
    nodatavalue: nodatavalue
    
    Outputs:
    mask: mask of missing points
    '''
    
    if np.isnan(nodatavalue):
      mask = np.isnan(q)
    else:
      mask = q == nodatavalue
	
    ny, nx = np.shape(q)
    def neighbors(i, j):
        return [((i + di) % ny, (j + dj) % nx)
                for di in [-1, 0, 1] for dj in [-1, 0, 1]
                if not (di == 0 and dj == 0)]

    stack = [(ny - 1, nx - 1)]
    mask[ny - 1, nx - 1] = False
    while stack:
        i, j = stack.pop()
        for k, l in neighbors(i, j):
            if mask[k, l]:
                stack.append((k, l))
            mask[k, l] = False

    return mask


# -----------------------------------------------
def interpolate_missing_data(x, y, q, nodatavalue=float('NaN'), r=20):
    '''
    Fill in nodatavalue data in the interior of the domain, where it can be
    interpolated from surrouding points with data
    
    Inputs:
    x,y: grid points
    q: data array
    nodatavalue: no data value
    r: radius to pull points for interpolation
    
    Outputs:
    p: data array with filled in holes
    
    '''
    ny, nx = np.shape(q)
    dx = x[1] - x[0]

    d = int(r / dx)

    mask = interior_missing_points(q, nodatavalue)
    p = np.copy(q)

    for i in range(d, ny - d):
        for j in range(d, nx - d):
            if mask[i, j]:
                pij = 0.0
                weights = 0.0

                for k in range(i - d, i + d + 1):
                    for l in range(j - d, j + d + 1):
                        if not(np.isnan(q[k, l])):
                            w = d + 1 - max(abs(k - i), abs(l - j))
                            pij += w * q[k, l]
                            weights += w

                p[i, j] = pij / weights

    return p