import numpy as np
from scipy.sparse import coo_matrix

def LoadCoo(filename):
    rawMat = np.loadtxt(filename,
                        dtype=[("row", int), ("col", int), ("val", float)])
    cooMat = coo_matrix(
        (rawMat['val'], (rawMat['row'], rawMat['col']))
        )
    return cooMat

def LoadVec(filename):
    return np.loadtxt(filename)

