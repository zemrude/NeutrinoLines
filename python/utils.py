import numpy as np

def psi_f(RA,decl):
    return np.arccos(np.cos(np.pi/2.-(-29.*np.pi/180))*np.cos(np.pi/2.-decl)\
                      +np.sin(np.pi/2.-(-29.*np.pi/180))*np.sin(np.pi/2.-decl)*\
                       np.cos(RA-266.*np.pi/180))

