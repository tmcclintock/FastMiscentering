"""
This is a python wrapper that calls the calc_Sigma_angular() c function.
This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
import ctypes, os, inspect
from ctypes import c_double,c_int,POINTER,cdll
sopath = os.path.join(os.path.dirname(__file__),"_fast_miscentering.so")
cfmlib = cdll.LoadLibrary(sopath)

def calc_Sigma_miscentered(Rm,R,Sigma,N=5):
    cfm = cfmlib.calc_Sigma_misc
    cfm.restype = c_int

    """
    Arguments are:
    Rp,R,Sigma,
    Sigma_misc,NR,N
    """
    NR = len(R)
    if NR != len(Sigma):
        raise Exception("len(R)!=len(Sigma)")
    cfm.argtypes=[c_double,POINTER(c_double),POINTER(c_double),\
                  POINTER(c_double),c_int,c_int]
    R_in = R.ctypes.data_as(POINTER(c_double))
    S_in = Sigma.ctypes.data_as(POINTER(c_double))

    Sigma_misc = np.zeros(NR)
    Sm_in = Sigma_misc.ctypes.data_as(POINTER(c_double))

    result = cfm(Rm,R_in,S_in,Sm_in,NR,N)

    if result != 0:
        raise Exception("Error message recieved in fast_misc.py")
    return Sigma_misc
