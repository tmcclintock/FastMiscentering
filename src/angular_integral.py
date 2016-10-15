"""
This is a python wrapper that calls the calc_Sigma_angular() c function.
This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
import ctypes
from ctypes import c_double,c_int,POINTER,cdll

def calc_Sigma_angular(Rp,R,Sigma,N=1000):
    csalib = cdll.LoadLibrary("src/c_angular_integral.so")
    csa = csalib.calc_Sigma_angular
    csa.restype = c_int

    """
    Arguments are:
    Rp,R,Sigma,
    Sigma_angular,NR,N
    """
    NR = len(R)
    if NR != len(Sigma):
        raise Exception("len(R)!=len(Sigma)")
    csa.argtypes=[c_double,POINTER(c_double),POINTER(c_double),\
                  POINTER(c_double),c_int,c_int]
    R_in = R.ctypes.data_as(POINTER(c_double))
    S_in = Sigma.ctypes.data_as(POINTER(c_double))

    Sigma_angular = np.zeros(NR)
    Sa_in = Sigma_angular.ctypes.data_as(POINTER(c_double))

    result = csa(Rp,R_in,S_in,Sa_in,NR,N)

    if result != 0:
        raise Exception("Error message recieved in angular_integral.py")
    return Sigma_angular
