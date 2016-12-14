"""
An example of how to interface with 
fast_miscentering.

This example compares Sigma_mis from
the gauss-legendre method to
that using splines, called SM.
"""
import sys,time
import fast_miscentering as fm
import numpy as np

R = np.genfromtxt("../test_data/R.txt")
Sigma = np.genfromtxt("../test_data/sigma_r.txt")
Rm = 0.249697 #Mpc/h offset used in DeltaSigma code
SM = np.genfromtxt("../test_data/miscentered_sigma_r.txt")

start = time.time()
Sigma_mis = fm.calc_Sigma_miscentered(Rm,R,Sigma)
end = time.time()
print "Sigma_miscentered time:",end-start

print (Sigma_mis[:10]/SM[:10])

import matplotlib.pyplot as plt
plt.loglog(R,Sigma,label=r"$\Sigma_{\rm centered}(R)$")
plt.loglog(R,Sigma_mis,label=r"$\Sigma_{\rm Miscentered}(R)$")
plt.loglog(R,SM,'k--',label="old")
plt.legend()
plt.ylabel(R"$\Sigma(R)\ [h{\rm M_\odot/Mpc^2}]$",fontsize=24)
plt.xlabel(r"$R\ [h^{-1}{\rm Mpc}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.show()
