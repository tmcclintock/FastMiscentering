"""
An example of how to interface with 
fast_miscentering.

At the moment only the angular integral is
implemented.
"""
import sys,time
sys.path.insert(0,"src/")
import fast_misc
import numpy as np
import matplotlib.pyplot as plt

R = np.genfromtxt("test_data/R.txt")
Sigma = np.genfromtxt("test_data/sigma_r.txt")
MS = np.genfromtxt("test_data/miscentered_sigma_r.txt")
Rm = 0.249697 #Used in DeltaSigma code

start = time.time()
Sigma_mis = fast_misc.calc_Sigma_miscentered(Rm,R,Sigma)
end = time.time()
print "Sigma_miscentered time:",end-start

plt.loglog(R,Sigma)
plt.loglog(R,MS)
plt.loglog(R,Sigma_mis)
plt.show()
