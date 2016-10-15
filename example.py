"""
An example of how to interface with 
fast_miscentering.

At the moment only the angular integral is
implemented.
"""
import sys,time
sys.path.insert(0,"src/")
import angular_integral
import numpy as np
import matplotlib.pyplot as plt

R = np.genfromtxt("test_data/R.txt")
Sigma = np.genfromtxt("test_data/sigma_r.txt")
MS = np.genfromtxt("test_data/miscentered_sigma_r.txt")
Rp = 0.249697 #Used in DeltaSigma code

start = time.time()
Sigma_angular = angular_integral.calc_Sigma_angular(Rp,R,Sigma)
end = time.time()
print "My time:",end-start

plt.loglog(R,Sigma)
plt.loglog(R,MS)
plt.loglog(R,Sigma_angular)
plt.show()
