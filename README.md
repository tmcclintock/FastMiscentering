# FastMiscentering
This is a fast implementation of the miscentering 
integrals needed to calculate the miscentered cluster column density.

Currently this contains a calculation of the angular integral
for Sigma(R) miscentered, which requires two integrals overall.
The angular integral can be done using Chebyshev-Gauss quadrature
and so can be done very quickly.