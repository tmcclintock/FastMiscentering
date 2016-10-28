FastMiscentering
================
This is a fast implementation of the miscentering 
integrals needed to calculate the miscentered cluster column density.

Currently this contains a calculation of the angular integral
for Sigma(R) miscentered, which requires two integrals overall.
The angular integral can be done using Chebyshev-Gauss quadrature
and so can be done very quickly.

Dependencies
------------
* numpy
* GSL

You must have a path setup to gsl/include called GSLI
and a path to gsl/lib called GSLL.

Installation
------------
From the FastMiscentering directory, run
```
python setup.py install
```

And if you care about keeping the root directory clean
```
python setup.py clean
```

Usage
-------
Please look at examples/example.py for an example of how to run.

In order to calculate, \Sigma_{miscentered} the 
call signature is simply
```python
Sigma_mis = fast_miscentering.calc_Sigma_miscentered(Rmis,R,Sigma)
```
where Rmis is the miscentering length, R is an array of radial positions,
and Sigma is the centered Sigma(R).

Running the examples/example.py code produces the following:

![alt text](https://github.com/tmcclintock/FastMiscentering/blob/master/figures/fm_example.png)