![](toc.png?raw=true "Estimating complete basis set extrapolation error through random walk")

<h1>Estimating complete basis set extrapolation errorthrough random walk</h1>

This repository contains a PYTHON implementation of a method of estimating the uncertainty of a result obtained through extrapolation to the complete basis set limit.
The method is based on an ensemble of random walkswhich simulate all possible extrapolation outcomes that could have been obtained if results from larger basis sets had been available.
The results assembled from a large collection of random walks can be then analyzed statistically, providing a route for uncertainty predictionat a confidence level required in a particular application.
The method is free of empirical parameters and compatible with any extrapolation scheme. The proposed technique has been tested in a series of numerical trials by comparing the determined confidence intervals with reliable reference data.

<h2>Usage</h2>

The script can be used as it is or <b>Estimator</b> class can be imported in other work.
If run as a "main" script it requires a path to a file containing the data to be extrapolated.<br>
An example of input file is:<br>
             &emsp;2    -0.0400183970<br>
             &emsp;3    -0.0411736630<br>
             &emsp;4    -0.0415978080<br>
             &emsp;5    -0.0417856800<br>
             &emsp;6    -0.0418812960<br>
             &emsp;7    -0.0419349210<br>
             &emsp;inf  -0.042044381<br>
with each line containing a pair (<i>l,e</i>), where <i>l</i> is a cardinal number of basis and <i>e</i> is a property to be extrapolated. If exact value is known, the <i>l</i> should be set up as "inf". 

As default, the Estimator class supports Helgaker's extrapolation $\infty \approx \frac{E_X\,X^3 - E_{X-1}(X-1)^3}{X^3 - (X-1)^3}$[<sup>1</sup>](https://pubs.aip.org/aip/jcp/article/106/23/9639/293509/Basis-set-convergence-of-correlated-calculations) and Lesiuk's extrapolation $E_\infty \approx E_X +  X^4(E_X - E_{X-1})(\zeta(4) - \sum_{l=1}^Xl^{-4}$[<sup>2</sup>](https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00705)

<b>As a script</b>

        Parameters
        ----------
        i : (str) The filepath to the file containing pairs of cardinal numbers and values to be extrapolated. If 'inf' string in the place of cardinal number, 
                        the value is taken as an exact value.
        CPU      : (int) Number of CPU cores to use for the procedure. Default is half of all available cores
        N_walk   : (int) Number of random walks to estimated the extrapolation and its uncertainty. Default is 1 000 000.
        N        : (int) Number of iteration for each random walk. Default is 300.
        thrs     : (float) Threshold for stopping the iteration for each random walk. Default is 1e-8.
    """ 

Run as 

````
python3 uncertainty_estimator.py --i extr_data/He_polar --CPU 10 --N_walk 10000000 --thrs 1e-9
````
<b>As an imported class</b>
````{python3}
from uncertainty_estimator import Estimator

estim = Estimator(CPU=10, verbose=2, N_walk=2000000, thrs=1e-9)
    
with open("extr_data/He_polar",'r') as fin:
    Lx = []
    Ex = []

    for line in fin:
        line = line.split()
        if line[0].upper() == "INF":
            Lx.append(line[0])
        else:
            Lx.append(int(line[0]))
        Ex.append(float(line[1]))
    v = [Ex, Lx]
    res = estim.estimate(v,k)
````
It is also possible to use a custom extrapolation method. The callable function has to take two lists (a list of cardinal numbers of a basis and properties to be extrapolated) and has to return a list of extrapolated values
````{python3}
def custom_extrapol(l,e,**kwargs):
    ...
    return e_x
    
with open("extr_data/He_polar",'r') as fin:
    Lx = []
    Ex = []

    for line in fin:
        line = line.split()
        if line[0].upper() == "INF":
            Lx.append(line[0])
        else:
            Lx.append(int(line[0]))
        Ex.append(float(line[1]))
    v = [Ex, Lx]
    res = estim.estimate(v, k, extr_method = custom_extrapol,**kwargs)
````

<b>DOI:</b>
         
