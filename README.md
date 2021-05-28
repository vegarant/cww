# CWW - Continuous Walsh sampling and wavelet reconstruction 
Code related to the paper "*Recovering wavelet coefficients from binary
samples using fast transforms*", by V. Antun.

This software implements fast matrix-vector multiplication with an NxM section of the infinite change-of-basis matrix between a Walsh sampling basis and an orthonormal wavelet reconstruction basis in one and two dimensions. It is accompanied by several examples. 

## Dependencies
* [SPGL1](http://www.cs.ubc.ca/~mpf/spgl1/). A solver for large-scale sparse reconstruction
* [Fastwht](https://bitbucket.org/vegarant/fastwht/). A fast implementation of MATLAB's `fwht`-function and Walsh functions.
* [WL - General purpose wavelet library](https://github.com/oyvindry/wl/tree/new_interface). Wavelet library for boundary wavelet support. Note that the branch `new_interface` is used.
* [CIlib - A software library for compressive imaging](https://github.com/vegarant/cilib). This library is used for one of the two-dimensional sampling patterns. 

## Get started 
Start by installing the dependencies listed above. In the `examples` folder, you will find a set of examples used to produce the figures in the paper. To run these, it is important to set all the default parameters used for plotting. This is done by running the script `etc/cww_set_detaults.m`. This script will produce a file called `var/cww_defaults.mat`, which is read by all of the scripts in the examples folder. Finally, add the following folders to your MATLAB path

```
addpath('/path/to/cww/var');
addpath('/path/to/cww/utils');
```


## Parallel computations
It is possible to speed up the computations of the two-dimensional algorithm by using MATLAB's [Parallel Computing Toolbox](https://se.mathworks.com/products/parallel-computing.html). To do so, replace write `parfor`, instead of `for`, in the script `utils/cww_handle_2d.m`. The relevant line numbers are 33, 39, 51 and 57. 



