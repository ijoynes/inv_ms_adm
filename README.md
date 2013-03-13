# Inverse Micro-Scale Atmospheric Dispersion Model #

## Introduction ##

This document is a manual for use of the custom Inverse Mirco-Scale Atmospheric Dispersion Model developed for Ian M. Joynes's Masters thesis ("Inverse Micro-Scale Dispersion Modelling for Fugitive Emissions in the Upstream Oil and Gas Industry", 2013).

## Domain Setup ##

The dispersion model was developed to accommodate 2D unstructured linear triangular (3 node triangle) meshes.  Although, the code could be extended to support higher order triangular elements such as:  

- quadratic (6 nodes)  
- cubic (10 nodes)  
- quartic (15 nodes)  
- quintic (21 nodes)  

The code can also be extended to support the quadrilateral elements:

- quadratic (8 nodes)
- biquadric (9 nodes (1 internal node))
- cubic (12 nodes)
- bicubic (16 nodes (4 internal nodes))

Finally, the code can also be extended to support 3D domains with elements of various interpolation orders such as:

- Tetrahedral
- Hexahedral
- Triangular Prism



## Flow Setup ##

## Operator Computation ##



## Optimization ##

The source parameter estimation step is formed as an optimization problem to minimize an objective function based on the residual error between receptor observation and predicted receptor observations from candidate emission source profiles.  This objective function is of the form 

	f(x) = 1/2*(c(x)-c*)'*B^-1*(c(x)-c*) + theta/2*(s(x)-s*)'*R^-1*(s(x)-s*)

Minimization of this objective function is handled by the general purpose Quasi-Newton algorithm [L-BFGS-B distributed by Jorge Nocedal](http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html).  This code uses the version [Lbfgsb.2.1](http://users.eecs.northwestern.edu/~nocedal/Software/Lbfgsb.2.1.tar.gz).

## L-BFGS/L-BFGS-B Mex Wrapper ##
To call the Fortran minimization routines of L-BFGS or L-BFGS-B in Matlab, the original Fortran code must be compiled into a mex file.  [Liam Stewart](http://www.cs.toronto.edu/~liam/software.shtml) has provided the necessary helper functions to translate routines.f (the L-BFGS/-B routine) into a Matlab mex file.  The various versions routines.f can be found.

- [L-BFGS](http://www.netlib.org/opt/lbfgs_um.shar)
<!-- -[L-BFGS](http://www.netlib.org/opt/lbfgs_bcm.shar) -->
- [L-BFGS-B v2.1/v3.0](http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html)
- [L-BFGS-B v2.4](http://www.cs.toronto.edu/~liam/lbfgsb/routines.f)
