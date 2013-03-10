# Inverse Micro-Scale Atmospheric Dispersion Model #

## Introduction ##

This document is a manual for use of the Inverse Mirco-Scale Atmospheric Dispersion Model.

## Domain Setup ##

## Flow Setup ##

## Operator Computation ##



## Optimization ##

The source parameter estimation step is formed as an optimization problem to minimize an objective function based on the residual error between receptor observation and predicted receptor observations from candidate emission source profiles.  This objective function is of the form 

	f(x) = 1/2*(c(x)-c*)'*B^-1*(c(x)-c*) + theta/2*(s(x)-s*)'*R^-1*(s(x)-s*)

Minimization of this objective function is handled by the general purpose Quasi-Newton algorithm L-BFGS-B distributed by Jorge Nocedal http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html.  This code uses the version Lbfgsb.2.1.