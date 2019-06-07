---
layout: default
---

# OpenMKM Outputs

Upon successful execution, OpenMKM prints data corresponding the reacting fluid
 state (density, temperature, pressure, mole fractions) and the surface 
coverages of surface species (if present). For batch and CSTR reactors, these 
are presented as a function of time. For PFR, the data is function of axial 
distance from the inlet. In addition to the above data, OpenMKM also prints 
data on the species such as formation enthalpies and on reactions such as rate 
constants, equilibrium constants, reaction enthalpies, etc. 

The data is 
printed into different output files to aid the user in sorting the generated 
data. The output files have *.out* extension. In these files, the data points 
corresponding to a single variable is presented in a single column and the 
columns are specarated by white space. Output files containing the steady state
data end with *_ss.out*. Similarly, output files containing the transient state
data end with *_tr.out*. Elementary data such as species formation enthalpies and
reaction rate constants do not contain *_ss* or *_tr* before *.out*.
