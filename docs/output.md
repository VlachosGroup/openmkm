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


## Species data

1. **species.out**: Lists the species supplied in the input starting with species in gas phase, followed by species in bulk phase and surface phases.

2. **Hform.out**: Lists the dimensionless formatin enthalpies (H/RT) of the species at the input temperature. These are calculated from the NASA polynomials specified in the input. 

3. **Sform.out**: Lists the dimensionless formatin entropies (S/R) of the species. These are calculated from the NASA polynomials specified in the input. 

## Reaction data

1. **reactions.out**: Lists the reactions supplied in the input file.

2. **Hrxn.out**: Lists the dimensionless enthalpies (H/RT) of the reactions at the input temperature. These are calculated from the formatin enthalpies and the stoichometry coefficients of the species participating in the reactions. 

3. **Srxn.out**: Lists the dimensionless entropies (S/R) of the reactions at the input temperature. These are calculated from the formatin entropies and the stoichometry coefficients of the species participating in the reactions. 

4. **Grxn.out**: Lists the dimensionless Gibbs free energies (G/RT) of the reactions at the input temperature. These are calculated from the formatin enthalpies, entropies, and the stoichometry coefficients of the species participating in the reactions. 

5. **kc.out**: Lists the equilibrium constants, $K$,  of the reactions at the input temperature. These are calculated from the Gibbs free energies of the reactions.

6. **kf.out**: Lists the forward rate constants, $k_f$,  of the reactions at the input temperature. These are calculated from the activation energies, and pre-exponentials supplied. 

7. **kr.out**: Lists the reverse rate constants, $k_r$,  of the reactions at the input temperature. These are calculated from the forward and equilibrium rate constants. 

## Reactor State data

These files have either *_ss.out* or *_tr.out* . For clarity, the appendix is omitted.

1. **gas\_mole**: 
