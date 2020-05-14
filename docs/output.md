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

The data is printed into different output files to aid the user in sorting the
generated data. The output files have *.out* extension. In these files, the
data points corresponding to a single variable is presented in a single column
and the columns are separated by white space. Output files containing the
steady state data end with *_ss.out*. Similarly, output files containing the
transient state data end with *_tr.out*. Elementary data such as species
formation enthalpies and reaction rate constants do not contain *_ss* or *_tr*
before *.out*.


## Species data

1. **species.out**:
   - Lists the species supplied in the input starting with
species in gas phase, followed by species in bulk phase and surface phases.
   - Provides specie's name, phase, atomic composition, and fractional surface
coverage (if a surface species)

2. **Hform.out**:
   - List each species dimensionless formation enthalpies (H/RT)
   - Calculated from the NASA polynomials at the input temperature and 
input pressure.
   - **Note**: These values are calculated assuming no coverage effects.

3. **Sform.out**:
   - Lists each species dimensionless formation entropies (S/R).
   - Calculated from the NASA polynomials specified at the input temperature 
   and input pressure.
   - **Note**: These values are calculated assuming no coverage effects.

## Reaction data

1. **reactions.out**:
   - Lists the reactions supplied in the input file.

2. **Hrxn.out**:
   - Lists the dimensionless enthalpies (H/RT) of the reactions at
   the input temperature.
   - Calculated from the formation enthalpies listed in Hform.out and
   the stoichiometry coefficients of the species participating in the reactions. 

3. **Srxn.out**:
   - Lists the dimensionless entropies (S/R) of the reactions at
   the input temperature.
   - Calculated from the standard entropies listed in Sform.out and
   the stoichiometry coefficients of the species participating in the reactions. 

4. **Grxn.out**:
   - Lists the dimensionless Gibbs free energies (G/RT) of the
   reactions at the input temperature.
   - Calculated at the reaction inlet temperatures from the reaction enthalpies from Hrxn.out
   and the reaction entropies from Srxn.out. 

5. **kc.out**:
   - Lists the equilibrium constants, **K<sub>c</sub>**, of the reactions at the
   input temperature.
   - Calculated at the reactor outlet temperature and standard pressure of 1 bar for gaseous species.
   - Formulas:
      1. Adsorption reaction:
      2. Surface reaction:

6. **kf.out**:
   - Lists the forward rate constants, **k<sub>f</sub>**, of the reactions
   - Calculated at the reactor inlet temperature and input pressure or calculated at the
   reactor outlet temperature and pressure.
   - Formulas:
      1. Adsorption reaction:
      2. Surface reaction:

7. **kr.out**:
   - Lists the reverse rate constants, **k<sub>r</sub>**, of the reactions
   - Formula: k<sub>r</sub> = k<sub>f</sub>/K<sub>c</sub>

## Reactor State data

These file specified below have either *_ss.out* or *_tr.out* extensions. For
clarity, the extensions are omitted.

Files with _ss.out extension:
1. **Batch reactor**:
   - 1<sup>st</sup> column values are **BLANK**

2. **CSTR**:
   - 1<sup>st</sup> column values are **BLANK**

3. **PFR**:
   - 1<sup>st</sup> column values are volumes (units of m<sup>3</sup>) correseponding to the amount of
   volume through which the reaction has occurred at that point
   - Number of rows will equal the number of nodes for the calculation
   - Each volume can be considered an individual CSTR
   - The first change in volume is different from all the others because of **BLANK**
   - 1<sup>st</sup> column values are **BLANK**

Files with _tr.out extension:
1. **Batch reactor**:
   - 1<sup>st</sup> column values are **BLANK**

2. **CSTR**:
   - 1<sup>st</sup> column values are **BLANK**

3. **PFR, numerical**:
   - 1<sup>st</sup> column values are **BLANK**

4. **PFR, analytical**:
   - 1<sup>st</sup> column values are **BLANK**

Description of specific files referring to values in all but the 1<sup>st</sup> column

1. **gas\_mole**:
   - Mole fractions of gas species.

2. **gas\_mass**:
 - Mass fractions of gas species.

3. **gas_msdot**: 
   - Production rates (units of **BLANK**) of the gas species on the catalyst surface

4. **rctr\_state**:
   - Temperature (in K), pressure (in Pascals), density
   (in kg/m<sup>3</sup>) and either specific internal energy (units of J/kg) or specific enthalpy
   (units og **BLANK**) depending on the type of reactor.

5. **surf\_cov**:
   - Coverage fractions in the range of [0,1] of the surface species.

6. **rates**:
   - Forward, reverse, and net rates of progress in (units of mol/s) and partial
   equilibrium index of the reactions Specif reaction listed in rightmost column).
   - Partial equilibrium formula: (forward rate)/(forward rate + reverse rate)
   - Displayed quatities are the values at the reactor outlet

## Simulation data

These files contain messages showing the status of the simulation such as
execution time, and any warnings and error messages.

1. **general\_info.out**:
   - Recording of the information output to the console during the simulation run

2. **console or screen**:
   - Status of the simulation such as execution time, and
   any run time messages, warnings and error messages.
