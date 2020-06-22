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
generated data. The output files have either of *.out*, *.dat*, and *.csv* extensions. 

In the *.out* files, the data points corresponding to a single variable is 
presented in a single column and the columns are separated by white space.

## Species data

1. **species.out**:
   - Lists the species supplied in the input starting with
species in gas phase, followed by species in bulk phase and surface phases.
   - Provides specie's name, phase, atomic composition, and fractional surface
coverage (if a surface species)
   - This file acts as input to RenView.

2. **Hform.out**:
   - List each species dimensionless formation enthalpies (H/RT)
   - Calculated from the NASA polynomials at the input temperature and 
input pressure.
   - **Note**: For surface species, these values are calculated assuming no coverage effects.

3. **Sform.out**:
   - Lists each species dimensionless formation entropies (S/R).
   - Calculated from the NASA polynomials specified at the input temperature 
   and input pressure.

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
   - Calculated at the reactor inlet temperature and input pressure.
   - Formulas:
      1. Adsorption reaction:
      2. Surface reaction:

7. **kr.out**:
   - Lists the reverse rate constants, **k<sub>r</sub>**, of the reactions
   - Formula: k<sub>r</sub> = k<sub>f</sub>/K<sub>c</sub>

## Reactor State data

The files specified below have either a *.csv* or *.dat* extension depending on the output file format selected in *reactor.yaml*. Further they may be qualified with *_ss.* or *_tr.* suffix indicating steady state or transient state respecitively. For clarity, the extensions are omitted at some places, but are implied.

Files with _ss.csv extension:
1. **Batch reactor**:
   - 1<sup>st</sup> column values are time (units of seconds).

2. **CSTR**:
   - 1<sup>st</sup> column values are time (units of seconds).

3. **PFR, numerical**:
   - 1<sup>st</sup> column values are volumes (units of m<sup>3</sup>) correseponding to the amount of
   volume through which the reaction has occurred at that point
   - Number of rows will equal the number of nodes for the calculation
   - Each volume can be considered an individual CSTR
   - The first change in volume is different from all the others because half of volume is used (to represent center of CSTR).

4. **PFR, analytical**:
   - 1<sup>st</sup> column values are length (units of m) of
   the PFR at which the values are sampled


Files with _tr.csv extension:
1. **Batch reactor**:
   - 1<sup>st</sup> column values are time (units of seconds).

2. **CSTR**:
   - 1<sup>st</sup> column values are time (units of seconds).

3. **PFR, numerical**:
   - 1<sup>st</sup> column values are time (units of seconds).

Description of specific files referring to values in all but the 1<sup>st</sup> column

1. **gas\_mole**:
   - Mole fractions of gas species.

2. **gas\_mass**:
   - Mass fractions of gas species.

3. **gas\_sdot**: 
   - Production rates (units of kmol/s) of the gas species on the catalyst surface. Due to the possibility of multiple surfaces, per unit surface area is not used.

4. **rctr\_state**:
   - Temperature (in K), pressure (in Pascals), density
   (in kg/m<sup>3</sup>), mass (in kg), volume (in m<sup>3</sup>), either specific internal energy (units of J/kg) or specific enthalpy
   (units of J/kg) depending on the type of reactor, mass flow rate into the reactor (in kg/s), mass flow rate out of the reactor (kg/s), mass produced from the surface (in kg/s), net rate of change of mass (kg/s) in the reactor.
   - **Note**: Mass flow rate out of the reactor is printed as -ve value. Similarly, if mass produced from surface is -ve, it means surface absorbs more mass than the mass released.

5. **surf\_cov**:
   - Coverage fractions in the range of [0,1] of the surface species.

5. **surf\_sdot**:
   - Production rates (units of kmol/m<sup>2</sup>/s) of the surface species on the catalyst surface.

6. **rates\_ss**:
   - Forward, reverse, and net rates of progress in (units of kmol/s) and partial
   equilibrium index of the specific reaction listed in rightmost column).
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

## Sensitivity Analysis Data

How to run a sensitivity analysis:

1. Run and analyze the results of a FIM run (typically for large mechanisms):
   - Used as a quick, first look to identify the important reactions to be examined in the slower, LSA run
   - In output file containing “_sensitivity” in the filename, there are two columns: rxnid (the reaction) and FIM_Diag
   (the corresponding FIM value). Identify the reactions with the largest FIM-Diag values (user defined as to what is
    “large”), and consider only these reactions for the LSA run

2. Run and analyze the results of a LSA run:
   - Include reactions identified in the FIM run
   - In the output file containing “_sensitivity” in the filename, computed sensitivity coefficients are printed. 

Description of additional files printed out in FIM and LSA runs:

1.	Filename containing “_gas”:
   - Identical information to gas_mass_ss file

2.	Filename containing “_surface”:
   - Identical information to surf_cov_ss file

3.	Filename containing “_state”:
   - Please fill this in, it is not the same information as in rctr_state_ss from the files I saw.

4.	Filename containing “_sensitivity”:
   - FIM run: the reaction and corresponding Fischer information matrix 
diagonal value (FIM_Diag). Any reaction with a relatively (when compared to the FIM diagonal element with highest value) small value of FIM diagonal matrix element is not crucial. 
   - Local Sensitivity Analysis (LSA) run: the sensitivity coefficients of all 
the species with respect to the specified reactions and species are given. 
Each row corresponds to the reactions and species input for LSA. Each column indicates the sensitivity coefficients of one species. First row lists the species and first column lists the LSA input reactions and species. A positive value indicates that the concentration of the output species increases with an increase in rate constant of the input reaction, where as a negative value indicates opposite.

