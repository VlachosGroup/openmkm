---
layout: default
---

# Input Files

## Introduction
To run OpenMKM, two input files are required:

1. A YAML file that specifies parameters related to the reactor.
2. An XML file that specifies thermodynamic and kinetic parameters.

When running OpenMKM in the command line, the following syntax is used:

```bash
omkm reactor.yaml thermo.xml
```

where ``reactor.yaml`` and ``thermo.xml`` correspond to (1) and (2) above.

These files can be generated manually or by using software, such as
[pMuTT][pmutt_omkm_example].

## YAML Input File.

The YAML file, as its name implies, uses the [YAML](https://yaml.org) format,
which is a human-friendly data format. This file contains information related
to the reactor, such as:

- operating conditions,
- phases,
- composition of gas or surface at initial conditions (t = 0 or z = 0),
- solver parameters.

### Generating YAML Files

To generate a YAML file, one could:

- write it manually,
- modify one of [OpenMKM's YAML examples][yaml_examples],
- use [pMuTT's write_yaml function][pmutt_write_yaml]

[This YAML tutorial][yaml_tutorial] introduces the YAML-language syntax. Since
OpenMKM's YAML file uses only simple types, readers do not need to write passed
the 'Advanced Options' heading.

There are several validators available to check your YAML file for syntax
errors.[This site][yaml_validator] will show errors and can even reformat
your file.


### OpenMKM YAML Format

Quantities with units can be specified in the YAML file as a float or a string.
If a float is entered, OpenMKM assumes it is in SI units. If a string is 
entered, then the units can be explictly specified. For example, the 
reactor volume could either be specified as:

```yaml
volume: 20 # Units assumed to be m3
```

or

```yaml
volume: "20 cm3" # Explictly state the volume is in cm3
```

Below is an exhaustive list of available options supported by OpenMKM. The level
indicates whether this field is nested under another.

[//]: # To regenerate table, copy HTML code and paste in https://www.tablesgenerator.com/html_tables, modify the contents and regenerate the HTML code

<table>
<thead>
  <tr>
    <th>1st Level</th>
    <th>2nd Level</th>
    <th>3rd Level</th>
    <th>Type</th>
    <th>Required</th>
    <th>Description</th>
    <th>Default Units</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>reactor</td>
    <td></td>
    <td></td>
    <td>dictionary</td>
    <td>Y</td>
    <td>Reactor parameters</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>reactor_type</td>
    <td></td>
    <td>string</td>
    <td>Y</td>
    <td>Type of reactor. Supported&nbsp;&nbsp;&nbsp;options: <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'pfr' (plug flow reactor)<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'pfr_0d' (plug flow reactor modeled as a series of CSTRs)<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'cstr' (continuously stirred tank reactor)<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'batch' (batch reactor)</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>mode</td>
    <td></td>
    <td>string</td>
    <td>Y</td>
    <td>Operation of reactor. Supported&nbsp;&nbsp;&nbsp;options:<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'isothermal' (constant temperature)<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'adiabatic' (no heat flow)</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>nodes</td>
    <td></td>
    <td>integer</td>
    <td>N</td>
    <td>Number of CSTRs to model the PFR. Only applicable if   ``reactor_type='pfr_0d'``</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>volume</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Volume of reactor</td>
    <td>m^3</td>
  </tr>
  <tr>
    <td></td>
    <td>temperature</td>
    <td></td>
    <td>float</td>
    <td>Y</td>
    <td>Temperature of reactor</td>
    <td>K</td>
  </tr>
  <tr>
    <td></td>
    <td>pressure</td>
    <td></td>
    <td>float</td>
    <td>Y</td>
    <td>Pressure of reactor</td>
    <td>Pa</td>
  </tr>
  <tr>
    <td></td>
    <td>area</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Surface area of reactor. Only applicable if ``reactor_type='pfr'``</td>
    <td>m^2</td>
  </tr>
  <tr>
    <td></td>
    <td>length</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Length of the reactor. Only applicable if ``reactor_type='pfr'``</td>
    <td>m</td>
  </tr>
  <tr>
    <td></td>
    <td>cat_abyv</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Catalyst surface area to reactor volume ratio. Only required if a surface phase is specified</td>
    <td>m^-1</td>
  </tr>
  <tr>
    <td>inlet_gas</td>
    <td></td>
    <td></td>
    <td>dictionary</td>
    <td>N</td>
    <td>Inlet gas properties. Not applicable&nbsp;&nbsp;&nbsp;if ``reactor_type = 'batch'``</td>
    <td></td>
  </tr>
  <tr>
    <td></td>
    <td>flow_rate</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Volumetric flow rate of inlet stream. Not required if&nbsp;&nbsp;&nbsp;``inlet_gas.residence_time`` or ``inlet_gas.mass_flow_rate`` is specified</td>
    <td>m^3/s</td>
  </tr>
  <tr>
    <td></td>
    <td>residence_time</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Residence time of reactor.  Not&nbsp;&nbsp;&nbsp;required if ``inlet_gas.flow_rate`` or ``inlet_gas.mass_flow_rate`` is&nbsp;&nbsp;&nbsp;specified</td>
    <td>s</td>
  </tr>
  <tr>
    <td></td>
    <td>mass_flow_rate</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Mass flow rate of inlet stream.&nbsp;&nbsp;&nbsp;&nbsp;Not required if ``inlet_gas.residence_time`` or&nbsp;&nbsp;&nbsp;``inlet_gas.flow_rate`` is specified</td>
    <td>kg/s</td>
  </tr>
  <tr>
    <td>simulation</td>
    <td></td>
    <td></td>
    <td></td>
    <td>N</td>
    <td>Simulation options</td>
    <td></td>
  </tr>
  <tr>
    <td></td>
    <td>end_time</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Reactor simulation time. For continuous reactors, the system is assumed   to reach steady state by this time</td>
    <td>s</td>
  </tr>
  <tr>
    <td></td>
    <td>transient</td>
    <td></td>
    <td>boolean</td>
    <td>N</td>
    <td>If True, transient results written&nbsp;&nbsp;&nbsp;to output files. Otherwise, transient files are empty</td>
    <td></td>
  </tr>
  <tr>
    <td></td>
    <td>stepping</td>
    <td></td>
    <td>string</td>
    <td>N</td>
    <td>Type of time stepping for&nbsp;&nbsp;&nbsp;transient operation. Pairs with ``simulation.step_size``. Supported&nbsp;&nbsp;&nbsp;options:<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'logarithmic'<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'regular'</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>step_size</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Step size. If ``simulation.stepping = 'logarithmic'``, represents the&nbsp;&nbsp;&nbsp;ratio between the next step and the current step. If ``simulation.stepping =&nbsp;&nbsp;&nbsp;'regular'``, represents the time between the next step and the current step&nbsp;&nbsp;&nbsp;in units of time.</td>
    <td>s</td>
  </tr>
  <tr>
    <td></td>
    <td>init_step</td>
    <td></td>
    <td>float</td>
    <td>N</td>
    <td>Initial time step</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>output_format</td>
    <td></td>
    <td>string</td>
    <td>N</td>
    <td>Format for output files.&nbsp;&nbsp;&nbsp;Supported options:<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'csv'<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 'dat'</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>solver</td>
    <td></td>
    <td>dictionary</td>
    <td>N</td>
    <td>Solver options</td>
    <td></td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>atol</td>
    <td>float</td>
    <td>N</td>
    <td>Absolute tolerance of solver</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>rtol</td>
    <td>float</td>
    <td>N</td>
    <td>Relative tolerance of solver</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>multi_input</td>
    <td></td>
    <td>dictionary</td>
    <td>N</td>
    <td>Multiple runs where temperature, pressure, and flow rate can be varied</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>multi_T</td>
    <td>list of float</td>
    <td>N</td>
    <td>Multiple temperatures of reactor</td>
    <td>K</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>multi_P</td>
    <td>list of float</td>
    <td>N</td>
    <td>Multiple pressures of reactor</td>
    <td>Pa</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>multi_flow_rate</td>
    <td>list of float</td>
    <td>N</td>
    <td>Multiple volumetric flow rates</td>
    <td>m^3/s</td>
  </tr>
  <tr>
    <td></td>
    <td>sensitivity</td>
    <td></td>
    <td>dictionary</td>
    <td>N</td>
    <td>Sensitivity analysis options</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>full</td>
    <td>boolean</td>
    <td>N</td>
    <td>If True, runs sensitivity analysis using the Fisher Information Matrix (FIM)</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>reactions</td>
    <td>list of str</td>
    <td>N</td>
    <td>IDs of reactions to perform local sensitivity analysis (LSA)</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>species</td>
    <td>list of str</td>
    <td>N</td>
    <td>Names of species to perform local sensitivity analysis (LSA)<br></td>
    <td></td>
  </tr>
  <tr>
    <td>phases</td>
    <td></td>
    <td></td>
    <td></td>
    <td>Y<br></td>
    <td>Phase properties</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>bulk</td>
    <td></td>
    <td>dictionary</td>
    <td>N</td>
    <td>Bulk phase properties</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>name</td>
    <td>string</td>
    <td>N</td>
    <td>Name of bulk phase</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>gas</td>
    <td></td>
    <td></td>
    <td>N</td>
    <td>Gas properties</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>name</td>
    <td>string</td>
    <td>N</td>
    <td>Name of gas phase</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>initial_state</td>
    <td>string</td>
    <td>N</td>
    <td>Non-zero initial mole fractions for gas phase. Multiple species should be separated by commas. For example: "H2: 0.4,N2: 0.6"</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td>surfaces</td>
    <td></td>
    <td>list of dictionaries</td>
    <td>N</td>
    <td>Surface phase properties. Note that multiple surface can be specified.</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>name</td>
    <td></td>
    <td>N</td>
    <td>Name of surface phase</td>
    <td>-</td>
  </tr>
  <tr>
    <td></td>
    <td></td>
    <td>initial_state</td>
    <td></td>
    <td>N</td>
    <td>Non-zero initial coverages for surface phase. Multiple species should be separated by commas. For example: "H2: 0.4,N2: 0.6"</td>
    <td>-</td>
  </tr>
</tbody>
</table>

### Sample YAML file.

```yaml
inlet_gas:
    flow_rate: "2 cm3/s"
phases:
    bulk:
        name: bulk
    gas:
        initial_state: "CH3CH3:1.0"
        name: gas
    surfaces:
    -   initial_state: "M(S):1.0"
        name: terrace
reactor:
    cat_abyv: "200 /cm"
    mode: "isothermal"
    nodes: 50
    pressure: "1 atm"
    temperature: 973
    type: "pfr_0d"
    volume: "1 cm3"
simulation:
    end_time: "1000 s"
    init_step: 1.0e-10
    output_format: "csv"
    step_size: 10
    stepping: "logarithmic"
    transient: true
```

## OpenMKM XML Format

The XML file specifies physical, thermodynamic, and kinetic parameters of:

- phases
- species
- reactions
- interactions
- BEPs

However, the XML can be difficult to generate manually so we recommend creating
a CTI file, which has a very similar syntax to Python code. CTI files can be
converted to an XML file using [ctml_writer.py][ctml_writer]. Note there are
some differences between OpenMKM's ctml_writer.py and the script supplied by
Cantera.

Alternatively, Chemkin input files can be converted to Cantera format using
[ck2cti.py][ck2cti_script].

## OpenMKM CTI Format

The OpenMKM CTI format borrows heavily from
[Cantera's CTI format][cantera_cti_docs]. Their documentation has a lot of
useful information about specifying types and the syntax.

### Units

Default units for the CTI file can be specified using a ``units`` directive.
It takes the following parameters:

- ``length`` -- Optional
- ``time`` -- Optional
- ``quantity`` -- Optional
- ``energy`` -- Optional. Used for lateral interactions.
- ``act_energy`` -- Optional. Used for activation energies and BEP
  relationships.
- ``pressure`` -- Optional
- ``mass`` -- Optional

Below, we show a sample ``units`` directive.

```python
units(length="m", time="s", quantity="mol", energy="kcal",
      act_energy="kcal/mol", pressure="atm", mass="g")
```

[The complete list of supported units is available here.][supported_units]

If you wish to override the units for a specific quantity, you can use a tuple
(note that this is different than the YAML syntax which uses strings).
For example, the ``units`` directive above specifies ``quantity="mol"`` and
``length="m"`` but we would like to specify the site density as
1.e15 molecules/cm2. This can be done using:

```python
site_density = (1.e15, "molec/cm2")
```

### Phases

Phases can control transport, kinetic, and thermodynamic properties of the
species. Every simulation will require at least one.
[The comprehensive set of parameters can be found here.][cantera_phase_docs]
Below, we show the three most common for heterogeneous catalysis.

#### Ideal Gas

The ``ideal_gas`` directive models species in the gas phase.

```python
ideal_gas(name="gas",
          elements="H C", # Elements present in species
          species="H2 CH4 CHCH CH2CH2 CH3CH3",  # Gaseous species
          reactions=["0001", "0003 to 0005"]) # Reaction IDs
```

#### Stoichiometric Solid

The ``stoichiometric_solid`` directive is used to specify a bulk solid phase.
This is useful for balancing empty sites.

```python
stoichiometric_solid(name="bulk",
                     elements="Ru", # Elements present in species
                     species="M(B)", # Bulk species
                     density=(12.45, 'g/cm3')) # Ignored property for isothermal
```

#### Interactive Interface

The ``interacting_interface`` directive is used to specify a catalytic surface.
For heterogeneous catalysis, this is the most important phase. It links two 
phases (usually an ``ideal_gas`` to a ``stoichiometric_solid``), contains the
surface intermediates, their interactions, reactions, and BEP relationships.

```python
interacting_interface(name="terrace",
                      elements="Ru H C", # Elements present in species
                      species="""M(S) H(S) C(S) CH(S) CH2(S) CH3(S) CH4(S) CC(S)
                                 CCH(S) CCH2(S) CHCH(S) CCH3(S) CHCH2(S)
                                 CHCH3(S) CH2CH2(S) CH2CH3(S) CH3CH3(S) """,
                      phases="gas bulk", # Connected phases
                      site_density=(2.49e-09, 'mol/cm2'), # Site density
                      interactions=["0000 to 0195"], # Lateral interactions on this phase
                      reactions=["0000 to 0004"], # Reactions on this phase
                      beps="C-C C-H") # BEPs associated with reactions on this phase
```

### Species

OpenMKM 
Each species present in your mechanism should have its thermodynamic properties
specified using one of the following:

- NASA7 polynomial
- NASA9 polynomial
- Shomate polynomial

[More information can be found on Cantera's documentation.][cantera_thermo_docs]

#### NASA7 Polynomial

NASA7 polynomials use 7 coefficients to describe the heat capacity, enthalpy,
and entropy of a species. Common species are available in the
[Burcat database][burcat_db] or can be generated using
[pMuTT's Nasa class][pmutt_nasa].

$$\frac {Cp} {R} = a_{1} + a_{2} T + a_{3} T^{2} + a_{4} T^{3} + a_{5} T^{4}$$

$$\frac {H} {RT} = a_{1} + a_{2} \frac {T} {2} + a_{3} \frac {T^{2}} {3} 
+ a_{4} \frac {T^{3}} {4} + a_{5} \frac {T^{4}} {5} + a_{6} \frac {1} {T}$$

$$\frac {S} {R} = a_{1} \ln {T} + a_{2} T + a_{3} \frac {T^{2}} {2} 
+ a_{4} \frac {T^{3}} {3} + a_{5} \frac {T^{4}} {4} + a_{7}$$

```python
species(name="CH4(S)", # Label used to identify the species
        atoms="C:1 H:4", # Atomic composition of species separated with spaces
        size=1, # Number of sites occupied by species. Only needed for surface species
        thermo=(NASA([250.0, 480.0], # This interval is valid between 250 - 480 K
                     [ 6.37200798E+00, -4.66939392E-03,  6.26492184E-06,
                       4.28698448E-08, -5.33259424E-11, -1.23940270E+04,
                      -2.09817174E+01]),
                NASA([480.0, 1500.0], # This interval is valid between 480 - 1500 K
                     [ 3.73282076E+00,  7.75238626E-03,  1.47755181E-06,
                      -3.18500889E-09,  8.84603159E-13, -1.20467882E+04,
                      -9.12054966E+00])))
```

#### NASA9 Polynomial

NASA9 polynomials use 9 coefficients to describe the heat capacity, enthalpy,
and entropy of a species. Common species are available in the
[Burcat database][burcat_db] or can be generated using
[pMuTT's Nasa9 class][pmutt_nasa9].

$$\frac {Cp} {R} = a_{1} T^{-2} + a_{2} T^{-1} + a_{3} + a_{4} T
+ a_{5} T^{2} + a_{6} T^{3} + a_{7} T^{4}$$

$$\frac {H} {RT} = -a_{1} \frac {T^{-2}} {2} +
a_{2} \frac {ln {T}} {T} + a_{3} + a_{4} \frac {T} {2} + a_{5}
\frac {T^{2}} {3} + a_{6} \frac {T^{3}} {4} + a_{7} \frac {T^{4}} {5} +
a_{8} \frac {1} {T}$$

$$\frac {S} {R} = -a_{1}\frac {T^{-2}} {2} - a_2 \frac {1} {T} +
a_{3} \ln {T} + a_{4} T + a_{5} \frac {T^{2}} {2} + a_{6}
\frac {T^{3}} {3} + a_{7}\frac {T^{4}} {4} + a_{9}$$

```python
species(name="CH4(S)", # Label used to identify the species
        atoms="C:1 H:4", # Atomic composition of species separated with spaces
        size=1, # Number of sites occupied by species. Only needed for surface species
        thermo=(NASA9([200.0, 1000.0], # This interval is valid between 200 - 1000 K
                      [22103.71497, -381.846182, 6.08273836, -0.00853091441,
                       1.384646189e-05, -9.62579362e-09, 2.519705809e-12,
                       710.846086, -10.76003744]),
                NASA9([1000.0, 6000.0], # This interval is valid between 1000 - 6000 K
                      [587712.406, -2239.249073, 6.06694922, -0.00061396855,
                       1.491806679e-07, -1.923105485e-11, 1.061954386e-15,
                       12832.10415, -15.86640027]),
                NASA9([6000.0, 20000.0], # This interval is valid between 6000 - 20000 K
                      [831013916.0, -642073.354, 202.0264635, -0.03065092046,
                       2.486903333e-06, -9.70595411e-11, 1.437538881e-15,
                       4938707.04, -1672.09974])))
```

#### Shomate

Shomate polynomials use 9 coefficients to describe the heat capacity, enthalpy,
and entropy of a species. Common species are available in the
[NIST database][nist_db] or can be generated using
[pMuTT's Shomate class][pmutt_shomate].

$$\frac{c_P}{R}=\frac{1}{R}\bigg(A+Bt+Ct^2+Dt^3+\frac{E}{t^2}
\bigg)$$

$$\frac{H}{RT}=\frac{1}{RT}\bigg(At+B\frac{t^2}{2}+C\frac{t^3}{3}
+D\frac{t^4}{4}-\frac{E}{t}+F\bigg)$$

$$\frac{S}{R}=\frac{1}{R}\bigg(A\ln(t)+Bt+C\frac{t^2}{2}+D
\frac{t^3}{3}-\frac{E}{2t^2}+G\bigg)$$

where $$t=\frac{T}{1000}$$ in K

```python
species(name="CH4", # Label used to identify the species
        atoms="C:1 H:4", # Atomic composition of species separated with spaces
        size=1, # Number of sites occupied by species. Only needed for surface species
        thermo=Shomate([298, 1300], # This interval is valid between 298 - 1300 K
                       [-7.03029000E-01,  1.08477300E+02, -4.25215700E+01,
                         5.86278800E+00,  6.78565000E-01, -7.68437600E+01,
                         1.58716300E+02]))
```

### Reactions

Cantera (and consequently OpenMKM) support a wide variety of reactions. See
[Cantera's documentation][cantera_reaction_docs] to see supported types. 


#### Surface Reactions

For surface chemistry, the ``surface_reaction`` directive will be the most
common, which uses a modified Arrhenius expression.

$$ k = A T^\beta \exp(-\frac {E_a}{RT})$$

where $$k$$ is the rate constant, $$A$$ is the pre-exponential factor,
$T$ is the temperature, $$\beta$$ adds explicit temperature dependence to the pre-exponential parameter, $$E_a$$ is the activation energy, and $R$ is the molar.

For surface reactions, $$A$$ is typically calculated using:

$$ A^{SR} = \frac {k_B} {h \sigma^{m-1}} $$

where $$k_B$$ is the Boltzmann constant, $$h$$ is the Planck's constant,
$$\sigma$$ is the site density, and $$m$$ is the number of surface species
(including empty sites).

If $$A$$ is represented this way, we typically replace $$E_a$$ with the 
Gibbs energy of activation, $$\Delta G^{\ddag}$$, to include entropic effects
based on transition state theory.

$$\Delta G^{\ddag} = \Delta H^{\ddag} - T \Delta S^{\ddag}$$

where $$\Delta G^{\ddag}$$ is the Gibbs energy of activation,
$$\Delta H^{\ddag}$$ is the enthalpy of activation, $$T$$ is the temperature,
and $$\Delta S^{\ddag}$$ is the entropy of activation.

```python
surface_reaction("CH4(S) + M(S) <=> CH3(S) + H(S) + M(B)", # Reaction string
                 [ 8.36812e+18, 1,  5.55117e+00], # A, beta, and Ea respectively
                 id="0001") # ID referenced by phase BEP directives
```

#### Adsorption Reactions

Adsorption reactions also use the ``surface_reaction`` directive but with a
couple of changes. The ``stick`` keyword must be specified and the first
parameter represents the sticking coefficient, $$s$$, instead of the
pre-exponential factor, $$A$$.

$$A^{ads} = \frac {s}{\sigma^{m}}\sqrt \frac{RT}{2\pi M_i}$$

where $$A^{ads}$$ is the pre-exponential for an adsorption reaction,
$$s$$ is the sticking coefficient, $$\sigma$$ is the site density,
$$m$$ is the number of surface species including empty sites,
$$R$$ is the molar gas constant, $$T$$ is the temperature, and
$$M_i$$ is the molecular weight of the gas species.

```python
surface_reaction("CH2CH2 + M(S) <=> CH2CH2(S) + M(B)", # Reaction string
                 stick( 5.00000e-01, 0.0,  3.67653e+00), # Stick keyword indicates adsorption, values are s, beta, and Ea
                 id="0001")
```

### Lateral Interactions

The interactions among adsorbates can cause weaker or tighter binding energies.
These are handled in OpenMKM using lateral interactions. Currently, we support
piecewise linear interactions.

$$ H_{ij} = \omega_{ijk} \theta_{j} + b_{ijk} $$

where $$i$$ is the species affected by species $$j$$, $$H_{ij}$$ is the
enthalpic change due the coverage effect, $$\theta_{j}$$ is the coverage
of species $$j$$, $$\omega_{ijk}$$ and $$b_{ijk}$$ are the slope (strength)
and intercept of the coverage effect at interval $$k$$.

They can be specified in two ways using the ``lateral_interaction`` directive.
Either pairwise or using a matrix.

#### Pairwise Interactions

```python
lateral_interaction(species="H(S) C(S)", # Species i, Species j
                    coverage_thresholds=[0, 0.11, 1], # Intervals to change 
                    strengths=[0.0, -19.0], # Slope of lateral interaction
                    id="0001")
```

In the above example, the lateral interaction strength between H(S) and C(S)
is 0 between 0 - 0.11 ML of C(S), and -19 kcal/mol (specified using the energy
parameter in the ``units`` directive) between 0.11 ML - 1 ML of C(S).

#### Matrix Interactions

The lateral interaction matrix takes the following parameters:
- ``species`` -- species which have lateral interactions
- ``interaction_matrix`` -- a matrix denoting the lateral interaction
  strengths, and
- ``coverage_thresholds`` -- coverage thresholds above which lateral
  interactions modify enthalpy of species.

```python
lateral_interactions(
    species = 'N(S1) H(S1) NH3(S1) NH2(S1) NH(S1)',
    interaction_matrix = [[-47.0179, -17.7545, -25.1631, -20.7620, -48.7823],
                          [-17.7545,  -6.7043,  -9.5019,  -7.8400, -18.4208],
                          [-25.1631,  -9.5019, -13.4668, -11.1115, -26.1074],
                          [-20.762,   -7.8400, -11.1115,  -9.1681, -21.5412],
                          [-48.7823, -18.4208, -26.1074, -21.5412, -50.6129]],
    coverage_thresholds = [0, 0, 0, 0, 0])
```

### Bell-Evans-Polanyi (BEP) relationships

Reaction activation energies can be specified using BEP relationships. To
specify BEP relationships, use the ``bep`` directive, which takes the following:

- ``slope`` -- slope of the BEP linear relationship 
- ``intercept`` -- intercept of the BEP relationship 
- ``direction`` -- Direction of the BEP. Supported values are *cleavage* and
  ``synthesis``. 
- ``id`` -- Optional argument to distinguish the BEP relationship. If no value
  is given, a numerical 4 digit integer starting with ``0001`` is used to
  distinguish the BEP relationships
- ``cleavage_reactions`` -- List of associated reaction ids where a bond is broken
- ``synthsis_reactions`` -- List of associated reaction ids where a bond is formed

``` python
bep(id='C-H', # Used by phases to identify the BEPs present
    slope=1.02, # Slope
    intercept=(24.44, 'kcal/mol'), # Intercept
    direction='cleavage', # Direction (cleavage or synthesis)
    cleavage_reactions=["0003", "0005", "0007 to 0009", "0013 to 0017"],
    synthesis_reactions=[])
```

### Chemkin Users 

Use the conversion script,
``\<CANTERA\_ROOT\>/interfaces/cython/cantera/ck2cti.py``, which parses gas.inp,
surf.inp and thermdat files to convert Chemkin files into the Cantera CTI input
format. For more information, refer to
[Cantera documentation][cantera_ck2cti_docs] on input file format.

The Chemkin input files are not sometimes parsed by ck2cti.py script due to the
bulk phase definition in ``surf.inp``. Remove the bulk phase definition and
retry. If it works, add the missing bulk phase definition directly into the CTI file using the ``stoichiometric_solid`` keyword (see section above).

[examples]: https://github.com/VlachosGroup/openmkm/tree/master/examples
[ctml_writer]: https://github.com/VlachosGroup/openmkm/blob/master/scripts/ctml_writer.py
[yaml_examples]: https://github.com/VlachosGroup/openmkm/tree/master/examples/rctr_input_files
[pmutt_write_yaml]: https://vlachosgroup.github.io/pMuTT/api/kinetic_models/omkm/pmutt.io.omkm.write_yaml.html
[yaml_tutorial]: https://rollout.io/blog/yaml-tutorial-everything-you-need-get-started/
[yaml_validator]: https://jsonformatter.org/yaml-validator
[ck2cti_script]: https://github.com/VlachosGroup/openmkm/blob/master/scripts/ck2cti.py
[cantera_cti_docs]: https://cantera.org/tutorials/cti/cti-syntax.html
[supported_units]: https://cantera.org/tutorials/cti/cti-syntax.html#recognized-units
[cantera_thermo_docs]: https://cantera.org/science/science-species.html#thermodynamic-property-models
[burcat_db]: http://garfield.chem.elte.hu/Burcat/burcat.html
[pmutt_nasa]: https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa.html
[pmutt_nasa9]: https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa9.html#
[nist_db]: https://webbook.nist.gov/chemistry/
[pmutt_shomate]: https://vlachosgroup.github.io/pMuTT/api/empirical/shomate/pmutt.empirical.shomate.Shomate.html
[cantera_reaction_docs]: https://cantera.org/science/reactions.html
[pmutt_omkm_example]: https://vlachosgroup.github.io/pMuTT/examples_jupyter/examples.html#openmkm-io-example
[cantera_ck2cti_docs]: https://cantera.org/tutorials/ck2cti-tutorial.html
[cantera_phase_docs]: https://cantera.org/documentation/dev/sphinx/html/yaml/phases.html#phase-definitions