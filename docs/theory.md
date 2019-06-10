---
layout: default
---
# Background Theory
Microkinetic modeling connects the elementary reaction kinetics to thermodynamic 
quantities measurable at reactor scale such as temperature, pressure, species
mass, etc. Here, theoretical underpinnings of microkinetic modeling such as the  
governing equations of the reactor models are presented.

## Reactor Models and Governing Equations
A reactor model is used to specify the conditions under which the chemical
reactions take place. Mathematically, reactor model imposes constraints to
conserve fundamental physical quantities, and these constraints are often
called governing equations, which are nothing but ordinary differential
equations. Different reactor models result in different governing equations. 
The governing equations of the reactor models implemented in OpenMKM are
presented below. 

### Batch Reactor
A batch reactor can be thought of fixed sized tank which is filled with the
reactants and left to evolve under well mixed conditions. The reacting fluid
composition, temperature, and pressure change as a function of time, but at any
point of time, due to the well mixed condition, are uniform through out the
reactor.

#### Mass Balance
The mass of the reacting fluid in a batch reactor is fixed. As as result, 

$$ \frac{dm}{dt} = \sum_k \frac{dm_k}{dt} = 0,$$

where $$m$$ and $$m_k$$ represent the total mass and individual mass of kth
species respectively. Based on $$m$$ and $$m_k$$, the mass fraction of the kth
species is defined as $$Y_k = m_k/m$$. The change in mass fraction of the kth
species is given as

$$ \rho\frac{dY_k}{dt} = (\dot{\omega_k} + \dot{s_k} \frac{A_{cat}}{V}) W_k, $$

where $$\rho$$ is the density of the reacting fluid, $$V$$ is the volume of the 
reactor, $$Y_{k,0}$$ is the initial mass fraction of the kth species,
$$\dot{\omega_k}$$ and $$\dot{s_k}$$ are the production rates of the kth
species in the reacting fluid phase and on catalyst surface respectively, 
$$\frac{A_{cat}}{V}$$ is the ratio of catalyst surface area to reactor volume,
W_k is the molar weight of the kth species.

#### Energy Balance
For modeling the kinetics of heterogeneous catalysis, the reactor is typically
operated under isothermal conditions, where the temperature of the reactor is
fixed. If the temperature of the reactor is allowed to change either due to heat 
released/absorbed due to chemical reactions or due to heat supplied from
external sources, energy balance has to be satisfied. 

$$\frac{dE}{dt} = \frac{dQ}{dt} + \frac{dW}{dt},$$

where $$E$$ is the internal energy of the reactor, $$Q$$ is the external heat
supplied, and $$W$$ is the work done on the reactor. $$W$$ is 0 if the reactor
volume is fixed, otherwise it is $$PdV$$.

The heat flux, $$\frac{dQ}{dt}$$, supplied to the reactor through an outer wall 
with an area, $$A_{wall}$$, and heat transfer cofficient, $$\hat{h}$$, from an
external heat source at temperature, $$T_{ext}$$, is given as 

$$\frac{dQ}{dt} = \hat{h}A_{wall}(T_{ext} - T), $$

where T is the reactor temperature. The internal energy of the reactor could be 
written as $$E = m \sum_k Y_k e_k$$, where $$e_k$$ is the specific internal 
energy of kth species w.r.t. unit mass given in terms of $$J/kg$$ or $$erg/g$$ 
for SI and CGS units respectively. Assuming fixed reactor size, a few more 
definitions to know are $$de_k = c_{v,k}dT$$ and $$c_v = \sum_k{Y_kc_{v,k}}$$, 
where $$c_v$$ is the mass specific heat at constant volume. Plugging the 
various definitions into the energy balance equation results in 

$$\rho c_v \frac{dT}{dt} = -\sum_k{e_k W_k (\dot{\omega_k} + \dot{s_k} \frac{A_{cat}}{V})} +  \hat{h}\frac{A_{wall}}{V}(T_{ext} - T).$$ 

### CSTR
A CSTR can be thought of batch reactor with inlet and outlet that continously
pump in and out the reactaning fluid. The reacting fluid is pumped in at a
volumetric flow rate of $$r$$ and the reacting fluid is well mixed. This 
results in the reacting fluid spending an average time, $$\tau$$ (called 
residence time), inside the reactor, after which it get expelled from the 
reactor through the outlet. The volumetric flow rate, $$r$$ is also defined in 
terms of mass flow rate, $$\dot{m}_0$$, which is given as 
$$\dot{m}_0 = \rho_0 r$$. Similarly, flow rate, $$r$$ and residence time, 
$$\tau$$ are related as $$\tau = V/r$$. At the beginning, the reacting fluid
contains only reactants, but as time progresses, it contains both reactants,
products, and reaction intermediates.  The reactor is typically operated at
steady state conditions, where the composition of the reacting fluid inside the 
CSTR, which is typically different from the composition of the initial feed,
does not change.

#### Mass Balance

$$ \rho\frac{dY_k}{dt} =  \frac{\dot{m}_0}{V} (Y_{k,0}- Y_k) + (\dot{\omega_k} + \dot{s_k} \frac{A_{cat}}{V}) W_k. $$

#### Energy Balance

$$\rho c_p \frac{dT}{dt} = \frac{\dot{m}_0}{V} \sum_k{Y_{k,0}(h_k,0- h_k)} - \sum_k{h_k W_k (\dot{\omega_k} + \dot{s_k} \frac{A_{cat}}{V})} +  \hat{h}\frac{A_{wall}}{V}(T_{ext} - T),$$ 

where $$c_p$$ is the specific heat for unit mass at constant pressure,
$$W_k h_k$$ is the molar specific enthalpy of kth species, and $$h_k$$ is the
specific enthalpy of kth species at unit mass. $$Y_{k,0}$$ and $$h_{k,0}$$
represent the initial mass fractions and initial specific enthalpies (w.r.t.
unit mass) of kth species respectively.

### PFR
A PFR is a tubular reactor with cross sectional area $$A_c$$. Reacting fluid
enters from one side, flows along the axial direction, and exits from the other
side. The state of the reacting fluid is uniform along the radial direction,
but varies along the axial direction. Here the equations presented assume a
fixed cross-sectional area for the PFR. Assuming steady state operating
conditions, the governing equations are formulated for a small differential
volume, $$dV$$, given as $$dV = A_c * dz,$$ where $$dz$$ is differential
length along the axial direction. The axial flow rate of the fluid is
represented with $$u$$ m/s. The volumetric flow rate, $$r$$, defined earlier
is then given as $$r = A_c u$$.

#### Mass balance

[//]: # ($$\frac{dm}{dt} = \int_{SA} \sum_k{\Dot{s_k}W_k dA} +  \int_{CV} \sum_k{\Dot{\omega_k}W_k dV})

$$\frac{d(\rho u)}{d z} = \sum_k{(\dot{w_k}+\dot{s_k}\frac{A_{cat}}{V})W_k}$$

The equation can be further simplified because there can be no net source/sink 
of mass in the reacting fluid phase, which implies
$$\sum_k{\dot{w_k}W_k} = 0 $$. This results in 

$$\frac{d(\rho u)}{d z} = \frac{A_{cat}}{V}\sum_k{\dot{s_k}W_k}.$$

At steady state, the amount of mass adsorbed onto a catalyst surface has to be
equal to the amount of desorbed from the catalyst surface. This results in a
further simplication of 

$$\frac{d(\rho u)}{d z} = 0.$$

$$\rho\frac{du}{d z} + u\frac{d\rho}{d z}  = 0.$$

#### Individual Species Mass Balance

$$ \rho u \frac{dY_k}{dz} +  \frac{A_{cat}}{V}Y_k\sum_k{\dot{s_k}  W_k} =   (\dot{\omega_k} + \dot{s_k} \frac{A_{cat}}{V}) W_k. $$

#### Energy Balance

$$\rho u c_p \frac{dT}{dz} =  - \sum_k{h_k W_k (\dot{\omega_k} + \dot{s_k} \frac{A_{cat}}{V})} +  \hat{h}\frac{A_{wall}}{V}(T_{ext} - T).$$ 

### Surface Coverages on Catalyst Surface
The afore-mentioned conservation equations consider only the species in the
reacting fluid. On a catalytic surface, the site density is considered as fixed.
This results in an additional constraint involving the coverages of surface
species, $$\theta_k$$.

$$\sum_k^{K_S}{\theta_k} = 1.$$

While the energy and mass balance equations are differential in nature, surface
coverage conservation equation is algebraic.