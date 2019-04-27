Welcome to Hetero_ct!  

Hetero_ct is a multiphysics and multiscale software aimed at Chemical Engineers interested in 
modeling chemical kinetics for heterogeneous catalytic reactions. Hetero_ct is opensource software 
and is developed at Delaware Energy Institute, University of Delaware. Hetero_ct uses 
Cantera (www.cantera.org) as backend library for kinetics and thermodynamics processing.

Hetero_ct is written in C++ and executed from command line. 

## Installation
Hetero_ct uses scons as package builder.

1. Download [a modified version of cantera source from Github](https://github.com/mbkumar/cantera/tree/hetero_ct).
Please note that the official version and the modified implement coverage effects differently.
Then [Install cantera from source](https://cantera.org/install/compiling-install.html). 
Cantera uses scons as builder. To compile the source, here is a sample scons command used.
~~~
scons build python_package=full f90_interface=y doxygen_docs=yes system_eigen=y system_sundials=y 
sundials_include=<sundials headers location, usually /usr/include> 
sundials_libdir=<sundials library file location, usually /usr/lib or /usr/lib64 or /usr/local/lib>
python_cmd=<path to python command>
python_prefix=<path to python site-packages folder such as ~/anaconda3/envs/my_env/lib/python3.7/site-packages>
extra_inc_dirs="/usr/include/eigen3:/usr/include" -j 4 
~~~
Here *-j 4* is used to speed up the compilation process.

## Usage


## Credits
 Development of hetero_ct is funded by 
RAPID Manufacturing Institute (www.aiche.org/rapid).


## References


