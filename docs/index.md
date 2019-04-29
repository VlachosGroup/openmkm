Welcome to Hetero_ct!  

Hetero_ct is a multiphysics and multiscale software aimed at Chemical Engineers interested in 
modeling chemical kinetics for heterogeneous catalytic reactions. Hetero_ct is opensource software 
and is developed at Delaware Energy Institute, University of Delaware. 
Hetero_ct is written in C++ and is compiled and executed from command line. 

## Installation

1. Hetero_ct uses [Cantera](http://www.cantera.org) as backend library for kinetics and 
thermodynamics processing.Download a
[modified version of cantera source from Github](https://github.com/mbkumar/cantera/tree/hetero_ct).
Please note that the official version and the modified version differ in their implementation of 
coverage effects. Then [install](https://cantera.org/install/compiling-install.html)  Cantera 
from source. Cantera uses scons for package building. To compile the source, here is a sample scons 
command used.
~~~ bash
scons build python_package=full f90_interface=y doxygen_docs=yes \
system_eigen=y system_sundials=y  \
sundials_include=<sundials headers location, usually /usr/include>  \
sundials_libdir=<sundials library file location, usually /usr/lib or /usr/lib64 or /usr/local/lib> \
python_cmd=<path to python command>      \
python_prefix=<path to python site-packages folder such as ~/anaconda3/envs/my_env/lib/python3.7/site-packages> \
extra_inc_dirs="/usr/include/eigen3:/usr/include" -j 4 
~~~
Here *-j 4* is used to speed up the compilation process.

2.
Hetero_ct also uses scons as package builder. Go to *hetero_ct/src*. Edit the *SConstruct* file to specify the
location of cantera and the sundials libraries.
Run *scons -j 4* command to compile. This builds *hetero_ct* executable in the same folder.

## Usage

Hetero_ct requires two input files, which are specified as arguments in the command line. 
~~~ bash
./hetero_ct <rctr.yaml> <input.xml>
~~~

1. The first argument *<rctr.yaml>* is the name of yaml file specifying the reactor model parameters, 
operating conditions, and the names of thermodynamic phases (which are defined in the Cantera XML file 
supplied as second argument) and the starting composition and coverages of gas and surface phases 
respectively.  For syntax of the yaml file, refer to the [rctr_tmplt.yaml](rctr_tmplt.yaml) file, 
which provides extensive comments on the keywords required by hetero_ct in the supplied yaml file .

2. The second argument *<input.xml>* is Cantera input file in XML format which provides the 
definitions of species, reactions, interactions, and gas, solid and  catalyst surface phases. For 
more information on the Cantera XML file format, refer to 
[Cantera documentation on input file format](https://cantera.org/tutorials/input-files.html).

Examples of both files can be found in *hetero_ct/test_files* folder. 
## Credits
 Development of hetero_ct is funded by [RAPID Manufacturing Institute](www.aiche.org/rapid).


## Reading Material

1. The documentation on Cantera can be found at its [website](http://www.cantera.org).

2. The governing equations of the  reactor models implemented can be found [here](ReactorModels.pdf).

3. Coverage dependent surface lateral interactions and how to specify them in the Cantera XML file 
are explained [here](SurfaceInteractions.pdf).

