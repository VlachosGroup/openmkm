Welcome to OpenMKM!  

OpenMKM is a multiphysics and multiscale software aimed at Chemical Engineers interested in 
modeling chemical kinetics for heterogeneous catalytic reactions. OpenMKM is opensource software 
and is developed at Delaware Energy Institute, University of Delaware. 
OpenMKM is written in C++ and is compiled and executed from command line. 

## Installation

1. OpenMKM uses [Cantera](http://www.cantera.org) as backend library for kinetics and 
thermodynamics processing. Download a
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
OpenMKM also uses scons as package builder. Go to *openmkm/src*. Edit the *SConstruct* file to specify the
location of cantera and the sundials libraries.
Run *scons -j 4* command to compile. This builds *omkm* executable in the same folder.

## Usage

OpenMKM requires two input files, which are specified as arguments in the command line. 
~~~ bash
./omkm <rctr.yaml> <input.xml>
~~~

1. The first argument *rctr.yaml* is the name of yaml file specifying the reactor model parameters, 
operating conditions, and the names of thermodynamic phases (which are defined in the Cantera XML file 
supplied as second argument) and the starting composition and coverages of gas and surface phases 
respectively.  

2. For syntax of the yaml file, refer to the [rctr_tmplt.yaml](rctr_tmplt.yaml) file, 
which provides extensive comments on the keywords required by hetero_ct in the supplied yaml file.
Multiple yaml files for various reactor models and operating modes are in located in *<openMKM_root>/test_files* folder. 

3. The second argument *input.xml* is Cantera input file in XML format which provides the 
definitions of species, reactions, interactions, and gas, solid and  catalyst surface phases. 
Scripts are available to convert Chemkin input files to Cantera files. Use 
*<CANTERA_ROOT>/interfaces/cython/cantera/ck2cti.py*, which parses gas.inp, surf.inp and thermdat files
to generate the Cantera input file in CTI format. The CTI file needs to converted again into XML file
using *<CANTERA_ROOT>/interfaces/cython/cantera/ctml_writer.py*.
For more information on the Cantera input file formats, refer to 
[Cantera documentation on input file format](https://cantera.org/tutorials/input-files.html).

4. The chemkin input files are not sometimes parsed by ck2cti.py script. Refer to [step by step guide](ck_conversion.md) on how to overcome those issues. 

5. Coverage effects need to be incorporated into the xml file. Since it is easy to work with CTI files, users have the option of specifying coverage effects in CTI file and then convert the CTI file into XML file with the python script *ctml_writer.py*.


## Credits
 Development of OpenMKM is funded by [RAPID Manufacturing Institute](www.aiche.org/rapid).


## Reading Material

1. The documentation on Cantera can be found at its [website](http://www.cantera.org).

2. The governing equations of the  reactor models implemented can be found [here](ReactorModels.pdf).

3. Coverage dependent surface lateral interactions and how to specify them in the Cantera XML file 
are explained [here](SurfaceInteractions.pdf).

