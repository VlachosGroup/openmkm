---
layout: default
---

## Running OpenMKM
When OpenMKM is installed and compiled by following the [instructions](install), *omkm* executable
is generated. Copy the *omkm* executable file to the simulation folder.
*omkm* executable requires two input files, which are specified as arguments in the command line. 
The first argument is the name of yaml file specifying the reactor model parameters, and the second argument
is the name of Cantera XML input file, which provides the definitions of species, reactions, interactions. 
 Once the three files are in the simulation folder, run the simulation with the following command:
~~~ bash
./omkm <rctr.yaml> <input.xml>
~~~

Notes:
1. The first argument *rctr.yaml* is the name of yaml file specifying the reactor model parameters, 
operating conditions, and the names of thermodynamic phases (which are defined in the Cantera XML file 
supplied as second argument) and the starting composition and coverages of gas and surface phases 
respectively.  

2. For syntax of the yaml file, refer to the [rctr_tmplt.yaml](rctr_tmplt.yaml) file, 
which provides extensive comments on the keywords required by hetero_ct in the supplied yaml file.
Multiple yaml files for various reactor models and operating modes are in located in *OpeMKM_Root/test_files* folder. 

3. The second argument *input.xml* is Cantera input file in XML format which provides the 
definitions of species, reactions, interactions, and gas, solid and  catalyst surface phases. 
Scripts are available to convert Chemkin input files to Cantera files. Use 
*CANTERA_ROOT/interfaces/cython/cantera/ck2cti.py*, which parses gas.inp, surf.inp and thermdat files
to generate the Cantera input file in CTI format. The CTI file needs to converted again into XML file
using *CANTERA_ROOT/interfaces/cython/cantera/ctml_writer.py*.
For more information on the Cantera input file formats, refer to 
[Cantera documentation on input file format](https://cantera.org/tutorials/input-files.html).

4. The chemkin input files are not sometimes parsed by ck2cti.py script. Refer to 
[step by step guide](ck_conversion.md) on how to overcome those issues. 

5. Coverage effects need to be incorporated into the xml file. Since it is easy to work with CTI files, users have the option of specifying coverage effects in CTI file and then convert the CTI file into XML file with the python script *ctml_writer.py*.


## Credits
 Development of OpenMKM is funded by [RAPID Manufacturing Institute](www.aiche.org/rapid).


## Reading Material

1. The documentation on Cantera can be found at its [website](http://www.cantera.org).

2. The governing equations of the  reactor models implemented can be found [here](ReactorModels.pdf).

3. Coverage dependent surface lateral interactions and how to specify them in the Cantera XML file 
are explained [here](SurfaceInteractions.pdf).

