---
layout: default
---

# OpenMKM Installation Guide

## Windows
On Windows machines, a GUI OpenMKM installer is available. 
After installing the OpenMKM, open the "Start Menu", browse to OpenMKM, and click on "OpenMKM CmdLine". 
This will open a new command line windows with OpenMKM added to the path. To test, type and run "omkm" in the command line 
which should print description about OpenMKM. 

To test further, open the "Start Menu", browse to OpenMKM, and click on *Examples*. This will open a new explorer window containing folder *examples*.
Switch to *examples\model_simul* and copy the *batch* folder to some directory lets say *C:\Users\<your_user_name>\tmp*. Switch to 
*C:\Users\<your_user_name>\tmp\batch* in the previously opened command line prompt and run 
```bash
omkm batch.yaml grimech30.xml
```
which should create lot of new files.


## Linux

### Using conan package manager (The Easy Way)
To reduce the complications associated with compiling so many dependencies, OpenMKM can be installed with conan package manager. 
Conan package manger has many prebuilt binaries compatible with OS and compiler versions. 
If prebuilt binaries are not available for any of the dependencies, they will be compiled and built during OpenMKM installation.

#### Installing and configuring conan
Conan can be installed simply by running 
```bash 
pip install conan
```
However, you may want to use anaconda virtual environment to install conan without polluting your system python. 

After conan is installed, configure it by running 
```bash
conan profile new default --detect 
```
If you are using gcc compiler suite, and its version is >= 5.1 (from the output of gcc --version command), run
```bash
conan profile update settings.compiler.libcxx=libstdc++11 default
```

Configure the conan remotes
```bash
conan remote add vklab https://api.bintray.com/conan/dei/vklab 
```
Here vklab is the name given to the remote specified in the url. The name is arbitrary and you can use any name you like. Running 
```bash
conan remote list 
```
should now show
```bash
conan-center: https://conan.bintray.com [Verify SSL: True]
vklab: https://api.bintray.com/conan/dei/vklab [Verify SSL: True]
```

To install OpenMKM, first create a directory and change the working directory to the newly created directory.
```bash
mkdir openmkm; cd openmkm
```
Now install OpenMKM by running
```bash
conan install openmkm/0.4@dei/vklab -g virtualenv --build missing
```
Installation may take anywhere between a minute to 30 minutes or even longer depending on various factors. 
Once the step is completed, there will be two files activate.sh, deactivate.sh in the directory. Run
```bash
source activate.sh
```
Now typing 
```bash
which omkm
```
should show the location of OpenMKM executable. Notice the long path of the openmkm executable. 
The root location of OpenMKM is one folder above the executable location. 


### Manual compiling (the hard way)

#### Dependencies
OpenMKM source code depends on Cantera, SUNDIALS, Eigen, yaml-cpp, and Boost libraries. 
Additionally compiling OpenMKM requries  users to install a compiler (typically gcc on Linux machines), 
git, scons, and cmake software packages. Here are the
steps that can be used to install dependencies and OpenMKM.
SUNDIALS and Eigen can be downloaded and installed as part of cantera installation. 
Refer to cantera installation steps below.

#### Basics
First install gcc, git, scons, and cmake using OS package managers. This can be quickly accomplished using software patterns.

1. Ubuntu:
```bash
sudo apt-get install build-essential git-all cmake scons
```

2. Fedora: 
```bash
sudo dnf groupinstall "Development Tools"
sudo dnf install git-all cmake scons
```

3: OpenSUSE
```bash
sudo zypper in -t pattern devel_basis
sudo zypper in cmake scons
```

After installing the packages, check the version of gcc installed.
```bash
gcc --version
```


#### Boost, Eigen (optional), & YAML-CPP 
Boost and Eigen are C++ template libraries and often are distributed with Linux. 
Similarly yaml-cpp is available on most of the Linux distros.
Refer to distro package manager documentation on how to obtain them. 

#### SUNDIALS (optional) installation
[SUNDIALS][sundials_page] is a numerical solver suite from Lawrence Livermore
National Lab. **OpenMKM v0.6.0** is tested with SUNDIALS 4.1. For Linux, SUNDIALS is
often available from distro package managers. Check if the version of the distro 
supplied SUNDIALS matches with 4.1. If the required version is not
supplied with distro, one can download v4.1.1 of [SUNDIALS][sundials_download] and follow
the supplied instructions to install. 

### Cantera installation
To install OpenMKM, Cantera needs to be built from source. Download the source
code from Bharat's fork on GitHub. With git this can be accomplished by 
``` bash
git clone https://github.com/mbkumar/cantera.git
```
and then checking out  the *openmkm* branch which contains modifications 
to official version of Cantera to support coverage effects.
``` bash
cd cantera
git checkout openmkm
```
Alternatively, one could select the openmkm branch and down the source code as a zip file.
Both Cantera and OpenMKM use *scons*, a Python based build tool that is
functionally similar to *cmake*, for build purposes. The following command
could be used to compile Cantera on Linux machines, where SUNDIALS is pre-installed.
``` bash
scons build optimize=False python_package=n f90_interface=n doxygen_docs=n \
system_eigen=y system_sundials=y sundials_libdir=/usr/lib64 \
python_cmd=~/anaconda3/envs/<env_name>/bin/python \
python_prefix=~/anaconda3/envs/<env_name>/lib/python3.7/site-packages \
prefix=~/cantera_install use_rpath_linkage=False \
extra_inc_dirs="/usr/include/eigen3:/usr/include"
```

A few points on the above Cantera compilation command:
1. System supplied SUNDIALS and Eigen3 software packages are preinstalled.
   SUNDIALS is installed in /usr/lib64 folder, and Eigen3 is installed at
   /usr/include/eigen3. Note that Eigen3 is a header only package. 

2. Use python_cmd and python_prefix to specify the Python command and the
   Python library paths. In the above command, instead of system Python,
   anaconda Python is used. \<env_name\> represents the virtual environment
   name and python3.7 is the version of Python used. 

3. cantera is installed at *~/cantera_install* folder.

Another way of installation is to install SUNDIALS and Eigen while installing cantera.
``` bash
scons build optimize=False python_package=n f90_interface=n doxygen_docs=n \
system_eigen=n system_sundials=n \
python_cmd=~/anaconda3/envs/<env_name>/bin/python \
python_prefix=~/anaconda3/envs/<env_name>/lib/python3.7/site-packages 
prefix=~/cantera_install use_rpath_linkage=False \
```
Now SUNDIALS library is included within the cantera library. The headers of Eigen and SUNDIALS are now available at *~/cantera_install/include/cantera/ext*.

## OpenMKM installation

1. Download OpenMKM from https://github.com/VlachosGroup/openmkm.

2. Go to OpenMKM_ROOT/src folder, where OpenMKM_ROOT is the top level directory
   of OpenMKM package, and edit SConstruct file to specify the dependencies.

    * Edit CPPPATH variable to specify the location of SUNDIALS, eigen3, and 
      cantera headers.

    * Edit LIBPATH variable to specify the location of sundials and cantera
      libraries. If SUNDIALS is installed as part of cantera installation, no need to specify SUNDIALS location

    * Adjust the CCFLAGS and LINKFLAGS variables to make sure the compiler and
      linker options are consistent with those generated by scons during
      Cantera building stage. 

3. Run ```scons```  

[sundials_page]: https://computation.llnl.gov/projects/sundials/
[sundials_download]: https://computation.llnl.gov/projects/sundials/sundials-software
