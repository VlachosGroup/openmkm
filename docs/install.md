---
layout: default
---

# OpenMKM Installation Guide

## Windows Linux or Mac OS (The easy way)
### Windows 
On Windows machines, the easiest way to install and test OpenMKM is by using [Docker](https://www.docker.com/products/docker-desktop/) containers. The latest build of OpenMKM is avaialble at the [VlachosGroup DockerHub repository](https://hub.docker.com/u/vlachosgroup).

1) First [install](https://docs.docker.com/desktop/windows/install/) Docker Desktop app to access the Docker API. 
2) You could use Windows Command Prompt or Windows Powershell to run the Docker command line utilities. (Optional) Use Microsoft Terminal for a more seamless experience. 
3) Install the OpenMKM container in one-line
``` 
docker pull vlachosgroup/openmkm
``` 
4) In the terminal or command prompt change directory to the path where you want to run your simulation, let's say *C:\Users\_your_user_name\tmp*. Then Copy the reactor specification `reactor.yaml` and thermodynamic specification `thermo.xml` to the current working directory. You can find example input files in the *examples\model_simul*. 
5) Run the omkm executable from within the docker container at the location with input files to generate results in the current folder. This should create a lot of new files in the current folder. For help with `docker run` please read [documentation](https://docs.docker.com/engine/reference/commandline/run/). Example command for running from Windows Powershell. 
  - Note 1: This does not work in Windows Command Prompt. 
  - Note 2: In this command we are mounting the current directory `$(pwd)` containing the `reactor.yaml` and `thermo.xml` to be used within the Docker container at the location `/data`, when the omkm binary is executed. Because we using Docker `mount` the results of the simulation are written to the current directory. For more information see [docker mount](https://docs.docker.com/storage/bind-mounts/).
```
docker run --rm -it --mount type=bind,source="$(pwd)",target=/data --workdir="/data" openmkm /bin/bash -c "omkm reactor.yaml thermo.xml"
```

### Linux 
1) Use the appropriate instructions for your linux flavor to [install Docker Desktop](https://docs.docker.com/desktop/linux/install/).
2. Rest of the instructions are similar to Windows. You can simply use Linux Terminal.  

### Mac OS 
1) Use the appropriate instructions for your Mac OS (Intel or Apple silicon) to [install Docker Desktop](https://docs.docker.com/desktop/mac/install/).
2. Rest of the instructions are similar to Windows. You can simply use Mac OS Terminal.  

## Compiling the source-code and all the dependencies (The Hard Way)

### Dependencies

OpenMKM source code depends on Cantera, yaml-cpp, and Boost libraries. Cantera depends on a lot of other libraries which can be automatically installed using it's internal build tool called scons. Additionally compiling OpenMKM requries  users to install a compiler (typically gcc on Linux/Mac machines, or MSVC for Windows machines). In addition to these we need general development build tools such as git, scons, and cmake software packages. Refer to cantera installation steps below.

### Instructions for Linux
#### Install Dev-Tools  
First install gcc, git, scons, and cmake using OS package managers. This can be quickly accomplished using software patterns.
  - For Example in Ubuntu `sudo apt-get install build-essential git-all cmake scons wget python3-setuptools`
  
  - After installing the packages, check the version of gcc installed with: `gcc --version`
  - Fedora: 
  ```
  sudo dnf groupinstall "Development Tools"
  sudo dnf install git-all cmake scons
  ```
  - OpenSUSE
  ```
  sudo zypper in -t pattern devel_basis
  sudo zypper in cmake scons
  ```

#### Install Boost 
Boost is a C++ template library that is often distributed with Linux. Refer to distro package manager documentation on how to obtain them. Example instructions for Ubuntu: `sudo apt-get install libboost-all-dev`. If you install boost using apt-get, you need to figure out the path of the boost header files and boost libraries. Typically, they are installed at /usr/lib/x86_64-linux-gnu/ or /usr/lib. Please search online to figure out where where boost is installed on you machine.   

However, if you want to install boost headers and libraries to your local install folder, you could do something like this. 
- Download version 1.71.0 [Boost library](https://boostorg.jfrog.io/artifactory/main/release/1.71.0/source/boost_1_71_0.tar.gz) directly from Boost website. Or to look at a list of version [Boost Versions](https://sourceforge.net/projects/boost/files/boost/). 
- Make sure to download and install boost to a known location. 
- Assuming you are working in your `HOME` software directory. 
- Make sure you explicitly set the `prefix` directory which is where you install boost headers and libs. You need this path. 

``` 
cd $HOME/software
mkdir boost/
wget https://boostorg.jfrog.io/artifactory/main/release/1.71.0/source/boost_1_71_0.tar.gz
tar -xzf boost_1_71_0.tar.gz
cd boost_1_71_0
./bootstrap.sh -prefix=/home/software/boost/boost-install && \
./b2 && \
./b2 install && \
``` 

#### Cantera installation
To install OpenMKM, Cantera needs to be built from source. Download the source-code from Bharat's fork on GitHub. With git this can be accomplished by the first line. Then checking out  the *openmkm* branch which contains modifications to official version of Cantera to support coverage effects.
```
cd <path to your local software folder>
git clone https://github.com/mbkumar/cantera.git
cd cantera
git checkout openmkm
```
Both Cantera and OpenMKM use *scons*, a Python based build tool that is
functionally similar to *cmake*, for build purposes. The following command could be used to compile Cantera on Linux machines. Please edit the `prefix` to describe where to install cantera and `extra_inc_dirs` to point to your local or system specific boost include location. 
```
scons build optimize=False python_package=n f90_interface=n \ doxygen_docs=n \
prefix=/home/software/cantera-install
extra_inc_dirs="/home/software/boost/boost-install/include"
```

Another way of installation is to create a cantera.conf file and update it with the options used by `scons`. This file is created by scons automatically if you followed the previous step. If this file exists, you could simple do `scons build` and `scons install`. Example cantera.conf file based on previous options
```
python_package = 'n'
f90_interface = 'n'
boost_inc_dir = '/usr/local/include/'
```
#### OpenMKM installation

1. Download OpenMKM from https://github.com/VlachosGroup/openmkm.

2. Go to OpenMKM_ROOT/src folder, where OpenMKM_ROOT is the top level directory
   of OpenMKM package, and edit SConstruct file to specify the dependencies.

    - Edit CPPPATH variable to specify the location of SUNDIALS, eigen3, and 
      cantera headers.
    - Edit LIBPATH variable to specify the location of sundials and cantera
      libraries. If SUNDIALS is installed as part of cantera installation, no need to specify SUNDIALS location
    - Adjust the CCFLAGS and LINKFLAGS variables to make sure the compiler and
      linker options are consistent with those generated by scons during
      Cantera building stage.
3. Example `Sconstruct` file for OpenMKM
```
env = Environment(CPPPATH=["/usr/local/include", "/usr/local/include/cantera/ext"],
                  CCFLAGS=["-std=c++11", "-Wall", 
                           "-g",
                           "-fprofile-arcs", "-ftest-coverage", 
                           "-pthread"],
                  LINKFLAGS=["-fprofile-arcs", "-ftest-coverage", 
                             "-pthread", "-g", "-std=c++11"])
sources = ["main.cpp", "util.cpp", "zerodReactor.cpp", "onedReactor.cpp", 
           "io.cpp", "pfr1d.cpp", "pfr1d_solver.cpp", 
           "IdealGasTRampReactor.cpp", "reactor_parser.cpp",
           "NonLinearSolver.cpp", "KIN_Solver.cpp", "ReactorNetHybrid.cpp"
           ]

env.Program(target="omkm", source=sources, 
            LIBS=["cantera", 
                  "boost_filesystem", "boost_system"],
            LIBPATH=["/usr/local/lib"])
```  
4. Run ```scons```
### Instructions for MacOS
#### Install Dev-Tools  
First install gcc, git, scons, and cmake using [HomeBrew](https://brew.sh/) or MacPorts. 
  - For example using brew you can do like so 
  ```brew install gcc
  brew install cmake
  brew install git 
  brew intall scons
  ```
  - After installing the packages, check the version of gcc installed with: `gcc --version`
#### Install Boost 
Boost is a C++ template library that is often distributed with Linux. Refer to distro package manager documentation on how to obtain them. Example instructions for HomeBrew: `brew install boost`

#### Cantera installation
To install OpenMKM, Cantera needs to be built from source. Download the source-code from Bharat's fork on GitHub. With git this can be accomplished by the first line. Then checking out  the *openmkm* branch which contains modifications to official version of Cantera to support coverage effects.
```
cd <path to your local software folder>
git clone https://github.com/mbkumar/cantera.git
cd cantera
git checkout openmkm
```
Both Cantera and OpenMKM use *scons*, a Python based build tool that is
functionally similar to *cmake*, for build purposes. The following command could be used to compile Cantera on Linux machines. Please edit the `prefix` to describe where to install cantera and `extra_inc_dirs` to point to your local or system specific boost include location. 
```
scons prefix=/Users/<yourusername>/software/cantera \
python_package='n' \
f90_interface='n' \
doxygen_docs='n' \
sphinx_docs='n' \
system_eigen='n' \
system_sundials='n' \
build
scons install 
```

Another way of installation is to create a cantera.conf file and update it with the options used by `scons`. This file is created by scons automatically if you followed the previous step. If this file exists, you could simple do `scons build` and `scons install`. Example cantera.conf file based on previous options
```
prefix = '/Users/<yourusername>/software/cantera'
python_package = 'n'
f90_interface = 'n'
system_eigen = 'n'
system_sundials = 'n'
```
#### OpenMKM installation

1. Download OpenMKM from https://github.com/VlachosGroup/openmkm.

2. Go to OpenMKM_ROOT/src folder, where OpenMKM_ROOT is the top level directory
   of OpenMKM package, and use CMAKE to compile OpenMKM.
   ```
   CC=gcc CXX=g++ cmake -S . -B build \
   -DCMAKE_INSTALL_PREFIX=/Users/<yourusername>/software/openmkm \
   -DCANTERA_PREFIX=<cantera-install-path-from-previous-step> \
   -DBOOST_ROOT=/opt/homebrew/Cellar/boost@1.76/1.76.0_1   #location of boost
   make 
   make install
   ```
3. Check the installation by running the examples. 

### Instructions for Windows
MS Visual Studio based installation will be added in the future. 
