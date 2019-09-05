# Building OpenMKM from source

The following pre-requisites are necessary to build OpenMKM from source:

  - Boost 1.67.0 or newer
  - The modified [https://github.com/mbkumar/cantera.git][Cantera library] (*openmkm* branch)

The Cantera library has many dependencies itself:  eigen, sundials, and fmt to name a few.  While it is permissible to build each of those pre-requisites separately, the Cantera build system can also download each of them into the source tree and embed the necessary features from each in the Cantera libraries.  The build procedure outlined in this document for Cantera uses that paradigm to keep Cantera as self-contained as possible.

## Boost

The Boost C++ library is required by Cantera and OpenMKM.  It is best to ensure the same copy of Boost is used by both.  Cantera uses some functionality that is not present in older releases of Boost (e.g. the 1.53.0 release native to CentOS 7.6) so a newer release is necessary (though Cantera's SCons build system doesn't check for a compatible version).

The C++ compiler used to build Cantera and OpenMKM **must** be the same compiler used to build the Boost library.  Boost relies heavily on the features of the underlying C++ compiler when implementing its functionality, so using a Boost library compiled with GCC 4.8 with the GCC 9.1 compiler suite will likely produce difficult-to-debug runtime issues.

## Cantera

Download the modified Cantera source code and checkout the `openmkm` branch:

```
$ cd ${WORKDIR}/sw
$ mkdir cantera
$ cd cantera
$ git clone https://github.com/mbkumar/cantera.git src
$ cd src
$ git checkout openmkm
```

On the Caviness cluster at UD, the appropriate shell environment for the build was setup as:

```
$ vpkg_require python/3.6.5 intel/2019u2 boost/1.71.0:gcc9
$ vpkg_require scons
```

The `cantera.conf` file must be modified to match each distinct computer system.  For the Caviness setup described above, the `cantera.conf` file looked like:

```
import os

MKL_LIB_DIR = os.path.join(os.environ['MKLROOT'],'lib/intel64')
MKL_INC_DIR = os.path.join(os.environ['MKLROOT'],'include')
CANTERA_INSTALL_DIR = os.path.join(os.environ['WORKDIR'],'sw/cantera/20190903-gcc9')

CXX = 'g++'
CC = 'gcc'
prefix = CANTERA_INSTALL_DIR
python_cmd = 'python3'
f90_interface = 'y'
FORTRAN = 'gfortran'
FORTRANFLAGS = ''
blas_lapack_libs = 'mkl_intel_lp64,mkl_gnu_thread,mkl_core'
blas_lapack_dir = MKL_LIB_DIR
thread_flags = '-pthread -lpthread -lgomp'
boost_inc_dir = '/opt/shared/boost/1.71.0-gcc9/include'
system_eigen = 'n'
system_fmt = 'n'
system_sundials = 'n'
googletest = 'none'
env_vars = 'MKLROOT,INTEL_LICENSE_FILE,PATH,LD_LIBRARY_PATH,PYTHONPATH'
cxx_flags = '-std=gnu++11 -m64 -I' + MKL_INC_DIR
cc_flags = '-m64 -I' + MKL_INC_DIR
extra_inc_dirs = '/home/1001/src/cantera/src/ext/sundials/include'
use_rpath_linkage = False
```

Note that the `extra_inc_dirs` path equates to `<git-clone>/src/ext/sundials/include` where `<git-clone>` is the directory in which you cloned the modified Cantera source from github.  Likewise, the `boost_inc_dir` is likely to be different across various computer systems and environments.  And naturally, `CANTERA_INSTALL_DIR` may need to be adjusted.

The *kinsol* sub-component of sundials is not built by Cantera's build environment by default, so a minor patch to their SCons configuration is necessary:

```
$ cd src
$ cat > SCons-patch.diff <<EOT
diff --git a/ext/SConscript b/ext/SConscript
index e4cfd9d..e0bbe59 100644
--- a/ext/SConscript
+++ b/ext/SConscript
@@ -65,7 +65,7 @@ if env['system_sundials'] == 'n':
                                      ConfigBuilder(sundials_configh)))

     # Copy sundials header files into common include directory
-    for subdir in ('sundials', 'nvector', 'cvodes', 'ida', 'sunmatrix', 'sunlinsol'):
+    for subdir in ('sundials', 'kinsol', 'nvector', 'cvodes', 'ida', 'sunmatrix', 'sunlinsol'):
         for header in mglob(env, 'sundials/include/'+subdir, 'h'):
             build(copyenv.Command('#include/cantera/ext/%s/%s' % (subdir, header.name),
                                   '#ext/sundials/include/%s/%s' % (subdir, header.name),
@@ -74,7 +74,7 @@ if env['system_sundials'] == 'n':
     # Compile Sundials source files
     subdirs = ['sundials', 'nvec_ser', 'cvodes', 'ida', 'sunmat_band',
                'sunmat_dense', 'sunmat_sparse', 'sunlinsol_dense',
-               'sunlinsol_band','sunlinsol_spgmr']
+               'sunlinsol_band','sunlinsol_spgmr', 'kinsol']
     if env['use_lapack']:
         subdirs.extend(('sunlinsol_lapackdense', 'sunlinsol_lapackband'))

EOT

$ git apply SCons-patch.diff
```

With the above `cantera.conf` file in place and the SCons scripts altered, the code can be built:

```
$ scons clean
$ scons build
```

The `clean` command ensures no previous builds are present in the source tree.  The `build` sub-command will download the eigen, sundials, and fmt libraries and configure and build them before building the Cantera code itself.  If successful, text will be displayed suggesting that you now install the library into whatever `CANTERA_INSTALL_DIR` was set to:

```
$ scons install
```

Voila:  you now have a Cantera library usable by OpenMKM.

### VALET package definition

On the UD Caviness cluster, VALET is used to automate changes to the shell environment.  A package definition file matching with the build above would look like:

```
cantera:
    description: "an open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport processes"
    prefix: /work/it_nss/sw/cantera
    url: "https://cantera.org"

    actions:
        - variable: PYTHON_CMD
          value: python3
        - pkgconfigdir: ${VALET_PATH_PREFIX}/lib64/pkgconfig
        - variable: PYTHONPATH
          operator: prepend-path
          value: ${VALET_PATH_PREFIX}

    versions:
        "20190903:gcc9":
            description: build from source, 2019-09-03, openmkm branch, GCC 9
            url: "https://github.com/mbkumar/cantera.git"
            dependencies:
                - python/3.6.5
                - intel/2019u2
                - boost/1.71.0:gcc9

```

The file (`cantera.vpkg_yaml`) can be saved to your `~/.valet` directory if this is a personal copy of Cantera, or to your workgroup's `${WORKDIR}/sw/valet` directory.

From a clean login shell, adding this version of Cantera to the environment looks like:

```
$ vpkg_require cantera/20190903:gcc9
Adding dependency `python/3.6.5` to your environment
Adding dependency `intel/2019u2` to your environment
Adding dependency `icu4c/61.1` to your environment
Adding dependency `gcc/9.1.0` to your environment
Adding dependency `boost/1.71.0:gcc9` to your environment
Adding package `cantera/20190903:gcc9` to your environment

$ which ck2cti
/work/it_nss/sw/cantera/20190903-gcc9/bin/ck2cti
```

## OpenMKM

With Boost and Cantera present on the system, OpenMKM can be built.  An SCons build script is present in the `src` sub-directory, as well as a `CMakeLists.txt` file for the CMake build system.  Both should be usable, however the CMake method is the one covered herein and tends to be easier.  CMake 3.1 or newer is required.

There are a number of CMake variables that control the build:

| Variable | Default value | Description |
|-|-|-|
| `CANTERA_PREFIX` | | **MANDATORY** The directory into which Cantera was installed (the value of `CANTERA_INSTALL_DIR` in the preceding section) |
| `SUNDIALS_PREFIX` | | The directory in which Sundials was installed (if you didn't use the embedded option for Cantera) |
| `KINSOL_PREFIX` | | The directory in which kinsol was installed (if you didn't use the embedded option for Cantera) |
| `BOOST_ROOT` | | The directory in which Boost is installed |
| `BOOST_INCLUDEDIR` | | The directory in which Boost header files are installed |
| `BOOST_LIBRARYDIR` | | The directory in which Boost libraries are installed |
| `USE_EXTERNAL_YAML_CPP` | `FALSE` | Set to TRUE if your system has a copy of the yaml-cpp library present |
| `YAML_CPP_PREFIX` | | The directory in which yaml-cpp is installed (if `USE_EXTERNAL_YAML_CPP` is `TRUE`) |

If `USE_EXTERNAL_YAML_CPP` is `FALSE`, CMake will pull the latest yaml-cpp source from [https://github.com/jbeder/yaml-cpp.git][github] and build it as part of the OpenMKM project (similar to Cantera's building eigen, sundials, and fmt itself).

Starting from a clean shell, the build setup proceeds as follows:

```
$ vpkg_devrequire cantera/20190903:gcc9 cmake/3.15.3
Adding dependency `python/3.6.5` to your environment
Adding dependency `intel/2019u2` to your environment
Adding dependency `icu4c/61.1` to your environment
Adding dependency `gcc/9.1.0` to your environment
Adding dependency `boost/1.71.0:gcc9` to your environment
Adding package `cantera/20190903:gcc9` to your environment
Adding package `cmake/3.15.3` to your environment

$ cd ${WORKDIR}/sw
$ mkdir cantera
$ cd cantera
$ git clone https://github.com/VlachosGroup/openmkm.git src
$ cd src
$ mkdir build
$ cd build
```

To bootstrap a build with CMake, the paths to Boost and Cantera are necessary.  We will also ensure that CMake will use the GCC 9.1 compilers we have configured into the environment:

```
$ CC=$(which gcc) CXX=$(which g++) cmake \
  -DBOOST_ROOT="$BOOST_PREFIX" \
  -DCANTERA_PREFIX="$CANTERA_PREFIX" \
  -DCMAKE_INSTALL_PREFIX="/work/it_nss/sw/openmkm/20190903-gcc9" \
  ..

-- The CXX compiler identification is GNU 9.1.0
-- Check for working CXX compiler: /opt/shared/gcc/9.1.0/bin/g++
-- Check for working CXX compiler: /opt/shared/gcc/9.1.0/bin/g++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Found Boost: /opt/shared/boost/1.71.0-gcc9/include (found suitable version "1.71.0", minimum required is "1.67.0") found components:  filesystem system
-- The C compiler identification is GNU 9.1.0
-- Check for working C compiler: /opt/shared/gcc/9.1.0/bin/gcc
-- Check for working C compiler: /opt/shared/gcc/9.1.0/bin/gcc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Performing Test FLAG_WEXTRA
-- Performing Test FLAG_WEXTRA - Success
-- Looking for C++ include pthread.h
-- Looking for C++ include pthread.h - found
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Failed
-- Looking for pthread_create in pthreads
-- Looking for pthread_create in pthreads - not found
-- Looking for pthread_create in pthread
-- Looking for pthread_create in pthread - found
-- Found Threads: TRUE
-- Configuring done
-- Generating done
-- Build files have been written to: /home/1001/sw/openmkm/src/src/build
$ make
   :
[ 95%] Building CXX object CMakeFiles/omkm.dir/KIN_Solver.cpp.o
[ 97%] Building CXX object CMakeFiles/omkm.dir/ReactorNetHybrid.cpp.o
[100%] Linking CXX executable omkm
[100%] Built target omkm
$ make install
[ 70%] Built target yaml-cpp
[100%] Built target omkm
Install the project...
-- Install configuration: ""
-- Installing: /work/it_nss/sw/openmkm/20190903-gcc9/bin/omkm
```

The `omkm` executable has been built and installed.

### VALET package definition

Finally, the VALET package definition matching the above build would look like:

```
openmkm:
    description: "open source software to model heterogeneous catalytic reactions"
    prefix: /work/it_nss/sw/openmkm
    url: "https://github.com/VlachosGroup/openmkm"

    actions:
        - variable: PYTHONPATH
          operator: prepend-path
          value: ${VALET_PATH_PREFIX}/lib/python3.6/site-packages

    versions:
        "20190903:gcc9":
            description: build from source, 2019-09-03, GCC 9
            dependencies:
                - cantera/20190903:gcc9

```

Save this file (`openmkm.vpkg_yaml`) to the appropriate VALET package definition directory.

From a clean login shell, adding this version of Cantera to the environment looks like:

```
$ vpkg_require openmkm/20190903:gcc9
Adding dependency `python/3.6.5` to your environment
Adding dependency `intel/2019u2` to your environment
Adding dependency `icu4c/61.1` to your environment
Adding dependency `gcc/9.1.0` to your environment
Adding dependency `boost/1.71.0:gcc9` to your environment
Adding dependency `cantera/20190903:gcc9` to your environment
Adding package `openmkm/20190903:gcc9` to your environment

$ which omkm
/work/it_nss/sw/openmkm/20190903/bin/omkm
```
