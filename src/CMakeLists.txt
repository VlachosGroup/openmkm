#
# CMakeLists.txt for OpenMKM
#
# Building OpenMKM requires:
#
#   - Boost 1.67.0 or newer
#   - Cantera
#   - sundials (if Cantera not built with bundled sundials)
#   - kinsol (if Cantera not built with bundled sundials)
#   - yaml-cpp (or using embedded internal build from source)
#
CMAKE_MINIMUM_REQUIRED (VERSION 3.11)
PROJECT (OpenMKM VERSION 0.25.0 LANGUAGES CXX)
SET (CMAKE_CXX_STANDARD 11)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
INCLUDE(GNUInstallDirs)

#
# Locate the Cantera library:
#
SET(CANTERA_PREFIX "" CACHE PATH "Cantera is installed under this directory")
FIND_PATH(Cantera_INCLUDE_DIRS "cantera/Interface.h" HINTS "${CANTERA_PREFIX}/include")
IF ( NOT Cantera_INCLUDE_DIRS )
    MESSAGE(FATAL_ERROR "Could not find the Cantera headers")
ENDIF ()
FIND_LIBRARY(Cantera_LIBRARIES "cantera_shared" HINTS "${CANTERA_PREFIX}/lib" "${CANTERA_PREFIX}/lib64" ENV LD_LIBRARY_PATH)
IF ( NOT Cantera_LIBRARIES )
    FIND_LIBRARY(Cantera_LIBRARIES "cantera" HINTS "${CANTERA_PREFIX}/lib" "${CANTERA_PREFIX}/lib64" ENV LD_LIBRARY_PATH)
    IF ( NOT Cantera_LIBRARIES )
        MESSAGE(FATAL_ERROR "Could not find the Cantera library")
    ENDIF ()
ENDIF ()

#
# Locate the sundials library:
#
SET(SUNDIALS_PREFIX "" CACHE PATH "sundials is installed under this directory")
FIND_PATH(sundials_INCLUDE_DIRS "sundials/sundials_config.h" HINTS "${SUNDIALS_PREFIX}/include")
IF ( NOT sundials_INCLUDE_DIRS )
    #
    # See if we can find sundials as part of Cantera:
    #
    FIND_PATH(sundials_INCLUDE_DIRS "sundials/sundials_config.h" HINTS "${CANTERA_PREFIX}/include/cantera/ext")
    IF ( NOT sundials_INCLUDE_DIRS )
        MESSAGE(FATAL_ERROR "Could not find the sundials header files")
    ENDIF ()
    SET(SUNDIALS_SKIP_LIBRARY TRUE CACHE BOOL "do not look for sundials libraries (assume they are part of Cantera)")
ENDIF ()
IF ( NOT SUNDIALS_SKIP_LIBRARY )
    FIND_LIBRARY(sundials_LIBRARIES
            NAMES "sundials_cvcodes" "sundials_ida" "sundials_kinsol" "sundials_nvecserial"
            HINTS "${SUNDIALS_PREFIX}/lib" "${SUNDIALS_PREFIX}/lib64" ENV LD_LIBRARY_PATH
        )
    IF ( NOT sundials_LIBRARIES )
        MESSAGE(FATAL_ERROR "Could not find the sundials library")
    ENDIF ()
ENDIF ()

#
# Locate the kinsol library:
#
SET(KINSOL_PREFIX "" CACHE PATH "kinsol is installed under this directory")
FIND_PATH(kinsol_INCLUDE_DIRS "kinsol/kinsol.h" HINTS "${KINSOL_PREFIX}/include")
MESSAGE("KINSOL_PREFIX=${KINSOL_PREFIX}")
MESSAGE("kinsol_INCLUDE_DIRS=${kinsol_INCLUDE_DIRS}")
IF ( NOT kinsol_INCLUDE_DIRS )
    #
    # See if we can find kinsol as part of Cantera:
    #
    FIND_PATH(kinsol_INCLUDE_DIRS "kinsol/kinsol.h" HINTS "${CANTERA_PREFIX}/include/cantera/ext")
    MESSAGE("KINSOL_INCLUDE_DIRS=${KINSOL_INCLUDE_DIRS}")
    IF ( NOT kinsol_INCLUDE_DIRS )
        MESSAGE(FATAL_ERROR "Could not find the sundials header files")
    ENDIF ()
    SET(KINSOL_SKIP_LIBRARY TRUE CACHE BOOL "do not look for kinsol libraries (assume they are part of Cantera)")
ENDIF ()
IF ( NOT KINSOL_SKIP_LIBRARY )
    FIND_LIBRARY(kinsol_LIBRARIES "kinsol" "sundials_kinsol" HINTS "${KINSOL_PREFIX}/lib" "${KINSOL_PREFIX}/lib64" ENV LD_LIBRARY_PATH)
    IF ( NOT kinsol_LIBRARIES )
        MESSAGE(FATAL_ERROR "Could not find the kinsol library")
    ENDIF ()
ENDIF ()

#
# Locate the Eigen library:
#
SET(EIGEN_PREFIX "" CACHE PATH "Eigen is installed under this directory")
FIND_PATH(eigen_INCLUDE_DIRS "Eigen/Eigen" HINTS "${EIGEN_PREFIX}")
IF ( NOT eigen_INCLUDE_DIRS )
    #
    # See if we can find eigen as part of Cantera:
    #
    FIND_PATH(eigen_INCLUDE_DIRS "Eigen/Eigen" HINTS "${CANTERA_PREFIX}/include/cantera/ext")
    IF ( NOT eigen_INCLUDE_DIRS )
        MESSAGE(FATAL_ERROR "Could not find the eigen header files")
    ENDIF ()
    SET(EIGEN_SKIP_LIBRARY TRUE CACHE BOOL "do not look for eigen libraries (assume they are part of Cantera)")
ENDIF ()
IF ( NOT EIGEN_SKIP_LIBRARY )
    FIND_LIBRARY(eigen_LIBRARIES
            NAMES "eigen_cvcodes" "eigen_ida" "eigen_kinsol" "eigen_nvecserial"
            HINTS "${EIGEN_PREFIX}/lib" "${EIGEN_PREFIX}/lib64" ENV LD_LIBRARY_PATH
        )
    IF ( NOT eigen_LIBRARIES )
        MESSAGE(FATAL_ERROR "Could not find the eigen library")
    ENDIF ()
ENDIF ()


#
# Locate the Boost libraries:
#
SET(BOOST_ROOT "" CACHE PATH "Boost is installed under this directory")
SET(BOOST_INCLUDEDIR "" CACHE PATH "Boost header files are found in this directory")
SET(BOOST_LIBRARYDIR "" CACHE PATH "Boost libraries are found in this directory")
FIND_PACKAGE(Boost 1.67.0 REQUIRED COMPONENTS filesystem system)

#
# Locate the yaml-cpp library _or_ build a copy of it:
#
OPTION(USE_EXTERNAL_YAML_CPP "Do not use the embedded yaml-cpp library" FALSE)
IF ( USE_EXTERNAL_YAML_CPP )
    SET(YAML_CPP_PREFIX "/usr/local" CACHE PATH "yaml-cpp is installed under this directory")
ENDIF ()
IF ( USE_EXTERNAL_YAML_CPP )
    FIND_PATH(yaml_cpp_INCLUDE_DIRS "yaml-cpp/yaml.h" HINTS "${YAML_CPP_PREFIX}/include")
    IF ( NOT yaml_cpp_INCLUDE_DIRS )
        MESSAGE(FATAL_ERROR "Could not find the yaml-cpp headers")
    ENDIF ()
    FIND_LIBRARY(yaml_cpp_LIBRARIES "yaml-cpp" HINTS "${YAML_CPP_PREFIX}/lib" "${YAML_CPP_PREFIX}/lib" ENV LD_LIBRARY_PATH)
    IF ( NOT yaml_cpp_LIBRARIES )
        MESSAGE(FATAL_ERROR "Could not find the yaml-cpp library")
    ENDIF ()
ELSE ()
    #execute_process(COMMAND git submodule update --init -- yaml-cpp 
    #                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    #SET(yaml_cpp_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp/include
    #CACHE PATH "yaml-cpp include directory")
    ADD_SUBDIRECTORY(yaml-cpp)
ENDIF ()

#
# Pthreads, please:
#
FIND_PACKAGE (Threads REQUIRED)

#
# Source files for the omkm target:
#
SET(omkm_SOURCES
        "main.cpp"
        "util.cpp"
        "zerodReactor.cpp"
        "onedReactor.cpp"
        "io.cpp"
        "pfr1d.cpp"
        "pfr1d_solver.cpp"
        "IdealGasTRampReactor.cpp"
        "reactor_parser.cpp"
        "NonLinearSolver.cpp"
        "KIN_Solver.cpp"
        "ReactorNetHybrid.cpp"
    )

#
# Define the omkm target:
#
ADD_EXECUTABLE(omkm ${omkm_SOURCES})
IF ( NOT USE_EXTERNAL_YAML_CPP )
    ADD_DEPENDENCIES(omkm yaml-cpp)
ENDIF ()
IF ( THREADS_HAVE_PTHREAD_ARG )
    TARGET_COMPILE_OPTIONS(omkm PUBLIC "-pthread")
ENDIF()
TARGET_INCLUDE_DIRECTORIES(omkm PUBLIC
        ${Cantera_INCLUDE_DIRS}
        ${sundials_INCLUDE_DIRS}
        ${kinsol_INCLUDE_DIRS}
        ${eigen_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${yaml_cpp_INCLUDE_DIRS}
    )
TARGET_LINK_LIBRARIES(omkm
        ${Cantera_LIBRARIES}
        ${sundials_LIBRARIES}
        ${kinsol_LIBRARIES}
        #${Boost_LIBRARIES}
        Boost::filesystem Boost::system
        ${yaml_cpp_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
    )

#
# What does "make install" do?
#
INSTALL(TARGETS omkm 
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT applications)
