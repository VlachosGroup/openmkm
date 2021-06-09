#
# Build this docker with
#
#     docker build -t omkm .
#
# To run this docker in an interactive bash shell, do
#
#     docker run -it omkm
# 
# If you want to mount a volume run
#     docker run -v /path/to/folder:/target/path -it omkm

# MAINTAINER FLB

##########################################
# Stage: install dependencies
##########################################
FROM ubuntu:20.04 AS dependencies

RUN apt-get -qq update --fix-missing
RUN apt-get install -qq -y --fix-missing python3-dev python3-pip git bash-completion
RUN pip3 install scons==4.1.0
RUN apt-get install -qq -y libyaml-cpp-dev libboost-dev libboost-filesystem-dev libboost-system-dev vim
##########################################
# Stage: build
##########################################
FROM dependencies AS build


RUN git clone https://github.com/SINTEF/cantera.git
RUN cd cantera && git checkout openmkm && scons build optimize=False python_package=n f90_interface=n doxygen_docs=n system_eigen=n system_sundials=n prefix=/cantera_install use_rpath_linkage=False && scons install

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN git clone https://github.com/SINTEF/openmkm.git

RUN /bin/bash -c "source /cantera_install/bin/setup_cantera; cd openmkm/src; scons"

RUN ln -s /openmkm/src/omkm /usr/bin/omkm

# Create and become a normal user
RUN useradd -ms /bin/bash user
USER user




