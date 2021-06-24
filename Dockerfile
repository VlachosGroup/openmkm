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
#
# To run omkm on files in current directory
#    docker run --rm -i --user="$(id -u):$(id -g)" -v "$PWD":/temp -w /temp omkm omkm file1.yaml file2.xml
# 
# To make omkm alias for direct exececution of omkm
#    alias omkm='docker run --rm -i --user="$(id -u):$(id -g)" -v "$PWD":/temp -w /temp omkm omkm "$@"'


# MAINTAINER Francesca L. Bleken, francescalb

##########################################
# Stage: install dependencies
##########################################
FROM ubuntu:20.04 AS dependencies

RUN apt-get -qq update --fix-missing
RUN apt-get install -qq -y --fix-missing python3-dev python3-pip git bash-completion
RUN pip3 install scons==4.1.0
RUN apt-get install -qq -y libyaml-cpp-dev libboost-dev libboost-filesystem-dev libboost-system-dev vim

# Link python to python3
RUN ln -s /usr/bin/python3 /usr/bin/python

##########################################
# Stage: build
##########################################
FROM dependencies AS build

# Create and become normal user
# RUN useradd -ms /bin/bash user
# USER user

# Install Cantera, openmkm-version
RUN mkdir -p /sw
RUN cd /sw && git clone https://github.com/SINTEF/cantera.git
RUN cd /sw/cantera && git checkout openmkm && scons build optimize=False python_package=n f90_interface=n doxygen_docs=n system_eigen=n system_sundials=n prefix=/sw/cantera_install use_rpath_linkage=False && scons install

# Set up OpenMKM

RUN mkdir -p /sw/openmkm

COPY docs /sw/openmkm/docs
COPY examples /sw/openmkm/examples
COPY hetero_ct /sw/openmkm/hetero_ct
COPY scripts /sw/openmkm/scripts
COPY src /sw/openmkm/src
COPY test_files /sw/openmkm/test_files

RUN /bin/bash -c "source /sw/cantera_install/bin/setup_cantera; cd /sw/openmkm/src; scons"

ENV PATH="/sw/openmkm/src:${PATH}"

#WORKDIR /home

# Default command

#CMD ["/bin/bash"]

