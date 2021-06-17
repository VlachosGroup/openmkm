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
RUN useradd -ms /bin/bash user
USER user

# Install Cantera, openmkm-version
RUN mkdir -p /home/user/sw
RUN cd /home/user/sw && git clone https://github.com/SINTEF/cantera.git
RUN cd /home/user/sw/cantera && git checkout openmkm && scons build optimize=False python_package=n f90_interface=n doxygen_docs=n system_eigen=n system_sundials=n prefix=/home/user/sw/cantera_install use_rpath_linkage=False && scons install

# Set up OpenMKM

RUN mkdir -p /home/user/sw/openmkm

COPY --chown=user:user docs /home/user/sw/openmkm/docs
COPY --chown=user:user examples /home/user/sw/openmkm/examples
COPY --chown=user:user hetero_ct /home/user/sw/openmkm/hetero_ct
COPY --chown=user:user scripts /home/user/sw/openmkm/scripts
COPY --chown=user:user src /home/user/sw/openmkm/src
COPY --chown=user:user test_files /home/user/sw/openmkm/test_files

RUN /bin/bash -c "source /home/user/sw/cantera_install/bin/setup_cantera; cd /home/user/sw/openmkm/src; scons"

ENV PATH="/home/user/sw/openmkm/src:${PATH}"

WORKDIR /home/user

# Default command

CMD ["/bin/bash"]

