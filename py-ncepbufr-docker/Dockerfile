# Base image is CentOS 7
FROM centos:7

MAINTAINER Jim Fluke <james.fluke@colostate.edu>

# Get the development tools including gfortran
RUN yum -y --setopt=tsflags=nodocs update && \
    yum clean all && \
    yum -y group install "Development Tools" && \
    yum -y install which vim

# Install miniconda
RUN curl -LO \
  http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /usr/local/miniconda -b
ENV PATH=/usr/local/miniconda/bin:${PATH}

# Update miniconda and install packages
RUN conda update -y conda \
  && conda install -y numpy pytz h5py \
  && rm -rf /usr/local/miniconda/pkgs/* && conda clean -afy \
  && find /usr/local/miniconda/ -follow -type f -name '*.a' -delete \
  && find /usr/local/miniconda/ -follow -type f -name '*.pyc' -delete \
  && find /usr/local/miniconda/ -follow -type f -name '*.js.map' -delete


WORKDIR /root

# Do these next steps all at once to make the clean up steps most effective??
# Install py-ncepbufr
# Use the py-ncepbufr setup.py to build the BUFR library and do the install
# This will install py-ncepbufr in /usr/local/miniconda/ to provide universal
# access
RUN curl -L -o py-ncepbufr-master.zip \
      https://github.com/JCSDA/py-ncepbufr/archive/master.zip \
  && unzip py-ncepbufr-master.zip && rm -v py-ncepbufr-master.zip \
  && cd py-ncepbufr-master/ && python setup.py install \
  && find . -name '*.o' | xargs rm -v

# Install the convenience script for compiling BUFR fortran programs
COPY buildf /usr/local/bin/

# Allow for a custom user environment
# The user running the container should use the docker run -e (--env) option
# to take advantage of this. In bash:
# CUSTOM_USER_ENVIRONMENT=<custom env settings>
# docker run -t -i --rm -e CUSTOM_USER_ENVIRONMENT <other options>
RUN echo '[ -n "$CUSTOM_USER_ENVIRONMENT" ] && eval "$CUSTOM_USER_ENVIRONMENT"' >> /root/.bashrc

ENTRYPOINT ["/bin/bash"]
