# Define base image
#FROM pymesh/pymesh
FROM --platform=linux/amd64 pymesh/pymesh

# install necessary dependencies
RUN apt-get update && \
	apt-get install -y wget git unzip cmake vim libgl1-mesa-glx dssp
	
# DOWNLOAD/INSTALL APBS
RUN mkdir /install
WORKDIR /install
RUN git clone https://github.com/Electrostatics/apbs-pdb2pqr
WORKDIR /install/apbs-pdb2pqr
RUN ls
RUN git checkout b3bfeec
RUN git submodule init
RUN git submodule update
RUN ls
RUN cmake -DGET_MSMS=ON apbs
RUN make
RUN make install
RUN cp -r /install/apbs-pdb2pqr/apbs/externals/mesh_routines/msms/msms_i86_64Linux2_2.6.1 /root/msms/
RUN curl https://bootstrap.pypa.io/pip/3.6/get-pip.py -o get-pip.py
RUN python get-pip.py

# INSTALL PDB2PQR
RUN pip install pdb2pqr

# Setup environment variables 
ENV MSMS_BIN /usr/local/bin/msms
ENV APBS_BIN /usr/local/bin/apbs
ENV MULTIVALUE_BIN /usr/local/share/apbs/tools/bin/multivalue
ENV PDB2PQR_BIN /usr/local/bin/pdb2pqr30

# DOWNLOAD reduce (for protonation)
WORKDIR /install
RUN git clone https://github.com/rlabduke/reduce.git
WORKDIR /install/reduce
RUN make install
RUN mkdir -p /install/reduce/build/reduce
WORKDIR /install/reduce/build/reduce
RUN cmake /install/reduce/reduce_src
WORKDIR /install/reduce/reduce_src
RUN make
RUN make install

# Install python libraries
RUN pip3 install matplotlib 
RUN pip3 install ipython Biopython sklearn tensorflow==1.12 networkx open3d==0.8.0.0 dask==1.2.2 packaging
#RUN pip install StrBioInfo

# Install small molecule dependencies
RUN pip3 install rdkit-pypi==2021.9.4
RUN pip3 install ProDy==2.0
RUN pip install openbabel-wheel==3.1.1.7

# Get an updated het_atm dictionary for reduce
WORKDIR /install/reduce
RUN python2 update_het_dict.py
ENV REDUCE_HET_DICT /install/reduce/reduce_wwPDB_het_dict.txt

# Change working directory
#WORKDIR /
WORKDIR /home

# We need to define the command to launch when we are going to run the image.
# We use the keyword 'CMD' to do that.
CMD [ "bash" ]
