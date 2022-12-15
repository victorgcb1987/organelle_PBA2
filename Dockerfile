# syntax=docker/dockerfile:1
FROM python:3.9-bullseye
RUN apt-get update --assume-yes && apt-get upgrade --assume-yes
RUN apt install curl --assume-yes
RUN apt install bzip2 --assume-yes
RUN apt install cmake --assume-yes
RUN apt install default-jdk --assume-yes
RUN apt install install default-jre
COPY . /organelle_pba2
WORKDIR /organelle_pba2
RUN chmod +x /organelle_pba2/orgpba2/assemble_organelle.py
RUN chmod +x /organelle_pba2/orgpba2/calculate_heteroplasmy.py
RUN python setup.py install
RUN mkdir dependencies
WORKDIR /dependencies/
RUN apt install ncbi-blast+ --assume-yes
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 --output minimap2.tar.bz2
RUN tar -jxvf minimap2.tar.bz2 
RUN rm minimap2.tar.bz2 
WORKDIR /dependencies
RUN git clone https://github.com/lh3/seqtk.git
WORKDIR /dependencies/seqtk
RUN make
WORKDIR /dependencies/
RUN git clone https://github.com/rrwick/Filtlong.git
WORKDIR /dependencies/Filtlong
RUN make -j
WORKDIR /dependencies/
RUN curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz --output canu-2.2.Linux.tar.xz
RUN tar -xJf canu-2.2.Linux.tar.xz
RUN rm canu-2.2.Linux.tar.xz
RUN git clone --recursive https://github.com/lbcb-sci/racon.git racon
WORKDIR /dependencies/racon
RUN mkdir build
WORKDIR /dependencies/racon/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j
WORKDIR /
RUN pip install --upgrade pip
RUN pip install biopython
ENV MINIMAP2_PATH="/dependencies/minimap2-2.24_x64-linux"
ENV CANU_PATH="/dependencies/canu-2.2/bin"
ENV SEQTK_PATH="/dependencies/seqtk"
ENV FILTLONG_PATH="/dependencies/Filtlong/bin"
ENV RACON_PATH="/dependencies/racon/build/bin"
ENV PATH="${PATH}:/organelle_pba2/orgpba2"