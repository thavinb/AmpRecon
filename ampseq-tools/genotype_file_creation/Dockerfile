FROM ubuntu:bionic-20191029

# 1. Install Python3.8
RUN apt-get update -y && \
    apt-get install -y python3.8
RUN apt-get install -y python3-pip

# 2. Install Python packages 
RUN pip3 install pandas \
    pip install PyVcf

# 3. Install bftools dependencies
RUN apt-get install gcc -y
RUN apt-get install make -y
RUN apt-get install libbz2-dev -y
RUN apt-get install zlib1g-dev -y
RUN apt-get install libncurses5-dev -y
RUN apt-get install libncursesw5-dev -y
RUN apt-get install liblzma-dev -y
RUN apt-get install -y tar gzip wget xz-utils

# 4. Install bcftools 1.8
RUN cd /tmp && wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2 && tar -C /tmp/  -vxjf /tmp/bcftools-1.8.tar.bz2 && cd bcftools-1.8 && make && cp -r * /usr/local/bin && cp -r * /usr/local/lib && cd /

RUN cd / && exec /bin/bash "$@"

RUN cd / && exec /bin/bash "$@"
