FROM ubuntu:bionic-20191029

# 1. Install Python3.8
RUN apt-get update -y && \
    apt-get install -y python3.8
RUN apt-get install -y python3-pip

# 2. Install Python packages 
RUN pip install PyVcf
