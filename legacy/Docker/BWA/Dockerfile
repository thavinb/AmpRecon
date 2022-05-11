FROM bashell/alpine-bash:latest

#Get GCC
RUN apk update \
 && apk upgrade \
 && apk add g++ \
 && apk add make

#Get git
RUN apk update \
 && apk upgrade \
 && apk add git

RUN apk update \
 && apk upgrade \
 && apk add libxml2-dev \
 && apk add gd-dev \
 && apk add gzip

WORKDIR /

RUN git clone https://github.com/lh3/bwa.git

WORKDIR /bwa

RUN git checkout tags/v0.7.17

RUN make

ENV PATH=$PATH:/bwa
