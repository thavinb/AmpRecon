FROM bashell/alpine-bash:latest

#Get perl
RUN apk update \
 && apk upgrade \
 && apk add perl
 
#Get GCC
RUN apk update \
 && apk upgrade \
 && apk add g++ \
 && apk add make

#Get git
RUN apk update \
 && apk upgrade \
 && apk add git
 
#Get htslib requirements
RUN apk update \
 && apk upgrade \
 && apk add autoconf \
 && apk add zlib \ 
 && apk add zlib-dev  \   
 && apk add libbz2 \  
 && apk add xz-dev \ 
 && apk add libcurl  \
 && apk add bzip2-dev \
 && apk add libtool \
 && apk add intltool \
 && apk add autoconf \
 && apk add libc-dev \
 && apk add curl-dev \ 
 && apk add automake

#Build htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
 && tar -C /tmp/ -xjf htslib-1.9.tar.bz2 \
 && rm htslib-1.9.tar.bz2

WORKDIR /tmp/htslib-1.9

RUN autoheader \
 && autoconf \
 && ./configure --prefix="/usr/local/bin/htslib/1.9/" \
 && make \
 && make install prefix="/usr/local/bin/htslib/1.9/"
 
#Get bambi requirements
RUN apk update \
 && apk upgrade \
 && apk add libxml2-dev \
 && apk add gd-dev
 
#build bambi

WORKDIR /tmp/

RUN git clone https://github.com/wtsi-npg/bambi.git 

WORKDIR /tmp/bambi

RUN git checkout tags/0.12.1 

RUN autoreconf -i 

RUN ./configure --prefix="/usr/local/bin/bambi" --with-htslib="/usr/local/bin/htslib/1.9/" LDFLAGS=-lz \
 && make install prefix="/usr/local/bin/bambi"
 
#RUN make check <-- requires samtools

#clean up
RUN rm -r /tmp/htslib-1.9

RUN rm -r /tmp/bambi
 
ENV PATH=$PATH:/usr/local/bin/bambi/bin

ENTRYPOINT ["bambi"]


 


