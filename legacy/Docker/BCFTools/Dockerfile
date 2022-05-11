FROM bashell/alpine-bash:latest

#Get git
RUN apk update \
 && apk upgrade \
 && apk add git

 #Get GCC
RUN apk update \
 && apk upgrade \
 && apk add g++ \
 && apk add make

 #Get bcf tool requirements
 RUN apk update \
  && apk upgrade \
  && apk add zlib \
  && apk add zlib-dev \
  && apk add bzip2-dev \
  && apk add xz-dev \
  && apk add curl-dev

#Build htslib
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
 && tar -C /tmp/ -xjf bcftools-1.10.2.tar.bz2 \
 && rm bcftools-1.10.2.tar.bz2

WORKDIR /tmp/bcftools-1.10.2

#build bcftools
RUN ./configure --prefix="/usr/local/bin/bcftools"  \
 && make \
 && make install

ENV PATH=$PATH:/usr/local/bin/bcftools/bin

ENTRYPOINT ["bcftools"]
