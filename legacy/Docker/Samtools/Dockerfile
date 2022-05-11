FROM bashell/alpine-bash:latest

ENV SAMTOOLS_VERSION 1.9

RUN apk update \
 && apk upgrade \
 && apk add --no-cache \
 build-base \
 zlib-dev \
 bzip2-dev \
 xz-dev \
 ncurses-dev \
 ca-certificates \
 wget

RUN wget -q https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2

RUN mkdir /samtools-1.9

RUN tar -C / -xjvf samtools-1.9.tar.bz2

WORKDIR /samtools-1.9/

RUN make

RUN mv /samtools-1.9/samtools /bin/ \
 && rm -rf /samtools*

RUN apk del build-base zlib-dev ca-certificates wget

RUN apk update \
 && apk upgrade \
 && apk add --no-cache python3

ENTRYPOINT samtools
