FROM bashell/alpine-bash:latest

WORKDIR /tmp

RUN apk update \
 && apk upgrade \
 && apk add --no-cache \
 g++ \
 bzip2-dev \
 xz-dev \
 libbz2 \
 zlib-dev \
make

RUN wget -q https://github.com/jkbonfield/io_lib/releases/download/io_lib-1-14-11/io_lib-1.14.11.tar.gz

RUN tar -xvzf io_lib-1.14.11.tar.gz \
&& rm /tmp/io_lib-1.14.11.tar.gz

WORKDIR /tmp/io_lib-1.14.11

RUN ./configure

RUN make && make install

WORKDIR /

RUN rm -r /tmp/io_lib-1.14.11
