#!/bin/bash

exit 0

############################

I ended up making a static binary, so I did not have to worry about
install issues. I would have to start from a fresh system to know
for sure if I included all the packages needed.

sudo apt-get install build-essential
sudo apt-get install yasm
sudo apt-get install autoconf autotools-dev libtool texinfo

git clone git://github.com/wbhart/mpir.git
cd mpir
./autogen.sh
./configure --enable-cxx --build=skylakeavx-pc-linux-gnu
cd ..

gcc -static -Wall -O3 -I./mpir ./mpir/.libs/libmpir.a -o run.sh vdf.c


