#!/bin/sh
export DEBIAN_FRONTEND=noninteractive &&
sudo apt-get install -y texinfo &&
sudo apt-get install -y cmake &&
sudo apt-get install -y m4 &&
cd ./gmp-source &&
sudo ./configure --prefix=/usr --enable-cxx && 
sudo make &&
cp ./.libs/libgmp.a  ./../  &&
cp ./.libs/libgmpxx.a  ./../  &&
cd ./../  &&
mkdir build
cd build/ &&
cmake -DCMAKE_BUILD_TYPE=Release ./../ &&
cmake --build . &&
cp ./bin/a.out ./../ &&
cd ./../ 
 
