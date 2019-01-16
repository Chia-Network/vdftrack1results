#!/bin/sh
sudo apt-get install libgmp3-dev
g++ -std=c++17 -O3 vdf4.cpp -lgmpxx -lgmp
