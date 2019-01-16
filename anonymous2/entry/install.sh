#!/bin/bash
sudo apt-get install build-essential;
gcc vdf.c -w -O3 -I../mpir libmpir.a --static -o vdf.exe;
