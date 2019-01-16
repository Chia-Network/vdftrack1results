#!/bin/sh
###############################################################################
# Installation script using GNU g++ compiler
# g++ version: g++ (SUSE Linux) 8.2.1 20181108 [gcc-8-branch revision 265914]
#
# Dependencies: libgmp
#
# g++ --version
# g++ -march=native -Q --help=target | egrep -- '-march=|-mtune' | cut -f3
###############################################################################

# Build-1: Basic set of compiler options
###############################################################################
#sudo apt-get install libgmp3-dev
g++ -o a.out -no-pie -pthread -pipe -Wall -O3 -march=native -mtune=native \
-fprefetch-loop-arrays --param max-inline-insns-single=2000 -std=c++1y \
main.cpp -Wl,-Bstatic -L./lib -l:libgmp.a -Wl,-Bdynamic
