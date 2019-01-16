#!/bin/sh
echo "
Install
======
Not required as this is a staticly linked binary.


Instructions to Build
====================

PARI scripting, compiled with MPIR optimized library.

1. Configured MPIR with "./configure --with-pic --build=skylakeavx-pc-linux-gnu"
2. Compiled MPIR with "make"
3. Tuned MPIR with -f 1000000
4. Recompile MPIR and copy MpirRoot/.libs/libmpir.a into GP root PariRoot and PariRoot/Olinux-x86_64/
5. Statically linked MPIR to GP by setting "GMPLIBS      = libmpir.a" in PariRoot/Olinux-x86_64/Makefile
6. Compile GP with MPIR using "make all"
"
