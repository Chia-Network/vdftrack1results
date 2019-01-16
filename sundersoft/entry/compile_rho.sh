#!/bin/bash
nvcc rho.cu -o rho --maxrregcount 64 -O3 -g -w --std c++14 --gpu-architecture sm_61 --resource-usage -lgmp -lgmpxx
