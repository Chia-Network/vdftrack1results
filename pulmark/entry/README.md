# Entry for Chia Network VDF Competition

## Introduction
The program implements VDF based on reference implementation. More information about the competition and reference samle code is available at Github: [Chia Network VDF Competition](https://github.com/Chia-Network/vdf-competition).

To build the program, run install.sh script that installs all dependencies and then compiles the code. 

The run.sh script executes the CLI program and it takes following mandatory arguments to calculate VDF:
* discriminent (hex format)
* number of iterations (decimal)

The script should output the result of the VDF, encoded as a, b of the final classgroup element.

## Dependencies
The program depends on [GNU Multi Precision library](https://gmplib.org) and uses it for big integer arithmetic. GMP is distributed under the dual licenses, GNU LGPL v3 and GNU GPL v2. These licenses make the library free to use, share, and improve, and allows user to pass on the result.

A custom binary for GMP library was created by using [Mercurial GMP repository](https://gmplib.org/repo/) and [build options](https://gmplib.org/manual/Build-Options.html#Build-Options). This library is optimized for competition by supporting the CPU instruction set and [Coffee Lake microarchitecture](https://en.wikichip.org/wiki/intel/microarchitectures/coffee_lake#Compiler_support) of Chia Network's test computer (CPU Intel i5-8400).

***WARNING***: Because of these CPU specific optimizations running this entry might crash on computers having different CPU than test computer.

## Description
The implementation of VDF is written in C/C++. 

The classgroup element is implemented as header-only in header file form. Header file defines C++ structures and classes for VDF functionality.

The run script executes the main program that parses input arguments from console terminal, inits the initial element and then uses methods implemented in header file for VDF calculation. After VDF calculation, results are displayed in console terminal.

## Approach
I'm just programmer (M.Sc. Electrical Engineering, 1989), started programming by using Fortran, PL/M and assemby languages at [University of Oulu](http://www.oulu.fi/university/) in 1982. Therefore I focused only into the reference sample code and tried to improve its performance.

The implementation follows the guidelines and principles presented in article [Optimizing Software in C++](https://www.agner.org/optimize/optimizing_cpp.pdf) written by Agner Fog.

#### Overview
* Fast big integer libraries exist already, so no need to reinvent the wheel.
* By using same library for big integer arithmetic as in the reference implementation allows rapid prototyping, verification of code changes and testing.
* Try to minimize use of big integer multiply and divide arithmetic operations that takes most time to compute.
* Optimize CPU memory usage:
    * try to avoid memory reallocations,
    * try to avoid CPU cache misses by keeping variables and functions close together.
* Optimize code for CPU parallel computing.

#### Studies, Improvements, Tests
* A code review of reference implementation and study of GMP API has been carried out to identify issues that could improve VDF performance.
* Replaced global variables and methods by defining C++ structures in separate header file.
* Running the VDF repeated squaring in a separate thread.
* Instrumented the code by using Linux timestamp counters followed by performance analysis to find the critical sections of code where program spends most time.
* Analyzed the GMP API calls utilized in reference implementation and replaced several API calls with more efficient ones.  
* Minimized temporary operations and variables.
* Changed the initialization of big integer variables to avoid memory reallocations.
* Static linking of GMP library for better performance.
* Simple optimizations into reduction algorithm, that uses around 86% of total execution time:
    * less normalize function calls,
    * optimized code in the inner loop of reduction algorithm by using less GMP API function calls compared to reference implementation

Compared to GMP library, following other external libraries where tested for big integer arithmetic:
* [Boost multi precision library](https://www.boost.org/doc/libs/1_66_0/libs/multiprecision/doc/html/index.html) for big integers: The performance was much worse than using GMP.
* [MPIR](http://mpir.org) for big integers: Using GMP API creates faster executable than MPIR.
* [asmlib](http://www.agner.org/optimize/asmlib.zip) by Agner Fog for basic memory and string operations: No performance improvements were noticed.

#### Performance
Currently the entry time is 113s for the test run in Chia Network's reference machine. The advantage compared to reference sample code is around 50-55sec or 32% improvement. The performance calculated for 10000x squaring iterations is 539ms. Compared to ASIC 2048bit RSA squaring implementation there is some discussion here at [ASIC latency performance for 2048-bit RSA squaring](https://github.com/Chia-Network/vdf-competition/issues/8).

During development, Linux [perf](https://perf.wiki.kernel.org) kernel tool was utilized to verify how code changes affected to CPU and memory usage. Outputs of performance counter stats for a test run are presented below for sample code entry and this entry using equal input arguments in my development machine for comparison.

**Chia Network's sample code entry performance stats**:


         240353.61 msec task-clock:u              #    1.000 CPUs utilized          
                 0      context-switches:u        #    0.000 K/sec                  
                 0      cpu-migrations:u          #    0.000 K/sec                  
               335      page-faults:u             #    1.394 M/sec                  
      809796241038      cycles:u                  # 3369195.479 GHz                 
      177404488196      stalled-cycles-frontend:u #   21.91% frontend cycles idle   
       98548189343      stalled-cycles-backend:u  #   12.17% backend cycles idle    
     1659894323306      instructions:u            #    2.05  insn per cycle         
                                                  #    0.11  stalled cycles per insn
      203780990636      branches:u                # 847840429.019 M/sec             
        5580050275      branch-misses:u           #    2.74% of all branches        

     240.406184670 seconds time elapsed

     240.340402000 seconds user
       0.015996000 seconds sys

**This entry performance stats**:


         165896.19 msec task-clock:u              #    0.999 CPUs utilized          
                 0      context-switches:u        #    0.000 K/sec                  
                 0      cpu-migrations:u          #    0.000 K/sec                  
               609      page-faults:u             #    3.671 M/sec                  
      561277808252      cycles:u                  # 3383311.281 GHz                 
      106533309706      stalled-cycles-frontend:u #   18.98% frontend cycles idle   
       67275130580      stalled-cycles-backend:u  #   11.99% backend cycles idle    
     1182594025347      instructions:u            #    2.11  insn per cycle         
                                                  #    0.09  stalled cycles per insn
      138719594017      branches:u                # 836184079.285 M/sec             
        3545242946      branch-misses:u           #    2.56% of all branches        

     165.979492988 seconds time elapsed

     165.878932000 seconds user
       0.017370000 seconds sys
       
The instructions per cycle(IPC) count: 2.11 is good. The stalled cycles means that CPU has to wait something, the front end stalls are more important, code cannot run fast if the instruction stream is not being kept up. Low numbers for stalls means good performance.

#### Further Improvements
Use of GPU for arithmetic operations should provide far better performance than CPU based approach. Due to tight competition schedule, it was not possible to use this option for competition. GPU option is very interesting, so considering to do some prototyping and testing with this option.

Off-loading reduction algorithm into GPU would be a good exercise. Reduction uses add, shift, substract, multiply, divide, compare and swap operations for big integer arithmetic. Unless GPU doesn't have build-in type for big integers or there is no tool to cross-compile current C++ code to be run on GPU, a big integer implementation is needed to share data and memory between CPU and GPU.

The following project report: [Large Integer Arithmetic in GPU for Cryptography](http://eprints.utar.edu.my/2494/1/CS-2017-1401837-1.pdf) describes efficient algorithms for big integer GPU implementation.

## License
Copyright &copy; [Markku Pulkkinen](https://keybase.io/pulmark). Released under the [Apache Version 2.0 License](https://www.apache.org/licenses/LICENSE-2.0).
