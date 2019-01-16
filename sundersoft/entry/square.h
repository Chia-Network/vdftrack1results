.intel_syntax noprefix

.global square
square:

#first argument is in RDI; points to an aligned array
#function result is in RAX
#uses linux calling convention

PUSH RBX
PUSH RBP
PUSH R12
PUSH R13
PUSH R14
PUSH R15

RDTSC
SHL RDX, 32
OR RAX, RDX
PUSH RAX

MOV R9, RDI;
ADD R9, 0x1020+0x68;

MOV R10, RDI;
ADD R10, 0x3820+0x68;

MOV R11, RDI;
ADD R11, 0x1020+0x68+0x10*0x8*0x20

MOV R8, 50000
bench_loop:

#YMM0 - temporary
#YMM1 - accumulator
#YMM14 - low value to multiply by
#YMM15 - high value to mulitply by
#RDI+0x1400 - start of output vector
#RDI+0x0400 - start of input vector

#no immediates: 4 bytes
#1 byte immediate: 5 bytes
#4 byte immediate: 8 bytes

#8-15: scalars
#0-3: vector cache
#4-7: temporary



#MOV RAX, RDI;
#MOV RBX, 0

MOV RAX, R9
MOV RBX, R10

loop:
.macro fma limit, offset=0, index=0
VPMULUDQ YMM4, YMM0, YMM8;
VMOVDQU YMM0, [RAX-0x68-0x00+0x20*(\index+\offset)];

VPMULUDQ YMM5, YMM1, YMM9;
VMOVDQU YMM1, [RAX-0x68-0x08+0x20*(\index+\offset)];

VPMULUDQ YMM6, YMM2, YMM10;
VMOVDQU YMM2, [RAX-0x68-0x10+0x20*(\index+\offset)];

VPMULUDQ YMM7, YMM3, YMM11;
VMOVDQU YMM3, [RAX-0x68-0x18+0x20*(\index+\offset)];

VPADDQ YMM4, YMM4, YMM5;
VPADDQ YMM6, YMM6, YMM7;
VPADDQ YMM4, YMM4, YMM6;

VPMULUDQ YMM6, YMM0, YMM12;
VPMULUDQ YMM7, YMM1, YMM13;
VPADDQ YMM5, YMM6, YMM7;

VPMULUDQ YMM6, YMM2, YMM14;
VPMULUDQ YMM7, YMM3, YMM15;
VPADDQ YMM6, YMM6, YMM7;

VPADDQ YMM5, YMM5, YMM6;
VPADDQ YMM4, YMM4, YMM5;

VPADDQ YMM4, YMM4,      [RBX-0x68+0x20*(\index+\offset)];
VMOVDQU                 [RBX-0x68+0x20*(\index+\offset)], YMM4;
.if \index-\limit+1
fma \limit, \offset, "(\index+1)"
.endif
.endm

#fma 1

fma 0x1, 0x0000
#fma 0x40, 0x0040
#fma 0x40, 0x0080
#fma 0x40, 0x00c0

ADD RAX, 0x1*0x20
ADD RBX, 0x1*0x20

CMP RAX, R11
JNE loop

DEC R8
JNE bench_loop

RDTSC
SHL RDX, 32
OR RAX, RDX
PUSH RAX

POP R11
POP R10
MOV RAX, R11
SUB RAX, R10

POP R15
POP R14
POP R13
POP R12
POP RBP
POP RBX

RET


#MOVQ  RAX, RDI  # copy first argument to %rax
#    IMULQ RAX, RDI  # multiply it by itself
#    INC RAX
                          # result is already in %rax
#RET               # return to caller
