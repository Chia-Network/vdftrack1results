#include "include.h"
#include "integer.h"
//#include "vdf_new.h"

#include "gpu_integer.h"
#include "gpu_integer_divide.h"

#include <cuda_runtime.h>

const int gpu_integer_size=32;
typedef gpu_integer<gpu_integer_size, false> gpu_integer_type;

__global__ void foo(gpu_integer_type* a, gpu_integer_type* b) {
    //int thread_number=blockIdx.x*blockDim.x + threadIdx.x;
    //a[pos]=pos;

    gpu_integer_type a_copy=*(gpu_integer_type*)a;
    gpu_integer_type b_copy=*(gpu_integer_type*)b;

    //decltype(&foo1) foos[]={&foo1, &foo2, &foo3, &foo4};

    gpu_integer_type res_a;
    res_a=a_copy/b_copy;

    a_copy=res_a;
    //b_copy=res_b;

    /*#pragma unroll 1
    for (int x=0;x<100;++x) {
        int num_limbs=0;
        #pragma unroll
        0
        a[0]=num_limbs;
    }*/

    *(gpu_integer_type*)a=a_copy;
}

/*__global__ void foo2(uint64* a, uint64* b) {
    uint64 c=(*a) / *b;
    uint64 d=(*a) % *b;
    
    *a=c;
    *b=d;

    //uint64 v=*(uint64*)a;
    //++v;
    //double w=*(double*)&v;

    //asm volatile( "rcp.rm.f64 %0, %0;" : "+d"(w) );
    //*a=w;

    //*a=1/w;
}**/

int main() {
    int* foo_buffer;
    const int foo_buffer_size=1024;
    error_check(cudaMallocManaged(&foo_buffer, foo_buffer_size*sizeof(foo_buffer[0])));

    gpu_integer_type int_a;
    gpu_integer_type int_b;

    foo<<<1, 1>>>(&int_a, &int_b); error_check(cudaGetLastError());
    //foo2<<<1, 1>>>((uint64*)int_a.data, (uint64*)int_b.data); error_check(cudaGetLastError());

    for (int x=0;x<foo_buffer_size;++x) {
        print(foo_buffer[x]);
    }
}