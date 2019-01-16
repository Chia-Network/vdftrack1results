const int gcd_head_size=12; //also have an extra 8 limbs
const int gcd_num_matrix_buffer_spills=0;
const int gcd_background_work_per_call=3;
const int gcd_num_iterations_background_work=10;
const int gcd_num_iterations=60;

const int reduce_head_size=12; //same
const int reduce_num_matrix_buffer_spills=2;
const int reduce_background_work_per_call=6;
const int reduce_num_iterations_background_work=5;
const int reduce_num_iterations=30;

const int num_spilled_registers=256;
const int num_immediates=1024;
const int num_asm_ints=128;
const int asm_int_size_limbs=256; //2048 bytes

//only half of the entries are actually read since the most significant bit is almost always 1. each entry is 8 bytes
const int divide_table_index_bits=12;

const int gcd_num_quotient_bits=31; //excludes sign bit
const int reduce_num_quotient_bits=15;

const bool asm_output_common_case_only=false; //for testing if everything fits in the uop cache / i-cache

const bool debug_reduce=false;
const bool test_asm_funcs=true;
const bool asm_output_tags=true;

const bool solve_linear_congruence_test_gcd=false;

#include "bit_manipulation.h"
#include "double_utility.h"
#include "integer.h"

#include "simd_integer.h"

#include "asm_types.h"
#include "asm_vm.h"
#include "simd_integer_fma_asm.h"
#include "simd_integer_gcd_asm.h"
#include "simd_integer_reduce_asm.h"

#include "asm_main.h"