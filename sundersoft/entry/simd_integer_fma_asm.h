namespace simd_integer_namespace {


struct simd_integer_asm {
    int start=-1; //multiple of 4, times 8 bytes
    int size=-1; //multiple of 4
    reg_scalar memory_base_reg; //32-aligned

    void assert_valid() {
        assert(start>=0 && size>0 && start%4==0 && size%4==0);
    }

    bool aliased_with(simd_integer_asm& t) {
        assert_valid();

        int s1=start;
        int e1=start+size-1;
        int s2=t.start;
        int e2=t.start+t.size-1;

        return
            (s1>=s2 && s1<=e2) || (e1>=s2 && e1<=e2) ||
            (s2>=s1 && s2<=e1) || (e2>=s1 && e2<=e1)
        ;
    }

    int current_size() {
        return size;
    }

    string operator[](int index) {
        assert(index>=-3 && index<size);
        return asm_memory(memory_base_reg, (start+index)*8);
    }

    //regs: 1x vector
    void is_negative(reg_alloc regs, string branch_negative, string branch_positive) {
        EXPAND_MACROS_SCOPE;

        reg_vector data=regs.bind_vector(m, "data");

        APPEND_M(str( "VMOVDQU `data, #", (*this)[size-4] ));
        APPEND_M(str( "VPTEST `data, #", asm_immediate({0, 0, 0, 1ull<<(data_size-1)}) )); //ZF=1 if positive

        assert(branch_negative.empty()!=branch_positive.empty());

        if (!branch_negative.empty()) {
            APPEND_M(str( "JNZ #", branch_negative ));
        }

        if (!branch_positive.empty()) {
            APPEND_M(str( "JZ #", branch_positive ));
        }
    }

    //regs: 2x vector, 2x vector argument
    void calculate_carry_loop(
        reg_alloc regs, reg_vector carry_mask_reg, reg_vector accumulator,
        int start_index, bool single_pass,
        int range_start, int range_end, bool preserve_accumulator
    ) {
        EXPAND_MACROS_SCOPE;

        reg_vector data=regs.bind_vector(m, "data");
        reg_vector carry=regs.bind_vector(m, "carry");

        if (!single_pass) {
            m.bind(accumulator, "accumulator");
        }

        m.bind(carry_mask_reg, "carry_mask");

        assert_valid();
        assert(start_index>=0 && start_index%4==0);
        assert(start_index<size);

        assert(range_start>=start_index && range_start%4==0);
        assert(range_end>range_start && range_end<=size && range_end%4==0);

        //need to go msb first since carries are cleared at the same time the output is written
        for (int index=range_end-4;index>=range_start;index-=4) {
            bool first=(index==size-4);

            APPEND_M(str( "VMOVDQU `carry, #", (*this)[index-1] ));
            APPEND_M(str( "VPSRAD `data, `carry, #", to_hex(data_size) )); //data used as temporary
            APPEND_M(str( "VPSRLQ `carry, `carry, #", to_hex(data_size) ));
            APPEND_M(str( "VPBLENDD `carry, `carry, `data, #", vpblendd_mask_8({0,1, 0,1, 0,1, 0,1}) )); //arithmetic right shift

            if (index==start_index) {
                //need to zero out the lsb carry so that it won't get added
                asm_immediate.assign(data, 0); //temporary
                APPEND_M(str( "VPBLENDD `carry, `data, `carry, #", vpblendd_mask_4({0, 1, 1, 1}) ));
            }

            //negates carry_mask to get data_mask
            APPEND_M(str( "VPANDN `data, `carry_mask, #", (*this)[index] ));

            string out_reg=(!single_pass && first && !preserve_accumulator)? "`accumulator" : "`data";

            APPEND_M(str( "VPADDQ #, `data, `carry", out_reg ));

            APPEND_M(str( "VMOVDQU #, #", (*this)[index], out_reg ));

            if (!single_pass && (!first || preserve_accumulator)) {
                APPEND_M(str( "VPOR `accumulator, `accumulator, #", out_reg ));
            }
        }
    }

    //regs: 3/4x vector (3 if single_pass==true)
    void calculate_carry(
        reg_alloc regs, int start_index=0, bool single_pass=false
    ) {
        EXPAND_MACROS_SCOPE;

        reg_vector carry_mask_reg=regs.bind_vector(m, "carry_mask");
        reg_vector accumulator;

        if (!single_pass) {
            accumulator=regs.bind_vector(m, "accumulator");
        }

        asm_immediate.assign(carry_mask_reg, carry_mask);

        string label_name;
        if (!single_pass) {
            label_name=m.alloc_label();
            APPEND_M(str( "#:", label_name ));
        }

        calculate_carry_loop(
            regs, carry_mask_reg, accumulator,
            start_index, single_pass,
            start_index, size, false
        );

        if (!single_pass) {
            APPEND_M(str( "VPTEST `carry_mask, `accumulator" )); //ZF = (carry_mask & accumulator)==0
            APPEND_M(str( "JNZ #", label_name ));
        }
    }

    /*void calculate_carry_lsb(reg_alloc regs) {
        EXPAND_MACROS_SCOPE;

        assert_valid();

        reg_scalar imcoming_carry=regs.bind_scalar(m, "incoming_carry");
        reg_scalar tmp=regs.bind_scalar(m, "tmp");

        asm_immediate.assign(incoming_carry, 0);

        for (int x=0;x<memory.size();++x) {
            //uint64& c=memory.at(x);
            APPEND_M(str( "MOV `tmp, #", (*this)[x] ));

            //c+=incoming_carry;
            APPEND_M(str( "ADD `tmp, `incoming_carry" ));

            //incoming_carry=sign_extend_carry(c);
            //return int64(v) >> (64-carry_size);
            APPEND_M(str( "MOV `incoming_carry, `tmp" ));
            APPEND_M(str( "SAR `incoming_carry, #", to_hex(64-carry_size) ));

            //c&=data_mask;
            APPEND_M(str( "AND `tmp, #", to_hex(data_mask) ));
            APPEND_M(str( "MOV #, `tmp", (*this)[x] ));
        }
    }

    //v should be 1 or -1
    void multiply_one(reg_alloc regs, reg_scalar v, bool skip_add=false) {
        EXPAND_MACROS_SCOPE;

        assert_valid();

        m.bind(v, "v");
        reg_scalar tmp=regs.bind_scalar(m, "tmp");
        reg_scalar tmpb=regs.bind_scalar(m, "tmpb");

        reg_vector tmp_1=regs.bind_vector(m, "tmp_1");
        reg_vector tmp_2=regs.bind_vector(m, "tmp_2");

        //uint64 mask=(v>>1); //-1 if negative, 0 if positive
        APPEND_M(str( "MOV `tmp, `v" ));
        APPEND_M(str( "SAR `tmp, 0x1" ));

        //mask&=data_mask;
        APPEND_M(str( "AND `tmp, #", to_hex(data_mask) ));

        APPEND_M(str( "VMOVQ `tmp_1_128, `tmp" ));
        APPEND_M(str( "VPBROADCASTQ `tmp_1, `tmp_1_128" ));

        //flip all of the bits if negative
        for (int x=0;x<size;x+=4) {
            //memory[x]^=mask;
            APPEND_M(str( "VPXOR `tmp_2, `tmp_1, #", (*this)[x] ));
            APPEND_M(str( "VMOVDQU #, `tmp_2", (*this)[x] ));
        }

        if (!skip_add) {
            //memory[0]+=uint64(v) >> 63; //1 if negative, 0 if positive
            APPEND_M(str( "MOV `tmp, `v" ));
            APPEND_M(str( "SAR `tmp, #", to_hex(63) ));

            APPEND_M(str( "MOV `tmpb, #", (*this)[0] ));
            APPEND_M(str( "ADD `tmpb, `tmp" ));
            APPEND_M(str( "MOV #, `tmpb", (*this)[0] ));

            //if ((memory[0] & carry_mask)!=0) {
            string skip_carry_label=m.alloc_label();

            APPEND_M(str( "TEST #, `tmpb", asm_immediate(carry_mask) ));
            APPEND_M(str( "JZ", skip_carry_label ));

            calculate_carry(regs);

            APPEND_M(str( "#:", skip_carry_label ));
        }
    }*/

    //regs: 1x vector
    void copy(reg_alloc regs, simd_integer_asm c) {
        EXPAND_MACROS_SCOPE;

        reg_vector data=regs.bind_vector(m, "data");

        assert_valid();
        c.assert_valid();
        assert(size==c.size);
        assert(!aliased_with(c));

        for (int x=0;x<size;x+=4) {
            APPEND_M(str( "VMOVDQU `data, #", (*this)[x] ));
            APPEND_M(str( "VMOVDQU #, `data", (*this)[x] ));
        }
    }

    //regs: 2x vector
    void clear(reg_alloc regs, int start_index=0) {
        EXPAND_MACROS_SCOPE;

        reg_vector data=regs.bind_vector(m, "data");
        reg_vector zero=regs.bind_vector(m, "zero");

        assert_valid();
        assert(start_index>=0 && start_index<size);

        asm_immediate.assign(zero, 0);

        if (start_index%4!=0) {
            int i=round_down(start_index, 4);
            APPEND_M(str( "VMOVDQU `data, #", (*this)[i] ));

            array<int, 4> mask;
            for (int x=0;x<4;++x) {
                mask[x]=(x+i>=start_index)? 0 : 1;
            }
            assert((mask!=array<int, 4>({0, 0, 0, 0})));
            assert((mask!=array<int, 4>({1, 1, 1, 1})));

            APPEND_M(str( "VPBLENDD `data, `zero, `data, #", vpblendd_mask_4(mask) ));
            APPEND_M(str( "VMOVDQU #, `data", (*this)[i] ));
        }

        for (int x=round_up(start_index, 4);x<size;x+=4) {
            APPEND_M(str( "VMOVDQU #, `zero", (*this)[x] ));
        }
    }

    //regs: 3x vector
    //should have no carry at the limb before partial_index and every limb after. values at and after partial_index must also be known
    void check_partial_carry_valid(reg_alloc regs, int partial_index) {
        EXPAND_MACROS_SCOPE;

        reg_vector data=regs.bind_vector(m, "data");
        reg_vector cmp_1=regs.bind_vector(m, "cmp_1");
        reg_vector cmp_2=regs.bind_vector(m, "cmp_2");

        assert_valid();

        uint64 min_v=1ull<<(carry_size-data_size-1);
        uint64 max_v=data_mask-(min_v+1);

        //uint64 v=memory.at(carry_start+1);
        //valid: v>=min_v && v<=max_v

        assert(partial_index>=0);
        int load_index=round_down(partial_index, 4);
        array<uint64, 4> mask;
        for (int x=0;x<4;++x) {
            mask[x]=(x+load_index==partial_index)? 1 : 0;
        }

        APPEND_M(str( "VMOVDQU `data, #", (*this)[load_index] ));
        APPEND_M(str( "VPCMPGTQ `cmp_1, `data, #", asm_immediate(min_v-1) )); // data>min_v-1 (data>=min_v)
        APPEND_M(str( "VPCMPGTQ `cmp_2, `data, #", asm_immediate(max_v) )); // data>max_v (!data<=max_v)
        APPEND_M(str( "VPANDN `cmp_1, `cmp_2, `cmp_1" )); // !(data>max_v) && (data>min_v-1)
        APPEND_M(str( "VPTEST `cmp_1, #", asm_immediate(mask) ));

        APPEND_M(str( "JZ #", m.alloc_error_label() ));
    }

    //regs: 4x vector + 1x vector argument
    void logical_shift_left(reg_alloc regs, simd_integer_asm c, reg_vector bits) {
        EXPAND_MACROS_SCOPE;

        assert(c.size==size);

        m.bind(bits, "bits");

        reg_vector this_data=regs.bind_vector(m, "this_data");
        reg_vector previous_data=regs.bind_vector(m, "previous_data");
        reg_vector data_mask_reg=regs.bind_vector(m, "data_mask");
        reg_vector data_size_minus_bits=regs.bind_vector(m, "data_size_minus_bits");

        asm_immediate.assign(data_mask_reg, data_mask);

        asm_immediate.assign(data_size_minus_bits, data_size);
        APPEND_M(str( "VPSUBQ `data_size_minus_bits, `data_size_minus_bits, `bits" ));

        assert_valid();

        //msb first
        for (int index=size-4;index>=0;index-=4) {
            APPEND_M(str( "VMOVDQU `previous_data, #", c[index-1] ));

            if (index==0) {
                asm_immediate.assign(this_data, 0); //temporary
                APPEND_M(str( "VPBLENDD `previous_data, `this_data, `previous_data, #", vpblendd_mask_4({0, 1, 1, 1}) ));
            }

            APPEND_M(str( "VMOVDQU `this_data, #", c[index] ));

            APPEND_M(str( "VPSLLQ `this_data, `this_data, `bits_128" ));
            APPEND_M(str( "VPAND `this_data, `this_data, `data_mask" ));

            APPEND_M(str( "VPSRLQ `previous_data, `previous_data, `data_size_minus_bits_128" ));

            APPEND_M(str( "VPOR `this_data, `this_data, `previous_data" ));

            APPEND_M(str( "VMOVDQU #, `this_data", (*this)[index] ));
        }
    }

    //regs: 4x vector + 1x vector argument
    void logical_shift_right(reg_alloc regs, reg_vector bits) {
        EXPAND_MACROS_SCOPE;

        m.bind(bits, "bits");

        reg_vector this_data=regs.bind_vector(m, "this_data");
        reg_vector next_data=regs.bind_vector(m, "next_data");
        reg_vector data_mask_reg=regs.bind_vector(m, "data_mask");
        reg_vector data_size_minus_bits=regs.bind_vector(m, "data_size_minus_bits");

        asm_immediate.assign(data_mask_reg, data_mask);

        asm_immediate.assign(data_size_minus_bits, data_size);
        APPEND_M(str( "VPSUBQ `data_size_minus_bits, `data_size_minus_bits, `bits" ));

        assert_valid();

        //lsb first
        for (int index=0;index<size;index+=4) {
            APPEND_M(str( "VMOVDQU `next_data, #", (*this)[index+1] ));

            if (index==size-4) {
                asm_immediate.assign(this_data, 0); //temporary
                APPEND_M(str( "VPBLENDD `next_data, `this_data, `next_data, #", vpblendd_mask_4({1, 1, 1, 0}) ));
            }

            APPEND_M(str( "VMOVDQU `this_data, #", (*this)[index] ));

            APPEND_M(str( "VPSRLQ `this_data, `this_data, `bits_128" ));

            APPEND_M(str( "VPSLLQ `next_data, `next_data, `data_size_minus_bits_128" ));
            APPEND_M(str( "VPAND `next_data, `next_data, `data_mask" ));

            APPEND_M(str( "VPOR `this_data, `this_data, `next_data" ));

            APPEND_M(str( "VMOVDQU #, `this_data", (*this)[index] ));
        }
    }

    //regs: 2x vector + 1x vector argument
    void fma(
        reg_alloc regs, simd_integer_asm c, simd_integer_asm a, reg_vector b, int shift_amount
    ) {
        EXPAND_MACROS_SCOPE;

        m.bind(b, "b");

        reg_vector a_reg=regs.bind_vector(m, "a");
        reg_vector zero=regs.bind_vector(m, "zero");

        assert_valid();
        a.assert_valid();
        if (c.size!=-1) {
            assert(c.size==size);
            c.assert_valid();

            if (c.start!=start) {
                assert(shift_amount==0 && a.size==size);
            }
        }

        assert(shift_amount>=0);
        if (shift_amount!=0) {
            assert(!aliased_with(a));
        }

        int output_start=round_down(shift_amount, 4);
        int output_end=round_up(a.size+shift_amount, 4);

        if (output_end>size) {
            output_end=size;
        }

        bool assigned_zero=false;

        //this is slow because it needs to process 8 b values at a time to get close to the theoretical performance
        //(only vector multiplications can actually use multiple b values at once though)
        for (int index=output_end-4;index>=output_start;index-=4) {
            int a_pos=index-shift_amount;
            assert(a_pos>=-3 && a_pos<a.size+3);

            APPEND_M(str( "VPMULDQ `a, `b, #", a[a_pos] ));

            array<int, 4> mask;
            for (int x=0;x<4;++x) {
                int i=a_pos+x;
                mask[x]=(i>=0 && i<a.size)? 1 : 0;
            }

            if (mask!=array<int, 4>({1, 1, 1, 1})) {
                if (!assigned_zero) {
                    asm_immediate.assign(zero, 0);
                    assigned_zero=true;
                }
                APPEND_M(str( "VPBLENDD `a, `zero, `a, #", vpblendd_mask_4(mask) ));
            }

            if (c.size!=-1) {
                APPEND_M(str( "VPADDQ `a, `a, #", c[index] ));
            }

            APPEND_M(str( "VMOVDQU #, `a", (*this)[index] ));
        }
    }

    //regs: 1x vector
    void zero_extend_data(reg_alloc regs) {
        EXPAND_MACROS_SCOPE;

        reg_vector data=regs.bind_vector(m, "data");

        assert_valid();

        APPEND_M(str( "VMOVDQU `data, #", (*this)[size-4] ));
        APPEND_M(str( "VPAND `data, `data, #", asm_immediate(data_mask) ));
        APPEND_M(str( "VMOVDQU #, `data", (*this)[size-4] ));
    }

    //regs: 2x vector
    void sign_extend_data(reg_alloc regs) {
        EXPAND_MACROS_SCOPE;

        //assumes the carry is zero

        reg_vector data=regs.bind_vector(m, "data");
        reg_vector mask_reg=regs.bind_vector(m, "mask");

        assert_valid();

        APPEND_M(str( "VMOVDQU `data, #", (*this)[size-4] ));
        
        APPEND_M(str( "VPCMPGTQ `mask, `data, #", asm_immediate({data_sign_mask-1}) ));
        
        APPEND_M(str( "VPAND `mask, `mask, #", asm_immediate({0, 0, 0, carry_mask}) ));
        APPEND_M(str( "VPOR `data, `data, `mask" ));

        APPEND_M(str( "VMOVDQU #, `data", (*this)[size-4] ));
    }
};
    
template<class type> void bind_spillable(type regs, string name) {
    for (int x=0;x<regs.size();++x) {
        assert((regs[x].first.value==-1) != (regs[x].second.value==-1));

        if (regs[x].first.value==-1) {
            m.bind(regs[x].second, str( "#_#", name, x ));
        } else {
            m.bind(regs[x].first, str( "#_#", name, x ));
        }
    }
}

//regs (size==3): 5x vector + 9x vector/spill arguments
//nothing spilled: 3x unaligned loads, 9x mul+add, 3x store ;                      3.0 load cycles / 6 alu cycles / 3 store cycles
//everything spilled: 3x unaligned loads, 9x aligned loads, 9x mul+add, 3x store ; 7.5 load cycles / 6 alu cycles / 3 store cycles
template<int size> void matrix_vector_multiply(
    reg_alloc regs,
    array<pair<reg_vector, reg_spill>, size*size> c_matrix,
    array<simd_integer_asm, size> out_vector, array<simd_integer_asm, size> mul_vector, array<simd_integer_asm, size> add_vector,
    int integer_start, int integer_end, int shift, int truncate=-1
) {
    EXPAND_MACROS_SCOPE;

    bind_spillable(c_matrix, "c_matrix");

    reg_vector accumulator=regs.bind_vector(m, "accumulator");
    reg_vector temp=regs.bind_vector(m, "temp");

    array<reg_vector, size> inputs;
    for (int x=0;x<size;++x) {
        inputs[x]=regs.get_vector();
    }
    m.bind(inputs, "inputs");

    assert(size>=2 && size<=3);

    int integer_size=mul_vector[0].size;

    assert(integer_start>=0 && integer_end>0 && integer_start%4==0 && integer_end%4==0 && integer_end>integer_start);

    assert(shift>=0 && shift<=3);

    for (int x=0;x<size;++x) {
        assert(add_vector[x].size==-1 || add_vector[x].size==integer_size);
        assert(mul_vector[x].size==integer_size);
        assert(out_vector[x].size==integer_size);
    }

    for (int x=integer_end-4;x>=integer_start;x-=4) {
        bool accumulator_is_zero=false;

        for (int y=0;y<size;++y) {
            int pos=x-shift;
            assert(pos>=-3 && pos<=integer_size-4); //can be less than -3 if shift is >=4

            APPEND_M(str( "VMOVDQU `inputs_#, #", y, mul_vector[y][pos] ));

            array<int, 4> mask;
            for (int i=0;i<4;++i) {
                mask[i]=(pos+i>=0 && (truncate==-1 || pos+i<truncate))? 1 : 0;
            }

            if (mask!=array<int, 4>({1, 1, 1, 1})) {
                if (!accumulator_is_zero) {
                    asm_immediate.assign(accumulator, 0);
                    accumulator_is_zero=true;
                }

                APPEND_M(str( "VPBLENDD `inputs_#, `accumulator, `inputs_#, #", y, y, vpblendd_mask_4(mask) ));
            }
        }

        for (int y=0;y<size;++y) {
            for (int z=0;z<size;++z) {
                if (z==0 && add_vector[y].size==-1) {
                    APPEND_M(str( "VPMULDQ `accumulator, `inputs_#, `c_matrix_#", z, y*size+z ));
                } else {
                    APPEND_M(str( "VPMULDQ `temp, `inputs_#, `c_matrix_#", z, y*size+z ));
                    if (z==0) {
                        APPEND_M(str( "VPADDQ `accumulator, `temp, #", add_vector[y][x] ));
                    }  else {
                        APPEND_M(str( "VPADDQ `accumulator, `temp, `accumulator" ));
                    }
                }
            }

            APPEND_M(str( "VMOVDQU #, `accumulator", out_vector[y][x] ));
        }
    }
}

//regs: 4x vector
void integer_multiply(reg_alloc regs, simd_integer_asm c, simd_integer_asm a, simd_integer_asm b) {
    EXPAND_MACROS_SCOPE;

    a.sign_extend_data(regs);
    b.sign_extend_data(regs);

    c.clear(regs);
    int current_num_adds=0;

    //slow
    for (int x=0;x<b.size;++x) {
        if (current_num_adds+1>max_num_adds) {
            c.calculate_carry(regs, 0, true);
            current_num_adds=0;
        }

        {
            EXPAND_MACROS_SCOPE;
            reg_alloc regs_copy=regs;

            reg_vector b_reg=regs_copy.bind_vector(m, "b");

            APPEND_M(str( "VPBROADCASTQ `b, #", b[x] ));

            c.fma(regs_copy, c, a, b_reg, x);
            ++current_num_adds;
        }
    }

    c.calculate_carry(regs);

    a.zero_extend_data(regs);
    b.zero_extend_data(regs);
}

//regs: 3x vector, 1x vector arguments
//result is a scalar and is max of index+1
template<unsigned long int num> void extract_bits_shifted_scan(
    reg_alloc regs, reg_vector max_index, array<simd_integer_asm, num> ints
) {
    EXPAND_MACROS_SCOPE;

    m.bind(max_index, "max_index");

    reg_vector index_reg=regs.bind_vector(m, "index");
    reg_vector current_value=regs.bind_vector(m, "current_value");
    reg_vector current_mask=regs.bind_vector(m, "current_mask");

    for (int x=1;x<num;++x) {
        assert(ints[x].size==ints[0].size);
        assert(ints[x].memory_base_reg.value==ints[0].memory_base_reg.value);
    }

    asm_immediate.assign(max_index, 0);
    asm_immediate.assign(index_reg, {1, 2, 3, 4});

    //incoming index is in `index
    auto accumulate=[&]() {
        APPEND_M(str( "VPCMPEQQ `current_mask, `current_value, #", asm_immediate(0) )); //true if current_value==0
        APPEND_M(str( "VPANDN `current_mask, `current_mask, `index" )); //0 if current_value==0, else index+1

        APPEND_M(str( "VPMAXUD `max_index, `max_index, `current_mask" )); //max_index=max(max_index, index+1)
    };

    for (int index=0;index<ints[0].size;index+=4) {
        for (int x=0;x<num;++x) {
            string command=(x==0)? "VMOVDQU `current_value" : "VPOR `current_value, `current_value";
            APPEND_M(str( "#, #", command, ints[x][index] ));
        }

        accumulate();

        APPEND_M(str( "VPADDQ `index, `index, #", asm_immediate(4) ));
    }

    //max reduce

    //max of lanes 0/1 in lanes 0 and 1 ; max of lanes 2/3 in lanes 2 and 3
    APPEND_M(str( "VPERMQ `index, `max_index, #", vpermq_mask({1, 0, 3, 2}) ));
    APPEND_M(str( "VPMAXUD `max_index, `max_index, `index" ));

    //max of 0/1 and 2/3 in lanes 0 and 1 ; max of 2/3 and 0/1 in lanes 2 and 3 (i.e. max of all lanes in each lane)
    APPEND_M(str( "VPERMQ `index, `max_index, #", vpermq_mask({2, 3, 0, 1}) ));
    APPEND_M(str( "VPMAXUD `max_index, `max_index, `index" ));
}

template<int num> void extract_bits_shifted(
    reg_alloc regs, array<reg_scalar, num> res, reg_scalar end_index, array<simd_integer_asm, num> ints
) {
    EXPAND_MACROS_SCOPE;

    m.bind(res, "res");
    m.bind(end_index, "end_index");
    m.bind(ints[0].memory_base_reg, "c_memory_base");

    reg_vector max_index_vector=regs.bind_vector(m, "max_index_vector");
    extract_bits_shifted_scan(regs, max_index_vector, ints);

    reg_scalar c_end_index=regs.bind_scalar(m, "c_end_index");
    reg_scalar mask=regs.bind_scalar(m, "mask");
    reg_scalar tmp=regs.bind_scalar(m, "tmp");

    assert(ints[0].size>=3);

    //todo //can probably use simd for this with scatter reads

    APPEND_M(str( "VMOVQ `end_index, `max_index_vector_128" ));
    APPEND_M(str( "DEC `end_index" ));

    APPEND_M(str( "MOV `c_end_index, `end_index" ));

    //need c_end_index to be at least 2 to avoid reading out of bounds
    string label_max_skip=m.alloc_label();
    APPEND_M(str( "CMP `c_end_index, 2" )); //c_end_index=max(c_end_index, 2)
    APPEND_M(str( "JGE #", label_max_skip )); //almost always taken
    APPEND_M(str( "MOV `c_end_index, 2" ));
    APPEND_M(str( "#:", label_max_skip ));

    //make c_end_index into a memory_base register so that it can be used as an offset
    APPEND_M(str( "SHL `c_end_index, 0x3" ));
    APPEND_M(str( "ADD `c_end_index, `c_memory_base" ));

    //load the msb values from each int (at c_end_index). also calculate the mask
    //if end_index<c_end_index, the mask will be 0
    for (int x=0;x<num;++x) {
        APPEND_M(str( "MOV `res_#, #", x, asm_memory(c_end_index, (ints[x].start)*8) ));
        APPEND_M(str( "# `mask, `res_#", (x==0)? "MOV" : "OR", x ));
    }

    //calculate num_leading_bits
    APPEND_M(str( "LZCNT `mask, `mask" ));
    APPEND_M(str( "SUB `mask, #", to_hex(carry_size) )); //max value is 64-carry_size which is <= data_size

    //calculate data_size+num_leading_bits
    APPEND_M(str( "MOV `tmp, `mask" ));
    APPEND_M(str( "ADD `tmp, #", data_size ));

    for (int x=0;x<num;++x) {
        //highest limb of result (nonnegative)
        APPEND_M(str( "SHLX `res_#, `res_#, `tmp", x, x )); //data_size+num_leading_bits
    }

    for (int x=0;x<num;++x) {
        //second highest limb
        APPEND_M(str( "SHLX `tmp, #, `mask", asm_memory(c_end_index, (ints[x].start-1)*8) ));
        APPEND_M(str( "OR `res_#, `tmp", x ));
    }

    //calculate data_size-num_leading_bits
    APPEND_M(str( "NEG `mask" ));
    APPEND_M(str( "ADD `mask, #", to_hex(data_size) ));

    for (int x=0;x<num;++x) {
        //lsb limb
        APPEND_M(str( "SHRX `tmp, #, `mask", asm_memory(c_end_index, (ints[x].start-2)*8) ));
        APPEND_M(str( "OR `res_#, `tmp", x ));
    }
}

//this doesn't actually resize a
/*void shrink(reg_alloc regs, simd_integer_asm a, int expected_size) {
    EXPAND_MACROS_SCOPE;

    assert(a.memory_base_reg.value==asm_memory.memory_base.value);

    assert(expected_size>=1 && a.size>=expected_size);

    reg_vector max_index_vector=regs.bind_vector(m, "max_index_vector");
    reg_scalar tmp=regs.bind_scalar(m, "tmp");

    //calculates max index + 1 which is also the size
    extract_bits_shifted_scan<1>(regs, max_index_vector, {a});

    APPEND_M(str( "VMOVQ `tmp, `max_index_vector_128" )); //tmp is nonzero size

    APPEND_M(str( "CMP `tmp, #", to_hex(expected_size) ));
    APPEND_M(str( "JA #", m.alloc_error_label() ));

    APPEND_M(str( "MOV `tmp, #", a[expected_size] ));
    APPEND_M(str( "TEST `tmp, #", to_hex(data_sign_mask) ));
    APPEND_M(str( "JNZ #", m.alloc_error_label() )); //jump taken if last data bit is set in a[expected_size]
} */

struct batched_bit_shifts_entry {
    map<int, int /*shift amount*/> sources;
};

struct batched_bit_shifts_entries {
    vector<batched_bit_shifts_entry> entries;

    batched_bit_shifts_entries(int t_output_size) {
        assert(t_output_size%4==0);
        entries.resize(t_output_size);
    }

    void add(int output_index, int input_index, int amount) {
        assert(amount>-64 && amount<64);
        auto& v=entries.at(output_index).sources;
        assert(!v.count(input_index));
        v[input_index]=amount;
    }
};

//if the input or output is a gmp integer, the start is set to 0 and the register is set to the address register for the gmp integer
//the entire output is buffered in registers, so should call this multiple times with different offsets to apply this to a large output
void batched_bit_shifts_impl(
    reg_alloc regs, simd_integer_asm input, simd_integer_asm output, vector<batched_bit_shifts_entry> entries, bool zero_carry
) {
    EXPAND_MACROS_SCOPE;

    vector<reg_vector> output_buffers;
    for (int x=0;x<output.size/4;++x) {
        output_buffers.push_back(regs.get_vector());
    }
    m.bind(output_buffers, "output_buffers");

    reg_vector input_buffer=regs.bind_vector(m, "input_buffer");
    reg_vector tmp=regs.bind_vector(m, "tmp");

    for (int x=0;x<output.size/4;++x) {
        asm_immediate.assign(output_buffers[x], 0);
    }

    assert(entries.size()==output.size);
    assert(entries.size()%4==0);
    for (int x=0;x<input.size;++x) {
        bool did_broadcast=false;

        for (int is_right=0;is_right<2;++is_right) {
            for (int y=0;y<entries.size();y+=4) {
                //array<int, 4> blend_mask={0, 0, 0, 0};
                array<uint64, 4> shift_amounts={~uint64(0), ~uint64(0), ~uint64(0), ~uint64(0)};

                bool has_shift=false;
                for (int z=0;z<4;++z) {
                    auto i=entries.at(y+z).sources.find(x);
                    if (i==entries.at(y+z).sources.end()) {
                        continue;
                    }

                    bool i_is_right=(i->second)<0;
                    if (i_is_right!=is_right) {
                        continue;
                    }

                    int amount=(i_is_right)? -i->second : i->second;
                    assert(amount>=0 && amount<64);

                    //blend_mask[z]=1;
                    shift_amounts[z]=amount;
                    has_shift=true;
                }

                if (!has_shift) {
                    continue;
                }

                if (!did_broadcast) {
                    APPEND_M(str( "VPBROADCASTQ `input_buffer, #", input[x] ));
                    did_broadcast=true;
                }

                if (is_right) {
                    APPEND_M(str( "VPSRLVQ `tmp, `input_buffer, #", asm_immediate(shift_amounts) ));
                } else {
                    APPEND_M(str( "VPSLLVQ `tmp, `input_buffer, #", asm_immediate(shift_amounts) ));
                }

                APPEND_M(str( "VPOR `output_buffers_#, `tmp, `output_buffers_#", y/4, y/4 ));
                //APPEND_M(str( "VPBLENDD `output_buffers_#, `output_buffers_#, `tmp, #", y/4, y/4, vpblendd_mask_4(blend_mask) ));
            }
        }
    }

    if (zero_carry) {
        asm_immediate.assign(tmp, data_mask);
    }

    for (int x=0;x<output.size/4;++x) {
        if (zero_carry) {
            APPEND_M(str( "VPAND `output_buffers_#, `output_buffers_#, `tmp", x, x ));
        }

        APPEND_M(str( "VMOVDQU #, `output_buffers_#", output[x*4], x ));
    }
}

void batched_bit_shifts(
    reg_alloc regs, simd_integer_asm input, simd_integer_asm output, batched_bit_shifts_entries entries, bool zero_carry
) {
    const int batch_size=14*4;

    vector<batched_bit_shifts_entry> batch_entries;

    assert(output.size==entries.entries.size());
    for (int x=0;x<output.size;x+=batch_size) {
        int end_x=x+batch_size;
        if (end_x>output.size) {
            end_x=output.size;
        }

        batch_entries.clear();
        for (int y=x;y<end_x;++y) {
            batch_entries.push_back(entries.entries.at(y));
        }

        simd_integer_asm c_output=output;
        c_output.start+=x;
        c_output.size=end_x-x;
        batched_bit_shifts_impl(regs, input, c_output, batch_entries, zero_carry);
    }
}

void simd_to_mpz(
    reg_alloc regs, simd_integer_asm input, simd_integer_asm output
) {
    input.assert_valid();
    output.assert_valid();

    batched_bit_shifts_entries entries(output.size);

    for (int input_index=0;input_index<input.size;++input_index) {
        for (int output_index=0;output_index<output.size;++output_index) {
            int input_start=data_size*input_index;
            int output_start=64*output_index;
            int amount=input_start-output_start;

            if (amount>-data_size && amount<64) {
                entries.add(output_index, input_index, amount);
            }
        }
    }

    batched_bit_shifts(regs, input, output, entries, false);
}

void mpz_to_simd(
    reg_alloc regs, simd_integer_asm input, simd_integer_asm output
) {
    input.assert_valid();
    output.assert_valid();

    batched_bit_shifts_entries entries(output.size);

    for (int input_index=0;input_index<input.size;++input_index) {
        for (int output_index=0;output_index<output.size;++output_index) {
            int input_start=64*input_index;
            int output_start=data_size*output_index;
            int amount=input_start-output_start;

            if (amount>-64 && amount<data_size) {
                entries.add(output_index, input_index, amount);
            }
        }
    }

    batched_bit_shifts(regs, input, output, entries, true);
}


}