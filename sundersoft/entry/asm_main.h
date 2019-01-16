namespace simd_integer_namespace {


generic_stats track_cycles_reduce;
generic_stats track_cycles_gcd;
array<generic_stats, 100> track_cycles_test;

array<simd_integer_asm, num_asm_ints> global_asm_ints;

//need to assign values when compiling asm code. don't need to during runtime
reg_alloc asm_init(bool assign_values) {
    assert(asm_int_size_limbs%4==0);
    assert(asm_int_size_limbs*8<=4096-32*2); //has padding before and after the integer

    reg_alloc regs;
    regs.init(true);

    asm_memory.memory_base=regs.bind_scalar(m, "memory_base", reg_rdi);
    reg_spill::memory_base=asm_memory.memory_base;

    generate_divide_table(assign_values);

    //todo //can tune the seed for this
    mt19937 c_rand;

    for (int x=0;x<num_asm_ints;++x) {
        //want the integers to be within a single page since they use unaligned reads and those are slow if they cross a page boundary
        //will store one integer per page
        int pos=asm_memory.alloc(4096, 4096);

        const int min_offset=32;
        const int max_offset=4096-32-asm_int_size_limbs*8;
        assert(max_offset>=min_offset);
        assert(min_offset%32==0);
        assert(max_offset%32==0);

        //to avoid false dependencies, want to have the integer position be random within a page
        int offset=uniform_int_distribution<int>(min_offset/32, max_offset/32)(c_rand)*32;

        assert((pos+offset)%32==0);
        assert(offset>=min_offset && offset<=max_offset);

        simd_integer_asm c;
        c.memory_base_reg=asm_memory.memory_base;
        c.start=(pos+offset)/8;
        c.size=asm_int_size_limbs;

        global_asm_ints[x]=c;
    }

    return regs;
}

uint64* get_address(simd_integer_asm targ) {
    return &(((uint64*)asm_memory.memory_buffer_base)[targ.start]);
}

void transfer(simd_integer_asm targ, simd_integer& source) {
    assert(targ.size>=source.current_size());
    assert(targ.memory_base_reg.value==asm_memory.memory_base.value);

    for (int x=0;x<source.current_size();++x) {
        ((uint64*)asm_memory.memory_buffer_base)[targ.start+x]=source.memory.at(x);
    }
}

void transfer(simd_integer& targ, simd_integer_asm source) {
    assert(source.size>=targ.current_size());
    assert(source.memory_base_reg.value==asm_memory.memory_base.value);

    for (int x=0;x<targ.current_size();++x) {
        targ.memory.at(x)=((uint64*)asm_memory.memory_buffer_base)[source.start+x];
    }
}

void transfer(reg_spill targ, array<uint64, 4> source) {
    assert(targ.value>=0 && targ.value<num_spilled_registers);

    for (int x=0;x<4;++x) {
        ((uint64*)(asm_memory.memory_buffer_base-32*(targ.value+1)))[x]=source[x];
    }
}

void transfer(array<uint64, 4>& targ, reg_spill source) {
    assert(source.value>=0 && source.value<num_spilled_registers);

    for (int x=0;x<4;++x) {
        targ[x]=((uint64*)(asm_memory.memory_buffer_base-32*(source.value+1)))[x];
    }
}

#ifdef COMPILE_ASM
    #define ASM_FUNC(NAME, ...)\
    void compile_asm_ ## NAME(reg_alloc regs) {\
        EXPAND_MACROS_SCOPE;\
        asm_function c_func( #NAME );\
        __VA_ARGS__\
    }
#else
    #define ASM_FUNC(NAME, ...)\
    extern "C" { int asm_func_ ## NAME(void* memory_base); }
#endif

ASM_FUNC(init_memory,
    asm_memory.assign_initial_data(regs);
)
#ifndef COMPILE_ASM
    void init_memory_asm() {
        asm_func_init_memory(asm_memory.memory_buffer_base);
    }
#endif

//
//

void transfer_impl(
    reg_alloc& regs, int input_size, int output_size, simd_integer_asm& input, simd_integer_asm& output, bool is_padded
) {
    simd_integer_asm args=global_asm_ints[0];

    reg_scalar input_address=regs.bind_scalar(m, "input_address");
    APPEND_M(str( "MOV `input_address, #", args[0] ));

    reg_scalar output_address=regs.bind_scalar(m, "output_address");
    APPEND_M(str( "MOV `output_address, #", args[1] ));

    input.start=0;
    input.size=input_size;
    input.memory_base_reg=input_address;

    output.start=0;
    output.size=output_size;
    output.memory_base_reg=output_address;

    if (is_padded) {
        output.clear(regs);
        output.start+=output.size;
    }
}

ASM_FUNC(transfer_padded_9,
    simd_integer_asm input;
    simd_integer_asm output;
    transfer_impl(regs, 20, 9*4, input, output, true);
    mpz_to_simd(regs, input, output);
)
ASM_FUNC(transfer_padded_18,
    simd_integer_asm input;
    simd_integer_asm output;
    transfer_impl(regs, 36, 18*4, input, output, true);
    mpz_to_simd(regs, input, output);
)
ASM_FUNC(transfer_9,
    simd_integer_asm input;
    simd_integer_asm output;
    transfer_impl(regs, 9*4, 20, input, output, false);
    simd_to_mpz(regs, input, output);
)
ASM_FUNC(transfer_18,
    simd_integer_asm input;
    simd_integer_asm output;
    transfer_impl(regs, 18*4, 36, input, output, false);
    simd_to_mpz(regs, input, output);
)
#ifndef COMPILE_ASM
    void transfer_args(void* input, void* output) {
        *(void**)get_address(global_asm_ints[0])=input;
        *(void**)(get_address(global_asm_ints[0])+1)=output;
    }

    //num_limbs: 9*4, 18*4
    simd_integer_asm transfer_padded_slow(simd_integer_asm targ, mpz_struct* a, int num_limbs) {
        targ.size=2*num_limbs;
        memset(get_address(targ), 0, 8*2*num_limbs);

        targ.start+=num_limbs;
        targ.size=num_limbs;

        mpz_export(get_address(targ), nullptr, -1, 8, -1, carry_size, a);

        return targ;
    }

    //num_limbs: 9*4, 18*4
    void transfer_slow(mpz_struct* a, simd_integer_asm targ, int num_limbs) {
        mpz_import(a, targ.size, -1, 8, -1, carry_size, get_address(targ));
    }

    simd_integer_asm transfer_padded_fast(simd_integer_asm targ, mpz_struct* a, int num_limbs) {
        assert(num_limbs==9*4 || num_limbs==18*4);
        int a_num_limbs=(num_limbs==9*4)? 20 : 36;

        mp_limb_t* limbs=mpz_limbs_modify(a, a_num_limbs); //this calls mpz_realloc

        //mpz doesn't zero out the padding
        for (int x=mpz_size(a);x<a_num_limbs;++x) {
            limbs[x]=0;
        }

        transfer_args(limbs, get_address(targ));

        ((num_limbs==9*4)? asm_func_transfer_padded_9 : asm_func_transfer_padded_18)(asm_memory.memory_buffer_base);

        targ.start+=num_limbs;
        targ.size=num_limbs;

        return targ;
    }

    void transfer_fast(mpz_struct* a, simd_integer_asm targ, int num_limbs) {
        assert(num_limbs==9*4 || num_limbs==18*4);
        int a_num_limbs=(num_limbs==9*4)? 20 : 36;

        mp_limb_t* limbs=mpz_limbs_write(a, a_num_limbs);
        transfer_args(get_address(targ), limbs);

        ((num_limbs==9*4)? asm_func_transfer_9 : asm_func_transfer_18)(asm_memory.memory_buffer_base);

        mpz_limbs_finish(a, a_num_limbs);
    }

    simd_integer_asm transfer_padded(simd_integer_asm targ, mpz_struct* a, int num_limbs) {
        simd_integer_asm res=transfer_padded_fast(targ, a, num_limbs);

        /*{
            todo
            simd_integer_asm res2=transfer_padded_slow(global_asm_ints[num_asm_ints-1], a, num_limbs);

            assert(res.size==res2.size);

            simd_integer_asm res1=res;
            res1.start-=res1.size;
            res1.size*=2;
            res2.start-=res2.size;
            res2.size*=2;

            for (int x=0;x<res1.size;++x) {
                assert(get_address(res1)[x] == get_address(res2)[x]);
            }
        }*/

        return res;
    }

    void transfer(mpz_struct* a, simd_integer_asm targ, int num_limbs) {
        transfer_fast(a, targ, num_limbs);

        /*{
            todo
            integer a_copy;
            transfer_slow(a_copy.impl, targ, num_limbs);

            integer a_copy_2;
            mpz_set(a_copy_2.impl, a);

            assert(a_copy==a_copy_2);
        }*/
    }
#endif

//
//

ASM_FUNC(integer_gcd,
    simd_integer_asm args=global_asm_ints[0];
    reg_spill v0_sign=regs.bind_spill(m, "v0_sign");

    simd_integer_asm v0=global_asm_ints[10];
    simd_integer_asm a_padded=global_asm_ints[1];
    simd_integer_asm b_padded=global_asm_ints[2];

    v0.size=9*4;
    a_padded.size=2*v0.size;
    b_padded.size=2*v0.size;

    array<simd_integer_asm, 3> large_buffers={global_asm_ints[3], global_asm_ints[4], global_asm_ints[5]};
    array<simd_integer_asm, 4> small_buffers={global_asm_ints[6], global_asm_ints[7], global_asm_ints[8], global_asm_ints[9]};

    for (simd_integer_asm& c : large_buffers) {
        c.size=v0.size;
    }

    for (simd_integer_asm& c : small_buffers) {
        c.size=12;
    }

    gcd(regs, v0_sign, a_padded, b_padded, v0, large_buffers, small_buffers, false, true);

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(v0_sign, "v0_sign");

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

        APPEND_M(str( "MOV `tmp, `v0_sign" ));
        APPEND_M(str( "MOV #, `tmp", args[0] ));
    }
)
#ifndef COMPILE_ASM
    //b>=0 ; b>=|a| (a and b are backwards) ; assumes gcd is 1 and won't assign g
    void integer_gcd_asm(mpz_struct* g, mpz_struct* s, mpz_struct* t, mpz_struct* a, mpz_struct* b) {
        assert(t==nullptr);

        const int max_size=1030;
        const int num_limbs=9*4;

        bool success=false;

        //this works if the two inputs are 1024 bits or less with a GMP limb size of 64 bits
        int a_bits=mpz_num_bits_upper_bound(a); //mpz_sizeinbase(a, 2)
        int b_bits=mpz_num_bits_upper_bound(b); //mpz_sizeinbase(b, 2)

        if (a_bits<=max_size && b_bits<=max_size) {
            //arguments are backwards because the first input is smaller than the second
            simd_integer_asm a_int=transfer_padded(global_asm_ints[2], a, num_limbs);
            simd_integer_asm b_int=transfer_padded(global_asm_ints[1], b, num_limbs);

            int res;
            {
                track_cycles c_track_cycles(track_cycles_gcd); // 25878
                res=asm_func_integer_gcd(asm_memory.memory_buffer_base);

                if (res!=0) {
                    c_track_cycles.abort();
                }
            }

            if (res==0) {
                simd_integer_asm v0_int=global_asm_ints[10];
                v0_int.start+=num_limbs;

                //need to do this before transfer because it also uses global_asm_ints[0]
                int64 v0_sign=*(int64*)get_address(global_asm_ints[0]);

                transfer(s, global_asm_ints[10], num_limbs);

                //b<0:
                // (-b)s===1 mod a
                // b(-s)===1 mod a
                int a_sign=(mpz_sgn(a)==-1)? -1 : 1;

                if (v0_sign*a_sign==-1) {
                    mpz_neg(s, s);
                }

                if (asm_inject_random_errors_rate!=0 && (rand()%asm_inject_random_errors_rate)==0) {
                    int bit_pos=rand()%1024;
                    mpz_combit(s, bit_pos);
                }

                //todo cout << "G\n";
                success=true;
            }
        }

        if (!success) {
            //todo cout << "g\n";
            mpz_gcdext(g, s, t, a, b);
        }
    }
#endif

//
//

ASM_FUNC(integer_reduce,
    simd_integer_asm args=global_asm_ints[0];
    simd_integer_asm a_padded=global_asm_ints[1];
    simd_integer_asm b_padded=global_asm_ints[2];
    simd_integer_asm c_padded=global_asm_ints[3];

    int size=18*4;
    a_padded.size=2*size;
    b_padded.size=2*size;
    c_padded.size=2*size;

    array<simd_integer_asm, 3> large_buffers={global_asm_ints[4], global_asm_ints[5], global_asm_ints[6]};
    array<simd_integer_asm, 6> small_buffers={
        global_asm_ints[7], global_asm_ints[8], global_asm_ints[9],
        global_asm_ints[10], global_asm_ints[11], global_asm_ints[12]
    };

    for (simd_integer_asm& c : large_buffers) {
        c.size=size;
    }

    for (simd_integer_asm& c : small_buffers) {
        c.size=12;
    }

    reg_spill b_sign=regs.bind_spill(m, "b_sign");

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(b_sign, "b_sign");

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

        APPEND_M(str( "MOV `tmp, #", args[0] ));
        APPEND_M(str( "MOV `b_sign, `tmp" ));
    }

    reduce(regs, a_padded, b_padded, c_padded, b_sign, large_buffers, small_buffers);

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(b_sign, "b_sign");

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

        APPEND_M(str( "MOV `tmp, `b_sign" ));
        APPEND_M(str( "MOV #, `tmp", args[0] ));
    }
)

#ifndef COMPILE_ASM
    //if this returns false then it shouldn't be called again with these inputs
    //result might not be completely finished. if so, do one iteration then call this again
    bool integer_reduce_asm(mpz_struct* a, mpz_struct* b, mpz_struct* c) {
        const int max_size=2060;
        const int num_limbs=18*4;

        bool success=false;

        int a_bits;
        int b_bits;
        int c_bits;

        {
            track_cycles c_track_cycles(track_cycles_test[0]); // 28
            //this works if the number of bits is 2048 or less
            a_bits=mpz_num_bits_upper_bound(a); //mpz_sizeinbase(a, 2)
            b_bits=mpz_num_bits_upper_bound(b); //mpz_sizeinbase(b, 2)
            c_bits=mpz_num_bits_upper_bound(c); //mpz_sizeinbase(c, 2)
        }

        if (
            a_bits<=max_size &&
            b_bits<=max_size &&
            c_bits<=max_size
        ) {
            simd_integer_asm a_int;
            simd_integer_asm b_int;
            simd_integer_asm c_int;

            //arguments are backwards because the first input is smaller than the second
            {
                track_cycles c_track_cycles(track_cycles_test[1]); // 474
                a_int=transfer_padded(global_asm_ints[1], a, num_limbs);
                b_int=transfer_padded(global_asm_ints[2], b, num_limbs);
                c_int=transfer_padded(global_asm_ints[3], c, num_limbs);
            }

            //both this and transfer_padded use global_asm_ints[0], so need to run this after
            int64& b_sign=*(int64*)get_address(global_asm_ints[0]);

            {
                track_cycles c_track_cycles(track_cycles_test[2]); // 24
                b_sign=(mpz_sgn(b)==-1)? -1 : 1;
            }

            int res;
            {
                //original:                             34762
                //no background work first iteration:   35012 (got rid of this)
                //skip multiplying identity matricies:  34078
                //6*(|q|+1) instead of 2*(|q|+1):       34154
                //cutoffs:                              33294
                track_cycles c_track_cycles(track_cycles_reduce);
                res=asm_func_integer_reduce(asm_memory.memory_buffer_base);
                if (res!=0) {
                    c_track_cycles.abort();
                }
            }
            if (res==0) {
                //both this and transfer use global_asm_ints[0], so need to run this before
                int64 b_sign_copy=b_sign;

                {
                    track_cycles c_track_cycles(track_cycles_test[3]); // 610
                    transfer(a, a_int, num_limbs);
                    transfer(b, b_int, num_limbs);
                    transfer(c, c_int, num_limbs);
                }

                if (b_sign_copy==-1) {
                    track_cycles c_track_cycles(track_cycles_test[4]); // 24
                    mpz_neg(b, b);
                }

                //errors in the sign of b won't be detected
                if (asm_inject_random_errors_rate!=0 && (rand()%asm_inject_random_errors_rate)==0) {
                    int which=rand()%3;
                    mpz_struct* target;

                    if (which==0) {
                        target=a;
                    } else
                    if (which==1) {
                        target=b;
                    } else {
                        target=c;
                    }

                    int bit_pos=rand()%1024;
                    mpz_combit(target, bit_pos);
                }

                success=true;
            }
        }

        return success;
    }
#endif

//
//

#ifndef COMPILE_ASM
    bool integer_gcd_asm(simd_integer& a, simd_integer& b, simd_integer& v0, int64& v0_sign) {
        simd_integer args;
        args.memory.resize(4);

        simd_integer_asm args_int=global_asm_ints[0];
        args_int.size=4;

        v0.memory.resize(9*4);
        assert(a.current_size()==9*4);
        assert(b.current_size()==9*4);

        simd_integer_asm a_int=global_asm_ints[1];
        a_int.start+=a.current_size();

        simd_integer_asm b_int=global_asm_ints[2];
        b_int.start+=b.current_size();

        transfer(args_int, args);
        transfer(a_int, a);
        transfer(b_int, b);

        int res=asm_func_integer_gcd(asm_memory.memory_buffer_base);

        if (res!=0) {
            //todo print(str( "label_error_#", res ));  assert(false);
            return false;
        }

        transfer(args, args_int);
        transfer(a, a_int); //this works no matter how big the result is
        transfer(v0, global_asm_ints[10]);

        v0_sign=args.memory.at(0);

        return true;
    }
#endif

ASM_FUNC(integer_multiply,
    simd_integer_asm c=global_asm_ints[0];
    simd_integer_asm a=global_asm_ints[1];
    simd_integer_asm b=global_asm_ints[2];

    c.size=2*9*4;
    a.size=9*4;
    b.size=9*4;

    integer_multiply(regs, c, a, b);
)
#ifndef COMPILE_ASM
    void integer_multiply_asm(simd_integer& c, simd_integer& a, simd_integer& b) {
        c.memory.resize(2*9*4);
        assert(a.current_size()==9*4);
        assert(b.current_size()==9*4);

        transfer(global_asm_ints[1], a);
        transfer(global_asm_ints[2], b);

        asm_func_integer_multiply(asm_memory.memory_buffer_base);

        transfer(c, global_asm_ints[0]);
    }
#endif

ASM_FUNC(integer_fma,
    simd_integer_asm c=global_asm_ints[0];
    simd_integer_asm a=global_asm_ints[1];
    simd_integer_asm b=global_asm_ints[2];

    c.size=9*4;
    a.size=9*4;
    b.size=1;

    reg_vector v=regs.bind_vector(m, "v");
    APPEND_M(str( "VPBROADCASTQ `v, #", b[0] ));

    c.fma(regs, simd_integer_asm(), a, v, 0);

    c.calculate_carry(regs);
)
#ifndef COMPILE_ASM
    void integer_fma_asm(simd_integer& c, simd_integer& a, simd_integer& b) {
        c.memory.resize(9*4);
        assert(a.current_size()==9*4);
        //assert(b.current_size()==1);
        
        transfer(global_asm_ints[1], a);
        transfer(global_asm_ints[2], b);
        
        asm_func_integer_fma(asm_memory.memory_buffer_base);
        
        transfer(c, global_asm_ints[0]);
    }
#endif

ASM_FUNC(divide_table,
    simd_integer_asm args=global_asm_ints[0];

    reg_scalar a=regs.bind_scalar(m, "a");
    reg_scalar b=regs.bind_scalar(m, "b");
    reg_scalar q=regs.bind_scalar(m, "q");
    reg_scalar r=regs.bind_scalar(m, "r", reg_rdx);

    APPEND_M(str( "MOV `a, #", args[0] ));
    APPEND_M(str( "MOV `b, #", args[1] ));

    divide_table(regs, a, b, q, r);

    APPEND_M(str( "MOV #, `q", args[0] ));
    APPEND_M(str( "MOV #, `r", args[1] ));
)
#ifndef COMPILE_ASM
    int64 divide_table_asm(int64 a, int64 b, int64& q) {
        simd_integer args;
        args.memory.resize(2);
        args.memory.at(0)=uint64(a);
        args.memory.at(1)=uint64(b);

        simd_integer_asm args_int=global_asm_ints[0];
        transfer(args_int, args);

        asm_func_divide_table(asm_memory.memory_buffer_base);

        transfer(args, args_int);

        q=args.memory.at(0);
        return args.memory.at(1);
    }
#endif

ASM_FUNC(gcd_64_iteration,
    simd_integer_asm args=global_asm_ints[0];

    reg_spill c_gcd_mask=regs.bind_spill(m, "c_gcd_mask");
    reg_scalar a_0=regs.bind_scalar(m, "a_0");
    reg_scalar a_1=regs.bind_scalar(m, "a_1");
    reg_vector uv_0=regs.bind_vector(m, "uv_0");
    reg_vector uv_1=regs.bind_vector(m, "uv_1");
    reg_scalar valid=regs.bind_scalar(m, "valid");

    APPEND_M(str( "VMOVDQU `uv_0, #", args[0] ));
    APPEND_M(str( "VMOVDQU `c_gcd_mask, `uv_0" ));

    APPEND_M(str( "MOV `a_0, #", args[4] ));
    APPEND_M(str( "MOV `a_1, #", args[5] ));

    APPEND_M(str( "VMOVDQU `uv_0, #", args[6] ));
    APPEND_M(str( "VMOVDQU `uv_1, #", args[10] ));

    string early_exit_label=m.alloc_label();

    APPEND_M(str( "MOV `valid, 0" ));
    gcd_64_iteration(regs, c_gcd_mask, {a_0, a_1}, {uv_0, uv_1}, early_exit_label);
    APPEND_M(str( "MOV `valid, 1" ));
    APPEND_M(str( "#:", early_exit_label ));

    APPEND_M(str( "MOV #, `valid", args[14] ));
    APPEND_M(str( "MOV #, `a_0", args[4] ));
    APPEND_M(str( "MOV #, `a_1", args[5] ));

    APPEND_M(str( "VMOVDQU #, `uv_0", args[6] ));
    APPEND_M(str( "VMOVDQU #, `uv_1", args[10] ));
)
#ifndef COMPILE_ASM
    bool gcd_64_iteration_asm(vector2& a, matrix2& uv, bool approximate) {
        simd_integer args;
        args.memory.resize(15);

        for (int x=0;x<4;++x) {
            args.memory.at(x)=(approximate)? gcd_mask_approximate[x] : gcd_mask_exact[x];
        }

        args.memory.at(4)=a[0];
        args.memory.at(5)=a[1];

        args.memory.at( 6)=uv[2*0+0]; args.memory.at( 7)=uv[2*0+1]; args.memory.at( 8)=0; args.memory.at( 9)=0;
        args.memory.at(10)=uv[2*1+0]; args.memory.at(11)=uv[2*1+1]; args.memory.at(12)=0; args.memory.at(13)=0;

        simd_integer_asm args_int=global_asm_ints[0];
        transfer(args_int, args);

        asm_func_gcd_64_iteration(asm_memory.memory_buffer_base);

        transfer(args, args_int);

        a[0]=args.memory.at(4);
        a[1]=args.memory.at(5);

        uv[2*0+0]=args.memory.at( 6); uv[2*0+1]=args.memory.at( 7); assert(0==args.memory.at( 8)); assert(0==args.memory.at( 9));
        uv[2*1+0]=args.memory.at(10); uv[2*1+1]=args.memory.at(11); assert(0==args.memory.at(12)); assert(0==args.memory.at(13));

        return (args.memory.at(14)==1);
    }
#endif

ASM_FUNC(multiply_matrix_batch_2,
    simd_integer_asm args=global_asm_ints[0];

    const int size=2;

    reg_vector tmp=regs.bind_vector(m, "tmp");

    array<array<reg_spill, size>, 4> res;
    array<array<reg_spill, size>, 4> current_batch;

    vector<reg_spill> spills;

    for (int x=0;x<4;++x) {
        for (int y=0;y<size;++y) {
            spills.push_back(regs.get_spill());
            current_batch[x][y]=spills.back();
        }
    }

    for (int x=0;x<4;++x) {
        for (int y=0;y<size;++y) {
            spills.push_back(regs.get_spill());
            res[x][y]=spills.back();
        }
    }

    m.bind(spills, "spills");
    for (int x=0;x<spills.size();++x) {
        APPEND_M(str( "VMOVDQU `tmp, #", args[x*4] ));
        APPEND_M(str( "VMOVDQU `spills_#, `tmp", x ));
    }

    multiply_matrix_batch(regs, res, current_batch);

    for (int x=0;x<spills.size();++x) {
        APPEND_M(str( "VMOVDQU `tmp, `spills_#", x ));
        APPEND_M(str( "VMOVDQU #, `tmp", args[x*4] ));
    }
)
#ifndef COMPILE_ASM
    template<unsigned long int size> void multiply_matrix_batch_2_asm(
        array<array<int64, size*size>, 4>& res, array<array<int64, size*size>, 4> current_batch
    ) {
        if (size!=2) {
            return;
        }

        simd_integer args;

        for (int z=0;z<4;++z) {
            for (int y=0;y<size;++y) {
                for (int x=0;x<size;++x) {
                    args.memory.push_back(current_batch[z][y*size+x]);
                }
                for (int x=size;x<4;++x) {
                    args.memory.push_back(0);
                }
            }
        }

        simd_integer_asm args_int=global_asm_ints[0];
        transfer(args_int, args);

        asm_func_multiply_matrix_batch_2(asm_memory.memory_buffer_base);

        args_int.start+=4*size*4;
        transfer(args, args_int);

        int i=0;
        for (int z=0;z<4;++z) {
            for (int y=0;y<size;++y) {
                for (int x=0;x<size;++x) {
                    res[z][y*size+x]=args.memory.at(i);
                    ++i;
                }
                for (int x=size;x<4;++x) {
                    assert(args.memory.at(i)==0);
                    ++i;
                }
            }
        }
    }
#endif

#ifdef COMPILE_ASM
    void compile_asm_multiplier_gcd(reg_alloc regs, int asm_func_multiplier_gcd_mode) {
        EXPAND_MACROS_SCOPE;
        asm_function c_func(str( "multiplier_gcd_#", asm_func_multiplier_gcd_mode ));

        const int num_vectors=2;
        const int size=2;

        matrix_multiplier_asm<2, 2> c_multiplier;

        c_multiplier.int_size=9*4;

        simd_integer_asm args=global_asm_ints[0];

        c_multiplier.ints[0][0]=global_asm_ints[1];
        c_multiplier.ints[0][1]=global_asm_ints[2];

        c_multiplier.ints[0][0].start+=c_multiplier.int_size;
        c_multiplier.ints[0][1].start+=c_multiplier.int_size;

        c_multiplier.ints[1][0]=global_asm_ints[3];
        c_multiplier.ints[1][1]=global_asm_ints[4];
        c_multiplier.ints_buffer[0]=global_asm_ints[5];
        c_multiplier.ints_buffer[1]=global_asm_ints[6];
        c_multiplier.extra_ints[0]=global_asm_ints[7];
        c_multiplier.extra_ints[1]=global_asm_ints[8];
        c_multiplier.extra_ints_multiplied[0]=global_asm_ints[9];
        c_multiplier.extra_ints_multiplied[1]=global_asm_ints[10];

        c_multiplier.ints[0][0].size=c_multiplier.int_size;
        c_multiplier.ints[0][1].size=c_multiplier.int_size;
        c_multiplier.ints[1][0].size=c_multiplier.int_size;
        c_multiplier.ints[1][1].size=c_multiplier.int_size;
        c_multiplier.ints_buffer[0].size=c_multiplier.int_size;
        c_multiplier.ints_buffer[1].size=c_multiplier.int_size;
        c_multiplier.extra_ints[0].size=12;
        c_multiplier.extra_ints[1].size=12;
        c_multiplier.extra_ints_multiplied[0].size=12;
        c_multiplier.extra_ints_multiplied[1].size=12;

        reg_scalar c_memory_base=regs.bind_scalar(m, "c_memory_base");

        c_multiplier.ints[0][0].memory_base_reg=c_memory_base;
        c_multiplier.ints[0][1].memory_base_reg=c_memory_base;

        c_multiplier.head_size=gcd_head_size+c_multiplier.extra_head_size;
        c_multiplier.c_scheduler.num_work_per_call=gcd_background_work_per_call;

        c_multiplier.matrix_buffer[0]=make_pair(regs.get_vector(), reg_spill());
        c_multiplier.matrix_buffer[1]=make_pair(regs.get_vector(), reg_spill());
        c_multiplier.matrix_buffer[2]=make_pair(regs.get_vector(), reg_spill());
        c_multiplier.matrix_buffer[3]=make_pair(regs.get_vector(), reg_spill());

        c_multiplier.carry_accumulator=c_multiplier.matrix_buffer[0].first;
        c_multiplier.carry_mask_reg=c_multiplier.matrix_buffer[1].first;

        array<array<reg_spill, 2>, num_vectors> c_matrix_spill;

        array<reg_vector, 2> c_matrix={regs.get_vector(), regs.get_vector()};
        m.bind(c_matrix, "c_matrix");

        vector<reg_spill> spills;

        for (int z=0;z<2;++z) {
            for (int x=0;x<4;++x) {
                for (int y=0;y<2;++y) {
                    spills.push_back(regs.get_spill());
                    c_multiplier.current_batch[z][x][y]=spills.back();
                }
            }

            for (int x=0;x<4;++x) {
                for (int y=0;y<2;++y) {
                    spills.push_back(regs.get_spill());
                    c_multiplier.previous_batch[z][x][y]=spills.back();
                }
            }

            for (int y=0;y<2;++y) {
                spills.push_back(regs.get_spill());
                c_matrix_spill[z][y]=spills.back();
            }
        }

        reg_spill c_memory_base_spill=regs.bind_spill(m, "c_memory_base_spill");
        spills.push_back(c_memory_base_spill);

        m.bind(spills, "spills");
        m.bind(c_multiplier.carry_accumulator, "tmp");
        m.bind(c_matrix_spill, "c_matrix_spill");

        for (int x=0;x<spills.size();++x) {
            APPEND_M(str( "VMOVDQU `tmp, #", args[x*4] ));
            APPEND_M(str( "VMOVDQU `spills_#, `tmp", x ));
        }

        APPEND_M(str( "MOV `c_memory_base, `c_memory_base_spill" ));

        if (asm_func_multiplier_gcd_mode>=0 && asm_func_multiplier_gcd_mode<=3) {
            for (int z=0;z<2;++z) {
                APPEND_M(str( "VMOVDQU `c_matrix_0, `c_matrix_spill_#_0", z ));
                APPEND_M(str( "VMOVDQU `c_matrix_1, `c_matrix_spill_#_1", z ));

                if (z==0) {
                    c_multiplier.multiply(regs, c_matrix, reg_scalar(), asm_func_multiplier_gcd_mode);
                } else {
                    c_multiplier.multiply_prepare(regs, {c_matrix}, asm_func_multiplier_gcd_mode);
                }
            }
        } else {
            assert(asm_func_multiplier_gcd_mode==4);
            c_multiplier.generate_background_work();
            c_multiplier.advance(regs, true);
        }

        APPEND_M(str( "MOV `c_memory_base_spill, `c_memory_base" ));

        for (int x=0;x<spills.size();++x) {
            APPEND_M(str( "VMOVDQU `tmp, `spills_#", x ));
            APPEND_M(str( "VMOVDQU #, `tmp", args[x*4] ));
        }
    }
#else
    extern "C" { int asm_func_multiplier_gcd_0(void* memory_base); }
    extern "C" { int asm_func_multiplier_gcd_1(void* memory_base); }
    extern "C" { int asm_func_multiplier_gcd_2(void* memory_base); }
    extern "C" { int asm_func_multiplier_gcd_3(void* memory_base); }
    extern "C" { int asm_func_multiplier_gcd_4(void* memory_base); }

    void multiplier_gcd_asm(
        array<array<simd_integer, 2>, 2>& ints,
        array<array<array<int64, 2*2>, 4>, 2>& previous_batch,
        array<array<array<int64, 2*2>, 4>, 2>& current_batch,
        array<array<int64, 2*2>, 2> c_matrix,
        int& shift_amount,
        int mode
    ) {
        assert(ints[0][0].current_size()==9*4);

        simd_integer zero;
        zero.memory.resize(2*9*4, 0);
        transfer(global_asm_ints[1], zero);
        transfer(global_asm_ints[2], zero);

        simd_integer_asm out_0_0=global_asm_ints[1];
        simd_integer_asm out_0_1=global_asm_ints[2];

        assert(shift_amount>=0 && shift_amount<9*4);

        assert(shift_amount%4==0);
        out_0_0.start+=9*4-shift_amount; //right shift; undoes the left shift
        out_0_1.start+=9*4-shift_amount;

        transfer(out_0_0, ints[0][0]);
        transfer(out_0_1, ints[0][1]);

        transfer(global_asm_ints[3], ints[1][0]);
        transfer(global_asm_ints[4], ints[1][1]);

        simd_integer args;

        const int size=2;

        for (int vector_index=0;vector_index<2;++vector_index) {
            for (int z=0;z<4;++z) {
                for (int y=0;y<size;++y) {
                    for (int x=0;x<size;++x) {
                        args.memory.push_back(current_batch[vector_index][z][y*size+x]);
                    }
                    for (int x=size;x<4;++x) {
                        args.memory.push_back(0);
                    }
                }
            }

            for (int z=0;z<4;++z) {
                for (int y=0;y<size;++y) {
                    for (int x=0;x<size;++x) {
                        args.memory.push_back(previous_batch[vector_index][z][y*size+x]);
                    }
                    for (int x=size;x<4;++x) {
                        args.memory.push_back(0);
                    }
                }
            }

            for (int y=0;y<size;++y) {
                for (int x=0;x<size;++x) {
                    args.memory.push_back(c_matrix[vector_index][y*size+x]);
                }
                for (int x=size;x<4;++x) {
                    args.memory.push_back(0);
                }
            }
        }

        args.memory.push_back(uint64(asm_memory.memory_buffer_base-8*shift_amount));
        args.memory.push_back(0);
        args.memory.push_back(0);
        args.memory.push_back(0);

        transfer(global_asm_ints[0], args);

        int res;
        if (mode==0) {
            res=asm_func_multiplier_gcd_0(asm_memory.memory_buffer_base);
        } else
        if (mode==1) {
            res=asm_func_multiplier_gcd_1(asm_memory.memory_buffer_base);
        } else
        if (mode==2) {
            res=asm_func_multiplier_gcd_2(asm_memory.memory_buffer_base);
        } else
        if (mode==3) {
            res=asm_func_multiplier_gcd_3(asm_memory.memory_buffer_base);
        } else {
            assert(mode==4);
            res=asm_func_multiplier_gcd_4(asm_memory.memory_buffer_base);
        }
        assert(res==0);

        transfer(args, global_asm_ints[0]);

        int i=0;
        for (int vector_index=0;vector_index<2;++vector_index) {
            for (int z=0;z<4;++z) {
                for (int y=0;y<size;++y) {
                    for (int x=0;x<size;++x) {
                        current_batch[vector_index][z][y*size+x]=args.memory.at(i);
                        ++i;
                    }
                    for (int x=size;x<4;++x) {
                        assert(args.memory.at(i)==0);
                        ++i;
                    }
                }
            }

            for (int z=0;z<4;++z) {
                for (int y=0;y<size;++y) {
                    for (int x=0;x<size;++x) {
                        previous_batch[vector_index][z][y*size+x]=args.memory.at(i);
                        ++i;
                    }
                    for (int x=size;x<4;++x) {
                        assert(args.memory.at(i)==0);
                        ++i;
                    }
                }
            }

            for (int y=0;y<size;++y) {
                for (int x=0;x<size;++x) {
                    c_matrix[vector_index][y*size+x]=args.memory.at(i);
                    ++i;
                }
                for (int x=size;x<4;++x) {
                    assert(args.memory.at(i)==0);
                    ++i;
                }
            }
        }

        char* c_memory_base=(char*)args.memory.at(i);
        ++i;

        uint64 memory_base_delta=asm_memory.memory_buffer_base-c_memory_base;
        assert(memory_base_delta%32==0);
        shift_amount=memory_base_delta/8;

        assert(shift_amount>=0 && shift_amount<9*4);

        out_0_0=global_asm_ints[1];
        out_0_1=global_asm_ints[2];

        out_0_0.start+=9*4-shift_amount;
        out_0_1.start+=9*4-shift_amount;

        transfer(ints[0][0], out_0_0);
        transfer(ints[0][1], out_0_1);

        transfer(ints[1][0], global_asm_ints[3]);
        transfer(ints[1][1], global_asm_ints[4]);
    }
#endif

ASM_FUNC(reduce_64_iteration,
    simd_integer_asm args=global_asm_ints[0];

    reg_scalar a_0=regs.bind_scalar(m, "a_0");
    reg_scalar a_1=regs.bind_scalar(m, "a_1");
    reg_scalar a_2=regs.bind_scalar(m, "a_2");
    reg_vector uv_0=regs.bind_vector(m, "uv_0");
    reg_vector uv_1=regs.bind_vector(m, "uv_1");
    reg_vector uv_2=regs.bind_vector(m, "uv_2");
    reg_vector max_uv=regs.bind_vector(m, "max_uv");
    reg_scalar valid=regs.bind_scalar(m, "valid");

    APPEND_M(str( "MOV `a_0, #", args[0] ));
    APPEND_M(str( "MOV `a_1, #", args[1] ));
    APPEND_M(str( "MOV `a_2, #", args[2] ));

    APPEND_M(str( "VMOVDQU `uv_0, #", args[4] ));
    APPEND_M(str( "VMOVDQU `uv_1, #", args[8] ));
    APPEND_M(str( "VMOVDQU `uv_2, #", args[12] ));
    APPEND_M(str( "VMOVDQU `max_uv, #", args[16] ));

    string early_exit_label=m.alloc_label();

    APPEND_M(str( "MOV `valid, 0" ));
    reduce_64_iteration(regs, {a_0, a_1, a_2}, {uv_0, uv_1, uv_2}, max_uv, early_exit_label);
    APPEND_M(str( "MOV `valid, 1" ));
    APPEND_M(str( "#:", early_exit_label ));

    APPEND_M(str( "MOV #, `valid", args[3] ));
    APPEND_M(str( "MOV #, `a_0", args[0] ));
    APPEND_M(str( "MOV #, `a_1", args[1] ));
    APPEND_M(str( "MOV #, `a_2", args[2] ));

    APPEND_M(str( "VMOVDQU #, `uv_0", args[4] ));
    APPEND_M(str( "VMOVDQU #, `uv_1", args[8] ));
    APPEND_M(str( "VMOVDQU #, `uv_2", args[12] ));
    APPEND_M(str( "VMOVDQU #, `max_uv", args[16] ));

)
#ifndef COMPILE_ASM
    bool reduce_64_iteration_asm(vector3& a, matrix3& uv, int64& max_uv) {
        simd_integer args;
        args.memory.resize(20);

        args.memory.at(0)=a[0];
        args.memory.at(1)=a[1];
        args.memory.at(2)=a[2];

        for (int y=0;y<3;++y) {
            args.memory.at(4+y*4+0)=uv[y*3+0];
            args.memory.at(4+y*4+1)=uv[y*3+1];
            args.memory.at(4+y*4+2)=uv[y*3+2];
            args.memory.at(4+y*4+3)=0;
        }

        args.memory.at(16)=max_uv;
        args.memory.at(17)=max_uv;
        args.memory.at(18)=max_uv;
        args.memory.at(19)=max_uv;

        simd_integer_asm args_int=global_asm_ints[0];
        transfer(args_int, args);

        asm_func_reduce_64_iteration(asm_memory.memory_buffer_base);

        transfer(args, args_int);

        a[0]=args.memory.at(0);
        a[1]=args.memory.at(1);
        a[2]=args.memory.at(2);

        for (int y=0;y<3;++y) {
            uv[y*3+0]=args.memory.at(4+y*4+0);
            uv[y*3+1]=args.memory.at(4+y*4+1);
            uv[y*3+2]=args.memory.at(4+y*4+2);
            assert(args.memory.at(4+y*4+3)==0);
        }

        max_uv=uint32(args.memory.at(16)); //upper 32 bits are invalid
        assert(max_uv==uint32(args.memory.at(17)));
        assert(max_uv==uint32(args.memory.at(18)));
        assert(max_uv==uint32(args.memory.at(19)));

        return (args.memory.at(3)==1);
    }
#endif

#ifndef COMPILE_ASM
    bool integer_reduce_asm(simd_integer& a, simd_integer& b, simd_integer& c, int64& b_sign) {
        int size=18*4;
        assert(a.current_size()==size);
        assert(b.current_size()==size);
        assert(c.current_size()==size);
        assert(b_sign==1 || b_sign==-1);

        simd_integer args;
        args.memory.resize(4);
        args.memory.at(0)=b_sign;

        simd_integer_asm args_int=global_asm_ints[0];
        args_int.size=4;

        simd_integer_asm a_int=global_asm_ints[1];
        a_int.start+=a.current_size();

        simd_integer_asm b_int=global_asm_ints[2];
        b_int.start+=b.current_size();

        simd_integer_asm c_int=global_asm_ints[3];
        c_int.start+=c.current_size();

        transfer(args_int, args);
        transfer(a_int, a);
        transfer(b_int, b);
        transfer(c_int, c);

        int res=asm_func_integer_reduce(asm_memory.memory_buffer_base);

        if (res!=0) {
            //todo print(str( "label_error_#", res ));  assert(false);
            return false;
        }

        transfer(args, args_int);
        transfer(a, a_int);
        transfer(b, b_int);
        transfer(c, c_int);

        b_sign=args.memory.at(0);

        return true;
    }
#endif

//
//

#ifdef COMPILE_ASM
    void compile_asm() {
        EXPAND_MACROS_SCOPE_PUBLIC;

        reg_alloc regs=asm_init(true);

        compile_asm_transfer_padded_9(regs);
        compile_asm_transfer_padded_18(regs);
        compile_asm_transfer_9(regs);
        compile_asm_transfer_18(regs);

        compile_asm_integer_multiply(regs);
        compile_asm_integer_fma(regs);
        compile_asm_integer_gcd(regs);
        compile_asm_divide_table(regs);
        compile_asm_gcd_64_iteration(regs);
        compile_asm_multiply_matrix_batch_2(regs);

        for (int x=0;x<5;++x) {
            compile_asm_multiplier_gcd(regs, x);
        }

        compile_asm_reduce_64_iteration(regs);
        compile_asm_integer_reduce(regs);

        compile_asm_init_memory(regs); //last

        m.output_tags=asm_output_tags;

        ofstream out( "asm_compiled.s" );
        out << m.format_res_text();
    }
#else
    void init_asm_runtime() {
        EXPAND_MACROS_SCOPE_PUBLIC;

        asm_init(false);
        asm_memory.init_buffer();
        init_memory_asm();
    }
#endif


}