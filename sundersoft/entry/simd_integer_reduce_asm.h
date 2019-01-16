namespace simd_integer_namespace {


void reduce_64_iteration(
    reg_alloc regs, array<reg_scalar, 3> a, array<reg_vector, 3> uv, reg_vector max_uv,
    string early_exit_label
) {
    EXPAND_MACROS_SCOPE;

    m.bind(a, "a");
    m.bind(uv, "uv");
    m.bind(max_uv, "max_uv");

    reg_scalar q=regs.bind_scalar(m, "q");
    reg_scalar r=regs.bind_scalar(m, "r", reg_rdx);
    reg_scalar tmp_a=regs.bind_scalar(m, "tmp_a");
    reg_scalar tmp_b=regs.bind_scalar(m, "tmp_b");
    reg_scalar new_a_1=regs.bind_scalar(m, "new_a_1");
    reg_scalar new_a_2=regs.bind_scalar(m, "new_a_2");

    //new_uv_0 = uv[2]
    reg_vector new_uv_1=regs.bind_vector(m, "new_uv_1");
    reg_vector new_uv_2=regs.bind_vector(m, "new_uv_2");
    reg_vector tmp_1=regs.bind_vector(m, "tmp_1");
    reg_vector tmp_2=regs.bind_vector(m, "tmp_2");

    //bool ac_sign_wrong=!(new_a[0]>0 && new_a[2]>0); if (ac_sign_wrong) { break; }
    APPEND_M(str( "MOV `tmp_a, `a_0" ));
    APPEND_M(str( "OR `tmp_a, `a_2" ));
    APPEND_M(str( "JS #", early_exit_label )); //early exit if either sign bit is 1

    //if (a[2]>=a[0] || (a[2]<<1)==0) { break; }
    APPEND_M(str( "CMP `a_2, `a_0" ));
    APPEND_M(str( "JGE #", early_exit_label ));

    //int64 q;
    //int64 r=divide_table(a[1]+a[2], a[2]<<1, q);
    APPEND_M(str( "MOV `tmp_a, `a_1" ));
    APPEND_M(str( "ADD `tmp_a, `a_2" ));
    APPEND_M(str( "MOV `tmp_b, `a_2" ));

    APPEND_M(str( "SHL `tmp_b, 0x1" ));
    APPEND_M(str( "JZ #", early_exit_label ));
    divide_table(regs, tmp_a, tmp_b, q, r);

    APPEND_M(str( "MOV `tmp_a, `q" ));
    APPEND_M(str( "SHL `tmp_a, #", to_hex(63-reduce_num_quotient_bits) ));
    APPEND_M(str( "SAR `tmp_a, #", to_hex(63-reduce_num_quotient_bits) ));
    APPEND_M(str( "CMP `tmp_a, `q" ));
    APPEND_M(str( "JNE #", early_exit_label )); //quotient is too big

    //vector3 new_a;
    //new_a[0]=a[2];
    //new_a[1]=a[2]-r;
    APPEND_M(str( "MOV `new_a_1, `a_2" ));
    APPEND_M(str( "SUB `new_a_1, `r" ));

    //new_a[2]=(((new_a[1]-a[1])*q)>>1) + a[0];
    APPEND_M(str( "MOV `new_a_2, `new_a_1" ));
    APPEND_M(str( "SUB `new_a_2, `a_1" )); // new_a[1]-a[1]
    APPEND_M(str( "IMUL `new_a_2, `q" )); // (new_a[1]-a[1])*q
    APPEND_M(str( "SAR `new_a_2, 0x1" )); // ((new_a[1]-a[1])*q)>>1
    APPEND_M(str( "ADD `new_a_2, `a_0" ));

    //matrix3 new_uv;
    //for (int x=0;x<3;++x) {
        //new_uv[0*3+x]=       uv[2*3+x];
        //new_uv[1*3+x]=(q<<1)*uv[2*3+x] -   uv[1*3+x];
        //new_uv[2*3+x]= (q*q)*uv[2*3+x] - q*uv[1*3+x] + uv[0*3+x];
    //}

    APPEND_M(str( "VMOVQ `tmp_1_128, `q" ));
    APPEND_M(str( "VPBROADCASTQ `tmp_1, `tmp_1_128" )); // tmp_1 = q

    APPEND_M(str( "VPADDQ `tmp_2, `tmp_1, `tmp_1" )); // tmp_2 = 2*q
    APPEND_M(str( "VPMULDQ `new_uv_1, `tmp_2, `uv_2" )); // new_uv[1]=2*q*uv[2]
    APPEND_M(str( "VPSUBQ `new_uv_1, `new_uv_1, `uv_1" )); // new_uv[1]=2*q*uv[2] - uv[1]

    APPEND_M(str( "VPMULDQ `tmp_2, `tmp_1, `tmp_1" )); // tmp_2 = q*q
    APPEND_M(str( "VPMULDQ `new_uv_2, `tmp_2, `uv_2" )); // new_uv[2]=q*q*uv[2]
    APPEND_M(str( "VPMULDQ `tmp_2, `tmp_1, `uv_1" )); // tmp_2 = q*uv[1]
    APPEND_M(str( "VPSUBQ `new_uv_2, `new_uv_2, `tmp_2" )); // new_uv[2]=q*q*uv[2] - q*uv[1]
    APPEND_M(str( "VPADDQ `new_uv_2, `new_uv_2, `uv_0" )); // new_uv[2]=q*q*uv[2] - q*uv[1] + uv[0]

    //overflow checking (can't calculate max_uv otherwise)
    //this also does the "max_uv<=data_mask" check
    APPEND_M(str( "VPADDQ `tmp_2, `new_uv_1, #", asm_immediate(1ull<<data_size) ));
    APPEND_M(str( "VPTEST `tmp_2, #", asm_immediate(carry_mask & (~(1ull<<data_size))) ));
    APPEND_M(str( "JNZ #", early_exit_label ));

    APPEND_M(str( "VPADDQ `tmp_2, `new_uv_2, #", asm_immediate(1ull<<data_size) ));
    APPEND_M(str( "VPTEST `tmp_2, #", asm_immediate(carry_mask & (~(1ull<<data_size))) ));
    APPEND_M(str( "JNZ #", early_exit_label ));

    //for (int x=0;x<9;++x) {
        //max_uv=max(max_uv, abs_int(new_uv[x]));
    //}

    //already have new_uv[0] in the accumulator because it is uv[2]
    APPEND_M(str( "VPABSD `tmp_2, `new_uv_1" )); //unused lanes are 0 or 1
    APPEND_M(str( "VPMAXUD `max_uv, `max_uv, `tmp_2" ));
    APPEND_M(str( "VPABSD `tmp_2, `new_uv_2" ));
    APPEND_M(str( "VPMAXUD `max_uv, `max_uv, `tmp_2" ));

    //max reduce

    //max of lanes 0/1 in lanes 0 and 1 ; max of lanes 2/3 in lanes 2 and 3
    APPEND_M(str( "VPERMQ `tmp_2, `max_uv, #", vpermq_mask({1, 0, 3, 2}) ));
    APPEND_M(str( "VPMAXUD `max_uv, `max_uv, `tmp_2" ));

    //max of 0/1 and 2/3 in lanes 0 and 1 ; max of 2/3 and 0/1 in lanes 2 and 3 (i.e. max of all lanes in each lane)
    APPEND_M(str( "VPERMQ `tmp_2, `max_uv, #", vpermq_mask({2, 3, 0, 1}) ));
    APPEND_M(str( "VPMAXUD `max_uv, `max_uv, `tmp_2" ));

    APPEND_M(str( "VMOVD `tmp_a_32, `max_uv_128" )); //zero extended to 64 bits

    //bool valid=( 6*(abs_int(q)+1)*max_uv < new_a[0]-abs_int(new_a[1]) ) && max_uv<=data_mask;
    // 2*(abs_int(q)+1)*max_uv + abs_int(new_a[1]) - new_a[0] < 0

    //q=|q|
    APPEND_M(str( "MOV `tmp_b, `q" ));
    APPEND_M(str( "SAR `tmp_b, #", to_hex(63) )); // mask = q>>63
    APPEND_M(str( "ADD `q, `tmp_b" )); // q = (q + mask)
    APPEND_M(str( "XOR `q, `tmp_b" )); // q = (q + mask) ^ mask;

    APPEND_M(str( "INC `q" )); // |q|+1
    //APPEND_M(str( "SHL `q, 0x1" )); // 6*(|q|+1)
    APPEND_M(str( "IMUL `q, `q, #", to_hex(6) )); // 6*(|q|+1)
    APPEND_M(str( "IMUL `q, `tmp_a" )); // 6*(|q|+1)*max_uv
    APPEND_M(str( "SUB `q, `a_2" )); // 6*(|q|+1)*max_uv - new_a[0]

    //r=|new_a[1]|
    APPEND_M(str( "MOV `r, `new_a_1" ));
    APPEND_M(str( "MOV `tmp_b, `new_a_1" ));
    APPEND_M(str( "SAR `tmp_b, #", to_hex(63) ));
    APPEND_M(str( "ADD `r, `tmp_b" ));
    APPEND_M(str( "XOR `r, `tmp_b" ));

    APPEND_M(str( "ADD `q, `r" )); // 6*(|q|+1)*max_uv + abs_int(new_a[1]) - new_a[0]
    APPEND_M(str( "JNS #", early_exit_label )); // not valid if >= 0 (i.e. sign flag is 0)

    //uv=new_uv;
    //a=new_a;
    APPEND_M(str( "MOV `a_0, `a_2" ));
    APPEND_M(str( "MOV `a_1, `new_a_1" ));
    APPEND_M(str( "MOV `a_2, `new_a_2" ));

    APPEND_M(str( "VMOVDQU `uv_0, `uv_2" ));
    APPEND_M(str( "VMOVDQU `uv_1, `new_uv_1" ));
    APPEND_M(str( "VMOVDQU `uv_2, `new_uv_2" ));
}

void reduce_iteration(reg_alloc regs, reg_spill terminated, reg_spill b_sign, matrix_multiplier_asm<3, 1>& c_multiplier, int pass) {
    EXPAND_MACROS_SCOPE;

    m.bind(terminated, "terminated");
    m.bind(b_sign, "b_sign");

    reg_scalar head_bits_end_index=regs.bind_scalar(m, "head_bits_end_index");

    array<reg_scalar, 3> head_bits={
        regs.bind_scalar(m, "head_bits_0"),
        regs.bind_scalar(m, "head_bits_1"),
        regs.bind_scalar(m, "head_bits_2")
    };

    array<reg_vector, 3> reduce_64_res={
        regs.bind_vector(m, "reduce_64_res_0"),
        regs.bind_vector(m, "reduce_64_res_1"),
        regs.bind_vector(m, "reduce_64_res_2")
    };
    reg_vector max_uv=regs.bind_vector(m, "max_uv");

    reg_scalar reduce_64_num_iterations=regs.bind_scalar(m, "reduce_64_num_iterations");

    array<simd_integer_asm, 3> head=c_multiplier.calculate_head();

    extract_bits_shifted<3>(regs, head_bits, head_bits_end_index, head);

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(b_sign, "b_sign");
        m.bind(head_bits, "head_bits");

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

        APPEND_M(str( "MOV `tmp, `b_sign" ));
        APPEND_M(str( "IMUL `head_bits_1, `tmp" ));
    }

    asm_immediate.assign(reduce_64_res[0], {1, 0, 0, 0});
    asm_immediate.assign(reduce_64_res[1], {0, 1, 0, 0});
    asm_immediate.assign(reduce_64_res[2], {0, 0, 1, 0});
    asm_immediate.assign(max_uv, 1);

    schedule_early_exit(
        regs,
        [&](reg_alloc c_regs, string early_exit_label, int iteration) {
            reduce_64_iteration(c_regs, head_bits, reduce_64_res, max_uv, early_exit_label);
        },
        reduce_num_iterations,
        reduce_num_iterations_background_work,
        c_multiplier.c_scheduler,
        reduce_64_num_iterations
    );

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(terminated, "terminated");
        m.bind(b_sign, "b_sign");
        m.bind(head_bits, "head_bits");
        m.bind(reduce_64_res, "reduce_64_res");
        m.bind(reduce_64_num_iterations, "reduce_64_num_iterations");

        reg_scalar tmp_a=regs_copy.bind_scalar(m, "tmp_a");
        reg_scalar tmp_b=regs_copy.bind_scalar(m, "tmp_b");

        reg_vector tmp_1=regs_copy.bind_vector(m, "tmp_1");
        reg_vector tmp_2=regs_copy.bind_vector(m, "tmp_2");

        //need to negate second column if b_sign is negative
        //need to negate second row if b_sign_new is negative

        //multiply second column by b_sign
        APPEND_M(str( "VPBROADCASTQ `tmp_1, `b_sign" ));
        asm_immediate.assign(tmp_2, 1);
        APPEND_M(str( "VPBLENDD `tmp_1, `tmp_2, `tmp_1, #", vpblendd_mask_4({0, 1, 0, 0}) ));
        APPEND_M(str( "VPMULDQ `reduce_64_res_0, `reduce_64_res_0, `tmp_1" ));
        APPEND_M(str( "VPMULDQ `reduce_64_res_1, `reduce_64_res_1, `tmp_1" ));
        APPEND_M(str( "VPMULDQ `reduce_64_res_2, `reduce_64_res_2, `tmp_1" ));

        APPEND_M(str( "SAR `head_bits_1, #", to_hex(63) )); //-1 if negative, 0 if positive
        APPEND_M(str( "SHL `head_bits_1, 0x1" )); //-2 if negative, 0 if positive
        APPEND_M(str( "INC `head_bits_1" )); //-1 if negative, 1 if positive (b_sign_new)

        //multiply second row by b_sign_new
        APPEND_M(str( "VMOVQ `tmp_1_128, `head_bits_1" ));
        APPEND_M(str( "VPBROADCASTQ `tmp_1, `tmp_1_128" ));
        APPEND_M(str( "VPMULDQ `reduce_64_res_1, `reduce_64_res_1, `tmp_1" ));

        if (pass==3) {
            string skip_assign_terminated_label=m.alloc_label();

            APPEND_M(str( "CMP `reduce_64_num_iterations, 0" ));
            APPEND_M(str( "JNE #", skip_assign_terminated_label ));

            asm_immediate.assign(tmp_a, 1);
            APPEND_M(str( "MOV `terminated, `tmp_a" ));

            APPEND_M(str( "#:", skip_assign_terminated_label ));
        }

        c_multiplier.multiply(regs, {reduce_64_res}, head_bits[1], pass);

        APPEND_M(str( "MOV `b_sign, `head_bits_1" ));
    }
}

//large buffers: same size as a/b/c
//small buffers: 12 limbs
void reduce(
    reg_alloc& regs,
    simd_integer_asm a_padded, simd_integer_asm b_padded, simd_integer_asm c_padded, reg_spill b_sign,
    array<simd_integer_asm, 3> large_buffers, array<simd_integer_asm, 6> small_buffers
) {
    EXPAND_MACROS_SCOPE;

    reg_spill terminated=regs.bind_spill(m, "terminated");
    reg_spill loop_count=regs.bind_spill(m, "loop_count");

    simd_integer_asm a=remove_padding(regs, a_padded);
    simd_integer_asm b=remove_padding(regs, b_padded);
    simd_integer_asm c=remove_padding(regs, c_padded);

    assert(a.size==b.size && a.size==c.size);

    assert(large_buffers[0].size==a.size);
    assert(large_buffers[1].size==a.size);
    assert(large_buffers[2].size==a.size);
    assert(small_buffers[0].size==12);
    assert(small_buffers[1].size==12);
    assert(small_buffers[2].size==12);
    assert(small_buffers[3].size==12);
    assert(small_buffers[4].size==12);
    assert(small_buffers[5].size==12);

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(terminated, "terminated");
        m.bind(loop_count, "loop_count");

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

        asm_immediate.assign(tmp, 0);
        APPEND_M(str( "MOV `terminated, `tmp" ));
        APPEND_M(str( "MOV `loop_count, `tmp" ));
    }

    matrix_multiplier_asm<3, 1> c_multiplier;
    c_multiplier.ints[0]={a, b, c};
    c_multiplier.ints_buffer={large_buffers[0], large_buffers[1], large_buffers[2]};
    c_multiplier.extra_ints={small_buffers[0], small_buffers[1], small_buffers[2]};
    c_multiplier.extra_ints_multiplied={small_buffers[3], small_buffers[4], small_buffers[5]};

    c_multiplier.init(regs, reduce_head_size, reduce_num_matrix_buffer_spills, reduce_background_work_per_call);
    c_multiplier.do_initial_shrink(regs);

    print( "reduce main loop:" );

    const vector<int> cutoffs={0, 16}; //sorted

    string loop_label=m.alloc_label();
    string loop_exit_label=m.alloc_label();

    vector<string> cutoff_labels;
    for (int x=0;x<cutoffs.size();++x) {
        cutoff_labels.push_back(m.alloc_label());
    }

    APPEND_M(str( "JMP #", loop_label ));

    string jump_table=m.alloc_label();
    {
        APPEND_M(str( ".balign 8" ));
        APPEND_M(str( "#:", jump_table ));
        for (int cutoff=0;cutoff<c_multiplier.int_size;++cutoff) {

            int cutoff_index=-1;
            for (int x=0;x<cutoffs.size();++x) {
                //static cutoff can be less than or equal to the calculated value
                if (cutoffs[x]>cutoff) {
                    break;
                }
                cutoff_index=x;
            }
            assert(cutoff_index!=-1);

            APPEND_M(str( ".quad #", cutoff_labels.at(cutoff_index) ));
        }
    }

    APPEND_M(str( "#:", loop_label ));
    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        reg_scalar cutoff_8=regs_copy.bind_scalar(m, "cutoff_8");
        reg_scalar table_address=regs_copy.bind_scalar(m, "table_address");
        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

        c_multiplier.get_cutoff_8(regs_copy, cutoff_8);

        string rip_label=m.alloc_label();

        //if you try to do this in the LEA instruction then the assembler will add garbage instead of the actual value
        //trying to disable position independent code makes the assembly code not work anymore; calling conventions might be different
        // or something
        APPEND_M(str( "MOV `tmp, #-#", jump_table, rip_label ));

        APPEND_M(str( "LEA `table_address, [RIP]" ));
        APPEND_M(str( "#:", rip_label ));

        APPEND_M(str( "ADD `table_address, `tmp" ));

        APPEND_M(str( "JMP [`table_address+`cutoff_8]" ));
    }

    for (int y=0;y<cutoffs.size();++y) {
        APPEND_M(str( "#:", cutoff_labels[y] ));
        c_multiplier.cutoff=cutoffs[y];

        {
            EXPAND_MACROS_SCOPE;
            reg_alloc regs_copy=regs;

            m.bind(loop_count, "loop_count");

            reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

            APPEND_M(str( "MOV `tmp, `loop_count" ));
            APPEND_M(str( "INC `tmp" ));
            APPEND_M(str( "MOV `loop_count, `tmp" ));
            APPEND_M(str( "CMP `tmp, #", to_hex(1000) ));
            APPEND_M(str( "JG #", m.alloc_error_label() ));
        }

        c_multiplier.generate_background_work();
        for (int x=0;x<4;++x) {
            reduce_iteration(regs, terminated, b_sign, c_multiplier, x);
        }
        c_multiplier.advance(regs, true);

        {
            EXPAND_MACROS_SCOPE;
            reg_alloc regs_copy=regs;

            m.bind(terminated, "terminated");

            reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");

            APPEND_M(str( "MOV `tmp, `terminated" ));
            APPEND_M(str( "CMP `tmp, 0" ));
            APPEND_M(str( "JE #", loop_label ));
            APPEND_M(str( "JMP #", loop_exit_label ));
        }
    }
    APPEND_M(str( "#:", loop_exit_label ));

    print( "reduce finish:" );

    c_multiplier.cutoff=0;
    c_multiplier.generate_background_work();
    for (int pass=0;pass<4;++pass) {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        array<reg_vector, 3> c_matrix={
            regs_copy.bind_vector(m, "c_matrix_0"),
            regs_copy.bind_vector(m, "c_matrix_1"),
            regs_copy.bind_vector(m, "c_matrix_2")
        };
        asm_immediate.assign(c_matrix[0], {1, 0, 0, 0});
        asm_immediate.assign(c_matrix[1], {0, 1, 0, 0});
        asm_immediate.assign(c_matrix[2], {0, 0, 1, 0});

        c_multiplier.multiply(regs_copy, c_matrix, reg_scalar(), pass);
    }
    c_multiplier.advance(regs, false);
}


}