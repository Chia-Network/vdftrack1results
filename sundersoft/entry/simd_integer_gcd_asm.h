namespace simd_integer_namespace {


void generate_divide_table(bool assign_values) {
    const int expected_size=1<<divide_table_index_bits;

    if (assign_values) {
        vector<uint64> res;

        const int max_index=bit_sequence(0, divide_table_index_bits);
        assert(max_index>=1);

        res.push_back(0);
        for (int index=1;index<=max_index;++index) {
            uint128 v = (~uint128(0)) / uint128(index);
            v>>=64;
            res.push_back(v);
        }

        assert(res.size()==expected_size);
        assert(asm_memory.alloc(res.size()*8, 8, (char*)&res[0])==0);
    } else {
        assert(asm_memory.alloc(expected_size*8, 8)==0);
    }
}

//regs: 1x scalar (RAX) + 4x scalar arguments (r==RDX)
//todo //test hit rate
void divide_table(reg_alloc regs, reg_scalar a, reg_scalar b, reg_scalar q, reg_scalar r) {
    EXPAND_MACROS_SCOPE;

    expand_macros_tag c_tag(m, "divide");

    regs.get_scalar(reg_rax);

    m.bind(a, "a");
    m.bind(b, "b");
    m.bind(q, "q");
    assert(r.value==reg_rdx.value);

    string b_shift_label=m.alloc_label();
    APPEND_M(str( "BSR `q, `b" )); // b_shift = bsr(b)
    APPEND_M(str( "SUB `q, #", to_hex(divide_table_index_bits-1) )); // b_shift = bsr(b)-(divide_table_index_bits-1)
    APPEND_M(str( "JNB #", b_shift_label ));
    asm_immediate.assign(q, 0); // if (b_shift<0) b_shift=0
    APPEND_M(str( "#:", b_shift_label ));

    APPEND_M(str( "SARX RAX, `b, `q" )); // b_approx = b>>b_shift
    APPEND_M(str( "MOV RAX, [`memory_base+RAX*8]" )); // b_approx_inverse = memory_base[b_approx]

    APPEND_M(str( "IMUL `a" )); // q = (b_approx_inverse*a)>>64
    APPEND_M(str( "SARX `q, RDX, `q" )); // q = q>>b_shift

    string wrong_remainder_label=m.alloc_label();
    APPEND_M(str( "MOV RAX, `q" ));
    APPEND_M(str( "IMUL RAX, `b" )); // r = q*b
    APPEND_M(str( "JO #", wrong_remainder_label )); // overflow
    APPEND_M(str( "MOV RDX, `a" ));
    APPEND_M(str( "SUB RDX, RAX" )); // r = a-q*b
    APPEND_M(str( "JO #", wrong_remainder_label )); // overflow

    APPEND_M(str( "CMP RDX, `b" ));
    APPEND_M(str( "JAE #", wrong_remainder_label )); // !(r>=0 && r<b)

    string end_label=m.alloc_label();
    APPEND_M(str( "JMP #", end_label ));

    if (!asm_output_common_case_only) {
        APPEND_M(str( "#:", wrong_remainder_label ));

        APPEND_M(str( "MOV RDX, `a" ));
        APPEND_M(str( "SAR RDX, #", to_hex(63) )); //all 1s if negative, all 0s if nonnegative

        APPEND_M(str( "MOV RAX, `a" ));
        APPEND_M(str( "IDIV `b" )); // RAX=a/b ; RDX=r=a%b
        APPEND_M(str( "MOV `q, RAX" ));
        APPEND_M(str( "CMP RDX, 0" ));
        APPEND_M(str( "JGE #", end_label )); // r>=0
        APPEND_M(str( "ADD RDX, `b" )); // r+=b
        APPEND_M(str( "DEC `q" ));
    }

    APPEND_M(str( "#:", end_label ));
}

//regs: 1x scalar argument
void schedule_early_exit(
    reg_alloc regs, function<void(reg_alloc regs, string early_exit_label, int iteration)> do_iteration,
    int num_iterations, int num_iterations_background_work, scheduler& c_scheduler, reg_scalar num_iterations_reg
) {
    EXPAND_MACROS_SCOPE;

    assert(num_iterations>=1 && num_iterations_background_work<=num_iterations);

    m.bind(num_iterations_reg, "num_iterations");

    vector<string> early_exit_labels;

    vector<vector<vector<string>>> background_work_recordings;

    for (int x=0;x<num_iterations;++x) {
        string early_exit_label=m.alloc_label();
        early_exit_labels.push_back(early_exit_label);

        do_iteration(regs, early_exit_label, x);
        if (x<num_iterations_background_work) {
            expand_macros_recording c;
            m.begin_recording(c);
            c_scheduler(regs);
            background_work_recordings.push_back(m.end_recording(c));
        }
    }

    APPEND_M(str( "MOV `num_iterations, #", num_iterations ));

    string end_label=m.alloc_label();
    APPEND_M(str( "JMP #", end_label ));

    if (!asm_output_common_case_only) {
        vector<string> background_work_labels;
        for (int x=0;x<num_iterations;++x) {
            APPEND_M(str( "#:", early_exit_labels.at(x) ));
            APPEND_M(str( "MOV `num_iterations, #", x ));
            if (x<num_iterations_background_work) {
                background_work_labels.push_back(m.alloc_label());
                APPEND_M(str( "JMP #", background_work_labels.back() ));
            } else {
                APPEND_M(str( "JMP #", end_label ));
            }
        }

        assert(background_work_labels.size()==num_iterations_background_work);

        for (int x=0;x<num_iterations_background_work;++x) {
            APPEND_M(str( "#:", background_work_labels.at(x) ));
            m.append_recording(background_work_recordings.at(x));
        }
    }

    APPEND_M(str( "#:", end_label ));
}

const array<uint64, 4> gcd_mask_approximate={1ull<<63, 1ull<<63, 0, 0};
const array<uint64, 4> gcd_mask_exact={0, 0, 0, 0};

//regs: 3x scalar, 3x vector, 2x scalar argument, 2x vector argument
void gcd_64_iteration(
    reg_alloc regs, reg_spill c_gcd_mask, array<reg_scalar, 2> a, array<reg_vector, 2> uv,
    string early_exit_label
) {
    EXPAND_MACROS_SCOPE;

    expand_macros_tag c_tag(m, "gcd_64");

    m.bind(c_gcd_mask, "c_gcd_mask");
    m.bind(a, "a");
    m.bind(uv, "uv");

    reg_scalar q=regs.bind_scalar(m, "q");
    reg_scalar r=regs.bind_scalar(m, "r", reg_rdx);

    reg_scalar tmp_a=regs.bind_scalar(m, "tmp_a");

    //new_uv_0 = uv[1]
    reg_vector new_uv_1=regs.bind_vector(m, "new_uv_1");
    reg_vector tmp_1=regs.bind_vector(m, "tmp_1");
    reg_vector tmp_2=regs.bind_vector(m, "tmp_2");

    APPEND_M(str( "CMP `a_1, 0" ));
    APPEND_M(str( "JZ #", early_exit_label ));

    divide_table(regs, a[0], a[1], q, r);
    APPEND_M(str( "MOV `tmp_a, `q" ));
    APPEND_M(str( "SHL `tmp_a, #", to_hex(63-gcd_num_quotient_bits) ));
    APPEND_M(str( "SAR `tmp_a, #", to_hex(63-gcd_num_quotient_bits) ));
    APPEND_M(str( "CMP `tmp_a, `q" ));
    APPEND_M(str( "JNE #", early_exit_label )); //quotient is too big

    APPEND_M(str( "MOV `a_0, `a_1" ));
    APPEND_M(str( "MOV `a_1, `r" ));

    APPEND_M(str( "VMOVQ `new_uv_1_128, `q" ));
    APPEND_M(str( "VPBROADCASTQ `new_uv_1, `new_uv_1_128" )); // new_uv_1 = q

    APPEND_M(str( "VPMULDQ `new_uv_1, `new_uv_1, `uv_1" )); // new_uv_1 = q*uv[1]
    APPEND_M(str( "VPSUBQ `new_uv_1, `uv_0, `new_uv_1" )); // new_uv_1 = uv[0] - q*uv[1]

    //overflow checking:
    //-the carry_mask bits must be all 0s or all 1s for each 64-bit entry
    //-if 1<<data_size is added, the carry_mask bits must be all 0s (negative) or 1<<data_size with the rest 0 (nonnegative)
    //-can add 1<<data_size, then check the carry_mask except the last bit
    APPEND_M(str( "VPADDQ `tmp_1, `new_uv_1, #", asm_immediate(1ull<<data_size) ));
    APPEND_M(str( "VPTEST `tmp_1, #", asm_immediate(carry_mask & (~(1ull<<data_size))) ));
    APPEND_M(str( "JNZ #", early_exit_label ));

    {
        APPEND_M(str( "VMOVQ `tmp_1_128, `a_0" ));
        APPEND_M(str( "VPBROADCASTQ `tmp_1, `tmp_1_128" )); // tmp_1 = a[0]
        APPEND_M(str( "VPADDQ `tmp_1, `tmp_1, `uv_1" )); // tmp_1 = a[0]+new_uv[0]

        APPEND_M(str( "VMOVQ `tmp_2_128, `a_1" ));
        APPEND_M(str( "VPBROADCASTQ `tmp_2, `tmp_2_128" )); // tmp_2 = a[1]
        APPEND_M(str( "VPADDQ `tmp_2, `tmp_2, `new_uv_1" )); // tmp_2 = a[1]+new_uv[1]

        APPEND_M(str( "VPSUBQ `tmp_1, `tmp_1, `tmp_2" )); // tmp_1 = a[0]+new_uv[0]-(a[1]+new_uv[1])

        APPEND_M(str( "VPOR `tmp_1, `tmp_1, `tmp_2" )); // sign is 1 if tmp_1<0 or tmp_2<0

        //approximate: ZF set if both signs of tmp_1 are 0 (i.e tmp_1>=0 and tmp_2>=0 for both lanes)
        //exact: ZF set always
        APPEND_M(str( "VPTEST `tmp_1, `c_gcd_mask" ));

        APPEND_M(str( "JNZ #", early_exit_label )); //taken if ZF==0

        //int64 delta=new_a[0]-new_a[1];
        //if (new_a[1]<-new_uv[1]) goto early_exit_label
        //if (delta<new_uv[1]-new_uv[0]) goto early_exit_label
        //if (new_a[1]+new_uv[1]<0) goto early_exit_label
        //if (new_a[0]+new_uv[0]-(new_a[1]+new_uv[1])<0) goto early_exit_label
    }

    APPEND_M(str( "VMOVDQU `uv_0, `uv_1" ));
    APPEND_M(str( "VMOVDQU `uv_1, `new_uv_1" ));
}

template<unsigned long int size> void assign_matrix_buffer(
    reg_alloc regs,
    array<pair<reg_vector, reg_spill>, size*size> matrix_buffer,
    array<reg_spill, size> c_matrix
) {
    EXPAND_MACROS_SCOPE;

    bind_spillable(matrix_buffer, "matrix_buffer");
    m.bind(c_matrix, "c_matrix");

    reg_vector tmp;

    for (int y=0;y<size;++y) {
        for (int x=0;x<size;++x) {
            bool out_is_reg=(matrix_buffer[y*size+x].first.value!=-1);

            if (out_is_reg) {
                APPEND_M(str( "VPBROADCASTQ `matrix_buffer_#, `c_matrix_#_#", y*size+x, y, x ));
            } else {
                if (tmp.value==-1) {
                    tmp=regs.bind_vector(m, "tmp");
                }

                APPEND_M(str( "VPBROADCASTQ `tmp, `c_matrix_#_#", y, x ));
                APPEND_M(str( "VMOVDQU `matrix_buffer_#, `tmp", y*size+x ));
            }
        }
    }
}

template<unsigned long int size> void assign_identity(
    reg_alloc regs,
    array<array<reg_spill, size>, 4> res
) {
    EXPAND_MACROS_SCOPE;

    m.bind(res, "res");

    reg_vector tmp=regs.bind_vector(m, "tmp");

    for (int z=0;z<4;++z) {
        for (int y=0;y<size;++y) {
            array<uint64, 4> data={0, 0, 0, 0};
            if (z==0) {
                data.at(y)=1;
            }

            asm_immediate.assign(tmp, data);
            APPEND_M(str( "VMOVDQU `res_#_#, `tmp", z, y ));
        }
    }
}

template<unsigned long int size> void assign_current_batch(
    reg_alloc regs,
    array<reg_spill, size> res,
    array<reg_vector, size> c_matrix
) {
    EXPAND_MACROS_SCOPE;

    m.bind(res, "res");
    m.bind(c_matrix, "c_matrix");

    for (int x=0;x<size;++x) {
        APPEND_M(str( "VMOVDQU `res_#, `c_matrix_#", x, x ));
    }
}

template<unsigned long int size> void bitwise_or_ints(
    reg_alloc regs,
    reg_vector accumulator,
    array<simd_integer_asm, size> ints,
    int start
) {
    EXPAND_MACROS_SCOPE;

    assert(start>=0 && start<=ints[0].size-4);

    m.bind(accumulator, "accumulator");

    for (int x=0;x<size;++x) {
        string cmd=(x==0)? "VMOVDQU `accumulator, " : "VPOR `accumulator, `accumulator, ";
        APPEND_M(str( "# #", cmd, ints[x][start] ));
    }
}

template<unsigned long int size> void multiply_matrix_batch(
    reg_alloc regs,
    array<reg_spill, size> res, //can alias with current_batch
    array<array<reg_spill, size>, 2> current_batch
) {
    EXPAND_MACROS_SCOPE;

    m.bind(res, "res");
    m.bind(current_batch, "current_batch");

    array<reg_vector, size> buffer;
    for (int x=0;x<size;++x) {
        buffer[x]=regs.bind_vector(m, str( "buffer_#", x ));
    }

    reg_vector a=regs.bind_vector(m, "a");
    reg_vector b=regs.bind_vector(m, "b");

    //outer product of column z of current_batch[0] and row z of current_batch[1], accumulated into buffer
    for (int z=0;z<size;++z) {
        APPEND_M(str( "VMOVDQU `b, `current_batch_1_#", z )); //row z of current_batch[1]
        for (int y=0;y<size;++y) {
            APPEND_M(str( "VPBROADCASTQ `a, `current_batch_0_#_#", y,z )); //column z of current_batch[0]

            string out_reg=(z==0)? str( "`buffer_#", y ) : "`a";
            APPEND_M(str( "VPMULDQ #, `a, `b", out_reg )); //row y of the outer product

            if (z!=0) {
                APPEND_M(str( "VPADDQ `buffer_#, `buffer_#, `a", y, y )); //accumulated outer products
            }
        }
    }

    for (int x=0;x<size;++x) {
        APPEND_M(str( "VMOVDQU `res_#, `buffer_#", x, x ));
    }
}

template<unsigned long int size> void multiply_matrix_batch(
    reg_alloc regs,
    array<array<reg_spill, size>, 4> res,
    array<array<reg_spill, size>, 4> current_batch
) {
    EXPAND_MACROS_SCOPE;

    expand_macros_tag c_tag(m, "batch");

    multiply_matrix_batch(regs, current_batch[0], {current_batch[0], current_batch[1]});
    multiply_matrix_batch(regs, current_batch[1], {current_batch[2], current_batch[3]});

    m.bind(res, "res");
    m.bind(current_batch, "current_batch");

    reg_scalar accumulator_low=regs.bind_scalar(m, "accumulator_low");
    reg_scalar accumulator_high=regs.bind_scalar(m, "accumulator_high");
    regs.get_scalar(reg_rax);
    regs.get_scalar(reg_rdx);

    reg_scalar data_mask_reg=regs.bind_scalar(m, "data_mask");
    asm_immediate.assign(data_mask_reg, data_mask);

    array<reg_scalar, size> buffer;
    for (int x=0;x<size;++x) {
        buffer[x]=regs.bind_scalar(m, str( "buffer_#", x ));
    }

    //will load a row from current_batch[0] into scalar registers and then calculate its dot product with each column of
    //current_batch[1], using two accumulator registers (128 bits)
    //each result is 128 bits and is split into 4 pieces (each of data_size) and outputted. the last piece is sign-extended to 64 bits

    for (int y=0;y<size;++y) {
        for (int x=0;x<size;++x) {
            APPEND_M(str( "MOV `buffer_#, `current_batch_0_#_#", x, y,x ));
        }

        for (int x=0;x<size;++x) {
            for (int z=0;z<size;++z) {
                APPEND_M(str( "MOV RAX, `buffer_#", z ));
                APPEND_M(str( "IMUL QWORD PTR `current_batch_1_#_#", z,x ));

                if (z==0) {
                    APPEND_M(str( "MOV `accumulator_low, RAX" ));
                    APPEND_M(str( "MOV `accumulator_high, RDX" ));
                } else {
                    APPEND_M(str( "ADD `accumulator_low, RAX" ));
                    APPEND_M(str( "ADC `accumulator_high, RDX" ));
                }
            }

            //final result for matrix position [y,x] is in the accumulator
            //limb 0: accumulator_low & data_mask
            //limb 1: (accumulator_low >> data_size) & data_mask [logical right shift]
            //limb 2: (SHRD [accumulator_high:accumulator_low] >> 2*data_size) & data_mask [logical right shift]
            //limb_3: accumulator_high >> (3*data_size-64) [arithmetic right shift]

            //todo //can use simd for this if i can be bothered. might need to shuffle data around

            APPEND_M(str( "MOV RAX, `accumulator_low" ));
            APPEND_M(str( "AND RAX, `data_mask" ));
            APPEND_M(str( "MOV `res_0_#_#, RAX", y,x ));

            APPEND_M(str( "MOV RAX, `accumulator_low" ));
            APPEND_M(str( "SHR RAX, #", to_hex(data_size) ));
            APPEND_M(str( "AND RAX, `data_mask" ));
            APPEND_M(str( "MOV `res_1_#_#, RAX", y,x ));

            APPEND_M(str( "MOV RAX, `accumulator_low" ));
            APPEND_M(str( "SHRD RAX, `accumulator_high, #", to_hex(2*data_size) ));
            APPEND_M(str( "AND RAX, `data_mask" ));
            APPEND_M(str( "MOV `res_2_#_#, RAX", y,x ));

            APPEND_M(str( "MOV RAX, `accumulator_high" ));
            APPEND_M(str( "SAR RAX, #", to_hex(3*data_size-64) ));
            APPEND_M(str( "MOV `res_3_#_#, RAX", y,x ));
        }
    }
}

template<int size, int num_vectors> struct matrix_multiplier_asm {
    const int extra_head_size=8;
    const int carry_adjusted_head_size=12; //will add data_sign_mask to this many msb limbs of the head before carrying

    //first vector is the primary and is shifted and shrunk, and has a head
    //secondary vectors are not shifted or shrunk, and have no head
    array<array<simd_integer_asm, size>, num_vectors> ints;
    array<simd_integer_asm, size> ints_buffer; //not shifted

    array<simd_integer_asm, size> extra_ints; //12 limbs
    array<simd_integer_asm, size> extra_ints_multiplied; //12 limbs

    array<array<array<reg_spill, size>, 4>, num_vectors> previous_batch; //lsb first
    array<array<array<reg_spill, size>, 4>, num_vectors> current_batch;

    int matrix_buffer_state=-1; //-1: invalid, 0-3: previous_batch
    bool preserve_carry_mask_and_accumulator=false;

    array<pair<reg_vector, reg_spill>, size*size> matrix_buffer;
    reg_vector carry_accumulator; //aliases with matrix_buffer
    reg_vector carry_mask_reg; //aliases with matrix_buffer

    int int_size=-1;
    int head_size=-1; //last 4 entries are potentially invalid; rest are valid

    scheduler c_scheduler;

    int cutoff=0; //this many lsb limbs in the tail are 0 (only for the first vector)

    void bind() {
        m.bind(carry_accumulator, "carry_accumulator");
        m.bind(carry_mask_reg, "carry_mask");
        m.bind(ints[0][0].memory_base_reg, "c_memory_base");

        m.bind(previous_batch, "previous_batch");
        m.bind(current_batch, "current_batch");

        bind_spillable(matrix_buffer, "matrix_buffer");
    }

    //this is the number of shrunked off limbs multiplied by 8
    void get_cutoff_8(reg_alloc regs, reg_scalar res) {
        EXPAND_MACROS_SCOPE;

        bind();
        m.bind(res, "res");

        //for each limb that was shrunk off, c_memory_base is reduced by 8
        APPEND_M(str( "MOV `res, `memory_base" ));
        APPEND_M(str( "SUB `res, `c_memory_base" )); // tmp = memory_base - c_memory_base
        //APPEND_M(str( "SHR `res, #", to_hex(3) )); // tmp = (memory_base - c_memory_base)/8
    }

    void has_tail(reg_alloc regs, reg_scalar tmp, string branch_has_tail, string branch_no_tail, bool small_head=true) {
        EXPAND_MACROS_SCOPE;

        int c_head_size=(small_head)? head_size-extra_head_size : head_size;

        bind();
        m.bind(tmp, "tmp");

        //tmp is c_memory_base if there is no tail
        //however, there might be additional shrinks if small_head is false, so c_memory_base can also be less than tmp if there is no
        // tail
        APPEND_M(str( "MOV `tmp, `memory_base" ));
        APPEND_M(str( "SUB `tmp, #", to_hex((int_size-c_head_size)*8) ));
        APPEND_M(str( "CMP `c_memory_base, `tmp" ));

        if (!branch_has_tail.empty()) {
            assert(branch_no_tail.empty());
            APPEND_M(str( "JA #", branch_has_tail ));
        }

        if (!branch_no_tail.empty()) {
            assert(branch_has_tail.empty());
            APPEND_M(str( "JBE #", branch_no_tail ));
        }
    }

    //resulting mask is in carry_mask. invalidates carry_accumulator
    void calculate_shrink(reg_alloc regs) {
        EXPAND_MACROS_SCOPE;

        bind();

        reg_scalar tmp=regs.bind_scalar(m, "tmp");
        reg_scalar tmp2=regs.bind_scalar(m, "tmp2");

        //need the 5 msb limbs to all be 0 to do a shrink. the shrink will remove 4 msb limbs
        bitwise_or_ints<size>(regs, carry_accumulator, ints[0], int_size-4);
        bitwise_or_ints<size>(regs, carry_mask_reg, ints[0], int_size-8);
        APPEND_M(str( "VPAND `carry_mask, `carry_mask, #", asm_immediate({0, 0, 0, ~uint64(0)}) ));
        APPEND_M(str( "VPOR `carry_accumulator, `carry_accumulator, `carry_mask" ));

        //tmp=1 if 5 msb limbs are all 0, else 0
        asm_immediate.assign(tmp, 0);
        APPEND_M(str( "VPTEST `carry_accumulator, #", asm_immediate(~uint64(0)) ));
        APPEND_M(str( "SETZ `tmp_8" ));

        //can't shrink if the tail is empty
        //some of the entries in the head aren't read by the gcd/reduce code so they are treated as part of the tail
        //those entries will eventually get zeroed out
        //tmp2=1 if there is a tail, else 0
        has_tail(regs, tmp2, "", "", true);
        APPEND_M(str( "MOV `tmp2, 0" )); //can't change the flags
        APPEND_M(str( "SETA `tmp2_8" ));

        //tmp=1 if shrinking, else 0
        APPEND_M(str( "AND `tmp, `tmp2" ));

        //broadcast tmp to all 4 lanes of carry_mask (1 if shrinking, else 0)
        APPEND_M(str( "MOVQ `carry_mask_128, `tmp" ));
        APPEND_M(str( "VPBROADCASTQ `carry_mask, `carry_mask_128" ));

        //if input is 0, result is 0. if input is 1, result is all 1s
        //carry_mask is all 1s if shrinking, else all 0s
        asm_immediate.assign(carry_accumulator, 0);
        APPEND_M(str( "VPCMPGTQ `carry_mask, `carry_mask, `carry_accumulator" ));
    }

    //shrink mask is in carry_mask
    void do_shrink(reg_alloc regs) {
        EXPAND_MACROS_SCOPE;

        bind();

        reg_scalar tmp=regs.bind_scalar(m, "tmp");
        reg_scalar tmp2=regs.bind_scalar(m, "tmp2");

        asm_immediate.assign(tmp, 32);
        APPEND_M(str( "MOVQ `tmp2, `carry_mask_128" ));
        APPEND_M(str( "AND `tmp, `tmp2" )); //32 if shrinking, else 0
        APPEND_M(str( "SUB `c_memory_base, `tmp" ));
    }

    void do_initial_shrink(reg_alloc regs) {
        EXPAND_MACROS_SCOPE;

        bind();

        string start_label=m.alloc_label();
        string end_label=m.alloc_label();

        APPEND_M(str( "#:", start_label ));

        calculate_shrink(regs);

        APPEND_M(str( "VPTEST `carry_mask, #", asm_immediate(~uint64(0)) ));
        APPEND_M(str( "JZ #", end_label ));

        do_shrink(regs);

        APPEND_M(str( "JMP #", start_label ));

        APPEND_M(str( "#:", end_label ));
    }

    //to finish everything:
    //-call advance (need to call multiply 4 times before calling advance, so use identity matricies until that is true)
    //-call multiply 4 times with identity matricies
    //-call advance again
    //to reduce int_size:
    //-call advance. it will do a shrink
    //-branch based on the new value of ints[0].memory_base_reg (relative to asm_memory.memory_base)
    //initial values:
    //-set ints to the input values. each int is padded on the left so that the total size is doubled
    //-set previous_batch to the identity matrix
    void advance(reg_alloc regs, bool allow_shrink) {
        EXPAND_MACROS_SCOPE;

        expand_macros_tag c_tag(m, "advance");

        bind();

        int num_background_work_calls=0;
        while (c_scheduler.has_work()) {
            c_scheduler(regs);
            ++num_background_work_calls;
        }

        print( "num_background_work_calls:", num_background_work_calls );

        reg_scalar tmp=regs.bind_scalar(m, "tmp");
        reg_scalar tmp2=regs.bind_scalar(m, "tmp2");
        reg_spill shrink_mask_spill=regs.bind_spill(m, "shrink_mask_spill");

        {
            string skip_carry_label=m.alloc_label();
            APPEND_M(str( "VPTEST `carry_accumulator, `carry_mask" ));
            APPEND_M(str( "JZ #", skip_carry_label ));

            string continue_carry_label=m.alloc_label();

            APPEND_M(str( "#:", continue_carry_label ));

            preserve_carry_mask_and_accumulator=false;
            generate_all_carry_work();

            while (c_scheduler.has_work()) {
                c_scheduler(regs);
            }

            APPEND_M(str( "VPTEST `carry_accumulator, `carry_mask" ));
            APPEND_M(str( "JNZ #", continue_carry_label ));

            APPEND_M(str( "#:", skip_carry_label ));
        }

        matrix_buffer_state=-1;
        preserve_carry_mask_and_accumulator=false;

        //current_batch is done being generated and previous_batch is done being used, so generate new previous_batch
        //if is_identity is true, the output wouldn't get used anyway so it is not calculated
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            multiply_matrix_batch<size>(regs, previous_batch[vector_index], current_batch[vector_index]);
        }

        //carry_mask is all 1s if doing a shrink, else all 0s
        if (allow_shrink) {
            calculate_shrink(regs);
        } else {
            asm_immediate.assign(carry_mask_reg, 0);
        }

        //need to add the extra simd generated by matrix multiplying the tail, to the head
        //(e.g. if the tail is 16 limbs, the multiplied tail is 19 limbs with a carry at the msb)
        //since the head is 4 matricies ahead of the tail, need to also multiply the extra simd by the 4 matricies (which are now in
        // previous_batch) to get an 8-limb result
        //
        //if doing a shrink, will multiply 8 limbs to get a 12-limb result
        //if not doing a shrink, the lower 4 limbs are 0 so they don't do anything

        //need to sign-extend the extra 4 limbs in ints_buffer to 8 limbs. it is sign-extended into the carry
        int start=int_size-head_size-4;

        for (int index=0;index<size;++index) {
            //first 4 limbs in extra_ints are 0 if not shrinking, or the last 4 limbs in the tail if shrinking
            APPEND_M(str( "VMOVDQU `carry_accumulator, #", ints[0][index][start] ));
            APPEND_M(str( "VPAND `carry_accumulator, `carry_accumulator, `carry_mask" ));
            APPEND_M(str( "VMOVDQU #, `carry_accumulator", extra_ints[index][0] ));

            //next 4 limbs are the extra limbs in ints_buffer
            APPEND_M(str( "VMOVDQU `carry_accumulator, #", ints_buffer[index][start+4] ));
            APPEND_M(str( "VMOVDQU #, `carry_accumulator", extra_ints[index][4] ));

            //last 4 limbs are sign-extended from the last extra limb. the last limb is sign-extended into the carry
            //set tmp to all 1s if the last data bit in the last limb is set, else 0
            APPEND_M(str( "MOV `tmp, #", ints_buffer[index][start+7] ));
            APPEND_M(str( "SHL `tmp, #", to_hex(64-data_size) ));
            APPEND_M(str( "SAR `tmp, #", to_hex(63) ));

            //tmp2 has the carry bits cleared
            APPEND_M(str( "MOV `tmp2, `tmp" ));
            APPEND_M(str( "AND `tmp2, #", to_hex(data_mask) ));

            //sign extend but not into carry
            APPEND_M(str( "MOVQ `carry_accumulator_128, `tmp2" ));
            APPEND_M(str( "VPBROADCASTQ `carry_accumulator, `carry_accumulator_128" ));
            APPEND_M(str( "VMOVDQU #, `carry_accumulator", extra_ints[index][8] ));

            //sign extend into last carry
            APPEND_M(str( "MOV #, `tmp", extra_ints[index][11] ));
        }

        //need to save this because it is being invalidated
        APPEND_M(str( "VMOVDQU `shrink_mask_spill, `carry_mask" ));

        //can now multiply it by the merged current_batch (which is now stored in previous_batch)
        //can merge this with the shrink code to avoid loading the matricies multiple times
        for (int pass=0;pass<4;++pass) {
            expand_macros_tag c_tag_2(m, "advance_mul");

            assign_matrix_buffer<size>(regs, matrix_buffer, previous_batch[0][pass]); //invalidates carry_mask and carry_accumulator
            matrix_vector_multiply<size>(
                regs,
                matrix_buffer,
                extra_ints_multiplied,
                extra_ints,
                (pass==0)? array<simd_integer_asm, size>() : extra_ints_multiplied,
                0, extra_ints[0].size, pass
            );
        }

        APPEND_M(str( "VMOVDQU `carry_mask, `shrink_mask_spill" ));

        for (int index=0;index<size;++index) {
            simd_integer_asm c=extra_ints_multiplied[index];

            //do one carry iteration. not really sure why this is required
            {
                expand_macros_tag c_tag_2(m, "advance_carry");
                c.calculate_carry(regs, 0, true);
            }

            //need to sign extend the last limb into the carry. this seems to work even if only one carry iteration is done
            APPEND_M(str( "MOV `tmp, #", c[c.size-1] ));
            APPEND_M(str( "SHL `tmp, #", to_hex(64-data_size) ));
            APPEND_M(str( "SAR `tmp, #", to_hex(64-data_size) ));
            APPEND_M(str( "MOV #, `tmp", c[c.size-1] ));
        }

        //if shrinking, need to zero out the old values of the 4 limbs that got multiplied, since they are being replaced
        //if not shrinking, the replacement values will be 0 so they should just be added
        for (int index=0;index<size;++index) {
            APPEND_M(str( "VMOVDQU `carry_accumulator, #", ints[0][index][start] ));
            APPEND_M(str( "VPANDN `carry_accumulator, `carry_mask, `carry_accumulator" )); //zero out if shrinking
            APPEND_M(str( "VMOVDQU #, `carry_accumulator", ints[0][index][start] ));
        }

        //add the multiplication results to the tail and head
        for (int index=0;index<size;++index) {
            for (int x=0;x<extra_ints_multiplied[index].size;x+=4) {
                APPEND_M(str( "VMOVDQU `carry_accumulator, #", ints[0][index][start+x] ));
                APPEND_M(str( "VPADDQ `carry_accumulator, `carry_accumulator, #", extra_ints_multiplied[index][x] ));
                APPEND_M(str( "VMOVDQU #, `carry_accumulator", ints[0][index][start+x] ));
            }
        }

        //finally need to shift ints[0] left by one simd. this is done by changing the c_memory_base pointer
        //this pointer is only used for ints[0] and each int in ints[0] is supposed to be zero padded before the lsb. the zero padding
        // will be preserved by future operations
        //(only done if shrinking)
        do_shrink(regs);

        //have carry only in the new head
        for (int index=0;index<size;++index) {
            expand_macros_tag c_tag_2(m, "advance_carry");
            ints[0][index].calculate_carry(regs, int_size-head_size);
        }
    }

    //should do the extra work before doing this
    void generate_background_work_carry(bool assign_accumulator) {
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            int c_head_size=(vector_index==0)? head_size : 0;
            int start=(vector_index==0)? cutoff : 0;

            assert(start%4==0 && start<=int_size-c_head_size);

            if (c_head_size==int_size) {
                continue;
            }

            for (int index=int_size-c_head_size-4;index>=start;index-=4) {
                //15/18 alu ops (size==3; 15 no accumulator ; 18 with accumulator)
                c_scheduler+=[=](reg_alloc regs) {
                    expand_macros_tag c_tag(m, "bg_carry");

                    assert(matrix_buffer_state==-1);
                    assert(preserve_carry_mask_and_accumulator);

                    for (int x=0;x<size;++x) {
                        ints[vector_index][x].calculate_carry_loop(
                            regs, carry_mask_reg, carry_accumulator,
                            0, !assign_accumulator,
                            index, index+4, true
                        );
                    }
                };
            }
        }
    }

    void generate_background_work_carry_extra(bool is_extra_pass, bool assign_accumulator) {
        if (is_extra_pass) {
            assert(!assign_accumulator);
        }

        if (head_size==int_size) {
            return;
        }

        int index=int_size-head_size;

        c_scheduler+=[=](reg_alloc regs) {
            expand_macros_tag c_tag(m, "bg_carry_extra");

            if (!preserve_carry_mask_and_accumulator) {
                asm_immediate.assign(carry_mask_reg, carry_mask);
            }

            if (assign_accumulator) {
                asm_immediate.assign(carry_accumulator, 0);
            }

            matrix_buffer_state=-1;
            preserve_carry_mask_and_accumulator=true;

            for (int x=0;x<size;++x) {
                simd_integer_asm c_int=ints[0][x];

                c_int.size=c_int.size-head_size+4;

                {
                    EXPAND_MACROS_SCOPE;
                    reg_alloc regs_copy=regs;

                    regs_copy.bind_vector(m, "tmp");

                    //back up old value
                    APPEND_M(str( "VMOVDQU `tmp, #", c_int[c_int.size-4] ));
                    APPEND_M(str( "VMOVDQU #, `tmp", extra_ints[0][x] ));

                    //assign extra value from buffer
                    APPEND_M(str( "VMOVDQU `tmp, #", ints_buffer[x][index] ));
                    APPEND_M(str( "VMOVDQU #, `tmp", c_int[c_int.size-4] ));
                }

                c_int.calculate_carry_loop(
                    regs, carry_mask_reg, carry_accumulator,
                    (is_extra_pass)? index : 0, !assign_accumulator,
                    index, index+4, true
                );

                {
                    EXPAND_MACROS_SCOPE;
                    reg_alloc regs_copy=regs;

                    regs_copy.bind_vector(m, "tmp");

                    //assign new extra value to buffer
                    APPEND_M(str( "VMOVDQU `tmp, #", c_int[c_int.size-4] ));
                    APPEND_M(str( "VMOVDQU #, `tmp", ints_buffer[x][index] ));

                    //restore old value
                    APPEND_M(str( "VMOVDQU `tmp, #", extra_ints[0][x] ));
                    APPEND_M(str( "VMOVDQU #, `tmp", c_int[c_int.size-4] ));
                }
            }
        };
    }

    void generate_all_carry_work() {
        int scheduler_start=c_scheduler.all_work.size();

        assert(!preserve_carry_mask_and_accumulator); //used to decide when to assign carry_mask

        generate_background_work_carry_extra(false, false);
        generate_background_work_carry(false);

        generate_background_work_carry_extra(true, false);
        generate_background_work_carry_extra(true, false);

        generate_background_work_carry_extra(false, true); //zeros out accumulator
        generate_background_work_carry(true);

        if (scheduler_start==c_scheduler.all_work.size()) {
            c_scheduler+=[=](reg_alloc regs) {
                matrix_buffer_state=-1;
                preserve_carry_mask_and_accumulator=true;

                asm_immediate.assign(carry_mask_reg, carry_mask);
                asm_immediate.assign(carry_accumulator, 0);
            };
        }
    }

    //carry_accumulator and carry_mask_reg can alias with the matrix_buffer registers
    void generate_background_work() {
        int initial_scheduler_pos=c_scheduler.all_work.size();

        //the final ints_buffer contains an extra simd after the end of the primary vector tail
        //this is added to the head
        for (int vector_index=num_vectors-1;vector_index>=0;--vector_index) {
            int c_head_size=(vector_index==0)? head_size : 0;

            if (c_head_size==int_size) {
                continue;
            }

            int start=(vector_index==0)? cutoff : 0;

            assert(start%4==0 && start<=int_size-c_head_size);

            for (int pass=0;pass<4;++pass) {
                int last_index=(vector_index==0)? int_size-c_head_size : int_size-4;

                for (int index=last_index;index>=start;index-=4) {
                    bool is_extra=(vector_index==0 && index==last_index);

                    //18 alu ops (size==3)
                    c_scheduler+=[=](reg_alloc regs) {
                        expand_macros_tag c_tag(m, "bg_mul");

                        int t_state=pass + vector_index*4;

                        if (matrix_buffer_state!=t_state) {
                            assign_matrix_buffer<size>(regs, matrix_buffer, previous_batch[vector_index][pass]);
                            matrix_buffer_state=t_state;
                        }

                        matrix_vector_multiply<size>(
                            regs,
                            matrix_buffer,
                            (pass==3 && !is_extra)? ints[vector_index] : ints_buffer,
                            ints[vector_index],
                            (pass==0)? array<simd_integer_asm, size>() : ints_buffer,
                            index, index+4, pass, int_size-c_head_size
                        );
                    };
                }
            }
        }

        generate_all_carry_work();

        int final_scheduler_pos=c_scheduler.all_work.size();
        print( "generate_background_work:", final_scheduler_pos-initial_scheduler_pos );
    }

    array<simd_integer_asm, size> calculate_head(bool remove_invalid=true) {
        //the head is ahead of the tail by 4-8 matricies. the number of invalid entries is this plus 1
        //if it is ahead by 8 matricies, advance will be called before the user code reads the head, so the maximum number of invalid
        // entries when the head is read by the user code is 8

        array<simd_integer_asm, size> res;
        for (int x=0;x<size;++x) {
            res[x]=ints[0][x];
            res[x].start+=int_size-head_size;
            res[x].size=head_size;

            if (remove_invalid && head_size!=int_size) {
                res[x].start+=extra_head_size;
                res[x].size-=extra_head_size;
            }
        }
        return res;
    }

    //c_matrix is for every vector except the first one
    //can call this before or after multiply
    void multiply_prepare(reg_alloc regs, array<array<reg_vector, size>, num_vectors-1> c_matrix, int pass) {
        for (int x=0;x<num_vectors-1;++x) {
            assign_current_batch<size>(regs, current_batch[x+1][3-pass], c_matrix[x]);
        }
    }

    //for pass 3, should call advance before reading the head since the lsb limb is partially invalid
    void multiply(reg_alloc regs, array<reg_vector, size> c_matrix, reg_scalar b_sign, int pass) {
        EXPAND_MACROS_SCOPE;

        expand_macros_tag c_tag(m, "multiply");

        bind();

        m.bind(c_matrix, "c_matrix");

        if (b_sign.value!=-1) {
            m.bind(b_sign, "b_sign");
        }

        reg_scalar tmp=regs.bind_scalar(m, "tmp");

        array<simd_integer_asm, size> c_head=calculate_head(false);

        reg_spill accumulator_copy;
        if (preserve_carry_mask_and_accumulator) {
            accumulator_copy=regs.bind_spill(m, "carry_accumulator_copy");
            APPEND_M(str( "VMOVDQU `carry_accumulator_copy, `carry_accumulator" ));
        }

        //left multiplying, so first matrix is stored last in current_batch
        assign_current_batch<size>(regs, current_batch[0][3-pass], c_matrix);

        string skip_multiply_label=m.alloc_label();
        {
            asm_immediate.assign(carry_accumulator, ~uint64(0));
            for (int x=0;x<size;++x) {
                array<uint64, 4> identity_row={0, 0, 0, 0};
                identity_row[x]=1;

                APPEND_M(str( "VPCMPEQQ `carry_mask, `c_matrix_#, #", x, asm_immediate(identity_row) ));
                APPEND_M(str( "VPAND `carry_accumulator, `carry_accumulator, `carry_mask" ));
            }

            //if all of the bits in the accumulator are set, the input is the identity matrix
            APPEND_M(str( "VPTEST `carry_accumulator, #", asm_immediate(~uint64(0)) )); //sets CF if input is the identity matrix
            APPEND_M(str( "JC #", skip_multiply_label ));
        }

        assign_matrix_buffer<size>(regs, matrix_buffer, current_batch[0][3-pass]);
        matrix_buffer_state=-1;

        {
            expand_macros_tag c_tag_2(m, "multiply_mul");
            matrix_vector_multiply<size>(
                regs,
                matrix_buffer,
                c_head,
                c_head,
                array<simd_integer_asm, size>(),
                0, head_size, 0
            );
        }

        asm_immediate.assign(carry_mask_reg, data_sign_mask);
        for (int index=0;index<size;++index) {
            //the result is always positive (except in rare cases when doing a reduce)
            //however, there is 0 padding in the front so carrying can take quadratic time. adding data_sign_mask, doing the carry, and
            // subtracting it will reduce the probability of carrying taking quadratic time
            //e.g. there is an integer where every limb is 0 except for two that are 2^data_size followed by 2^64-1
            //if carrying lsb first, then the sign extended carry for the lsb limb is 1, which is added to 2^64-1 to get 0, so the
            // resulting integer is 0
            //if carrying msb first (which is required for simd), then the sign extended carry for 2^64-1 is also 2^64-1 which is added
            // to 0. the sign extended carry for 2^data_size is added to 2^64-1 masked by data_mask which is 2^data_size. the result of
            // a single carry iteration is to shift everything left one limb, so O(n) iterations are required to zero everything out
            //adding something to each limb then subtracting causes that pattern to cancel itself out
            //note: these patterns only happen if there are subtractions but the result is nonnegative
            for (int x=head_size-carry_adjusted_head_size;x<head_size;x+=4) {
                if (x<0) {
                    continue;
                }
                //add data_sign_mask to some of the limbs (the ones that are potentially 0)
                APPEND_M(str( "VPADDQ `carry_accumulator, `carry_mask, #", c_head[index][x] ));
                APPEND_M(str( "VMOVDQU #, `carry_accumulator", c_head[index][x] ));
            }

            //do a single carry round; this is usually sufficient to cancel out any patterns that could cause carrying to take quadratic time
            {
                expand_macros_tag c_tag_2(m, "multiply_carry");
                c_head[index].calculate_carry(regs, 0, true);
            }
        }

        asm_immediate.assign(carry_mask_reg, uint64(-int64(data_sign_mask)));
        for (int index=0;index<size;++index) {
            for (int x=head_size-carry_adjusted_head_size;x<head_size;x+=4) {
                if (x<0) {
                    continue;
                }
                //subtract out the data_sign_mask that was added
                APPEND_M(str( "VPADDQ `carry_accumulator, `carry_mask, #", c_head[index][x] ));
                APPEND_M(str( "VMOVDQU #, `carry_accumulator", c_head[index][x] ));
            }

            //do the rest of the carries (usually 1 is sufficient)
            {
                expand_macros_tag c_tag_2(m, "multiply_carry");
                c_head[index].calculate_carry(regs);
            }
        }

        string no_partial_carry_label=m.alloc_label();

        //don't need to partial carry if there is no tail. the full head size is used
        //for reduce this doesn't do anything
        has_tail(regs, tmp, "", no_partial_carry_label, false);

        for (int index=0;index<size;++index) {
            //at pass 0, the head is ahead of the tail by 4 matricies
            //this increases by 1 matrix until advance is called
            c_head[index].check_partial_carry_valid(regs, 4+pass+1);
        }

        APPEND_M(str( "#:", no_partial_carry_label ));

        //for reduce, b's new sign can be wrong which makes it become negative
        //if this happens, need to negate b, do a full carry, negate row 2 of c_matrix, and negate b_sign
        if (b_sign.value!=-1 && !asm_output_common_case_only) {
            string skip_negate_label=m.alloc_label();

            c_head[1].is_negative(regs, "", skip_negate_label);

            asm_immediate.assign(carry_accumulator, uint64(int64(-1)));
            c_head[1].fma(regs, simd_integer_asm(), c_head[1], carry_accumulator, 0);

            c_head[1].calculate_carry(regs);
            c_head[1].check_partial_carry_valid(regs, pass+1); //this is only used for reduce, where there is always a tail

            assert(num_vectors==1); //can't modify c_matrix if there are more vectors since it is used for multiply_prepare
            APPEND_M(str( "VPMULDQ `c_matrix_1, `c_matrix_1, #", asm_immediate(int64(-1)) ));
            APPEND_M(str( "NEG `b_sign" ));

            assign_current_batch<size>(regs, current_batch[0][3-pass], c_matrix);

            APPEND_M(str( "#:", skip_negate_label ));
        }

        APPEND_M(str( "#:", skip_multiply_label ));

        if (preserve_carry_mask_and_accumulator) {
            APPEND_M(str( "VMOVDQU `carry_accumulator, `carry_accumulator_copy" ));
            asm_immediate.assign(carry_mask_reg, carry_mask);
        }
    }

    void init(reg_alloc& regs, int t_head_size, int num_matrix_buffer_spills, int background_work_per_call) {
        EXPAND_MACROS_SCOPE;

        vector<reg_vector> matrix_buffer_regs;
        for (int x=0;x<matrix_buffer.size();++x) {
            if (x<num_matrix_buffer_spills) {
                matrix_buffer[x].second=regs.get_spill();
            } else {
                matrix_buffer[x].first=regs.get_vector();
                matrix_buffer_regs.push_back(matrix_buffer[x].first);
            }
        }

        carry_accumulator=matrix_buffer_regs.at(0);
        carry_mask_reg=matrix_buffer_regs.at(1);
        reg_scalar c_memory_base=regs.bind_scalar(m, "c_memory_base");

        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            for (int index=0;index<size;++index) {
                assert(ints[vector_index][index].memory_base_reg.value==asm_memory.memory_base.value);
            }
        }
        for (int index=0;index<size;++index) {
            assert(ints_buffer[index].memory_base_reg.value==asm_memory.memory_base.value);
        }

        for (int index=0;index<size;++index) {
            ints[0][index].memory_base_reg=c_memory_base;
        }

        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            for (int pass=0;pass<4;++pass) {
                for (int row=0;row<size;++row) {
                    previous_batch[vector_index][pass][row]=regs.get_spill();
                    current_batch[vector_index][pass][row]=regs.get_spill();
                }
            }
        }

        APPEND_M(str( "MOV `c_memory_base, `memory_base" ));
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            assign_identity(regs, previous_batch[vector_index]);
        }

        int_size=ints[0][0].size;
        head_size=t_head_size+extra_head_size;

        c_scheduler.num_work_per_call=background_work_per_call;
    }
};

//integer is padded on the lsb so that its total size is doubled. this also clears out the padding
simd_integer_asm remove_padding(reg_alloc regs, simd_integer_asm a_padded) {
    simd_integer_asm a=a_padded;
    a.size/=2;
    a.start+=a.size;
    assert(a.size%4==0);

    simd_integer_asm a_padding=a_padded;
    a_padding.size=a.size;
    a_padding.clear(regs);

    return a;
}

void gcd_iteration(
    reg_alloc regs, reg_spill total_num_iterations, reg_spill terminated, reg_spill c_gcd_mask,
    matrix_multiplier_asm<2, 2>& c_multiplier, int pass
) {
    EXPAND_MACROS_SCOPE;

    m.bind(total_num_iterations, "total_num_iterations");
    m.bind(terminated, "terminated");
    m.bind(c_gcd_mask, "c_gcd_mask");

    reg_scalar head_bits_end_index=regs.bind_scalar(m, "head_bits_end_index");

    array<reg_scalar, 2> head_bits={
        regs.bind_scalar(m, "head_bits_0"),
        regs.bind_scalar(m, "head_bits_1")
    };

    array<reg_vector, 2> gcd_64_res={
        regs.bind_vector(m, "gcd_64_res_0"),
        regs.bind_vector(m, "gcd_64_res_1")
    };

    reg_scalar gcd_64_num_iterations=regs.bind_scalar(m, "gcd_64_num_iterations");

    array<simd_integer_asm, 2> head=c_multiplier.calculate_head();

    extract_bits_shifted<2>(regs, head_bits, head_bits_end_index, head);

    asm_immediate.assign(gcd_64_res[0], {1, 0, 0, 0});
    asm_immediate.assign(gcd_64_res[1], {0, 1, 0, 0});

    {
        reg_alloc regs_copy=regs;
        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");
        reg_vector tmp_2=regs_copy.bind_vector(m, "tmp_2");

        string approximate_label=m.alloc_label();

        //for an exact gcd, size of a and b must be 2 limbs or less (i.e. the end index is -1, 0, or 1)
        APPEND_M(str( "CMP `head_bits_end_index, 1" ));
        APPEND_M(str( "JG #", approximate_label ));

        //if there is a tail then the gcd isn't exact
        c_multiplier.has_tail(regs_copy, head_bits_end_index, approximate_label, "");

        asm_immediate.assign(tmp_2, gcd_mask_exact);
        APPEND_M(str( "VMOVDQU `c_gcd_mask, `tmp_2" ));

        APPEND_M(str( "#:", approximate_label ));
    }
    //head_bits_end_index invalidated

    schedule_early_exit(
        regs,
        [&](reg_alloc c_regs, string early_exit_label, int iteration) {
            gcd_64_iteration(c_regs, c_gcd_mask, head_bits, gcd_64_res, early_exit_label);
        },
        gcd_num_iterations,
        gcd_num_iterations_background_work,
        c_multiplier.c_scheduler,
        gcd_64_num_iterations
    );

    //the u/v vector is multiplied by the adjusted matrix so they don't become negative
    //the a/b vector is multiplied by the result matrix

    c_multiplier.multiply(regs, gcd_64_res, reg_scalar(), pass);

    //both v0 and v1 are positive; need to negate columns in the matrix if they should be negative
    //the output values of v0/v1 may be negative. need to negate rows to make them positive
    //combined effect:
    // int delta_parity = parity*new_parity
    // new_parity=((gcd_64_num_iterations%2==0)? 1 : -1)*parity;
    // new_parity=d*parity ; delta_parity=d*parity^2=d ; d=((gcd_64_num_iterations%2==0)? 1 : -1)
    // c_matrix[0*2+0]*=delta_parity;
    // c_matrix[1*2+0]*=-delta_parity;
    // c_matrix[0*2+1]*=-delta_parity;
    // c_matrix[1*2+1]*=delta_parity;
    //this is the same as taking the absolute value of every entry in the matrix

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(total_num_iterations, "total_num_iterations");
        m.bind(terminated, "terminated");
        m.bind(gcd_64_res, "gcd_64_res");
        m.bind(gcd_64_num_iterations, "gcd_64_num_iterations");

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");
        reg_scalar tmp_b=regs_copy.bind_scalar(m, "tmp_b");

        reg_vector tmp_1=regs_copy.bind_vector(m, "tmp_1");
        reg_vector tmp_2=regs_copy.bind_vector(m, "tmp_2");

        //take the absolute value of each matrix entry to get the u/v adjusted matrix
        //need to zero out the high 32 bits of each entry since a 32-bit absolute value is being used
        asm_immediate.assign(tmp_1, 0);
        APPEND_M(str( "VPBLENDD `gcd_64_res_0, `tmp_1, `gcd_64_res_0, #", vpblendd_mask_8({1,0, 1,0, 1,0, 1,0}) ));
        APPEND_M(str( "VPBLENDD `gcd_64_res_1, `tmp_1, `gcd_64_res_1, #", vpblendd_mask_8({1,0, 1,0, 1,0, 1,0}) ));
        APPEND_M(str( "VPABSD `gcd_64_res_0, `gcd_64_res_0" ));
        APPEND_M(str( "VPABSD `gcd_64_res_1, `gcd_64_res_1" ));

        //add gcd_64_num_iterations to total_num_iterations
        APPEND_M(str( "MOV `tmp_b, `total_num_iterations" ));
        APPEND_M(str( "ADD `tmp_b, `gcd_64_num_iterations" ));
        APPEND_M(str( "MOV `total_num_iterations, `tmp_b" ));

        if (pass==3) {
            string skip_assign_terminated_label=m.alloc_label();

            APPEND_M(str( "CMP `gcd_64_num_iterations, 0" ));
            APPEND_M(str( "JNE #", skip_assign_terminated_label ));

            asm_immediate.assign(tmp, 1);
            APPEND_M(str( "MOV `terminated, `tmp" ));

            APPEND_M(str( "#:", skip_assign_terminated_label ));
        }
    }

    c_multiplier.multiply_prepare(regs, {gcd_64_res}, pass);
}

//large buffers: same size as a/b/v0
//small buffers: 12 limbs
void gcd(
    reg_alloc& regs, reg_spill v0_sign,
    simd_integer_asm a_padded, simd_integer_asm b_padded, simd_integer_asm v0,
    array<simd_integer_asm, 3> large_buffers, array<simd_integer_asm, 4> small_buffers,
    bool v0_is_1, bool gcd_is_1
) {
    EXPAND_MACROS_SCOPE;

    reg_spill total_num_iterations=regs.bind_spill(m, "total_num_iterations");
    reg_spill terminated=regs.bind_spill(m, "terminated");
    reg_spill loop_count=regs.bind_spill(m, "loop_count");
    reg_spill c_gcd_mask=regs.bind_spill(m, "c_gcd_mask");

    simd_integer_asm a=remove_padding(regs, a_padded);
    simd_integer_asm b=remove_padding(regs, b_padded);

    assert(a.size==b.size);

    assert(large_buffers[0].size==a.size);
    assert(large_buffers[1].size==a.size);
    assert(large_buffers[2].size==a.size);
    assert(small_buffers[0].size==12);
    assert(small_buffers[1].size==12);
    assert(small_buffers[2].size==12);
    assert(small_buffers[3].size==12);

    simd_integer_asm v1=large_buffers[0];

    assert(v0.size==a.size && v1.size==a.size);

    v0.clear(regs);
    v1.clear(regs);

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(total_num_iterations, "total_num_iterations");
        m.bind(terminated, "terminated");
        m.bind(loop_count, "loop_count");
        m.bind(c_gcd_mask, "c_gcd_mask");

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");
        reg_vector tmp_2=regs_copy.bind_vector(m, "tmp_2");

        asm_immediate.assign(tmp, 0);
        APPEND_M(str( "MOV `total_num_iterations, `tmp" ));
        APPEND_M(str( "MOV `terminated, `tmp" ));
        APPEND_M(str( "MOV `loop_count, `tmp" ));

        asm_immediate.assign(tmp_2, gcd_mask_approximate);
        APPEND_M(str( "VMOVDQU `c_gcd_mask, `tmp_2" ));

        asm_immediate.assign(tmp, 1);
        simd_integer_asm v=(v0_is_1)? v0 : v1;
        APPEND_M(str( "MOV #, `tmp", v[0] ));
    }

    matrix_multiplier_asm<2, 2> c_multiplier;
    c_multiplier.ints[0]={a, b};
    c_multiplier.ints[1]={v0, v1};
    c_multiplier.ints_buffer={large_buffers[1], large_buffers[2]};
    c_multiplier.extra_ints={small_buffers[0], small_buffers[1]};
    c_multiplier.extra_ints_multiplied={small_buffers[2], small_buffers[3]};

    c_multiplier.init(regs, gcd_head_size, gcd_num_matrix_buffer_spills, gcd_background_work_per_call);
    c_multiplier.do_initial_shrink(regs);

    print( "gcd main loop:" );

    {
        string loop_label=m.alloc_label();
        APPEND_M(str( "#:", loop_label ));

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
            gcd_iteration(regs, total_num_iterations, terminated, c_gcd_mask, c_multiplier, x);
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
        }
    }

    print( "gcd finish:" );

    c_multiplier.generate_background_work();
    for (int pass=0;pass<4;++pass) {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        array<reg_vector, 2> c_matrix={
            regs_copy.bind_vector(m, "c_matrix_0"),
            regs_copy.bind_vector(m, "c_matrix_1")
        };
        asm_immediate.assign(c_matrix[0], {1, 0, 0, 0});
        asm_immediate.assign(c_matrix[1], {0, 1, 0, 0});

        c_multiplier.multiply(regs_copy, c_matrix, reg_scalar(), pass);
        c_multiplier.multiply_prepare(regs_copy, {c_matrix}, pass);
    }
    c_multiplier.advance(regs, false);

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        reg_scalar tmp=regs_copy.bind_scalar(m, "tmp");
        reg_vector tmp_2=regs_copy.bind_vector(m, "tmp_2");

        //b has to be 0; if so then there probably won't be a tail
        //also, for a to be in a static location, there can't be a tail
        c_multiplier.has_tail(regs_copy, tmp, m.alloc_error_label(), "");

        array<simd_integer_asm, 2> head=c_multiplier.calculate_head();

        //make sure b is 0 by checking its head
        assert(head[1].size%4==0);
        for (int x=0;x<head[1].size;x+=4) {
            APPEND_M(str( "VMOVDQU `tmp_2, #", head[1][x] ));
            APPEND_M(str( "VPTEST `tmp_2, #", asm_immediate(~uint64(0)) ));
            APPEND_M(str( "JNZ #", m.alloc_error_label() ));
        }

        if (gcd_is_1) {
            APPEND_M(str( "VMOVDQU `tmp_2, #", head[0][0] ));
            APPEND_M(str( "VPCMPEQQ `tmp_2, `tmp_2, #", asm_immediate({1, 0, 0, 0}) )); //all 1s if test passes
            APPEND_M(str( "VPTEST `tmp_2, #", asm_immediate(~uint64(0)) )); //sets CF to 1 if tmp_2 is all 1s
            APPEND_M(str( "JNC #", m.alloc_error_label() )); //gcd is not 1 if CF is 0

            for (int x=4;x<head[0].size;x+=4) {
                APPEND_M(str( "VMOVDQU `tmp_2, #", head[0][x] ));
                APPEND_M(str( "VPTEST `tmp_2, #", asm_immediate(~uint64(0)) ));
                APPEND_M(str( "JNZ #", m.alloc_error_label() ));
            }
        }
    }

    {
        EXPAND_MACROS_SCOPE;
        reg_alloc regs_copy=regs;

        m.bind(total_num_iterations, "total_num_iterations");
        m.bind(v0_sign, "v0_sign");

        reg_scalar parity=regs_copy.bind_scalar(m, "parity");
        reg_vector tmp_2=regs_copy.bind_vector(m, "tmp_2");

        // parity=((total_num_iterations%2==0)? 1 : -1) ; 1 if even, -1 if odd
        APPEND_M(str( "MOV `parity, `total_num_iterations" ));
        APPEND_M(str( "ANDN `parity, `parity, #", asm_immediate(1) )); //1 if even, 0 if odd
        APPEND_M(str( "SHL `parity, 1" )); //2 if even, 0 if odd
        APPEND_M(str( "DEC `parity" )); //1 if even, -1 if odd

        if (!v0_is_1) {
            APPEND_M(str( "NEG `parity" ));
        }

        APPEND_M(str( "MOV `v0_sign, `parity" ));

        /*APPEND_M(str( "VMOVQ `tmp_2_128, `parity" ));
        APPEND_M(str( "VPBROADCASTQ `tmp_2, `tmp_2_128" ));

        //not calculating v1 because it's not used
        v0.fma(regs_copy, simd_integer_asm(), v0, tmp_2, 0);
        v0.calculate_carry(regs_copy);*/
    }
}


}