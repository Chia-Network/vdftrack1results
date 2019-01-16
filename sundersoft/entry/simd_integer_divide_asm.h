namespace simd_integer_namespace {


void normalize_divisor(
    reg_alloc regs, simd_integer_asm res, simd_integer_asm b, reg_scalar divisor_size, reg_scalar shift_amount
) {
    EXPAND_MACROS_SCOPE;

    assert(b.memory_base_reg.value==asm_memory.memory_base.value);

    regs.bind_scalar(divisor_size, "divisor_size");
    regs.bind_scalar(shift_amount, "shift_amount");

    reg_vector max_index_vector=regs.bind_vector(m, "max_index_vector");
    reg_scalar tmp=regs.bind_scalar(m, "tmp");
    reg_scalar tmp2=regs.bind_scalar(m, "tmp2");

    //calculates max index + 1 which is also the size
    extract_bits_shifted_scan<1>(regs, max_index_vector, {b});

    APPEND_M(str( "VMOVQ `divisor_size, `max_index_vector_128" ));

    APPEND_M(str( "MOV `shift_amount, `divisor_size" ));
    APPEND_M(str( "DEC `shift_amount" ));
    APPEND_M(str( "SHL `shift_amount, 3" ));
    APPEND_M(str( "ADD `shift_amount, `memory_base" ));
    APPEND_M(str( "MOV `shift_amount, [`shift_amount]" )); //msb value

    APPEND_M(str( "BSR `shift_amount, `shift_amount" )); //bit index of highest set bit
    APPEND_M(str( "NEG `shift_amount" ));
    APPEND_M(str( "ADD `shift_amount, #", to_hex(data_size-1) )); // data_size-1-[index of highest bit] ; amount to shift left

    APPEND_M(str( "VMOVQ `max_index_vector_128, `shift_amount" )); //zeroes out upper 128 bits

    res.logical_shift_left(regs, b, max_index_vector);
}

//result put in low input
void calculate_reciprocal(reg_alloc regs, reg_scalar high, reg_scalar low) {
    EXPAND_MACROS_SCOPE;

    m.bind(high, "high");
    m.bind(low, "low");

    reg_vector tmp=regs.bind_vector(m, "tmp");
    reg_vector tmp2=regs.bind_vector(m, "tmp2");

    //uint64 both_source=low | (high<<data_size);
    APPEND_M(str( "SHL `high, #", to_hex(data_size) ));
    APPEND_M(str( "OR `high, `low" ));

    //uint64 both=both_source;
    //both>>=2*data_size-53;
    APPEND_M(str( "MOV `low, `high" ));
    APPEND_M(str( "SHR `low, #", to_hex(2*data_size-53) ));

    //both&=~(1ull<<52);
    APPEND_M(str( "AND `low, #", to_hex(~(1ull<<52)) ));

    //assert(both>1)
    APPEND_M(str( "CMP `low, 0x1" ));
    APPEND_M(str( "JBE #", m.alloc_error_label() ));

    //--both;
    APPEND_M(str( "DEC `low" ));

    //double_bits bits; //bits>1 and bits<2
    //bits.fraction=both;
    //bits.set_exponent(0);
    APPEND_M(str( "OR `low, #", to_hex(1023ull<<52) ));
    APPEND_M(str( "VMOVQ `tmp_128, `low" ));

    //bits=double_bits(1.0/bits.to_double());
    double one_double=1.0;
    asm_immediate.assign(tmp2, *(uint64*)&one_double);
    APPEND_M(str( "VDIVSD `tmp_128, `tmp2_128, `tmp_128" ));
    APPEND_M(str( "VMOVQ `low, `tmp_128" ));

    //res=bits.fraction;
    APPEND_M(str( "AND `low, #", to_hex(bit_sequence(0, 52)) ));

    //++res;
    APPEND_M(str( "INC `low" ));

    //res|=1ull<<52;
    APPEND_M(str( "OR `low, #", to_hex(1ull<<52) ));

    //res<<=(62-52);
    APPEND_M(str( "SHL `low, #", to_hex(62-52) ));
}

//result put in low input
void calculate_quotient(reg_alloc regs, reg_scalar high, reg_scalar low, reg_scalar reciprocal) {
    EXPAND_MACROS_SCOPE;

    m.bind(high, "high");
    m.bind(low, "low");
    m.bind(reciprocal, "reciprocal");

    regs.get_scalar(reg_rax);
    regs.get_scalar(reg_rdx);

    //uint64 both=low | (high<<data_size);
    APPEND_M(str( "SHL `high, #", to_hex(data_size) ));
    APPEND_M(str( "OR `high, `low" ));

    //uint64 product_high=(uint128(both)*uint128(reciprocal))>>64;
    APPEND_M(str( "MOV RAX, `high" ));
    APPEND_M(str( "MUL `reciprocal" )); //RDX:RAX <- RAX * r/m64

    //++product_high;
    APPEND_M(str( "INC RDX" ));

    //uint64 res=product_high>>(data_size-2);
    APPEND_M(str( "SHR RDX, #", to_hex(data_size-2) ));

    //assert(res<1ull<<data_size)
    APPEND_M(str( "CMP RDX, #", to_hex(1ull<<data_size) ));
    APPEND_M(str( "JAE #", m.alloc_error_label() ));

    APPEND_M(str( "MOV `low, RDX" ));
}

//top 3 limbs of a should be 0
//b is normalized
void divide_integers_impl(
    reg_alloc regs, simd_integer_asm a_main, simd_integer_asm b, simd_integer_asm q, simd_integer_asm r,
    reg_scalar divisor_size, reg_scalar shift_amount
) {
    EXPAND_MACROS_SCOPE;

    m.bind(divisor_size, "divisor_size");
    m.bind(shift_amount, "shift_amount");

    reg_scalar reciprocal=regs.bind_scalar(m, "reciprocal");
    reg_scalar high=regs.bind_scalar(m, "high");
    reg_scalar low=regs.bind_scalar(m, "low");
    reg_vector tmp=regs.bind_vector(m, "tmp");

    assert(r.size==a.size);
    assert(q.size==a.size);
    q.clear(regs);

    simd_integer_asm a=r;
    a.logical_shift_left(regs, a_main, shift_amount);

    string jump_table=m.alloc_label();
    string jump_table_exit=m.alloc_label();

    APPEND_M(str( ".section .rodata" ));
    APPEND_M(str( ".balign 8" ));
    APPEND_M(str( "#:", jump_table ));
    vector<string> divisor_size_labels;
    for (int divisor_size_static=0;divisor_size_static<=b.size;++divisor_size_static) {
        divisor_size_labels.push_back(m.alloc_label());
        APPEND_M(str( ".quad #", divisor_size_labels.back() ));
    }
    APPEND_M(str( ".section .text" ));

    APPEND_M(str( "JMP [#+`divisor_size*8]", jump_table ));

    simd_integer_asm a_source=a;

    for (int divisor_size_static=0;divisor_size_static<=b.size;++divisor_size_static) {
        APPEND_M(str( "#:", divisor_size_labels.at(divisor_size_static) ));

        if (divisor_size_static==0) {
            APPEND_M(str( "JMP #", m.alloc_error_label() ));
            continue;
        }

        a=a_source;
        int a_size=a.size-3;
        assert(a_size>=0);
        assert(a.size>=a_size+1);

        int n=divisor_size_static;
        int c_m=(a_size+1)-divisor_size_static;

        APPEND_M(str( "MOV `high, #", b[divisor_size_static-1] ));

        if (divisor_size_static==1) {
            asm_immediate.assign(low, 0);
        } else {
            APPEND_M(str( "MOV `low, #", b[divisor_size_static-2] ));
        }

        calculate_reciprocal(regs, high, low);
        APPEND_M(str( "MOV `reciprocal, `low" ));

        int current_num_adds=0;
        for (int j=c_m-1;j>=0;--j) {
            assert(a_size==n+j);

            APPEND_M(str( "MOV `high, #", a[a_size] ));
            APPEND_M(str( "MOV `low, #", a[a_size-1] ));
            calculate_quotient(regs, high, low, reciprocal); //qj is in low

            if (current_num_adds+1>max_num_adds) {
                a.calculate_carry(regs);
                current_num_adds=0;
            }

            string skip_add_label=m.alloc_label();
            APPEND_M(str( "CMP `low, 0" ));
            APPEND_M(str( "JE #", skip_add_label )); //the partial carry check will fail if this isn't done

            APPEND_M(str( "MOV `high, `low" ));
            APPEND_M(str( "NEG `high" )); //-qj
            APPEND_M(str( "VMOVQ `tmp_128, `high" ));
            APPEND_M(str( "VPBROADCASTQ `tmp, `tmp_128" )); // tmp = -qj

            a.fma(regs, a, b, tmp, j);
            ++current_num_adds;

            int carry_start=a_size-4;
            if (carry_start<=0) {
                carry_start=0;
            }

            a.calculate_carry(regs, carry_start);
            if (carry_start!=0) {
                a.check_partial_carry_valid(regs, carry_start+1);
            }

            APPEND_M(str( "#:", skip_add_label ));

            --a_size;
            int new_a_size=round_up(a_size+2, 4);
            assert(new_a_size<=a.size && new_a_size>=4);
            a.size=new_a_size;

            APPEND_M(str( "MOV #, `low", q[j] ));
        }

        APPEND_M(str( "JMP #", jump_table_exit ));
    }

    APPEND_M(str( "#:", jump_table_exit ));

    a=a_source;

    a.calculate_carry();
    a.logical_shift_right(a, amount);
}

//a will be padded by one simd; will set up sizes of everything else automatically but the memory buffer should be at least:
// max(a.size, b.size)+4
//values of a and b are unchanged
//a_sign should be 1 or -1; if it is not bound then it is assumed to be 1
void divide_integers(
    reg_alloc regs, simd_integer_asm a, simd_integer_asm b, simd_integer_asm q, simd_integer_asm r,
    simd_integer_asm buffer,
    reg_scalar a_sign, int expected_q_size
) {
    EXPAND_MACROS_SCOPE;

    if (a_sign.value!=-1) {
        m.bind(a_sign, "a_sign");
    }

    reg_scalar divisor_size=regs.bind_scalar(m, "divisor_size");
    reg_scalar shift_amount=regs.bind_scalar(m, "shift_amount");
    reg_vector tmp=regs.bind_vector(m, "tmp");

    a.size+=4;
    a.clear(regs, a.size-4);

    buffer.size=b.size;
    normalize_divisor(regs, buffer, b, divisor_size, shift_amount);

    q.size=a.size;
    r.size=a.size;
    divide_integers_impl(regs, a, buffer, q, r);

    r.size=b.size;

    shrink(regs, q, expected_q_size);

    if (a_sign.value!=-1) {
        string positive_label=m.alloc_label();
        APPEND_M(str( "CMP `a_sign, 1" ));
        APPEND_M(str( "JE #", positive_label ));

        asm_immediate.assign(divisor_size, -1);
        q.multiply_one(regs, divisor_size, true); // q'=-q-1

        asm_immediate.assign(tmp, uint64(int64(-1)));
        r.fma(regs, b, r, tmp, 0); //r'=b-r (inputs and result are nonnegative)
        r.calculate_carry(regs);

        APPEND_M(str( "#:", positive_label ));
    }
}


}