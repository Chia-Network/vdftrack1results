namespace simd_integer_namespace {


vector3 mul_matrix_vector(matrix3 c_matrix, vector3 v) {
    vector3 res;
    for (int y=0;y<3;++y) {
        res[y]=0;
        for (int x=0;x<3;++x) {
            res[y]+=c_matrix[y*3+x]*v[x];
        }
    }
    return res;
}

void reduce_64(vector3 start_a, pair<matrix3, vector3>& res, int& num_iterations, int max_iterations) {
    matrix3 uv={
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    };
    vector3 a=start_a;

    num_iterations=0;

    int64 max_uv=1;

    int asm_num_iterations=0;
    matrix3 uv_asm=uv;
    vector3 a_asm=a;
    int64 max_uv_asm=0;

    while (true) {
        if (test_asm_funcs) {
            if (reduce_64_iteration_asm(a_asm, uv_asm, max_uv_asm)) {
                ++asm_num_iterations;
            }
        }

        if (a[2]>=a[0] || (a[2]<<1)==0) {
            break;
        }

        //can't do the divide if these don't hold; quotient is still correct
        //not using the quotient because the predicted sign of b has a high chance of being wrong
        bool ac_sign_wrong=!(a[0]>0 && a[2]>0);
        if (ac_sign_wrong) {
            break;
        }

        assert(a[1]+a[2]==int128(a[1])+int128(a[2]));
        assert(a[2]<<1==int128(a[2])<<1);

        int64 q;
        int64 r=divide_table(a[1]+a[2], a[2]<<1, q);
        {
            int shift_amount=63-reduce_num_quotient_bits;
            if ((q<<shift_amount)>>shift_amount!=q) {
                break;
            }
        }

        vector3 new_a;

        new_a[0]=a[2];
        new_a[1]=a[2]-r;
        new_a[2]=(((new_a[1]-a[1])*q)>>1) + a[0];

        vector3 new_a2;
        new_a2[0]=    a[2];
        new_a2[1]=2*q*a[2] -   a[1];
        new_a2[2]=q*q*a[2] - q*a[1] + a[0];

        matrix3 new_uv;

        //can use simd for this if everything is 32 bits. need to detect overflow
        //need to update the second row so that b will have a positive sign if it came out negative
        for (int x=0;x<3;++x) {
            new_uv[0*3+x]=       uv[2*3+x];
            new_uv[1*3+x]=(q<<1)*uv[2*3+x] -   uv[1*3+x];
            new_uv[2*3+x]= (q*q)*uv[2*3+x] - q*uv[1*3+x] + uv[0*3+x];
        }

        assert(mul_matrix_vector(new_uv, start_a)==new_a);

        //need to make this better
        //it is a lot better with simd and 32 bits. this also has overflow deteection since multiplier output is 64 bits and
        // adds are 64 bits
        //simd version will be good enough
        //the signs of everything in the matrix are determined by the iteration count since the quotient is nonnegative
        for (int x=0;x<9;++x) {
            max_uv=max(max_uv, abs_int(new_uv[x]));
        }

        //2(|q|+1)E < a'-|b'|
        // E = 3*max_uv
        //also max_uv has to not be too big
        bool valid=( 6*(abs_int(q)+1)*max_uv < new_a[0]-abs_int(new_a[1]) ) && max_uv<=data_mask;

        if (!valid) {
            break;
        }

        assert(new_a[1]==int128(a[2])-int128(r));
        assert(new_a[2]==(((int128(new_a[1])-int128(a[1]))*int128(q))>>1) + int128(a[0]));

        assert(new_a2[1]==2*q*int128(a[2]) -   int128(a[1]));
        assert(new_a2[2]==q*q*int128(a[2]) - q*int128(a[1]) + int128(a[0]));

        assert(new_a==new_a2);

        uv=new_uv;
        a=new_a;
        ++num_iterations;

        if (test_asm_funcs) {
            assert(uv==uv_asm);
            assert(a==a_asm);
            assert(num_iterations==asm_num_iterations);
            assert(max_uv_asm==max_uv);
        }

        if (debug_reduce) {
            cerr << to_hex(q) << ", ";
        }

        if (num_iterations>=max_iterations) {
            break;
        }
    }

    if (test_asm_funcs) {
        assert(uv==uv_asm);
        assert(a==a_asm); //reduce code does not update a if uv is not updated
        assert(num_iterations==asm_num_iterations);
        //assert(max_uv_asm==max_uv); will get assigned at the end so this is different
    }

    assert(mul_matrix_vector(uv, start_a)==a);

    res=make_pair(uv, a);
}

void reduce_iteration(
    bool& terminated, int64& b_sign, matrix_multiplier_both<3, 1>& c_multiplier, int pass
) {
    array<simd_integer, 3> head=c_multiplier.calculate_head();

    int head_bits_end_index;
    vector3 start_a=extract_bits_shifted<3>(get_ptr(head), head_bits_end_index);

    start_a[1]*=b_sign;

    pair<matrix3, vector3> reduce_64_res;
    int reduce_64_num_iterations;
    reduce_64(start_a, reduce_64_res, reduce_64_num_iterations, reduce_num_iterations);

    int new_b_sign=sign_int(reduce_64_res.second[1]);

    matrix3 m=reduce_64_res.first;

    //need to negate second column if b_sign is negative
    //need to negate second row if b_sign_new is negative
    //combined effect:
    m[0*3+1]*=b_sign;
    m[2*3+1]*=b_sign;
    m[1*3+0]*=new_b_sign;
    m[1*3+1]*=b_sign*new_b_sign;
    m[1*3+2]*=new_b_sign;

    b_sign=new_b_sign;

    if (reduce_64_num_iterations==0) {
        terminated=true;
    }

    c_multiplier.multiply({m}, &b_sign, pass);
}

bool reduce_impl(simd_integer& a, simd_integer& b, simd_integer& c, int64& b_sign, int expected_size) {
    assert(a.current_size()==b.current_size() && a.current_size()==c.current_size());

    bool terminated=false;

    matrix_multiplier_both<3, 1> c_multiplier;
    c_multiplier.assign_int(0, 0, a);
    c_multiplier.assign_int(0, 1, b);
    c_multiplier.assign_int(0, 2, c);

    c_multiplier.init(reduce_head_size);

    while (!terminated) {
        for (int x=0;x<4;++x) {
            reduce_iteration(terminated, b_sign, c_multiplier, x);
        }
        c_multiplier.advance(true);
    }

    for (int pass=0;pass<4;++pass) {
        c_multiplier.multiply(get_identity_matricies<3, 1>(), nullptr, pass);
    }
    c_multiplier.advance(false);

    bool converged;
    {
        array<simd_integer, 3> head=c_multiplier.calculate_head();

        int head_bits_end_index;
        vector3 start_a=extract_bits_shifted<3>(get_ptr(head), head_bits_end_index);

        converged=(start_a[2]>start_a[0]);
    }

    simd_integer a_shifted=c_multiplier.get_int(0, 0);
    simd_integer b_shifted=c_multiplier.get_int(0, 1);
    simd_integer c_shifted=c_multiplier.get_int(0, 2);

    for (int x=0;x<a.current_size();++x) {
        int pos=x+c_multiplier.c_fast.shift_amount;

        if (pos>=a.current_size()) {
            a.memory.at(x)=0;
            b.memory.at(x)=0;
            c.memory.at(x)=0;
        } else {
            a.memory.at(x)=a_shifted.memory.at(pos);
            b.memory.at(x)=b_shifted.memory.at(pos);
            c.memory.at(x)=c_shifted.memory.at(pos);
        }
    }

    if (converged) {
        a.shrink(expected_size);
        b.shrink(expected_size);
        c.shrink(expected_size);
    }

    return converged;
}

void reduce_slow_iteration(simd_integer& a, simd_integer& b, simd_integer& c, int64& b_sign) {
    simd_integer c_b;
    c_b.memory.resize(a.memory.size());
    c_b.fma(&c, b, {uint64(b_sign)}, 0, false, false);
    c_b.calculate_carry_lsb();

    int64 c_b_sign;
    to_sign_magnitude(c_b, c_b_sign);

    simd_integer c_2=c;
    c_2.logical_shift_left(c_2, 1);

    simd_integer s_padded;
    simd_integer r_unused;
    divide_integers(c_b, c_2, s_padded, r_unused, c_b_sign, 1);

    int64 s=sign_extend_data(s_padded.memory.at(0));

    //A=c
    simd_integer new_a=c;

    //B = ((s*c)<<1) - b
    simd_integer new_b=b;
    new_b.multiply_one(-b_sign);
    new_b.fma(&new_b, c, {uint64(2*s)}, 0, false, false);
    new_b.calculate_carry_lsb();

    //C = c*s*s - b*s+a;
    int64 s_s=s*s;
    uint64 s_s_low=s_s & data_mask;
    uint64 s_s_high=uint64(s_s >> data_size);

    simd_integer new_c=a;
    new_c.fma(&new_c, b, {uint64(int64(-s*b_sign))}, 0, false, false);
    new_c.fma(&new_c, c, {s_s_low}, 0, false, false);
    new_c.fma(&new_c, c, {s_s_high}, 1, false, false);
    new_c.calculate_carry();

    {
        integer a_int=to_integer(a);
        integer b_int=to_integer(b);
        integer c_int=to_integer(c);

        if (b_sign==-1) {
            b_int=-b_int;
        }

        integer s_int = (c_int+b_int)/(c_int<<1);

        integer new_a_int = c_int;
        integer new_b_int = ((s_int*c_int)<<1) - b_int;
        integer new_c_int = c_int*s_int*s_int - b_int*s_int+a_int;

        integer expected_new_a_int=to_integer(new_a);
        integer expected_new_b_int=to_integer(new_b);
        integer expected_new_c_int=to_integer(new_c);
        integer expected_s_int=to_integer(s_padded);

        assert(new_a_int==expected_new_a_int);
        assert(new_b_int==expected_new_b_int);
        assert(new_c_int==expected_new_c_int);
        assert(s_int==expected_s_int);
    }

    if (debug_reduce) {
        cerr << "SLOW " << to_hex(s) << ", ";
    }

    a=new_a;
    b=new_b;
    c=new_c;

    to_sign_magnitude(b, b_sign);
}

void reduce(simd_integer& a, simd_integer& b, simd_integer& c, int64& b_sign, int expected_size) {
    if (debug_reduce) {
        print( "=== reduce simd ===" );
    }

    while (true) {
        bool res=reduce_impl(a, b, c, b_sign, expected_size);
        if (res) {
            break;
        } else {
            cerr << "S";
            reduce_slow_iteration(a, b, c, b_sign);
        }
    }

    if (debug_reduce) {
        print( "" );
        print( "=== end reduce simd ===" );
    }
}


}