namespace simd_integer_namespace {


int divide_table_stats_calls=0;
int divide_table_stats_table=0;

generic_stats gcd_64_num_iterations;

//used for both gcd and reduce
int64 divide_table_lookup(int64 index) {
    assert(index>=0 && index<=bit_sequence(0, divide_table_index_bits));

    uint128 res = (~uint128(0)) / uint128(max(uint64(index), uint64(1)));
    res>>=64;

    return res;
}

int64 divide_table_64(int64 a, int64 b, int64& q) {
    assert(b>0);

    q=a/b;
    int64 r=a%b;

    if (r<0) {
        r+=b;
        --q;
    }

    assert(r>=0 && r<b && q*b+r==a);

    return r;
}

//note: this floors the quotient instead of truncating it like the div instruction
int64 divide_table(int64 a, int64 b, int64& q) {
    ++divide_table_stats_calls;
    
    assert(b>0);

    //b_shift=(64-divide_table_index_bits) - lzcnt(b)
    //bsr(b)=63-lzcnt(b)
    //63-bsr(b)=lzcnt(b)
    //b_shift=(64-divide_table_index_bits) - 63-bsr(b)
    //b_shift=64-divide_table_index_bits - 63 + bsr(b)
    //b_shift=1-divide_table_index_bits + bsr(b)
    //b_shift=bsr(b) - (divide_table_index_bits-1)

    int b_shift = (64-divide_table_index_bits) - __builtin_clzll(b);
    if (b_shift<0) { //almost never happens
        b_shift=0;
    }

    int64 b_approx = b >> b_shift;
    int64 b_approx_inverse = divide_table_lookup(b_approx);

    q = (int128(a)*int128(b_approx_inverse)) >> 64; //high part of product
    q >>= b_shift;

    int128 qb_128=int128(q)*int128(b);
    int64 qb_64=int64(qb_128);

    int128 r_128=int128(a)-int128(qb_64);
    int64 r_64=int64(r_128);

    //int128 r=int128(a)-int128(q)*int128(b);
    //if (uint128(r)>=b) {

    bool invalid_1=(int128(qb_64)!=qb_128 || int128(r_64)!=r_128 || uint64(r_64)>=b);

    int128 r_2=int128(a)-int128(q)*int128(b);

    bool invalid_2=(uint128(r_2)>=b);

    assert(invalid_1==invalid_2);
    if (!invalid_2) {
        assert(r_64==int64(r_2));
    }

    int64 r=r_2;
    if (invalid_2) {
        r=divide_table_64(a, b, q);
    } else {
        ++divide_table_stats_table;
    }

    int64 q_expected;
    int64 r_expected=divide_table_64(a, b, q_expected);

    assert(q==q_expected);
    assert(r==r_expected);

    if (test_asm_funcs) {
        int64 q_asm;
        int64 r_asm=divide_table_asm(a, b, q_asm);

        assert(q_asm==q_expected);
        assert(r_asm==r_expected);
    }

    return r;
}

void gcd_64(vector2 start_a, pair<matrix2, vector2>& res, int& num_iterations, bool approximate, int max_iterations) {
    matrix2 uv={1, 0, 0, 1};
    vector2 a=start_a;

    num_iterations=0;

    if (approximate && (start_a[0]==start_a[1] || start_a[1]==0)) {
        res=make_pair(uv, a);
        return;
    }

    int asm_num_iterations=0;
    matrix2 uv_asm=uv;
    vector2 a_asm=a;

    while (true) {
        if (test_asm_funcs) {
            if (gcd_64_iteration_asm(a_asm, uv_asm, approximate)) {
                ++asm_num_iterations;
            }
        }

        if (a[1]==0) {
            break;
        }

        assert(a[0]>a[1] && a[1]>0);

        int64 q;
        int64 r=divide_table(a[0], a[1], q);
        {
            int shift_amount=63-gcd_num_quotient_bits;
            if ((q<<shift_amount)>>shift_amount!=q) {
                break;
            }
        }


        vector2 new_a={a[1], r};

        matrix2 new_uv;
        for (int x=0;x<2;++x) {
            new_uv[0*2+x]=uv[1*2+x];
            new_uv[1*2+x]=uv[0*2+x] - q*uv[1*2+x];
        }

        bool valid=true;

        if (approximate) {
            assert(new_uv[1*2+0]!=0);
            bool is_even=(new_uv[1*2+0]<0);

            bool valid_exact;
            if (is_even) {
                valid_exact=(new_a[1]>=-new_uv[1*2+0] && new_a[0]-new_a[1]>=new_uv[1*2+1]-new_uv[0*2+1]);
            } else {
                valid_exact=(new_a[1]>=-new_uv[1*2+1] && new_a[0]-new_a[1]>=new_uv[1*2+0]-new_uv[0*2+0]);
            }

            //valid=valid_exact;
            valid=
                (new_a[1]>=-new_uv[1*2+0] && new_a[0]-new_a[1]>=new_uv[1*2+1]-new_uv[0*2+1]) &&
                (new_a[1]>=-new_uv[1*2+1] && new_a[0]-new_a[1]>=new_uv[1*2+0]-new_uv[0*2+0])
            ;

            assert(valid==valid_exact);

            if (valid) {
                assert(valid_exact);
            }
        }

        //have to do this even if approximate is false
        for (int x=0;x<4;++x) {
            if (abs_int(new_uv[x])>data_mask) {
                valid=false;
            }
        }

        if (!valid) {
            break;
        }

        uv=new_uv;
        a=new_a;
        ++num_iterations;

        if (test_asm_funcs) {
            assert(uv==uv_asm);
            assert(a==a_asm);
            assert(num_iterations==asm_num_iterations);
        }

        if (num_iterations>=max_iterations) {
            break;
        }
    }

    gcd_64_num_iterations.add(num_iterations);

    for (int x=0;x<4;++x) {
        assert(abs_int(uv[x])<=data_mask);
    }

    if (test_asm_funcs) {
        assert(uv==uv_asm);
        //assert(a==a_asm); the asm code will update a even if it becomes invalid; fine since it's not used
        assert(num_iterations==asm_num_iterations);
    }

    res=make_pair(uv, a);
}

template<long unsigned int num> array<int64, num> extract_bits_shifted(array<simd_integer*, num> i, int& end_index) {
    uint32 mask=0;
    end_index=0;

    //this is simd code with no branches
    for (int x=i[0]->current_size()-1;x>=2;--x) {
        uint32 t_mask=0;
        for (int y=0;y<num;++y) {
            t_mask|=i[y]->memory.at(x);
        }

        int32 t_mask_nonzero=(t_mask==0)? 0 : ~uint32(0);

        int t_end_index=(x+1) & t_mask_nonzero;
        end_index=max(end_index, t_end_index);

        if (end_index==(x+1)) {
            mask=t_mask;
        }
    }

    --end_index;

    array<int64, num> res;

    int c_end_index=max(end_index, 2);

    int num_leading_bits=(mask==0)? 64 : __builtin_clzll(mask); //LZCNT instruction
    num_leading_bits-=carry_size;
    assert(num_leading_bits<=data_size);
    //num_leading_bits=min(num_leading_bits, data_size); todo //doubt this does anything

    for (int x=0;x<num;++x) {
        uint64 a=i[x]->memory.at(c_end_index);
        uint64 b=i[x]->memory.at(c_end_index-1);
        uint64 c=i[x]->memory.at(c_end_index-2);

        res[x] = a<<(data_size+num_leading_bits) | b<<(num_leading_bits) | c>>(data_size-num_leading_bits);
    }

    return res;
}

template<long unsigned int size> array<simd_integer*, size> get_ptr(array<simd_integer, size>& ints) {
    array<simd_integer*, size> res;
    for (int x=0;x<size;++x) {
        res[x]=&ints[x];
    }
    return res;
}

template<long unsigned int size> array<simd_integer*, size> get_null_ptr() {
    array<simd_integer*, size> res;
    for (int x=0;x<size;++x) {
        res[x]=nullptr;
    }
    return res;
}

template<long unsigned int size> void assign_identity(
    array<array<int64, size*size>, 4>& res
) {
    for (int z=0;z<4;++z) {
        for (int y=0;y<size;++y) {
            for (int x=0;x<size;++x) {
                res[z][y*size+x]=(x==y && z==0)? 1 : 0;
            }
        }
    }
}

template<long unsigned int size> uint64 bitwise_or_ints(
    array<simd_integer, size>& ints,
    int start
) {
    uint64 res=0;
    for (int x=0;x<size;++x) {
        res|=ints[x].memory.at(start);
    }
    return res;
}

template<long unsigned int size> void multiply_matrix_batch(
    array<int64, size*size>& res,
    array<array<int64, size*size>, 2> current_batch
) {
    array<int64, size*size> buffer;
    for (int x=0;x<size*size;++x) {
        buffer[x]=0;
    }

    //outer product of column z of current_batch[0] and row z of current_batch[1], accumulated into buffer
    for (int z=0;z<size;++z) {
        for (int y=0;y<size;++y) {
            for (int x=0;x<size;++x) {
                buffer[y*size+x]+=current_batch[1][z*size+x] * current_batch[0][y*size+z];
            }
        }
    }

    res=buffer;
}

template<long unsigned int size, long unsigned int num> array<int128, size*size> multiply_matrix_batch_slow(
    array<array<int64, size*size>, num> current_batch
) {
    array<int128, size*size> res;

    for (int y=0;y<size;++y) {
        for (int x=0;x<size;++x) {
            res[y*size+x]=(x==y)? 1 : 0;
        }
    }

    for (int pass=0;pass<num;++pass) {
        array<int64, size*size> c=current_batch[num-1-pass];

        array<int128, size*size> new_res;
        for (int x=0;x<size*size;++x) {
            new_res[x]=0;
        }

        for (int y=0;y<size;++y) {
            for (int x=0;x<size;++x) {
                for (int z=0;z<size;++z) {
                    new_res[y*size+x]+=int128(c[y*size+z]) * res[z*size+x];
                }
            }
        }

        res=new_res;
    }

    return res;
}

template<long unsigned int size> void multiply_matrix_batch(
    array<array<int64, size*size>, 4>& res,
    array<array<int64, size*size>, 4> current_batch
) {
    array<int128, size*size> expected_0=multiply_matrix_batch_slow<size, 2>({current_batch[0], current_batch[1]});
    array<int128, size*size> expected_1=multiply_matrix_batch_slow<size, 2>({current_batch[2], current_batch[3]});
    array<int128, size*size> expected=multiply_matrix_batch_slow<size>(current_batch);

    multiply_matrix_batch<size>(current_batch[0], {current_batch[0], current_batch[1]});
    multiply_matrix_batch<size>(current_batch[1], {current_batch[2], current_batch[3]});

    for (int x=0;x<size*size;++x) {
        assert(expected_0[x]==current_batch[0][x]);
        assert(expected_1[x]==current_batch[1][x]);
    }

    for (int y=0;y<size;++y) {
        array<int64, size> buffer;
        for (int x=0;x<size;++x) {
            buffer[x]=current_batch[0][y*size+x];
        }

        for (int x=0;x<size;++x) {
            int128 accumulator=0;
            for (int z=0;z<size;++z) {
                accumulator+=int128(buffer[z])*current_batch[1][z*size+x];
            }

            //final result for matrix position [y,x] is in the accumulator

            assert(expected[y*size+x]==accumulator);

            uint64 accumulator_low=accumulator;
            uint64 accumulator_high=accumulator >> 64;

            res[0][y*size+x]=accumulator_low & data_mask;
            res[1][y*size+x]=(accumulator_low >> data_size) & data_mask;
            res[2][y*size+x]=(accumulator >> (2*data_size)) & data_mask;
            res[3][y*size+x]=int64(accumulator_high) >> (3*data_size-64);
        }
    }
}

template<int size> array<int64, size*size> get_identity_matrix() {
    array<int64, size*size> res;
    for (int y=0;y<size;++y) {
        for (int x=0;x<size;++x) {
            res[y*size+x]=(x==y)? 1 : 0;
        }
    }
    return res;
}

template<int size, int num_vectors> array<array<int64, size*size>, num_vectors> get_identity_matricies() {
    array<array<int64, size*size>, num_vectors> identity;
    for (int x=0;x<num_vectors;++x) {
        identity[x]=get_identity_matrix<size>();
    }
    return identity;
}

void add_ints_buffer(simd_integer& target, simd_integer& source, int start) {
    for (int x=0;x<source.current_size();++x) {
        target.memory.at(x+start)+=source.memory.at(x);
    }
}

template<int size, int num_vectors> struct matrix_multiplier_fast {
    const int extra_head_size=8;

    bool is_primary=false;

    array<array<simd_integer, size>, num_vectors> ints;

    array<array<array<int64, size*size>, 4>, num_vectors> previous_batch;
    array<array<array<int64, size*size>, 4>, num_vectors> current_batch;

    int int_size=-1;
    int head_size=-1; //last 4 entries are potentially invalid; rest are valid

    int shift_amount=0;

    bool has_tail(bool small_head=true) {
        int c_head_size=(small_head)? head_size-extra_head_size : head_size;
        return int_size-shift_amount>c_head_size;
    }

    void initial_shrink() {
        while (true) {
            //need the 5 msb limbs to all be 0 to do a shrink. the shrink will remove 4 msb limbs
            bool do_shrink=true;
            for (int x=0;x<5;++x) {
                if (bitwise_or_ints<size>(ints[0], int_size-x-1)!=0) {
                    do_shrink=false;
                }
            }

            //can't shrink if the tail is empty
            //the extra entires in the head aren't read by the gcd/reduce code so they are treated as part of the tail
            //those entries will eventually get zeroed out
            if (!has_tail(true)) {
                do_shrink=false;
            }

            if (!do_shrink) {
                break;
            }

            if (do_shrink) {
                for (int index=0;index<size;++index) {
                    ints[0][index].left_shift_limbs(4);
                }
                shift_amount+=4;
            }
        }
    }

    //should also call carry_work_extra before calling this
    void carry_work(uint64& carry_accumulator, bool assign_accumulator) {
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            int c_head_size=(vector_index==0)? head_size : 0;

            if (c_head_size==int_size) {
                continue;
            }

            for (int index=int_size-c_head_size-4;index>=0;index-=4) {
                for (int x=0;x<size;++x) {
                    ints[vector_index][x].calculate_carry_loop(
                        carry_accumulator,
                        0, !assign_accumulator, //false, false,
                        index, index+4, true
                    );
                }
            }
        }
    }

    void carry_work_extra(
        uint64& carry_accumulator, array<simd_integer, size>& ints_buffer, bool is_extra_pass, bool assign_accumulator
    ) {
        if (is_extra_pass) {
            assert(!assign_accumulator);
        }

        if (head_size==int_size) {
            return;
        }

        int index=int_size-head_size;

        for (int x=0;x<size;++x) {
            simd_integer& c_int=ints[0][x];

            simd_integer c_int_subset=c_int.subset(0, c_int.current_size()-head_size);
            for (int y=0;y<4;++y) {
                c_int_subset.memory.push_back(ints_buffer[x].memory.at(index+y));
            }

            c_int_subset.calculate_carry_loop(
                carry_accumulator,
                (is_extra_pass)? index : 0, !assign_accumulator, //false, false,
                index, index+4, true
            );

            for (int y=0;y<4;++y) {
                uint64& to=ints_buffer[x].memory.at(index+y);
                uint64 from=c_int_subset.memory[c_int_subset.memory.size()-4+y];
                to=from;
            }
            c_int_subset.memory.resize(c_int_subset.memory.size()-4);

            c_int.assign_subset(0, c_int_subset);
        }
    }

    void advance(bool allow_shrink) {
        array<simd_integer, size> ints_buffer=ints[0];

        //background work:
        //the final ints_buffer contains an extra simd after the end of the primary vector tail
        //this is added to the head
        for (int vector_index=num_vectors-1;vector_index>=0;--vector_index) {
            int c_head_size=(vector_index==0)? head_size : 0;

            if (c_head_size==int_size) {
                continue;
            }

            for (int pass=0;pass<4;++pass) {
                int last_index=(vector_index==0)? int_size-c_head_size : int_size-4;

                for (int index=last_index;index>=0;index-=4) {
                    bool is_extra=(vector_index==0 && index==last_index);

                    matrix_vector_multiply<size>(
                        previous_batch[vector_index][pass],
                        (pass==3 && !is_extra)? get_ptr(ints[vector_index]) : get_ptr(ints_buffer),
                        get_ptr(ints[vector_index]),
                        (pass==0)? get_null_ptr<size>() : get_ptr(ints_buffer),
                        index, index+4, pass, int_size-c_head_size
                    );
                }
            }
        }

        uint64 carry_accumulator;

        carry_work_extra(carry_accumulator, ints_buffer, false, false);
        carry_work(carry_accumulator, false);

        carry_work_extra(carry_accumulator, ints_buffer, true, false);
        carry_work_extra(carry_accumulator, ints_buffer, true, false);

        carry_accumulator=0;
        carry_work_extra(carry_accumulator, ints_buffer, false, true);
        carry_work(carry_accumulator, true);

        while ((carry_accumulator & carry_mask)!=0) {
            carry_accumulator=0;
            carry_work_extra(carry_accumulator, ints_buffer, false, true);
            carry_work(carry_accumulator, true);

            if (is_primary) {
                cerr << "!";
            }
        }

        //current_batch is done being generated and previous_batch is done being used, so generate new previous_batch
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            multiply_matrix_batch<size>(previous_batch[vector_index], current_batch[vector_index]);

            if (test_asm_funcs) {
                array<array<int64, size*size>, 4> previous_batch_asm;
                multiply_matrix_batch_2_asm<size>(previous_batch_asm, current_batch[vector_index]);

                if (size==2) {
                    assert(previous_batch_asm==previous_batch[vector_index]);
                }
            }
        }

        bool do_shrink=false;

        if (allow_shrink) {
            //need the 5 msb limbs to all be 0 to do a shrink. the shrink will remove 4 msb limbs
            do_shrink=true;
            for (int x=0;x<5;++x) {
                if (bitwise_or_ints<size>(ints[0], int_size-x-1)!=0) {
                    do_shrink=false;
                }
            }

            //can't shrink if the tail is empty
            //the extra entires in the head aren't read by the gcd/reduce code so they are treated as part of the tail
            //those entries will eventually get zeroed out
            if (!has_tail(true)) {
                do_shrink=false;
            }
        }

        //need to add the extra simd generated by matrix multiplying the tail, to the head
        //(e.g. if the tail is 16 limbs, the multiplied tail is 19 limbs with a carry at the msb)
        //since the head is 4 matricies ahead of the tail, need to also multiply the extra simd by the 4 matricies (which are now in
        // previous_batch) to get an 8-limb result
        //
        //if doing a shrink, will multiply 8 limbs to get a 12-limb result

        {
            array<simd_integer, size> extra_ints;
            array<simd_integer, size> extra_ints_multiplied;

            //need to sign-extend the extra 4 limbs in ints_buffer to 8 limbs. it is sign-extended into the carry
            int start=(do_shrink)? int_size-head_size-4 : int_size-head_size;
            int num=(do_shrink)? 8 : 4;

            for (int index=0;index<size;++index) {
                simd_integer c;
                for (int x=0;x<num;++x) {
                    int pos=(start+x);

                    simd_integer& c_buffer=(pos<int_size-head_size)? ints[0][index] : ints_buffer[index];
                    uint64 v=c_buffer.memory.at(pos);

                    c.memory.push_back(v);
                }

                uint64 sign_extend_c=(c.memory.back() & data_sign_mask)? ~uint64(0) : uint64(0);

                for (int x=0;x<3;++x) {
                    c.memory.push_back(sign_extend_c & data_mask);
                }
                c.memory.push_back(sign_extend_c);

                extra_ints[index]=c;

                c.memory.clear();
                c.memory.resize(num+4, 0);
                extra_ints_multiplied[index]=c;
            }

            //can now multiply it by the merged current_batch (which is now stored in previous_batch)
            //can merge this with the shrink code to avoid loading the matricies multiple times
            for (int pass=0;pass<4;++pass) {
                matrix_vector_multiply<size>(
                    previous_batch[0][pass],
                    get_ptr(extra_ints_multiplied),
                    get_ptr(extra_ints),
                    (pass==0)? get_null_ptr<size>() : get_ptr(extra_ints_multiplied),
                    0, num+4, pass
                );
            }

            for (int index=0;index<size;++index) {
                simd_integer& c=extra_ints_multiplied[index];

                c.calculate_carry(0, true);
                c.memory.back()=sign_extend_data(c.memory.back());
            }

            //need to zero out the old values of the 4 limbs that got multiplied, since they are being replaced
            if (do_shrink) {
                for (int index=0;index<size;++index) {
                    for (int x=0;x<4;++x) {
                        ints[0][index].memory.at(start+x)=0;
                    }
                }
            }

            for (int index=0;index<size;++index) {
                add_ints_buffer(ints[0][index], extra_ints_multiplied[index], start);
            }

            //have carry only in the head
            for (int index=0;index<size;++index) {
                ints[0][index].calculate_carry(start, false, 1, false);
            }
        }

        if (do_shrink) {
            //finally need to shift ints[0] left by one simd. this is done by changing the c_memory_base pointer
            //this pointer is only used for ints[0] and each int in ints[0] is supposed to be zero padded before the lsb.
            //the zero padding will be preserved by future operations
            for (int index=0;index<size;++index) {
                ints[0][index].left_shift_limbs(4);
            }
            shift_amount+=4;
        }
    }

    array<simd_integer, size> calculate_head(bool remove_invalid=true) {
        //the head can be ahead of the tail by up to 7 matricies

        int c_head_size=(remove_invalid && head_size!=int_size)? head_size-extra_head_size : head_size;

        array<simd_integer, size> res;
        for (int x=0;x<size;++x) {
            res[x]=ints[0][x].subset(int_size-c_head_size, c_head_size);
        }
        return res;
    }

    //for pass 3, should call advance before reading the head since the lsb limb is partially invalid
    void multiply(array<array<int64, size*size>, num_vectors> c_matrix, int64* b_sign, int pass) {
        array<simd_integer, size> c_head=calculate_head(false);

        //left multiplying, so first matrix is stored last in current_batch
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            current_batch[vector_index][3-pass]=c_matrix[vector_index];
        }

        matrix_vector_multiply<size>(
            c_matrix[0],
            get_ptr(c_head),
            get_ptr(c_head),
            get_null_ptr<size>(),
            0, head_size, 0
        );

        for (int index=0;index<size;++index) {
            for (int x=head_size-12;x<head_size;++x) {
                c_head[index].memory[x]+=data_sign_mask;
            }

            c_head[index].calculate_carry(0, true);

            for (int x=head_size-12;x<head_size;++x) {
                c_head[index].memory[x]-=data_sign_mask;
            }

            c_head[index].calculate_carry(0, false, (is_primary && b_sign==nullptr)? 1 : -1);
        }

        if (has_tail(false)) {
            for (int index=0;index<size;++index) {
                //at pass 0, the head is ahead of the tail by 4 matricies
                //this increases by 1 matrix until advance is called
                c_head[index].check_partial_carry_valid(4+pass+1);
            }
        }

        //for reduce, b's new sign can be wrong which makes it become negative
        //if this happens, need to negate b, do a full carry, negate row 2 of c_matrix, and negate b_sign
        if (b_sign!=nullptr && c_head[1].is_negative()) {
            c_head[1].fma(nullptr, c_head[1], {uint64(int64(-1))}, 0, false, false);
            c_head[1].calculate_carry();
            c_head[1].check_partial_carry_valid(4+pass+1); //this is only used for reduce, where there is always a tail

            for (int x=0;x<3;++x) {
                c_matrix[0][1*size+x]=-c_matrix[0][1*size+x];
            }

            *b_sign=-*b_sign;

            current_batch[0][3-pass]=c_matrix[0];
        }

        for (int index=0;index<size;++index) {
            ints[0][index].assign_subset(int_size-c_head[index].current_size(), c_head[index]);
        }
    }

    void init(int t_head_size) {
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            assign_identity<size>(previous_batch[vector_index]);
            //current_batch[vector_index]=get_identity_matricies<size, 4>();
        }

        int_size=ints[0][0].current_size();
        head_size=t_head_size+extra_head_size;

        initial_shrink();
    }
};

template<int size, int num_vectors> struct matrix_multiplier_basic {
    array<array<simd_integer, size>, num_vectors> ints;

    int int_size=-1;
    int head_size=-1;

    int shift_amount=0;

    bool has_tail() {
        return int_size-shift_amount>head_size;
    }

    void initial_shrink() {
        while (true) {
            //need the 5 msb limbs to all be 0 to do a shrink. the shrink will remove 4 msb limbs
            bool do_shrink=true;
            for (int x=0;x<5;++x) {
                if (bitwise_or_ints<size>(ints[0], int_size-x-1)!=0) {
                    do_shrink=false;
                }
            }

            //can't shrink if the tail is empty
            //the extra entires in the head aren't read by the gcd/reduce code so they are treated as part of the tail
            //those entries will eventually get zeroed out
            if (!has_tail()) {
                do_shrink=false;
            }

            if (!do_shrink) {
                break;
            }

            if (do_shrink) {
                for (int index=0;index<size;++index) {
                    ints[0][index].left_shift_limbs(4);
                }
                shift_amount+=4;
            }
        }
    }

    void advance(bool allow_shrink) {
        if (!allow_shrink) {
            return;
        }

        //need the 5 msb limbs to all be 0 to do a shrink. the shrink will remove 4 msb limbs
        bool do_shrink=true;
        for (int x=0;x<5;++x) {
            if (bitwise_or_ints<size>(ints[0], int_size-x-1)!=0) {
                do_shrink=false;
            }
        }

        //can't shrink if the tail is empty
        //4 of the entries in the head aren't read by the gcd/reduce code so they are treated as part of the tail
        //those entries will eventually get zeroed out
        if (!has_tail()) {
            do_shrink=false;
        }

        if (!do_shrink) {
            return;
        }

        //finally need to shift ints[0] left by one simd. this is done by changing the c_memory_base pointer
        //this pointer is only used for ints[0] and each int in ints[0] is supposed to be zero padded before the lsb.
        //the zero padding will be preserved by future operations
        for (int index=0;index<size;++index) {
            ints[0][index].left_shift_limbs(4);
        }
        shift_amount+=4;
    }

    array<simd_integer, size> calculate_head() {
        array<simd_integer, size> res;
        for (int x=0;x<size;++x) {
            res[x]=ints[0][x].subset(int_size-head_size, head_size);
        }
        return res;
    }

    void multiply(array<array<int64, size*size>, num_vectors> c_matrix, int64* b_sign, int pass) {
        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            matrix_vector_multiply<size>(
                c_matrix[vector_index],
                get_ptr(ints[vector_index]),
                get_ptr(ints[vector_index]),
                get_null_ptr<size>(),
                0, int_size, 0
            );
        }

        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            for (int index=0;index<size;++index) {
                ints[vector_index][index].calculate_carry();
            }
        }

        if (b_sign!=nullptr && ints[0][1].is_negative()) {
            ints[0][1].fma(nullptr, ints[0][1], {uint64(int64(-1))}, 0, false, false);
            ints[0][1].calculate_carry();

            *b_sign=-*b_sign;
        }
    }

    void init(int t_head_size) {
        int_size=ints[0][0].current_size();
        head_size=t_head_size;

        initial_shrink();
    }
};

template<int size, int num_vectors> struct matrix_multiplier_both {
    matrix_multiplier_basic<size, num_vectors> c_basic;
    matrix_multiplier_fast<size, num_vectors> c_fast;

    int next_pass=0;

    void check_consistency() {
        matrix_multiplier_fast<size, num_vectors> c_fast_copy=c_fast;
        c_fast_copy.is_primary=false;

        if (next_pass!=0) {
            for (int pass=next_pass;pass<4;++pass) {
                c_fast_copy.multiply(get_identity_matricies<size, num_vectors>(), nullptr, pass);
            }
            c_fast_copy.advance(false);
        }

        for (int pass=0;pass<4;++pass) {
            c_fast_copy.multiply(get_identity_matricies<size, num_vectors>(), nullptr, pass);
        }
        c_fast_copy.advance(false);

        assert(c_fast_copy.shift_amount==c_basic.shift_amount);

        for (int vector_index=0;vector_index<num_vectors;++vector_index) {
            for (int index=0;index<size;++index) {
                assert(
                    c_basic.ints[vector_index][index].memory==
                    c_fast_copy.ints[vector_index][index].memory
                );
            }
        }
    }

    void assign_int(int vector_index, int index, simd_integer& value) {
        c_basic.ints[vector_index][index]=value;
        c_fast.ints[vector_index][index]=value;
    }

    simd_integer& get_int(int vector_index, int index) {
        assert(
            c_basic.ints[vector_index][index].memory==
            c_fast.ints[vector_index][index].memory
        );

        check_consistency();

        return c_fast.ints[vector_index][index];
    }

    bool has_tail() {
        bool res_basic=c_basic.has_tail();
        bool res_fast=c_fast.has_tail();
        assert(res_basic==res_fast);

        check_consistency();

        return res_fast;
    }

    void advance(bool allow_shrink) {
        //todo allow_shrink=false;

        c_basic.advance(allow_shrink);
        c_fast.advance(allow_shrink);
        next_pass=0;

        check_consistency();
    }

    array<simd_integer, size> calculate_head() {
        array<simd_integer, size> res_basic=c_basic.calculate_head();
        array<simd_integer, size> res_fast=c_fast.calculate_head();

        for (int index=0;index<size;++index) {
            assert(res_basic[index].memory==res_fast[index].memory);
        }

        check_consistency();

        return res_fast;
    }

    void multiply(array<array<int64, size*size>, num_vectors> c_matrix, int64* b_sign, int pass) {
        assert(pass==next_pass && pass>=0 && pass<4);
        ++next_pass;

        int64 b_sign_basic=(b_sign==nullptr)? 0 : *b_sign;
        int64 b_sign_fast=b_sign_basic;

        c_basic.multiply(c_matrix, (b_sign==nullptr)? nullptr : &b_sign_basic, pass);
        c_fast.multiply(c_matrix, (b_sign==nullptr)? nullptr : &b_sign_fast, pass);

        if (b_sign!=nullptr) {
            assert(b_sign_basic==b_sign_fast);
            *b_sign=b_sign_fast;
        }

        check_consistency();
    }

    void init(int t_head_size) {
        c_fast.is_primary=true;

        c_basic.init(t_head_size);
        c_fast.init(t_head_size);

        check_consistency();
    }
};

void gcd_iteration(
    bool& terminated, int64& total_num_iterations, matrix_multiplier_both<2, 2>& c_multiplier, int pass
) {
    array<simd_integer, 2> head=c_multiplier.calculate_head();

    int head_bits_end_index;
    vector2 start_a=extract_bits_shifted<2>(get_ptr(head), head_bits_end_index);

    pair<matrix2, vector2> gcd_64_res;
    int gcd_64_num_iterations;
    gcd_64(start_a, gcd_64_res, gcd_64_num_iterations, c_multiplier.has_tail() || head_bits_end_index>1, gcd_num_iterations);

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
    //multiplying by: [d, -d] ; [-d d]

    //implement this in asm code (multiplying by two different matricies)
    matrix2 m=gcd_64_res.first;

    int64 d=(((~uint64(gcd_64_num_iterations)) & 1ull)<<1)-1;
    assert(d==(gcd_64_num_iterations%2==0)? -1 : 1);

    //this appears to be the same as just taking the absolute value of each matrix entry which is 2 instructions
    gcd_64_res.first[0*2+0]*=d;
    gcd_64_res.first[1*2+0]*=-d;
    gcd_64_res.first[0*2+1]*=-d;
    gcd_64_res.first[1*2+1]*=d;

    for (int x=0;x<4;++x) {
        //can use the absolute value simd instruction if this is true
        assert(gcd_64_res.first[x]>=0);
    }

    //parity*=d;
    total_num_iterations+=gcd_64_num_iterations;

    if (gcd_64_num_iterations==0) {
        terminated=true;
    }

    matrix_multiplier_fast<2, 2> old_state=c_multiplier.c_fast;

    c_multiplier.multiply({m, gcd_64_res.first}, nullptr, pass);

    matrix_multiplier_fast<2, 2> new_state=c_multiplier.c_fast;

    if (test_asm_funcs) {
        multiplier_gcd_asm(
            old_state.ints, old_state.previous_batch, old_state.current_batch, {m, gcd_64_res.first}, old_state.shift_amount, pass
        );

        //todo old_state.ints[0][0].calculate_carry();
        //todo old_state.ints[0][1].calculate_carry();

        assert(old_state.ints==new_state.ints);
        assert(old_state.previous_batch==new_state.previous_batch);

        for (int vector_index=0;vector_index<2;++vector_index) {
            for (int x=0;x<=pass;++x) {
                assert(old_state.current_batch[vector_index][pass]==new_state.current_batch[vector_index][pass]);
            }
        }

        assert(old_state.shift_amount==new_state.shift_amount);
    }
}

void gcd(simd_integer& a, simd_integer& b, simd_integer& v0, int64& v0_sign, bool v0_is_1) {
    assert(a.current_size()==b.current_size());

    v0.memory.clear();
    v0.memory.resize(a.current_size(), 0);

    simd_integer v1;
    v1.memory.resize(a.current_size(), 0);

    bool terminated=false;

    ((v0_is_1)? v0 : v1).memory.at(0)=1;

    matrix_multiplier_both<2, 2> c_multiplier;
    c_multiplier.assign_int(0, 0, a);
    c_multiplier.assign_int(0, 1, b);
    c_multiplier.assign_int(1, 0, v0);
    c_multiplier.assign_int(1, 1, v1);

    c_multiplier.init(gcd_head_size);

    int64 total_num_iterations=0;

    while (!terminated) {
        for (int x=0;x<4;++x) {
            gcd_iteration(terminated, total_num_iterations, c_multiplier, x);
        }

        matrix_multiplier_fast<2, 2> old_state=c_multiplier.c_fast;

        c_multiplier.advance(true);

        matrix_multiplier_fast<2, 2> new_state=c_multiplier.c_fast;

        if (test_asm_funcs) {
            multiplier_gcd_asm(
                old_state.ints, old_state.previous_batch, old_state.current_batch,
                get_identity_matricies<2, 2>(), old_state.shift_amount, 4
            );
            assert(old_state.ints==new_state.ints);
            assert(old_state.previous_batch==new_state.previous_batch);
            assert(old_state.shift_amount==new_state.shift_amount);
        }
    }

    for (int pass=0;pass<4;++pass) {
        matrix_multiplier_fast<2, 2> old_state=c_multiplier.c_fast;

        c_multiplier.multiply(get_identity_matricies<2, 2>(), nullptr, pass);

        matrix_multiplier_fast<2, 2> new_state=c_multiplier.c_fast;

        if (test_asm_funcs) {
            multiplier_gcd_asm(
                old_state.ints, old_state.previous_batch, old_state.current_batch,
                get_identity_matricies<2, 2>(), old_state.shift_amount, pass
            );
            assert(old_state.ints==new_state.ints);
            assert(old_state.previous_batch==new_state.previous_batch);

            for (int vector_index=0;vector_index<2;++vector_index) {
                for (int x=0;x<=pass;++x) {
                    assert(old_state.current_batch[vector_index][pass]==new_state.current_batch[vector_index][pass]);
                }
            }

            assert(old_state.shift_amount==new_state.shift_amount);
        }
    }

    c_multiplier.advance(false);

    assert(!c_multiplier.has_tail());
    {
        array<simd_integer, 2> head=c_multiplier.calculate_head();
        for (int x=0;x<head[1].current_size();++x) {
            assert(head[1].memory[x]==0);
        }
    }

    a=c_multiplier.get_int(0, 0).subset(a.current_size()-gcd_head_size, gcd_head_size);
    b=c_multiplier.get_int(0, 1).subset(b.current_size()-gcd_head_size, gcd_head_size);
    v0=c_multiplier.get_int(1, 0);

    int64 parity=(v0_is_1)? 1 : -1;

    {
        int64 d=(((~uint64(total_num_iterations)) & 1ull)<<1)-1;
        assert(d==(total_num_iterations%2==0)? -1 : 1);

        parity*=d;
    }

    v0_sign=parity;

    //v0.fma(nullptr, v0, {uint64(parity)}, 0, false, false);
    //v0.calculate_carry(0, false, 2, false);
}


}