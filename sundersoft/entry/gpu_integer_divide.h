template<int size> GFUNC void normalize_divisor(gpu_integer<size, false>& b, int& shift_limbs, int& shift_bits) {
    shift_limbs=0;
    #pragma unroll
    for (int x=0;x<size;++x) {
        if (b[size-1]==0) {
            ++shift_limbs;
            b.left_shift_limbs(1);
        } else {
            break;
        }
    }

    shift_bits=clz(b[size-1]);
    b<<=shift_bits;
}

//ptxas does not like this and will waste half of the registers if done
//result is >= the actual reciprocal; max result is 2^63
GFUNC uint64 calculate_reciprocal(uint32 high, uint32 low) {
    uint64 both_source=uint64(low) | (uint64(high)<<32);

    uint64 both=both_source;
    both>>=2*32-53;

    both&=~(1ull<<52);

    uint64 res;

    if (both<=1) {
        res=1ull<<63;
    } else {
        --both;
        //assert(both>=1); because i forgot to comment this out, the compiler has wasted half of the registers

        uint64 bits=both;
        bits|=1023ull<<52;

        double bits_double=*(double*)&bits;
        bits_double=rcp(bits_double);
        bits=*(uint64*)&bits_double;

        bits&=(1ull<<52)-1;

        res=bits;
        ++res;

        res|=1ull<<52;
        res<<=(62-52);
    }

    return res;
}

//result is >= the actual quotient
GFUNC uint32 calculate_quotient(uint32 high, uint32 low, uint64 reciprocal) {
    uint64 both=uint64(low) | (uint64(high)<<32);

    uint64 product_high=mul_uint64_high(both, reciprocal);
    ++product_high;

    uint64 res=product_high>>(32-2);

    if (res>=1ull<<32) {
        res=(1ull<<32)-1;
    }

    return uint32(res);
}

//should pad a by 1 limb then left shift it by num_bits
template<int size_a, int size_b>
GFUNC void divide_integers_impl(
    gpu_integer<size_a, false> a, gpu_integer<size_b, false> b, int b_shift_limbs,
    gpu_integer<size_a-1, false>& q, gpu_integer<size_b, false>& r
) {
    const int max_quotient_size=size_a-1;
    gpu_integer<max_quotient_size, false> res;

    uint64 reciprocal;
    if (size_b>=2) {
        reciprocal=calculate_reciprocal(b[size_b-1], b[size_b-2]);
    } else {
        reciprocal=calculate_reciprocal(b[size_b-1], 0);
    }

    gpu_integer<size_a, false> b_shifted;
    b_shifted=b;
    b_shifted.left_shift_limbs(size_a-size_b-1); //it is already left shifted by b_shift_limbs

    int quotient_size=size_a-(size_b-b_shift_limbs);

    #pragma unroll
    for (int x=0;x<max_quotient_size;++x) {
        //this is more efficient than having an if statement without a break because of the compiler
        if (x>=quotient_size) {
            break;
        }
        {
            uint32 qj=calculate_quotient(a[size_a-1-x], a[size_a-2-x], reciprocal);

            //this is slower than using the doubles even though the doubles waste half the registers
            //ptxas generates horrible code which isn't scheduled properly
            //uint64 qj_64=((uint64(a[size_a-1-x])<<32) | uint64(a[size_a-2-x])) / uint64(b[size_b-1]);
            //uint32 qj=uint32(min( qj_64, uint64(~uint32(0)) ));

            auto b_shifted_qj=b_shifted;
            b_shifted_qj*=qj;
            a-=b_shifted_qj;

            #pragma unroll 1
            //while ((a[size_a-1-x] & 1u<<31) != 0) {
            while ((a[size_a-1] & 1u<<31) != 0) {
                --qj;
                a+=b_shifted;
            }

            //make the compiler optimize stuff
            //this doesn't do anything; will just leave it in anyway
            /*#pragma unroll
            for (int y=size_a-1-x;y<size_a;++y) {
                a[y]=0;
            }*/

            b_shifted.right_shift_limbs(1);

            res[max_quotient_size-1-x]=qj;
        }
    }

    #pragma unroll 1
    for (int x=0;x<max_quotient_size;++x) {
        if (quotient_size>=max_quotient_size) {
            break;
        }

        res.right_shift_limbs(1);
        ++quotient_size;
    }
    
    q=res;
    r=a;
}

template<int size_a, int size_b, bool a_is_signed>
GFUNC void divide_integers(
    gpu_integer<size_a, a_is_signed> a, gpu_integer<size_b, false> b,
    gpu_integer<size_a, a_is_signed>& q, gpu_integer<size_b, false>& r
) {
    int shift_limbs;
    int shift_bits;

    auto b_normalized=b;
    normalize_divisor(b_normalized, shift_limbs, shift_bits);

    gpu_integer<size_a+1, false> a_shifted;
    a_shifted=a;
    a_shifted<<=shift_bits;

    gpu_integer<size_a, false> q_unsigned;
    divide_integers_impl(a_shifted, b_normalized, shift_limbs, q_unsigned, r);

    r>>=shift_bits;

    if (a_is_signed && a.negative) {
        if (r==gpu_integer<size_b, false>(0u)) {
            q=q_unsigned;
            q=-q;
        } else {
            q=q_unsigned+gpu_integer<size_a, false>(1u);
            q=-q; //q'=-q-1
            r=b-r;
        }
    } else {
        q=q_unsigned;
    }

    #ifdef TEST_GPU_CODE
    {
        integer a_int=integer(a);
        integer b_int=integer(b);
        integer q_expected=a_int/b_int;
        integer r_expected=a_int%b_int;
        integer q_actual=q;
        integer r_actual=r;
        assert(q_expected==q_actual);
        assert(r_expected==r_actual);
    }
    #endif
}

//the compiler is too stupid to merge these into a single call to divide_integers if both are used on the same inputs

template<int size_a, int size_b, bool a_is_signed>
GFUNC gpu_integer<size_a, a_is_signed> operator/(
    gpu_integer<size_a, a_is_signed> a, gpu_integer<size_b, false> b
) {
    gpu_integer<size_a, a_is_signed> q;
    gpu_integer<size_b, false> r;
    divide_integers(a, b, q, r);
    return q;
}

template<int size_a, int size_b, bool a_is_signed>
GFUNC gpu_integer<size_b, false> operator%(
    gpu_integer<size_a, a_is_signed> a, gpu_integer<size_b, false> b
) {
    gpu_integer<size_a, a_is_signed> q;
    gpu_integer<size_b, false> r;
    divide_integers(a, b, q, r);
    return r;
}