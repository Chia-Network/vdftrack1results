namespace simd_integer_namespace {


//input should have no carry
//the most significant bit of b is set after this is done so it appears negative, but it is unsigned
//this can shrink b but not grow it
//
//can assume the highest two limbs of b are nonzero
int normalize_divisor(simd_integer& b) {
    assert(!b.is_negative());

    int amount=0;
    while (true) {
        uint64 b_msb=b.memory.at(b.current_size()-1);
        if (b_msb==0) {
            assert(b.current_size()>=2); //can't divide by 0
            b.memory.resize(b.current_size()-1);
            continue;
        }

        amount=data_size - (64-__builtin_clzll(uint32(b_msb)));
        break;
    }

    assert(amount<data_size);

    b.logical_shift_left(b, amount);
    return amount;
}

//result is >= the actual reciprocal; max result is 2^63
uint64 calculate_reciprocal(uint64 high, uint64 low) {
    assert( ( high & (1ull<<(data_size-1)) ) != 0 );
    assert(2*data_size>=53 && 2*data_size<=64);

    uint64 both_source=low | (high<<data_size);

    uint64 both=both_source;
    both>>=2*data_size-53;

    both&=~(1ull<<52);

    uint64 res;

    if (both<=1) {
        res=1ull<<63;
    } else {
        --both;
        assert(both>=1);

        double_bits bits; //bits>1 and bits<2
        bits.fraction=both;
        bits.set_exponent(0);
        
        bits=double_bits(1.0/bits.to_double()); //bits>=1/2 and bits<1
        assert(bits.exponent==1022);

        res=bits.fraction;
        ++res;

        res|=1ull<<52;
        res<<=(62-52);
    }

    //testing
    //res|=bit_sequence(0, 40);

    return res;
}

//result is >= the actual quotient
uint64 calculate_quotient(uint64 high, uint64 low, uint64 reciprocal) {
    uint64 both=low | (high<<data_size);

    uint64 product_high=(uint128(both)*uint128(reciprocal))>>64;
    ++product_high;

    uint64 res=product_high>>(data_size-2);

    if (res>=1ull<<data_size) {
        res=(1ull<<data_size)-1;
    }

    return res;
}

//grows a by 1 limb
//size of q should be size of a + 1
void divide_integers_impl(simd_integer& a, simd_integer& b, simd_integer* q, int amount, bool calculate_remainder) {
    assert(!a.is_negative());

    int n=b.current_size();
    int c_m=(a.current_size()+1)-b.current_size();

    if (q!=nullptr) {
        assert(q->current_size()==a.current_size()-1);
        q->clear();
    }

    if (c_m<=0) {
        //quotient is 0 and remainder is a
        return;
    }

    a.memory.resize(a.current_size()+2, 0);
    a.logical_shift_left(a, amount);

    assert(b.is_negative()); //the most significant bit has to be set (not actually negativel unsigned)
    uint64 reciprocal=calculate_reciprocal(b.memory.at(b.current_size()-1), (b.current_size()==1)? 0 : b.memory.at(b.current_size()-2));
    //uint64 reciprocal=calculate_reciprocal(b.memory_msb(0), 0);

    int current_num_adds=0;
    for (int j=c_m-1;j>=0;--j) {
        assert(a.current_size()==n+j+2);
        assert(a.memory.at(a.current_size()-1)==0);

        //uint64 qj_slow=(a.memory_msb(1)<<data_size | a.memory_msb(2))/b.memory_msb(0);

        uint64 qj=calculate_quotient(a.memory.at(a.current_size()-2), a.memory.at(a.current_size()-3), reciprocal);

        //this assert doesn't work if b.memory_msb(1)!=0
        //assert(qj>=qj_slow);

        if (current_num_adds+1>max_num_adds) {
            a.calculate_carry();
            current_num_adds=0;
        }

        if (qj!=0) {
            a.fma(&a, b, {-qj}, j, false, false);
            ++current_num_adds;

            int carry_start=a.calculate_carry_msb_start(4);
            a.calculate_carry(carry_start);
            a.check_partial_carry_valid(carry_start+1);

            //if (!a.calculate_carry_msb(carry_start)) {
                //can just return false if this happens; can add a threshold to calculate_carry_msb to make it less likely
                //a.calculate_carry();
            //}

            int overshoot=0;

            //this is so rare that division can just fail if it happens
            while (a.is_negative()) {
                --qj;
                a.fma(&a, b, {1}, j, false, false);
                a.calculate_carry();

                ++overshoot;
            }
        }

        //print(overshoot);

        assert(a.memory.at(a.current_size()-1)==0 && a.memory.at(a.current_size()-2)==0);
        a.memory.resize(a.current_size()-1);
        assert(!a.is_negative());

        if (q!=nullptr) {
            q->memory.at(j)=qj;
        }
    }

    a.calculate_carry();
    a.logical_shift_right(a, amount);
}

void divide_integers(simd_integer& a, simd_integer& b, simd_integer& q, simd_integer& r, int64 a_sign, int expected_q_size) {
    simd_integer tmp=b;
    int amount=normalize_divisor(tmp);

    q.memory.resize(a.current_size()-1, 0);

    r=a;
    divide_integers_impl(r, tmp, &q, amount, true);
    assert(r.memory.size()<=b.memory.size() && !r.is_negative());
    r.memory.resize(b.memory.size(), 0);

    q.shrink(expected_q_size);

    //a, b positive
    // a/b: a=r+qb ; 0<=r<b
    // (-a)/b: -a = -r-qb = -r-(q+1)b+b = (b-r) - (q+1)b
    //need to add 1 to q and negate it
    //-this is the same as negating q then subtracting 1
    //-negating q is the same as inverting the bits and adding 1
    //-so, this is the same as inverting all of the bits in q
    //also need to set the remainder to b-r
    //-this negates r then adds b to it. same as inverting r, then adding b+1

    if (a_sign==-1) {
        q.multiply_one(-1, true); // q'=-q-1
        r.fma(&b, r, {uint64(int64(-1))}, 0, false, false); //r'=b-r (inputs and result are nonnegative)
        r.calculate_carry();
    }
}


}