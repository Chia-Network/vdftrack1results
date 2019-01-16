namespace simd_integer_namespace {


typedef array<int64, 3> vector3;
typedef array<int64, 9> matrix3; //row major
typedef array<int64, 2> vector2;
typedef array<int64, 4> matrix2; //row major

const int data_size=29;
const int carry_size=64-data_size; //35
const uint64 data_mask=bit_sequence(0, data_size);
const uint64 carry_mask=bit_sequence(data_size, carry_size);
const uint64 data_sign_mask=1ull<<(data_size-1);

const int calculate_carry_msb_num=2;

//multiplication output has twice as many bits, and there is one sign bit
//after a single carry iteration, the carry will have some bits (carry_bits-data_bits+1 = 7), so need to use up one add
const int max_num_adds=(1<<(64-1-2*data_size))-2;

int ceil_div(int a, int b) {
    return (a+b-1)/b;
}

int round_down(int a, int b) {
    return a/b*b;
}

int round_up(int a, int b) {
    return ceil_div(a, b)*b;
}

uint64 sign_extend_data(uint64 v) {
    return int64(v<<(64-data_size)) >> (64-data_size);
}

uint64 zero_extend_data(uint64 v) {
    return v & data_mask;
}

uint64 sign_extend_carry(uint64 v) {
    return int64(v) >> (64-carry_size);
}

uint64 zero_extend_carry(uint64 v) {
    return v >> (64-carry_size);
}

struct simd_integer {
    vector<uint64> memory;

    int current_size() {
        return memory.size();
    }

    bool is_negative() {
        return int64(sign_extend_data(memory.back()))<0;
    }

    simd_integer subset(int start, int size) {
        simd_integer res;
        for (int x=0;x<size;++x) {
            res.memory.push_back(memory.at(start+x));
        }
        return res;
    }

    void assign_subset(int start, simd_integer& data) {
        for (int x=0;x<data.current_size();++x) {
            memory.at(start+x)=data.memory.at(x);
        }
    }

    bool operator==(const simd_integer& c) const {
        return memory==c.memory;
    }

    USED simd_integer compare(simd_integer& b) {
        simd_integer res;
        for (int x=0;x<current_size();++x) {
            uint64 v=memory[x]-((x>=b.current_size())? 0 : b.memory[x]);
            res.memory.push_back(v);
        }
        for (int x=res.current_size();x<b.current_size();++x) {
            res.memory.push_back(-b.memory[x]);
        }
        return res;
    }

    void left_shift_limbs(int num) {
        assert(!is_negative());

        for (int x=memory.size()-num;x<memory.size();++x) {
            assert(memory.at(x)==0);
        }

        vector<uint64> new_memory;
        for (int x=0;x<num;++x) {
            new_memory.push_back(0);
        }

        for (int x=0;x<memory.size()-num;++x) {
            new_memory.push_back(memory.at(x));
        }

        memory=new_memory;
        assert(!is_negative());
    }

    //
    //

    void shrink(int expected_size) {
        assert(expected_size>=1);
        while (current_size()>expected_size) {
            assert(memory.back()==0);
            memory.pop_back();
            assert(!is_negative());
        }
    }

    void calculate_carry_loop(
        uint64& accumulator,
        int start_index, bool single_pass,
        int range_start, int range_end, bool preserve_accumulator
    ) {
        //need to go msb first since carries are cleared at the same time the output is written
        for (int x=range_end-1;x>=range_start;--x) {
            bool first=(x==current_size()-1);

            assert(x>=start_index);

            uint64 data=0;
            uint64 carry=0;

            if (x!=start_index) {
                carry=sign_extend_carry(memory.at(x-1));

                if (x==start_index) {
                    //zero out previous 4 carries
                    for (int d=0;d<3;++d) {
                        memory.at(x-1-d)=memory.at(x-1-d) & data_mask;
                    }
                }
            }

            data=zero_extend_data(memory.at(x));

            uint64 out=data+carry;
            memory.at(x)=out;

            if (!single_pass && first && !preserve_accumulator) {
                accumulator=out;
            }

            if (!single_pass && (!first || preserve_accumulator)) {
                accumulator|=out;
            }
        }
    }

    void calculate_carry(
        int start_index=0, bool single_pass=false, int max_iterations=-1, bool debug=false
    ) {
        vector<simd_integer> states;
        int iter=0;

        if (debug) {
            states.push_back(*this);
        }

        while (true) {
            uint64 accumulator=0;
            calculate_carry_loop(
                accumulator,
                start_index, single_pass,
                start_index, current_size(), false
            );

            if (debug) {
                states.push_back(*this);
            }
            ++iter;

            if (single_pass || (accumulator & carry_mask)==0) {
                break;
            }
        }

        if (debug) {
            print( "calculate_carry:", iter );
        }

        if (max_iterations!=-1) {
            assert(iter<=max_iterations);
        }
    }

    void calculate_carry_lsb() {
        uint64 incoming_carry=0;
        for (int x=0;x<memory.size();++x) {
            uint64& c=memory.at(x);

            c+=incoming_carry;
            incoming_carry=sign_extend_carry(c);
            c&=data_mask;
        }
    }

    void multiply_one(int64 v, bool skip_add=false) {
        assert(v==1 || v==-1);

        uint64 mask=(v>>1); //-1 if negative, 0 if positive
        mask&=data_mask;

        //flip all of the bits if negative
        for (int x=0;x<memory.size();++x) {
            memory[x]^=mask;
        }

        if (!skip_add) {
            memory[0]+=uint64(v) >> 63; //1 if negative, 0 if positive

            //very unlikely to happen
            if ((memory[0] & carry_mask)!=0) {
                calculate_carry();
            }
        }
    }

    //can multiply by 1
    void copy(simd_integer& c) {
        assert(current_size()==c.current_size());
        memory=c.memory;
    }

    //can multiply by 0
    void clear() {
        for (int x=0;x<current_size();++x) {
            memory.at(x)=0;
        }
    }

    int calculate_carry_msb_start(int size) {
        int carry_start=current_size()-size-calculate_carry_msb_num;
        if (carry_start<=0) {
            carry_start=0;
        }

        return carry_start;
    }

    //should have no carry at the limb before partial_index and every limb after. values at and after partial_index must also be known
    void check_partial_carry_valid(int partial_index) {
        //memory[carry_start-1] can contain a signed integer with carry_size bits
        //when this is carried to memory[carry_start], it will invalidate the data and produece a carry
        // of up to:
        //data = 0x1FFFFFFF ; carry = 0x3FFFFFFFF ; data+carry>>data_size = 0x20
        //data = 0x0 ; carry = 0xFFFFFFFBFFFFFFFF ; data+carry>>data_size = 0xFFFFFFFFFFFFFFDF = -0x21
        //the carry will invalidate the data of memory[carry_start+1] but will only produce a carry if the previous data was
        //close to the min or max range

        if (partial_index==1) {
            return;
        }

        uint64 v=memory.at(partial_index);
        assert(zero_extend_carry(v)==0);

        assert(carry_size>data_size);
        uint64 threshold=1ull<<(carry_size-data_size-1); //0x20
        assert(!(v<threshold || v>data_mask-(threshold+1)));
    }

    //assumes this is carried
    //can just multiply by 1<<bits
    void logical_shift_left(simd_integer& c, int bits) {
        assert(bits>=0 && bits<=data_size);

        //msb first
        for (int x=current_size()-1;x>=0;--x) {
            uint64 this_data=memory.at(x);
            uint64 previous_data=0;
            if (x!=0) {
                previous_data=memory.at(x-1);
            }

            assert(zero_extend_carry(this_data)==0);
            assert(zero_extend_carry(previous_data)==0);

            uint64 out=(this_data << bits) & data_mask;
            out|=previous_data >> (data_size-bits);

            memory.at(x)=out;
        }
    }

    //assumes this is carried
    void logical_shift_right(simd_integer& c, int bits) {
        assert(bits>=0 && bits<=data_size);

        //lsb first
        for (int x=0;x<current_size();++x) {
            uint64 this_data=memory.at(x);
            uint64 next_data=0;
            if (x!=current_size()-1) {
                next_data=memory.at(x+1);
            }

            assert(zero_extend_carry(this_data)==0);
            assert(zero_extend_carry(next_data)==0);

            uint64 out=this_data >> bits;
            out|=(next_data << (data_size-bits)) & data_mask;

            memory.at(x)=out;
        }
    }

    //c can alias with this or be null
    //if sign_extend_b is false, can negate each entry in b to do a subtraction
    //can be in-place for the multiplication (c==nullptr, &a=this) only if sizeof(b)==1
    //can always be in-place for the addition
    void fma(
        simd_integer* c, simd_integer& a, vector<uint64> b, int shift_amount, bool sign_extend_a, bool sign_extend_b
    ) {
        assert(c==nullptr || c->current_size()==current_size());

        assert(b.size()>=1);
        if (sign_extend_b) {
            b.back()=sign_extend_data(b.back());
        }

        for (int b_index=0;b_index<b.size();++b_index) {
            uint64 b_value=b[b_index];

            int output_start=shift_amount+b_index;
            int output_end=a.current_size()+shift_amount+b_index;

            if (output_end>current_size()) {
                output_end=current_size();
            }

            for (int x=output_end-1;x>=output_start;--x) {
                int a_pos=x-(shift_amount+b_index);

                uint64 a_value=a.memory.at(a_pos);

                if (a_pos==a.current_size()-1 && sign_extend_a) {
                    a_value=sign_extend_data(a_value);
                }

                //the multiply is always signed even if nothing is sign extended; sign bits will be 0 in that case
                uint64 ab_value=int64(a_value)*int64(b_value);

                if (c==nullptr) {
                    memory.at(x)=ab_value;
                } else {
                    memory.at(x)=c->memory.at(x)+ab_value;
                }
            }
        }
    }
};

simd_integer from_integer(const integer& t, int size) {
    assert(size>=ceil_div(t.num_bits()+1, data_size));

    simd_integer res;
    res.memory.resize(size, 0);

    mpz_export(&res.memory[0], nullptr, -1, 8, -1, carry_size, t.impl);

    if (t<0) {
        res.fma(nullptr, res, {uint64(int64(-1))}, 0, true, false);
        res.calculate_carry();
    }

    return res;
}

integer to_integer(simd_integer t, bool is_signed=true) {
    t.calculate_carry();

    bool is_negative=false;
    if (is_signed && t.is_negative()) {
        t.fma(nullptr, t, {uint64(int64(-1))}, 0, true, false);
        t.calculate_carry();
        is_negative=true;
    }

    integer res;
    mpz_import(res.impl, t.current_size(), -1, 8, -1, carry_size, &t.memory[0]);

    if (is_negative) {
        res=-res;
    }

    return res;
}

int sign_int(int64 i) {
    return 1 | (i>>63);
}

int64 abs_int(int64 v) {
    int64 mask = v>>63;
    return (v + mask) ^ mask;
}

void to_sign_magnitude(simd_integer& b, int64& b_sign) {
    b_sign=sign_int(b.memory.at(b.current_size()-1) << carry_size);
    b.multiply_one(b_sign);
}

void from_sign_magnitude(simd_integer& b, int64 b_sign) {
    b.multiply_one(b_sign);
}


};