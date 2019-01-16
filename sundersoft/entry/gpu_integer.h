#ifdef __CUDA_ARCH__
    #define GFUNC __host__ __device__ __forceinline__
#else
    #define GFUNC
#endif

static uint32 global_carry_flag=0;

constexpr GFUNC int max_constexpr(int a, int b) {
    if (a>b) {
        return a;
    } else {
        return b;
    }
}

#ifdef __CUDA_ARCH__
    GFUNC uint32 addc_cc(uint32 a, uint32 b) {
        uint32 res;
        asm volatile( "addc.cc.u32 %0, %1, %2;" : "=r"(res) : "r"(a), "r"(b) );
        return res;
    }

    GFUNC uint32 subc_cc(uint32 a, uint32 b) {
        uint32 res;
        asm volatile( "subc.cc.u32 %0, %1, %2;" : "=r"(res) : "r"(a), "r"(b) );
        return res;
    }

    GFUNC uint32 add_cc(uint32 a, uint32 b) {
        uint32 res;
        asm volatile( "add.cc.u32 %0, %1, %2;" : "=r"(res) : "r"(a), "r"(b) );
        return res;
    }

    GFUNC uint32 sub_cc(uint32 a, uint32 b) {
        uint32 res;
        asm volatile( "sub.cc.u32 %0, %1, %2;" : "=r"(res) : "r"(a), "r"(b) );
        return res;
    }

    GFUNC uint32 addc(uint32 a, uint32 b) {
        uint32 res;
        asm volatile( "addc.u32 %0, %1, %2;" : "=r"(res) : "r"(a), "r"(b) );
        return res;
    }

    GFUNC uint32 subc(uint32 a, uint32 b) {
        uint32 res;
        asm volatile( "subc.u32 %0, %1, %2;" : "=r"(res) : "r"(a), "r"(b) );
        return res;
    }

    GFUNC uint32 clz(uint32 a) {
        uint32 res;
        asm( "clz.b32 %0, %1;" : "=r"(res) : "r"(a) );
        return res;
    }

    GFUNC double rcp(double a) {
        double res;
        asm( "rcp.rz.f64 %0, %1;" : "=d"(res) : "d"(a) );
        return res;
    }

    GFUNC uint64 mul_uint64_high(uint64 a, uint64 b) {
        uint64 res;
        asm( "mul.hi.u64 %0, %1, %2;" : "=l"(res) : "l"(a), "l"(b) );
        return res;
    }
#else
    GFUNC uint32 addc_cc(uint32 a, uint32 b) {
        uint64 res=uint64(a) + uint64(b) + uint64(global_carry_flag);
        global_carry_flag=uint32(res>>32);
        return uint32(res);
    }

    GFUNC uint32 subc_cc(uint32 a, uint32 b) {
        uint64 res=uint64(a) - uint64(b) - uint64(global_carry_flag);
        global_carry_flag=uint32(res>>32) & 1;
        return uint32(res);
    }

    GFUNC uint32 add_cc(uint32 a, uint32 b) {
        uint64 res=uint64(a) + uint64(b);
        global_carry_flag=uint32(res>>32);
        return uint32(res);
    }

    GFUNC uint32 sub_cc(uint32 a, uint32 b) {
        uint64 res=uint64(a) - uint64(b);
        global_carry_flag=uint32(res>>32) & 1;
        return uint32(res);
    }

    GFUNC uint32 addc(uint32 a, uint32 b) {
        uint64 res=uint64(a) + uint64(b) + uint64(global_carry_flag);
        return uint32(res);
    }

    GFUNC uint32 subc(uint32 a, uint32 b) {
        uint64 res=uint64(a) - uint64(b) - uint64(global_carry_flag);
        return uint32(res);
    }

    GFUNC uint32 clz(uint32 a) {
        if (a==0) {
            return 32;
        } else {
            return __builtin_clz(a);
        }
    }

    GFUNC double rcp(double a) {
        //assumes the round-to-zero rounding mode is used
        return 1/a;
    }

    GFUNC uint64 mul_uint64_high(uint64 a, uint64 b) {
        return uint64((uint128(a)*uint128(b))>>64);
    }
#endif

//all "=" operators truncate ; all operators that return a separate result will pad the result as necessary
template<int size, bool is_signed> struct gpu_integer {
    uint32 data[size]; //little endian
    bool negative=false;

    GFUNC gpu_integer() {
        #ifndef USE_GPU
            assert(size>=1);
        #endif

        #pragma unroll
        for (int x=0;x<size;++x) {
            data[x]=0;
        }
        negative=false;
    }

    gpu_integer(const integer& i) : gpu_integer() {
        assert(i.num_bits()<=size*32);
        if (i<0) {
            assert(is_signed);
            negative=true;
        }

        mpz_export(data, nullptr, -1, 4, -1, 0, i.impl);
    }

    operator integer() const {
        integer res;
        mpz_import(res.impl, size, -1, 4, -1, 0, data);

        if (is_signed && negative) {
            res=-res;
        }

        return res;
    }

    GFUNC gpu_integer(uint32 v) : gpu_integer() {
        #ifndef USE_GPU
            assert(size>=1);
        #endif
        data[0]=v;
    }

    GFUNC gpu_integer(int32 v) : gpu_integer() {
        #ifndef USE_GPU
            assert(size>=1);
            assert(is_signed);
        #endif
        data[0]=(v<0)? -v : v;
        negative=(v<0);
    }

    GFUNC gpu_integer(uint64 v) : gpu_integer() {
        #ifndef USE_GPU
            assert(size>=2);
        #endif
        data[0]=v;
        data[1]=v>>32;
    }

    GFUNC gpu_integer(int64 v) : gpu_integer() {
        #ifndef USE_GPU
            assert(is_signed);
            assert(size>=2);
        #endif
        data[0]=(v<0)? -v : v;
        data[1]=((v<0)? -v : v) >> 32;
        negative=(v<0);
    }

    //truncation
    template<int t_size, bool t_is_signed> GFUNC explicit gpu_integer(gpu_integer<t_size, t_is_signed> t) {
        if (is_signed) {
            negative=(t_is_signed)? t.negative : false;
        }

        #pragma unroll
        for (int x=0;x<size;++x) {
            data[x]=(x<t_size)? t[x] : 0;
        }
    }

    GFUNC gpu_integer& operator=(uint32 v) { return *this=gpu_integer(v); }
    GFUNC gpu_integer& operator=(int32 v) { return *this=gpu_integer(v); }
    GFUNC gpu_integer& operator=(uint64 v) { return *this=gpu_integer(v); }
    GFUNC gpu_integer& operator=(int64 v) { return *this=gpu_integer(v); }
    gpu_integer& operator=(const integer& v) { return *this=gpu_integer(v); }
    template<int t_size, bool t_is_signed> GFUNC gpu_integer& operator=(gpu_integer<t_size, t_is_signed> t) { return *this=gpu_integer(t); }

    GFUNC uint32& operator[](int pos) {
        #ifndef USE_GPU
            assert(pos>=0 && pos<size);
        #endif
        return data[pos];
    }

    GFUNC const uint32& operator[](int pos) const {
        #ifndef USE_GPU
            assert(pos>=0 && pos<size);
        #endif
        return data[pos];
    }

    //nonnegative: unchanged
    //negative: result is -v-1 ; still need to add 1
    GFUNC void invert_if_negative() {
        if (!is_signed) {
            return;
        }

        #pragma unroll
        for (int x=0;x<size;++x) {
            data[x]^=(negative)? ~0u : 0u;
        }
    }

    GFUNC void add_small(uint32 extra) {
        uint64 v=uint64(data[0])+extra;
        data[0]=uint32(v);

        if (v>=1ull<<32) {
            #pragma unroll
            for (int x=1;x<size;++x) {
                if (x==1) {
                    (*this)[x]=add_cc((*this)[x], 1);
                } else {
                    (*this)[x]=addc_cc((*this)[x], 0);
                }
            }
        }
    }

    GFUNC void restore_negative() {
        if (!is_signed) {
            return;
        }

        negative=int32(data[size-1])<0;

        invert_if_negative();
        add_small((negative)? 1 : 0);
    }

    GFUNC gpu_integer operator-() const {
        #ifndef USE_GPU
            assert(is_signed);
        #endif
        gpu_integer res=*this;
        res.negative=!res.negative;
        return res;
    }

    GFUNC void add_twos_compliment(gpu_integer b) {
        if (is_signed) {
            invert_if_negative();
            b.invert_if_negative();

            uint32 extra=(negative)? 1 : 0;
            extra+=(b.negative)? 1 : 0;
            add_small(extra);
        }

        #pragma unroll
        for (int x=0;x<size;++x) {
            if (x==0) {
                (*this)[x]=add_cc((*this)[x], b[x]);
            } else {
                (*this)[x]=addc_cc((*this)[x], b[x]);
            }
        }
    }

    GFUNC void operator+=(gpu_integer b) {
        add_twos_compliment(b);
        restore_negative();
    }

    GFUNC void subtract_twos_compliment(gpu_integer b) {
        if (is_signed) {
            b.negative=!b.negative;
            add_twos_compliment(b);
        } else {
            #pragma unroll
            for (int x=0;x<size;++x) {
                if (x==0) {
                    (*this)[x]=sub_cc((*this)[x], b[x]);
                } else {
                    (*this)[x]=subc_cc((*this)[x], b[x]);
                }
            }
        }
    }

    GFUNC void operator-=(gpu_integer b) {
        subtract_twos_compliment(b);
        restore_negative();
    }

    GFUNC void operator*=(uint32 v) {
        gpu_integer low;
        gpu_integer high_shifted;

        uint32 previous_res_high=0;

        #pragma unroll
        for (int x=0;x<size;++x) {
            uint64 res=uint64(data[x])*uint64(v);
            uint32 res_low=uint32(res);
            uint32 res_high=uint32(res>>32);

            low[x]=res_low;
            high_shifted[x]=previous_res_high;
            previous_res_high=res_high;
        }

        *this=low;
        *this+=high_shifted;
    }

    GFUNC void operator*=(int32 v) {
        #ifndef USE_GPU
            assert(is_signed);
        #endif
        *this *= uint32(v) & ~(1u<<31);
        negative^=v & (1u<<31);
    }

    //start and end must be statically known at compile time to not spill registers
    //compiler is broken so this is static
    template<int t_size, bool t_is_signed, int this_size, bool this_is_signed>
    GFUNC static gpu_integer<t_size, t_is_signed> subset(
        gpu_integer<this_size, this_is_signed> this_v, int start
    ) {
        const int end=start+t_size;

        gpu_integer<t_size, t_is_signed> res;
        if (t_is_signed && this_is_signed) {
            res.negative=this_v.negative;
        }

        #pragma unroll
        for (int x=start;x<end;++x) {
            int pos=x-start;
            res[x]=(pos>=0 && pos<this_size)? this_v[x] : 0;
        }

        return res;
    }

    //amount must be statically known at compile time to not spill registers
    GFUNC void left_shift_limbs(uint32 amount) {
        #pragma unroll
        for (int x=size-1;x>=0;--x) {
            int pos=x-amount;
            (*this)[x] = (pos>=0 && pos<size)? (*this)[pos] : 0;
        }
    }

    //amount must be statically known at compile time to not spill registers
    GFUNC void right_shift_limbs(uint32 amount) {
        #pragma unroll
        for (int x=0;x<size;++x) {
            int pos=x+amount;
            (*this)[x] = (pos>=0 && pos<size)? (*this)[pos] : 0;
        }
    }

    //must be <= 32 bits
    GFUNC void operator<<=(uint32 amount) {
        //shift by 32 is undefined on intel; returns 0 on nvidia
        #ifndef __CUDA_ARCH__
            if (amount==0) {
                return;
            }
        #endif

        #pragma unroll
        for (int x=size-1;x>=0;--x) {
            uint32 previous=(x==0)? 0 : (*this)[x-1];
            (*this)[x] = ((*this)[x]<<amount) | (previous>>(32-amount));
        }
    }

    //must be <= 32 bits
    GFUNC void operator>>=(uint32 amount) {
        #ifndef __CUDA_ARCH__
            if (amount==0) {
                return;
            }
        #endif

        #pragma unroll
        for (int x=0;x<size;++x) {
            uint32 next=(x==size-1)? 0 : (*this)[x+1];
            (*this)[x] = ((*this)[x]>>amount) | (next<<(32-amount));
        }
    }

    template<int b_size, bool b_is_signed>
    GFUNC gpu_integer<size+b_size, is_signed || b_is_signed> operator*(
        gpu_integer<b_size, b_is_signed> b
    ) const {
        const int output_size=size+b_size;
        const bool output_is_signed=is_signed || b_is_signed;
        gpu_integer<output_size, false> res;

        #pragma unroll
        for (int x=0;x<b_size;++x) {
            auto r=subset<output_size, false>(*this, 0);
            r*=b[x];
            r.left_shift_limbs(x);
            res+=r;
        }

        auto res_signed=subset<output_size, output_is_signed>(res, 0);

        bool a_negative=(is_signed)? negative : false;
        bool b_negative=(b_is_signed)? b.negative : false;
        if (output_is_signed) {
            res_signed.negative=(a_negative!=b_negative);
        }

        return res_signed;
    }

    template<int b_size, bool b_is_signed>
    GFUNC gpu_integer<max_constexpr(size, b_size)+1, is_signed || b_is_signed> operator+(
        gpu_integer<b_size, b_is_signed> b
    ) const {
        const int output_size=max_constexpr(size, b_size)+1;
        const bool output_is_signed=is_signed || b_is_signed;

        auto res=subset<output_size, output_is_signed>(*this, 0);
        auto add=subset<output_size, output_is_signed>(b, 0);

        res+=add;
        return res;
    }

    GFUNC gpu_integer<size+1, is_signed> operator<<(int num) const {
        auto res=subset<size+1, is_signed>(*this, 0);
        res<<=num;
        return res;
    }

    //this rounds to 0 so it is different from division unless the input is divisible by 2^num
    GFUNC gpu_integer<size, is_signed> operator>>(int num) const {
        auto res=subset<size, is_signed>(*this, 0);
        res>>=num;
        return res;
    }

    template<int b_size, bool b_is_signed>
    GFUNC gpu_integer<max_constexpr(size, b_size)+1, is_signed || b_is_signed> operator-(
        gpu_integer<b_size, b_is_signed> b
    ) const {
        const int output_size=max_constexpr(size, b_size)+1;
        const bool output_is_signed=is_signed || b_is_signed;

        auto res=subset<output_size, output_is_signed>(*this, 0);
        auto sub=subset<output_size, output_is_signed>(b, 0);

        res-=sub;
        return res;
    }

    template<int b_size, bool b_is_signed>
    GFUNC bool operator>=(gpu_integer<b_size, b_is_signed> b) const {
        const int output_size=max_constexpr(size, b_size)+1;
        const bool output_is_signed=is_signed || b_is_signed;

        auto res=subset<output_size, output_is_signed>(*this, 0);
        auto sub=subset<output_size, output_is_signed>(b, 0);

        res.subtract_twos_compliment(sub);
        return (res[output_size-1] & 1u<<31) == 0;
    }

    template<int b_size, bool b_is_signed>
    GFUNC bool operator==(gpu_integer<b_size, b_is_signed> b) const {
        const int output_size=max_constexpr(size, b_size);

        bool this_negative=(is_signed)? negative : false;
        bool b_negative=(b_is_signed)? b.negative : false;

        bool equal=true;
        auto this_padded=subset<output_size, false>(*this, 0);
        auto b_padded=subset<output_size, false>(b, 0);

        #pragma unroll
        for (int x=0;x<output_size;++x) {
            equal=equal && (this_padded[x]==b_padded[x]);
        }

        return (this_negative==b_negative) && equal;
    }

    template<int b_size, bool b_is_signed>
    GFUNC bool operator<(gpu_integer<b_size, b_is_signed> b) const {
        return !(*this >= b);
    }

    template<int b_size, bool b_is_signed>
    GFUNC bool operator<=(gpu_integer<b_size, b_is_signed> b) const {
        return b >= *this;
    }

    template<int b_size, bool b_is_signed>
    GFUNC bool operator>(gpu_integer<b_size, b_is_signed> b) const {
        return !(b >= *this);
    }

    template<int b_size, bool b_is_signed>
    GFUNC bool operator!=(gpu_integer<b_size, b_is_signed> b) const {
        return !(b == *this);
    }
};

template<int size, bool is_signed> gpu_integer<size, is_signed> abs(gpu_integer<size, is_signed> v) {
    #ifndef USE_GPU
        assert(is_signed);
    #endif
    v.negative=false;
    return v;
}