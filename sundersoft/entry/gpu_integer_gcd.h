template<int size> struct gpu_gcd_res {
    gpu_integer<size, false> gcd;
    gpu_integer<size, true> s;
};

template<int size> GFUNC gpu_gcd_res<size> gcd(gpu_integer<size, true> a_signed, gpu_integer<size, true> b_signed) {
    bool a_negative=a_signed.negative;
    bool b_negative=b_signed.negative;

    gpu_integer<size, false> a;
    a=a_signed;

    gpu_integer<size, false> b;
    b=b_signed;

    gpu_integer<size, false> u0;
    gpu_integer<size, false> u1;
    int parity;

    if (a<b) {
        auto a_copy=a;
        a=b;
        b=a_copy;

        u0=0u;
        u1=1u;
        parity=-1;
    }  else {
        u0=1u;
        u1=0u;
        parity=1;
    }

    #pragma unroll 1
    while (b!=gpu_integer<size, false>(0u)) {
        gpu_integer<size, false> q;
        gpu_integer<size, false> r;

        divide_integers(a, b, q, r);

        a=b;
        b=r;

        auto u1_new=u0 + q*u1;

        u0=u1;
        u1=u1_new;
        parity=-parity;
    }

    // sa+bt=g ; all nonnegative
    // (-s)(-a)+bt=g
    // sa+(-b)(-t)=g
    // (-s)(-a)+(-b)(-t)=g
    // sign of each cofactor is the sign of the input

    gpu_gcd_res<size> res;
    res.gcd=a;
    res.s=u0;
    res.s.negative=a_negative != (parity==-1);

    #ifdef TEST_GPU_CODE
    {
        auto expected_gcd_res=gcd(integer(a_signed), integer(b_signed));
        assert(expected_gcd_res.gcd==integer(res.gcd));
        assert(expected_gcd_res.s==integer(res.s));
    }
    #endif

    return res;
}