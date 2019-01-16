constexpr GFUNC int bits_to_words(int bits) {
    if (bits%32==0) {
        return bits/32;
    } else {
        return (bits/32)+1;
    }
}

template<int half_d_bits> struct gpu_form {
    static const int size_1=bits_to_words(half_d_bits*1);
    static const int size_2=bits_to_words(half_d_bits*2);
    static const int size_3=bits_to_words(half_d_bits*3);
    static const int size_4=bits_to_words(half_d_bits*4);

    gpu_integer<size_1, false> a;
    gpu_integer<size_1, true> b;
    gpu_integer<size_2, false> c;

    GFUNC gpu_form() {}
    gpu_form(const form& f) {
        a=f.a;
        b=f.b;
        c=f.c;
    }

    operator form() const {
        form res;
        res.a=a;
        res.b=b;
        res.c=c;
        return res;
    }

    gpu_form& operator=(const form& f) { return *this=gpu_form(f); }

    GFUNC int hash() const {
        uint64 res=uint64(c[0]) | (uint64(c[1])<<32);
        return int((res>>4) & ((1ull<<31)-1));
    }
};

template<int size> GFUNC gpu_integer<size, false> three_gcd(
    gpu_integer<size, true> a, gpu_integer<size, true> b, gpu_integer<size, true> c
) {
    auto res1=gcd(a, b);
    auto res2=gcd(gpu_integer<size, true>(res1.gcd), c);
    return res2.gcd;
}

//a and b are N bits and m is M bits; outputs are M bits
//a and b are signed and m is unsigned
//mu and v are unsigned
template<bool calculate_v, int size_ab, int size_m>
GFUNC void solve_linear_congruence(
    gpu_integer<size_ab, true> a, gpu_integer<size_ab, true> b, gpu_integer<size_m, false> m,
    gpu_integer<size_m, false>& mu, gpu_integer<size_m, false>& v
) {
    #ifndef USE_GPU
        assert(size_ab>=size_m);
    #endif

    auto gcd_res=gcd(a, gpu_integer<size_ab, true>(m));
    auto g=gpu_integer<size_m, false>( gcd_res.gcd );
    auto d=gcd_res.s;

    auto q=b/g;

    mu=(q*d)%m;
    if (calculate_v) {
        v=m/g;
    }
}

template<int half_d_bits> gpu_form<half_d_bits> GFUNC multiply(gpu_form<half_d_bits> f1, gpu_form<half_d_bits> f2) {
    const int size_1=gpu_form<half_d_bits>::size_1;
    const int size_2=gpu_form<half_d_bits>::size_2;
    const int size_3=gpu_form<half_d_bits>::size_3;
    const int size_4=gpu_form<half_d_bits>::size_4;

    auto g = gpu_integer<size_1, true>( (f2.b + f1.b) >> 1 );
    auto h = gpu_integer<size_1, true>( (f2.b - f1.b) >> 1 );

    auto w = gpu_integer<size_1, false>( three_gcd(gpu_integer<size_1, true>(f1.a), gpu_integer<size_1, true>(f2.a), g) );

    auto s = gpu_integer<size_1, false>( f1.a / w );
    auto t = gpu_integer<size_1, false>( f2.a / w );
    auto u = gpu_integer<size_1, true>( g / w );

    auto tu = gpu_integer<size_2, true>( t*u );
    auto f1c_s = gpu_integer<size_2, false>( f1.c*s );
    auto hu = gpu_integer<size_2, true>( h*u );
    auto st = gpu_integer<size_2, false>(s * t);

    gpu_integer<size_2, false> k_temp;
    gpu_integer<size_2, false> constant_factor_copy;
    solve_linear_congruence<true>(
        tu,
        gpu_integer<size_2, true>(hu + f1c_s),
        gpu_integer<size_2, false>(st),
        k_temp,
        constant_factor_copy
    );

    auto constant_factor = gpu_integer<size_1, false>( constant_factor_copy );

    auto constant_factor_t = gpu_integer<size_2, false>( constant_factor*t );
    auto k_temp_t = gpu_integer<size_3, false>( k_temp*t );

    gpu_integer<size_1, false> n;
    gpu_integer<size_1, false> constant_factor_2_unused;
    solve_linear_congruence<false>(
        gpu_integer<size_3, true>(constant_factor_t),
        gpu_integer<size_3, true>(h - k_temp_t),
        s,
        n,
        constant_factor_2_unused
    );

    auto constant_factor_n = gpu_integer<size_2, false>( constant_factor * n );
    auto k = gpu_integer<size_2, false>(k_temp + constant_factor_n);

    auto kt = gpu_integer<size_3, false>( k*t );
    auto l = gpu_integer<size_2, true>( gpu_integer<size_3, true>(kt - h) / s );

    auto ktu = gpu_integer<size_4, true>( k*tu );

    auto m = gpu_integer<size_2, true>( gpu_integer<size_4, true>(ktu - hu - f1c_s) / (st) );

    auto uw = gpu_integer<size_1, true>( u*w );

    auto ls = gpu_integer<size_3, true>( l*s );

    auto kl = gpu_integer<size_4, true>( k*l );

    auto mw = gpu_integer<size_2, true>( m*w );

    auto f3_a = gpu_integer<size_4, false>( st ); // padded 2 -> 4
    auto f3_b = gpu_integer<size_4, true>( uw - kt - ls ); // padded 3 -> 4
    auto f3_c = gpu_integer<size_4, false>( kl - mw );

    normalize(f3_a, f3_b, f3_c);

    auto f3_a_2 = gpu_integer<size_2, false>( f3_a );
    auto f3_b_2 = gpu_integer<size_2, true>( f3_b );
    auto f3_c_2 = gpu_integer<size_2, false>( f3_c );

    reduce(f3_a_2, f3_b_2, f3_c_2);

    gpu_form<half_d_bits> f3;
    f3.a=f3_a_2;
    f3.b=f3_b_2;
    f3.c=f3_c_2;

    #ifdef TEST_GPU_CODE
    {
        form f1_int=f1;
        form f2_int=f2;
        form f3_int=f1*f2;
        assert(f3_int.a==integer(f3.a));
        assert(f3_int.b==integer(f3.b));
        assert(f3_int.c==integer(f3.c));
    }
    #endif

    return f3;
}