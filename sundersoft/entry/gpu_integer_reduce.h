//assumes the msb bits are 0
template<int size_a, int size_b, int size_c> GFUNC void normalize(
    gpu_integer<size_a, false>& a, gpu_integer<size_b, true>& b, gpu_integer<size_c, false>& c
) {
    #ifndef USE_GPU
        assert(size_c>=size_b && size_b>=size_a);
    #endif

    gpu_integer<size_b, true> q;
    gpu_integer<size_a, false> r;

    {
        auto a2=gpu_integer<size_a, false>(a);
        a2<<=1;

        auto ab=gpu_integer<size_b, true>(a-b);

        divide_integers(ab, a2, q, r);
    }

    auto A = gpu_integer<size_a, false>(a);
    auto B = gpu_integer<size_a, true>( gpu_integer<size_a, true>(a) - gpu_integer<size_a, true>(r) );

    auto Bb = gpu_integer<size_a, true>( B+b );
    auto Bbq = gpu_integer<size_a+size_b, true>( Bb*q );
    auto C = gpu_integer<size_c, false>( (Bbq>>1) + c );

    a=A;
    b=gpu_integer<size_b, true>( B );
    c=C;
}

//assumes the msb bits are 0
//must call normalize first before calling this
template<int size> GFUNC void reduce(
    gpu_integer<size, false>& a, gpu_integer<size, true>& b, gpu_integer<size, false>& c
) {
    #ifdef TEST_GPU_CODE
        integer a_int=a;
        integer b_int=b;
        integer c_int=c;
    #endif

    #pragma unroll 1
    while (a>c || (a==c && b<gpu_integer<size, true>(0u))) {
        auto new_a=c;
        auto new_b=-b;
        auto new_c=a;

        normalize(new_a, new_b, new_c);

        a=new_a;
        b=new_b;
        c=new_c;
    }

    normalize(a, b, c);

    #ifdef TEST_GPU_CODE
        reduce(a_int, b_int, c_int);
        assert(integer(a)==a_int);
        assert(integer(b)==b_int);
        assert(integer(c)==c_int);
    #endif
}