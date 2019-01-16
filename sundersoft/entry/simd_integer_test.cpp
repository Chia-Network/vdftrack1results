#include "include.h"

#include "headers.h"

#include "simd_integer_fma.h"
#include "simd_integer_divide.h"
#include "simd_integer_gcd.h"
#include "simd_integer_reduce.h"

#include "vdf_new.h"

int main(int argc, char** argv) {
    using namespace simd_integer_namespace;

    set_rounding_mode();

    init_asm_runtime();

    parse_args(argc, argv);

    integer a;
    integer b;
    integer c;
    generator_for_discriminant(arg_discriminant, a, b, c);

    int num_failed=0;
    int num_good=0;

    for (int x=0;x<arg_iterations;++x) {
        print(x);

        square(a, b, c);

        integer a_copy=a;
        integer b_copy=b;
        integer c_copy=c;

        reduce(a, b, c);

        if (x>20) {

            normalize(a_copy, b_copy, c_copy);

            const int simd_size=18*4;

            simd_integer a_simd=from_integer(a_copy, simd_size);
            simd_integer b_simd=from_integer(b_copy, simd_size);
            simd_integer c_simd=from_integer(c_copy, simd_size);

            int64 b_sign;
            to_sign_magnitude(b_simd, b_sign);

            //reduce(a_simd, b_simd, c_simd, b_sign, 9*4);
            if (!integer_reduce_asm(a_simd, b_simd, c_simd, b_sign)) {
                cerr << "!!!!";
                continue;
            }

            from_sign_magnitude(b_simd, b_sign);

            integer a_reduced=to_integer(a_simd);
            integer b_reduced=to_integer(b_simd);
            integer c_reduced=to_integer(c_simd);

            if (a_reduced<c_reduced) {
                cerr << "G";
                normalize(a_reduced, b_reduced, c_reduced);

                assert(a_reduced==a);
                assert(b_reduced==b);
                assert(c_reduced==c);
            }
        }
    }

    /*int a_num_bits=1024;
    int b_num_bits=1024;

    int start_iter=0;
    int end_iter=-1;

    const bool do_print=false;

    if (argc>=2) {
        a_num_bits=from_string<int>(argv[1]);
        //b_num_bits=from_string<int>(argv[2]);
        b_num_bits=a_num_bits;
    }

    if (argc>=3) {
        start_iter=from_string<int>(argv[2]);
    }

    const int a_simd_size=ceil_div(a_num_bits+1, data_size*4)*4;
    const int b_simd_size=ceil_div(b_num_bits+1, data_size*4)*4;

    int failures=0;

    int iter=0;
    while (end_iter==-1 || iter<end_iter) {
        integer a=rand_integer(a_num_bits);
        integer b=rand_integer(b_num_bits);

        if (iter<start_iter) {
            ++iter;
            continue;
        }

        if (iter!=0 && iter%100==0) {
            print(iter, failures);
        }

        if (a<=b) {
            swap(a, b);
            if (a==b) {
                a=a+integer(1);
            }
        }

        if (do_print) {
            print(a.to_string());
            print( "" );

            print(b.to_string());
            print( "" );
        }

        simd_integer a_simd=from_integer(a, a_simd_size);
        simd_integer b_simd=from_integer(b, b_simd_size);

        simd_integer a_v0_simd;
        //simd_integer b_v0_simd;
        
        //simd_integer a_simd_copy=a_simd;
        //simd_integer b_simd_copy=b_simd;

        //int64 b_v0_simd_sign;
        //gcd(a_simd, b_simd, b_v0_simd, b_v0_simd_sign, false);
        //from_sign_magnitude(b_v0_simd, b_v0_simd_sign);

        //int64 a_v0_simd_sign;
        //gcd(a_simd_copy, b_simd_copy, a_v0_simd, a_v0_simd_sign, true);
        //from_sign_magnitude(a_v0_simd, a_v0_simd_sign);

        if (!integer_gcd_asm(a_simd, b_simd, a_v0_simd)) {
        ... need to use v0_sign, etc
        ... also forced the result to be 1 and calculated the other cofactor
            ++failures;
            ++iter;
            continue;
        }

        integer g=to_integer(a_simd);
        integer s=to_integer(a_v0_simd);
        //integer t=to_integer(b_v0_simd);

        gcd_res expected_gcd=gcd(a, b);

        if (do_print) {
            print( "g", g.to_string() );
            print( "s", s.to_string() );
            //print( "t", t.to_string() );
            //print( "a*s+b*t", (a*s+b*t).to_string() );

            print( "gcd_res-expected.gcd", (g-expected_gcd.gcd).to_string() );
            print( "s-expected.s", (s-expected_gcd.s).to_string() );
            //print( "t-expected.t", (t-expected_gcd.t).to_string() );
        } else {
            assert(g==expected_gcd.gcd);
            assert(s==expected_gcd.s);
            //assert(t==expected_gcd.t);
        }

        ++iter;
    }*/

    /*const int num_bits=1024;
    const int simd_size=9*4;

    //integer a=rand_integer(1000); //num_bits);
    //integer b=rand_integer(24); //num_bits);
    integer a=rand_integer(num_bits);
    integer b=rand_integer(num_bits);

    //a=-a;
    //b=-b;

    simd_integer a_simd=from_integer(a, simd_size);
    simd_integer b_simd=from_integer(b, simd_size);

    simd_integer c_simd;
    //c_simd.memory.resize(a_simd.current_size()); //a_simd.current_size()+b_simd.current_size());
    c_simd.memory.resize(a_simd.current_size()+b_simd.current_size());

    //integer_multiply(c_simd, a_simd, b_simd);
    integer_multiply_asm(c_simd, a_simd, b_simd);

    //c_simd.fma(nullptr, a_simd, {b_simd.memory.at(0)}, 0, false, false);
    //integer_fma_asm(c_simd, a_simd, b_simd);

    print( "a", to_integer(a_simd).to_string());
    print( "" );

    print( "b", to_integer(b_simd).to_string());
    print( "" );

    print( "c", to_integer(c_simd).to_string());
    print( "" );

    print( "a*b", (a*b).to_string());

    print( "c-a*b", (to_integer(c_simd)-a*b).to_string());
    print( "" );*/

    //
    //

    /*int amount=normalize_divisor(b_simd);

    simd_integer q_simd;
    q_simd.memory.resize(c_simd.current_size()+1, 0);
    divide_integers(c_simd, b_simd, &q_simd, amount, true);

    print( "q", to_integer(q_simd).to_string());
    print( "" );
    
    print( "r", to_integer(c_simd).to_string());
    print( "" );

    print( "q-a", (to_integer(q_simd)-a).to_string());
    print( "" );*/
}