#include "include.h"

//used for the final submission and correctness testing
#define VDF_MODE 0

//used for performance or other testing
//#define VDF_MODE 1

//
//

#if VDF_MODE==0
    const bool enable_track_cycles=false;
    const bool test_correctness=false;
    const bool assert_on_rollback=false;
    const bool debug_rollback=false;
    const int asm_inject_random_errors_rate=0;
    const int repeated_square_checkpoint_interval=1<<10; //should be a power of 2
#endif

#if VDF_MODE==1
    const bool enable_track_cycles=true;
    const int test_correctness=true;
    const bool assert_on_rollback=true;
    const bool debug_rollback=false;
    const int asm_inject_random_errors_rate=0;
    const int repeated_square_checkpoint_interval=1<<10;
#endif

const int reduce_max_iterations=10000;

#include "headers.h"
#include "vdf_original.h"

#include "vdf_new.h"

using namespace std;
using simd_integer_namespace::track_cycles_test;

generic_stats track_cycles_total;

USED string to_string(mpz_struct* t) {
    integer t_int;
    mpz_set(t_int.impl, t);
    return t_int.to_string();
}

//setting is_reduce to true will negate b but not swap a or c; caller needs to do that
bool normalize_fast(integer& a_memory, integer& b_memory, integer& c_memory, bool is_reduce) {
    track_cycles c_track_cycles(track_cycles_test[5]); // 334

    mpz_struct* a=a_memory.impl;
    mpz_struct* b=b_memory.impl;
    mpz_struct* c=c_memory.impl;

    static mpz_t a2, a4, ab, ab_sum, q, r, B, Bb, Bbq, half_Bbq, C;

    static bool is_init=false;
    if (!is_init) {
        mpz_inits(a2, a4, ab, ab_sum, q, r, B, Bb, Bbq, half_Bbq, C, NULL);
        is_init=true;
    }

    if (is_reduce) {
        track_cycles c_track_cycles(track_cycles_test[6]); // 142
        mpz_add(ab, a, b); // ab = a-(-b)
    } else {
        track_cycles c_track_cycles(track_cycles_test[7]); // 106
        mpz_sub(ab, a, b); // ab = a-b
    }

    if (!is_reduce) {
        track_cycles c_track_cycles(track_cycles_test[8]); // 126
        mpz_add(ab_sum, a, b); // ab_sum = a+b

        // q=0:
        // A=a
        // B=b
        // C=c

        // -a<b<=a
        // -a<b && b<=a
        // 0<a+b && b-a<=0
        // a+b>0 && a-b>=0
        bool already_normalized=(mpz_sgn(ab_sum)>0 && mpz_sgn(ab)>=0); // a+b>0 && a-b>=0
        if (already_normalized) {
            //for normalize this is true about 80% of the time
            return true;
        }
    }

    {
        track_cycles c_track_cycles(track_cycles_test[9]); // 110
        mpz_mul_2exp(a2, a, 1); // a2 = 2a
    }

    if (!is_reduce) {
        // q=1:
        // a-b>=2a && a-b<4a
        // A=a
        // B=b+2a
        // C=a+b+c

        track_cycles c_track_cycles(track_cycles_test[10]); // 358
        mpz_mul_2exp(a4, a2, 1); // a2 = 2a

        bool q_is_1=mpz_cmp(ab, a2)>=0 && mpz_cmp(ab, a4)<0;

        if (q_is_1) {
            mpz_add(C, ab_sum, c); // C = a+b+c
            mpz_add(B, b, a2); // B = b+2a

            // a <- A (already done)
            mpz_swap(b, B); // b <- B
            mpz_swap(c, C); // c <- C

            //q is always 0 or 1 for normalize experimentally
            return true;
        }
    }

    if (mpz_sgn(a2)==0) {
        return false;
    }

    {
        track_cycles c_track_cycles(track_cycles_test[11]); // 506
        mpz_fdiv_qr(q, r, ab, a2); // q = (a-b)/(2a) ; r = (a-b)%(2a)
    }

    // A = a
    mpz_sub(B, a, r); // B = a-r

    if (is_reduce) {
        mpz_sub(Bb, B, b); // Bb = B+(-b)
    } else {
        mpz_add(Bb, B, b); // Bb = B+b
    }

    {
        track_cycles c_track_cycles(track_cycles_test[12]); // 158
        mpz_mul(Bbq, Bb, q); // Bbq = (B+b)*q
    }
    {
        track_cycles c_track_cycles(track_cycles_test[13]); // 170
        mpz_fdiv_q_2exp(half_Bbq, Bbq, 1); // half_Bbq = ((B+b)*q)/2
    }
    mpz_add(C, half_Bbq, c); // C = ((B+b)*q)/2 + c

    {
        track_cycles c_track_cycles(track_cycles_test[14]); // 70
        // a <- A (already done)
        mpz_swap(b, B); // b <- B
        mpz_swap(c, C); // c <- C
    }

    return true;
}

bool reduce_fast(integer& a_memory, integer& b_memory, integer& c_memory) {
    if (!normalize_fast(a_memory, b_memory, c_memory, false)) {
        return false;
    }

    //todo bool use_asm=false;
    bool use_asm=true;

    int iter=0;

    while (true) {
        if (use_asm) {
            //this will do the reduce all the way 80% of the time; need to call it multiple times the other 20%
            if (!simd_integer_namespace::integer_reduce_asm(a_memory.impl, b_memory.impl, c_memory.impl)) {
                //the asm code could not make any progress so it is going to fail again if it is called
                //this is very rare
                use_asm=false;
            }
        }

        int a_c_cmp;
        {
            track_cycles c_track_cycles(track_cycles_test[15]); // 60
            a_c_cmp=mpz_cmp(a_memory.impl, c_memory.impl); // this can be any integer, not just -1, 0, or 1
        }
        bool a_greater_than_c=(a_c_cmp>0);
        bool a_equals_c=(a_c_cmp==0);

        bool keep_going=a_greater_than_c || (a_equals_c && mpz_sgn(b_memory.impl)==-1);

        if (!keep_going) {
            break;
        }

        {
            track_cycles c_track_cycles(track_cycles_test[16]); // 1826
            if (!normalize_fast(c_memory, b_memory, a_memory, true)) {
                return false;
            }
        }
        mpz_swap(a_memory.impl, c_memory.impl);

        ++iter;
        if (iter>reduce_max_iterations) {
            return false;
        }
    }

    if (!normalize_fast(a_memory, b_memory, c_memory, false)) {
        return false;
    }

    return true;
}

bool square_fast(integer& a_memory, integer& b_memory, integer& c_memory) {
    mpz_struct* a=a_memory.impl;
    mpz_struct* b=b_memory.impl;
    mpz_struct* c=c_memory.impl;

    track_cycles c_track_cycles(track_cycles_test[17]); // 36444

    static mpz_t g, s, cs, u, A, au, au2, B, bu, buc, buca, uu, C;

    static bool is_init=false;
    if (!is_init) {
        mpz_inits(g, s, cs, u, A, au, au2, B, bu, buc, buca, uu, C, NULL);
        is_init=true;
    }

    //this assumes the gcd is 1, which it is when doing squaring; it does not assign the gcd output
    //this also assumes that the second input is >= the first input in magnitude, which is true
    //also assumes that the second input is nonnegative
    //todo
    simd_integer_namespace::integer_gcd_asm(g, s, nullptr, b, a); // s = modular inverse of b wrt. a
    //mpz_gcdext(g, s, nullptr, b, a);

    // b*u-c === x mod a
    // b*u === x+c mod a
    // b*((c*s)%a) === x+c mod a ; c*s = q*a + r ; 0<=r<a ; r = c*s-q*a
    // b*((q*a+r)%a) === x+c mod a
    // b*r === x+c mod a
    // b*(c*s - q*a) === x+c mod a
    // b*c*s - b*q*a === x+c mod a
    // b*c*s === x+c mod a
    // b*c*s-c === x mod a
    // (b*s-1)*c === x mod a
    // b*s-1 === y mod a
    // 1-1 === y mod a ; s is the modular inverse of b wrt. a, so b*s===1 mod a
    // 0 === y mod a
    // (b*s-1)*c === 0*c === 0 mod a
    // b*u-c === 0 mod a
    //remainder of (b*u-c)/a is 0
    //
    //   (b*u-c)/a
    // = (b*r - c)/a
    // = (b*(c*s-q*a) - c)/a
    // = (b*c*s - b*q*a - c)/a
    // = -b*q + (b*c*s - c)/a
    // = -b*q + ((b*s-1)*c)/a
    // = -b*q + c * (b*s-1)/a ; remainder is 0 so can remove c
    // b*s + a*t = 1
    // = -b*q + c * (1 - a*t - 1)/a
    // = -b*q + c * (- a*t)/a
    // = -(b*q + c*t)
    //didn't calculate the other cofactor for the gcd so this part isn't useful

    if (mpz_sgn(a)==0) {
        return false;
    }

    {
        track_cycles c_track_cycles(track_cycles_test[18]); // 600
        mpz_mul(cs, c, s); // cs = c*s
    }
    {
        track_cycles c_track_cycles(track_cycles_test[19]); // 1224
        mpz_mod(u, cs, a); // u = (c*s) % a
    }

    {
        track_cycles c_track_cycles(track_cycles_test[20]); // 376
        mpz_mul(A, a, a); // A = a*a
    }
    {
        track_cycles c_track_cycles(track_cycles_test[21]); // 534
        mpz_mul(au, a, u); // au = a*u
    }
    {
        track_cycles c_track_cycles(track_cycles_test[22]); // 88
        mpz_mul_2exp(au2, au, 1); // au2 = 2*a*u
    }
    {
        track_cycles c_track_cycles(track_cycles_test[23]); // 202
        mpz_sub(B, b, au2); // B = b - 2*a*u
    }

    {
        track_cycles c_track_cycles(track_cycles_test[24]); // 532
        mpz_mul(bu, b, u); // bu = b*u
    }
    {
        track_cycles c_track_cycles(track_cycles_test[25]); // 146
        mpz_sub(buc, bu, c); // buc = b*u-c
    }
    {
        track_cycles c_track_cycles(track_cycles_test[26]); // 642
        mpz_divexact(buca, buc, a); // buca = (b*u-c)/a
    }
    {
        track_cycles c_track_cycles(track_cycles_test[27]); // 360
        mpz_mul(uu, u, u); // uu = u*u
    }
    {
        track_cycles c_track_cycles(track_cycles_test[28]); // 122
        mpz_sub(C, uu, buca); // C = u*u - (b*u-c)/a
    }

    {
        track_cycles c_track_cycles(track_cycles_test[29]); // 40
        mpz_swap(a, A); // a <- A
        mpz_swap(b, B); // b <- B
        mpz_swap(c, C); // c <- C
    }

    return true;
}

void square_original(form& f) {
    vdf_original::form f_in;
    f_in.a[0]=f.a.impl[0];
    f_in.b[0]=f.b.impl[0];
    f_in.c[0]=f.c.impl[0];

    vdf_original::form& f_res=*vdf_original::square(f_in);

    mpz_set(f.a.impl, f_res.a);
    mpz_set(f.b.impl, f_res.b);
    mpz_set(f.c.impl, f_res.c);
}

void output_error(form start, int location) {
    print( "=== error ===" );
    print(start.a.to_string());
    print(start.b.to_string());
    print(start.c.to_string());
    print(location);
    assert(false);
}

bool square_fast(form& f) {
    track_cycles c_track_cycles(track_cycles_total);

    static form f_copy;
    if (test_correctness) {
        f_copy=f;
    }

    const int d_bits=2048;
    const int extra_d_bits=256; //calculated_d_bits is an upper bound so it is too high

    if (!square_fast(f.a, f.b, f.c)) {
        if (test_correctness) {
            output_error(f_copy, 0);
        }
        return false;
    }

    //todo reduce(f.a, f.b, f.c);
    if (!reduce_fast(f.a, f.b, f.c)) {
        if (test_correctness) {
            output_error(f_copy, 1);
        }
        return false;
    }

    // d=b^2-4ac
    // log2(d) ~ max(2*log2(|b|), 2+log2(a)+log2(c))+1
    // log2(x) < number of bits in x

    track_cycles c_track_cycles_2(track_cycles_test[30]);
    int a_bits=mpz_num_bits_upper_bound(f.a.impl);
    int b_bits=mpz_num_bits_upper_bound(f.b.impl);
    int c_bits=mpz_num_bits_upper_bound(f.c.impl);

    int calculated_d_bits=max(2*b_bits, 2+a_bits+c_bits)+1;

    //want to make sure the values aren't growing without bound due to a bug
    bool res=calculated_d_bits<d_bits+extra_d_bits;

    //a and c must be positive (if a or c are 0 then the discriminant is positive which is wrong)
    res&=(mpz_sgn(f.a.impl)==1);
    res&=(mpz_sgn(f.c.impl)==1);

    if (test_correctness) {
        if (!res) {
            output_error(f_copy, 2);
        }

        form f_copy_2=f_copy;
        square_original(f_copy_2);
        if (!(f==f_copy_2)) {
            output_error(f_copy, -1);
        }
    }

    return res;
}

struct repeated_square {
    integer d;

    int64 checkpoint_iteration=0;
    form checkpoint;

    int64 current_iteration=0;
    form current;

    int64 num_iterations=0;

    bool error_mode=false;

    bool is_checkpoint() {
        return
            current_iteration==num_iterations ||
            (current_iteration & (repeated_square_checkpoint_interval-1)) == 0
        ;
    }

    void advance_fast(bool& did_rollback) {
        bool is_error=false;

        if (!square_fast(current)) {
            is_error=true;
        }

        if (!is_error) {
            ++current_iteration;
            if (is_checkpoint() && !current.check_valid(d)) {
                is_error=true;
            }
        }

        if (is_error) {
            if (debug_rollback) {
                print( "Rollback", current_iteration, " -> ", checkpoint_iteration );
            }

            current_iteration=checkpoint_iteration;
            current=checkpoint;
            error_mode=true;
            did_rollback=true;
            assert(!assert_on_rollback);
        }
    }

    void advance_error() {
        square_original(current);
        ++current_iteration;
    }

    void advance() {
        bool did_rollback=false;
        if (error_mode) {
            advance_error();
        } else {
            advance_fast(did_rollback);
        }

        if (!did_rollback && is_checkpoint()) {
            checkpoint_iteration=current_iteration;
            checkpoint=current;
            error_mode=false;
        }
    }

    repeated_square(integer t_d, form initial, int64 t_num_iterations) {
        d=t_d;
        checkpoint=initial;
        current=initial;
        num_iterations=t_num_iterations;

        while (current_iteration<num_iterations) {
            advance();
        }
    }
};

int main(int argc, char* argv[]) {
    #if VDF_MODE!=0
        print( "=== Test mode ===" );
    #endif

    set_rounding_mode();

    simd_integer_namespace::init_asm_runtime();
    vdf_original::init();

    integer d(argv[1]);
    int64 num_iterations=from_string<int64>(argv[2]);

    repeated_square c_square(d, form::generator(d), num_iterations);

    cout << c_square.current.a.impl << "\n";
    cout << c_square.current.b.impl;

    if (enable_track_cycles) {
        print( "" );
        print( "" );

        simd_integer_namespace::track_cycles_reduce.output( "track_cycles_reduce" );
        simd_integer_namespace::track_cycles_gcd.output( "track_cycles_gcd" );
        track_cycles_total.output( "track_cycles_total" );

        for (int x=0;x<track_cycles_test.size();++x) {
            if (track_cycles_test[x].entries.empty()) {
                continue;
            }
            track_cycles_test[x].output(str( "track_cycles_test_#", x ));
        }
    }
}