/**
Copyright 2018 Chia Network Inc
Modifications copyright (C) 2019 Akashnil Dutta

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
**/

#include <iostream>
#include <gmpxx.h>
#include <cassert>
#include <cmath>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <thread>

#define LOG2(X) (63 - __builtin_clzll((X)))

using namespace std;

struct form {
    // y = ax^2 + bxy + y^2
    mpz_t a;
    mpz_t b;
    mpz_t c;

    mpz_t d; // discriminant
 };

ostream& operator<<(ostream& os, const form& f) {
    return os << "a: " <<  f.a << endl << "b: " << f.b << endl << "c: " << f.c << endl;
}

mpz_t g, d, e, q, a2, mu, denom;
mpz_t r, r_, r_old, r_old_, s, s_, s_old, s_old_, t, t_, t_old, t_old_, quo;
mpz_t faa, fab, fac, fba, fbb, fbc, fca, fcb, fcc;

form f_;

// int64_t total_time, gcd_time, reduce_time;

inline int64_t get_time() {
    auto now = std::chrono::system_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);

    auto value = now_ms.time_since_epoch();
    return value.count();
}

const int64_t THRESH = 1UL<<31;
const int64_t DOUBLE_THRESH = 1UL<<31;
const int64_t EXP_THRESH = 31;
const int64_t MAXV = ((1UL<<63) - 1);

inline void mpz_addmul_si(mpz_t rop, const mpz_t op1, signed long int op2) {
    if (op2 >= 0) mpz_addmul_ui(rop, op1, op2);
    else mpz_submul_ui(rop, op1, -op2);
}

inline uint64_t signed_shift(uint64_t op, int shift) {
    if (shift > 0) return op << shift;
    if (shift <= -64) return 0;
    return op >> (-shift);
}

inline int mpz_sign(const mpz_t op) {
    return mpz_sgn(op);
}

// Return an approximation x of the large mpz_t op by an int64_t and the exponent e adjustment.
// We must have (x * 2^e) / op = constant approximately.
inline int64_t mpz_get_si_2exp (signed long int *exp, const mpz_t op) {
    uint64_t size = mpz_size(op);
    uint64_t last = mpz_getlimbn(op, size - 1);
    uint64_t ret;
    int lg2 = LOG2(last) + 1;
    *exp = lg2; ret = signed_shift(last, 63 - *exp);
    if (size > 1) {
        *exp += (size-1) * 64;
        uint64_t prev = mpz_getlimbn(op, size - 2);
        ret += signed_shift(prev, -1 - lg2);
    }
    if (mpz_sgn(op) < 0) return - ((int64_t)ret);
    return ret;
}

inline void normalize(form& f) {
    mpz_add(mu, f.b, f.c);
    mpz_mul_ui(a2, f.c, 2);
    mpz_fdiv_q(denom, mu, a2);

    mpz_set(f_.a, f.c);

    mpz_mul_ui(a2, denom, 2);
    mpz_neg(f_.b, f.b);
    mpz_addmul(f_.b, f.c, a2);

    mpz_set(f_.c, f.a);
    mpz_submul(f_.c, f.b, denom);
    mpz_mul(denom, denom, denom);
    mpz_addmul(f_.c, f.c, denom);

    mpz_set(f.a, f_.a);
    mpz_set(f.b, f_.b);
    mpz_set(f.c, f_.c);
}

inline void slow_reduce(form& f) {
    while (true) {
        normalize(f);
        int cmp = mpz_cmp(f.a, f.c);
        if (cmp < 0 || (cmp == 0 && mpz_sgn(f.b) >= 0)) break;
    }
}

// Test if f is reduced. If it almost is but a, c are swapped, 
// then just swap them to make it reduced.
inline bool test_reduction(form& f) {
    int a_b = mpz_cmpabs(f.a, f.b);
    int c_b = mpz_cmpabs(f.c, f.b);

    if (a_b < 0 || c_b < 0) return false;

    int a_c = mpz_cmp(f.a, f.c);

    if (a_c > 0) {
        mpz_swap(f.a, f.c); mpz_neg(f.b, f.b);
    }

    if (a_c == 0 && mpz_sgn(f.b) < 0) {
        mpz_neg(f.b, f.b);
    }
    
    return true;
}

// This is a replacement of the original slow_reduce procudure. Please refer to the readme file
// for the algorithm.
inline void fast_reduce(form& f) {

    int64_t u, v, w, x, u_, v_, w_, x_;
    int64_t delta, gamma, sgn;
    int64_t a, b, c, a_, b_, c_;
    int64_t aa, ab, ac, ba, bb, bc, ca, cb, cc;
    signed long int a_exp, b_exp, c_exp, max_exp, min_exp;

    while (!test_reduction(f)) {
        
        a = mpz_get_si_2exp(&a_exp, f.a);
        b = mpz_get_si_2exp(&b_exp, f.b);
        c = mpz_get_si_2exp(&c_exp, f.c);

        max_exp = a_exp;
        min_exp = a_exp;

        if (max_exp < b_exp) max_exp = b_exp;
        if (min_exp > b_exp) min_exp = b_exp;

        if (max_exp < c_exp) max_exp = c_exp;
        if (min_exp > c_exp) min_exp = c_exp;

        if (max_exp - min_exp > EXP_THRESH) {
            normalize(f); continue;
        }
        max_exp++; // for safety vs overflow

        // Ensure a, b, c are shifted so that a : b : c ratios are same as f.a : f.b : f.c
        // a, b, c will be used as approximations to f.a, f.b, f.c
        a >>= (max_exp - a_exp);
        b >>= (max_exp - b_exp);
        c >>= (max_exp - c_exp);

        u_ = 1; v_ = 0; w_ = 0; x_ = 1;

        // We must be very careful about overflow in the following steps
        do {
            u = u_; v = v_; w = w_; x = x_;
            // Ensure that delta = floor ((b+c) / 2c)
            delta = b >= 0 ? (b+c) / (c<<1) : - (-b+c) / (c<<1);
            a_ = c;
            c_ = c * delta;
            b_ = -b + (c_ << 1);
            gamma = b - c_;
            c_ = a - delta * gamma;

            a = a_; b = b_; c = c_;

            u_ = v;
            v_ = -u + delta * v;
            w_ = x;
            x_ = -w + delta * x;
        // The condition (abs(v_) | abs(x_)) <= THRESH protects against overflow
        } while ((abs(v_) | abs(x_)) <= THRESH && a > c && c > 0);

        if ((abs(v_) | abs(x_)) <= THRESH) {
            u = u_; v = v_; w = w_; x = x_;
        }

        aa = u * u; ab = u * w; ac = w * w;
        ba = u * v << 1; bb = u * x + v * w; bc = w * x << 1;
        ca = v * v; cb = v * x; cc = x * x;

        // The following operations take 40% of the overall runtime.

        mpz_mul_si(faa, f.a, aa);
        mpz_mul_si(fab, f.b, ab);
        mpz_mul_si(fac, f.c, ac);

        mpz_mul_si(fba, f.a, ba);
        mpz_mul_si(fbb, f.b, bb);
        mpz_mul_si(fbc, f.c, bc);

        mpz_mul_si(fca, f.a, ca);
        mpz_mul_si(fcb, f.b, cb);
        mpz_mul_si(fcc, f.c, cc);

        mpz_add(f.a, faa, fab);
        mpz_add(f.a, f.a, fac);

        mpz_add(f.b, fba, fbb);
        mpz_add(f.b, f.b, fbc);

        mpz_add(f.c, fca, fcb);
        mpz_add(f.c, f.c, fcc);
    }
}

inline form generator_for_discriminant(mpz_t* d) {
    form x;
    mpz_init_set_ui(x.a, 2);
    mpz_init_set_ui(x.b, 1);
    mpz_init(x.c);
    mpz_init_set(x.d, *d);

    // c = b*b - d
    mpz_mul(x.c, x.b, x.b);
    mpz_sub(x.c, x.c, x.d);

    // denom = 4a
    mpz_mul_ui(denom, x.a, 4);

    mpz_fdiv_q(x.c, x.c, denom);

    fast_reduce(x);
    return x;
}

/**
 * This algorithm is the same as the composition/multiply algorithm,
 * but simplified to where both inputs are equal (squaring). It also
 * assumes that the discriminant is a negative prime. Algorithm:
 *
 * 1. solve for mu: b(mu) = c mod a
 * 2.  A = a^2
 *     B = B - 2a * mu
 *     C = mu^2 - (b * mu - c)/a
 * 3. reduce f(A, B, C)
 **/
inline void square(form &f) {
    //gcd_time -= get_time();
    mpz_gcdext(g, e, d, f.a, f.b);
    //gcd_time += get_time();
    mpz_mul(mu, f.c, d);
    mpz_fdiv_qr(q, mu, mu, f.a);
    mpz_mul(f.c, e, f.c);
    mpz_addmul(f.c, q, f.b);
    mpz_addmul(f.c, mu, mu);

    // New b
    mpz_mul_ui(a2, f.a, 2);
    mpz_submul(f.b, a2, mu);

    // New a
    mpz_mul(f.a, f.a, f.a);

    //reduce_time -= get_time();
    fast_reduce(f);
    //reduce_time += get_time();
}

// Performs the VDF squaring iterations
inline void repeated_square(form *f, uint64_t iterations) {
    // total_time -= get_time();
    for (uint64_t i=0; i < iterations; i++) {
        square(*f);
    }
    // total_time += get_time();
}

int main(int argc, char* argv[]) {

    mpz_inits(denom, g, d, e, q, a2, mu, f_.a, f_.b, f_.c, f_.d, NULL);
    mpz_inits(r, r_, r_old, r_old_, s, s_, s_old, s_old_, t, t_, t_old, t_old_, quo, NULL);
    mpz_inits(faa, fab, fac, fba, fbb, fbc, fca, fcb, fcc, NULL);

    mpz_t discriminant;
    int ret = mpz_init_set_str(discriminant, argv[1], 0);
    uint64_t iterations = stoi(argv[2]);
    assert(ret == 0);
    form x = generator_for_discriminant(&discriminant);
    repeated_square(&x, iterations);

    // Outputs a and b of final element
    cout << x.a << endl << x.b;

    // cout << "gcd_time: " << (float) gcd_time / (float) total_time << endl;
    // cout << "reduce_time: " << (float) reduce_time / (float) total_time << endl;
}
