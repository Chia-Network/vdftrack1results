/**
Copyright 2018 Chia Network Inc

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

#include "mpir.h"

typedef struct form {
    // y = ax^2 + bxy + y^2
    mpz_t a;
    mpz_t b;
    mpz_t c;

    mpz_t d; // discriminant
} form;

mpz_t discriminant;
mpz_t negative_a, r, denom, old_b, ra, s, x, old_a, g, d, e, q, w, u, a, b, m, k, mu, v, sigma, lambda, h, t, l, j;
form f3;

inline void normalize(form *f) {
    mpz_neg(negative_a, f->a);
    
    if (mpz_cmp(f->b, negative_a) > 0 && mpz_cmp(f->b, f->a) <= 0) {
        // Already normalized
        return;
    }
    
    // r = (a - b) / 2a
    // a = a
    // b = b + 2ra
    // c = ar^2 + br + c
    mpz_sub(r, f->a, f->b);

    //mpz_mul_ui(denom, f->a, 2);
    mpz_mul_2exp(denom, f->a, 1);

    // r = (a-b) / 2a
    mpz_fdiv_q(r, r, denom);

    mpz_set(old_b, f->b);

    mpz_mul(ra, r, f->a);
    mpz_add(f->b, f->b, ra);
    mpz_add(f->b, f->b, ra);

    // c += ar^2
    mpz_mul(ra, ra, r);
    mpz_add(f->c, f->c, ra);

    // c += rb
    mpz_set(ra, r);
    mpz_mul(ra, ra, old_b);
    mpz_add(f->c, f->c, ra);
}

inline void reduce(form *f) {
    normalize(f);
    while ((mpz_cmp(f->a, f->c) > 0) ||
           (mpz_cmp(f->a, f->c) == 0 && mpz_cmp_si(f->b, 0) < 0)) {
        mpz_add(s, f->c, f->b);

        // x = 2c
        //mpz_mul_ui(x, f->c, 2);
        mpz_mul_2exp(x, f->c, 1);
        mpz_fdiv_q(s, s, x);

        mpz_set(old_a, f->a);
        mpz_set(old_b, f->b);

        // b = -b
        mpz_set(f->a, f->c);
        mpz_neg(f->b, f->b);

        // x = 2sc
        mpz_mul(x, s, f->c);
        //mpz_mul_ui(x, x, 2);
        mpz_mul_2exp(x, x, 1);


        // b += 2sc
        mpz_add(f->b, f->b, x);

        // c = cs^2
        mpz_mul(f->c, f->c, s);
        mpz_mul(f->c, f->c, s);

        // x = bs
        mpz_mul(x, old_b, s);

        // c -= bs
        mpz_sub(f->c, f->c, x);

        // c += a
        mpz_add(f->c, f->c, old_a);
    }
    normalize(f);
}


form xyz;
inline form *generator_for_discriminant(mpz_t* d) {

    mpz_init_set_ui(xyz.a, 2);
    mpz_init_set_ui(xyz.b, 1);
    mpz_init(xyz.c);
    mpz_init_set(xyz.d, *d);

    // c = b*b - d
    mpz_mul(xyz.c, xyz.b, xyz.b);
    mpz_sub(xyz.c, xyz.c, xyz.d);

    // denom = 4a
    //mpz_mul_ui(denom, xyz.a, 4);
    mpz_mul_2exp(denom, xyz.a, 2);

    mpz_fdiv_q(xyz.c, xyz.c, denom);
    //reduce(x); //already reduced
    return &xyz;
}

// Returns mu and v, solving for x:  ax = b mod m
// such that x = u + vn (n are all integers). Assumes that mu and v are initialized.
// Returns 0 on success, -1 on failure

// Faster version without check, and without returning v
inline int solve_linear_congruence(mpz_t *mu, mpz_t *a, mpz_t *b, mpz_t *m) {
    mpz_gcdext(g, d, e, *a, *m);
    mpz_fdiv_q(q, *b, g);
    mpz_mul(*mu, q, d);
    mpz_mod(*mu, *mu, *m);
    return 0;
}

// Takes the gcd of three numbers
inline void three_gcd(mpz_t *ret, mpz_t *a, mpz_t *b, mpz_t *c) {
    mpz_gcd(*ret, *a, *b);
    mpz_gcd(*ret, *ret, *c);
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
inline form* square(form *f1) { // was &
    int ret = solve_linear_congruence(&mu, &f1->b, &f1->c, &f1->a);
    //assert(ret == 0);

    mpz_mul(m, f1->b, mu);
    mpz_sub(m, m, f1->c);
    mpz_fdiv_q(m, m, f1->a);

    // New a
    mpz_set(old_a, f1->a);
    mpz_mul(f3.a, f1->a, f1->a);

    // New b
    mpz_mul(a, mu, old_a);
    mpz_mul_ui(a, a, 2);
    mpz_sub(f3.b, f1->b, a);

    // New c
    mpz_mul(f3.c, mu, mu);
    mpz_sub(f3.c, f3.c, m);
    mpz_set(f3.d, f1->d);
    reduce(&f3);
    return &f3;
}

// Performs the VDF squaring iterations
inline form *repeated_square(form *f, uint64_t iterations) {
    for (uint64_t i=0; i < iterations; i++) {
        f = square(f);
    }
    return f;
}

int main(int argc, char* argv[])
{
    //mpz_inits(negative_a, r, denom, old_a, old_b, ra, s, x, g, d, e, q, w, m,
    //          u, a, b, k, mu, v, sigma, lambda, f3.a, f3.b, f3.c, f3.d,
    //          NULL);

    mpz_init2(negative_a, 4096);
    mpz_init2(r, 4096);
    mpz_init2(denom, 4096);
    mpz_init2(old_a, 4096);
    mpz_init2(old_b, 4096);
    mpz_init2(ra, 4096);
    mpz_init2(s, 4096);
    mpz_init2(x, 4096);
    mpz_init2(g, 4096);
    mpz_init2(d, 4096);
    mpz_init2(e, 4096);
    mpz_init2(q, 4096);
    mpz_init2(w, 4096);
    mpz_init2(m, 4096);
    mpz_init2(u, 4096);
    mpz_init2(a, 4096);
    mpz_init2(b, 4096);
    mpz_init2(k, 4096);
    mpz_init2(mu, 4096);
    mpz_init2(v, 4096);
    mpz_init2(sigma, 4096);
    mpz_init2(lambda, 4096);
    mpz_init2(f3.a, 4096);
    mpz_init2(f3.b, 4096);
    mpz_init2(f3.c, 4096);
    mpz_init2(f3.d, 4096);


    int ret = mpz_init_set_str(discriminant, argv[1], 0);
    uint64_t iterations = atoll(argv[2]);
    
    //assert(ret == 0);

    //form *x = generator_for_discriminant(&discriminant);

    mpz_init_set_ui(xyz.a, 2);
    mpz_init_set_ui(xyz.b, 1);
    mpz_init(xyz.c);
    mpz_init_set(xyz.d, discriminant);

    // c = b*b - d
    mpz_mul(xyz.c, xyz.b, xyz.b);
    mpz_sub(xyz.c, xyz.c, xyz.d);

    // denom = 4a
    //mpz_mul_ui(denom, xyz.a, 4);
    mpz_mul_2exp(denom, xyz.a, 2);

    mpz_fdiv_q(xyz.c, xyz.c, denom);


    form *y = repeated_square(&xyz, iterations);

    gmp_printf("\n\n%Zd\n%Zd\n\n", y->a, y->b);
}
