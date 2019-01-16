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


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpir.h>

#define PARANOIA_CHECKS 0

typedef struct {
    // y = ax^2 + bxy + cy^2
    mpz_t a;
    mpz_t b;
    mpz_t c;
} form_t;

// globals
mpz_t discriminant;
form_t f;

// normalize "locals"
mpz_t negative_a, r, denom, old_b, ra;

// reduce "locals"
mpz_t s, x, old_a;

mpz_t L; // partial_extended_gcd global parameter
mpz_t etmp_q,etmp_r,etmp; // partial_extended_gcd "locals"

mpz_t c, u, v, v2, v3; // nudupl "locals"
mpz_t a2,b2,c2,d1;
mpz_t e,d,g;
mpz_t dtmp1;


void normalize(void)
{
    mpz_neg(negative_a, f.a);
    if (mpz_cmp(f.b, negative_a) > 0 && mpz_cmp(f.b, f.a) <= 0) {
        // Already normalized
        return;
    }
    // r = (a - b) / 2a
    // a = a
    // b = b + 2ra
    // c = ar^2 + br + c
    mpz_sub(r, f.a, f.b);

    //mpz_mul_ui(denom, f.a, 2);
    mpz_mul_2exp(denom, f.a, 1);

    // r = (a-b) / 2a
    mpz_fdiv_q(r, r, denom);

    mpz_set(old_b, f.b);

    mpz_mul(ra, r, f.a);
    mpz_add(f.b, f.b, ra);
    mpz_add(f.b, f.b, ra);

    // c += ar^2
    mpz_mul(ra, ra, r);
    mpz_add(f.c, f.c, ra);

    // c += rb
    mpz_set(ra, r);
    mpz_mul(ra, ra, old_b);
    mpz_add(f.c, f.c, ra);
}

void reduce(void)
{
    normalize();
    while ((mpz_cmp(f.a, f.c) > 0) ||
           (mpz_cmp(f.a, f.c) == 0 && mpz_cmp_si(f.b, 0) < 0)) {
        mpz_add(s, f.c, f.b);

        // x = 2c
        //mpz_mul_ui(x, f.c, 2);
        mpz_mul_2exp(x, f.c, 1);
        mpz_fdiv_q(s, s, x);

        mpz_set(old_a, f.a);
        mpz_set(old_b, f.b);

        // b = -b
        mpz_set(f.a, f.c);
        mpz_neg(f.b, f.b);

        // x = 2sc
        mpz_mul(x, s, f.c);
        //mpz_mul_ui(x, x, 2);
        mpz_mul_2exp(x, x, 1);


        // b += 2sc
        mpz_add(f.b, f.b, x);

        // c = cs^2
        mpz_mul(f.c, f.c, s);
        mpz_mul(f.c, f.c, s);

        // x = bs
        mpz_mul(x, old_b, s);

        // c -= bs
        mpz_sub(f.c, f.c, x);

        // c += a
        mpz_add(f.c, f.c, old_a);
    }
}

int partial_extended_gcd(mpz_t *a, mpz_t *b, mpz_t *u, mpz_t *v)
{
    if(mpz_cmpabs(*b,L) <= 0)
        return 0;

    int count=0;
    int run=1;
    mpz_set_ui(*u, 0);
    mpz_set_ui(*v, 1);

    while(1) {
        // etmp = u - v*(a/b)
        mpz_tdiv_qr(etmp_q, etmp_r, *a, *b);
        mpz_mul(etmp, *v, etmp_q);
        mpz_sub(etmp, *u, etmp);

        mpz_set(*a, *b);
        mpz_set(*b, etmp_r);

        mpz_set(*u, *v);
        mpz_set(*v, etmp);

        ++count;

        //NOTE: according to "A NOTE ON NUCOMP", Alfred J. Van Der Poorten
        //  It is beneficial to do one more step
        if(run==0)
            break;
        if(mpz_cmpabs(*b,L) <= 0)
            run=0;
    }
    return count;
}


// Shanks-Atkins NUDUPL
void nudupl(void) {
    // gcdext(g,u,v,a,b) --> g = au + bv
    mpz_gcdext(d1, u, dtmp1, f.b, f.a);
    if(PARANOIA_CHECKS && mpz_cmp_ui(d1, 1)) {
        fprintf(stderr,"ERR: GCD!=1\n");
    }

    // c =  - (u * f1.c) % f1.a
    mpz_mul(c,u,f.c);
    mpz_neg(c,c);
    mpz_mod(c,c,f.a);


    // if( |c| > |c-a| ) c -= a;
    mpz_sub(dtmp1,c,f.a);
    if (mpz_cmpabs(c,dtmp1) > 0)
        mpz_set(c, dtmp1);

    mpz_set(d,f.a);
    mpz_set(v3,c);
    int count = partial_extended_gcd(&d, &v3, &v, &v2);

    mpz_mul(a2,d,d);
    mpz_mul(c2,v3,v3);

    if(count==0){
        mpz_mul(g, v3, f.b);
        mpz_add(g, g, f.c);
        if(PARANOIA_CHECKS){
            mpz_tdiv_qr(g, dtmp1, g, d);
            if(mpz_cmp_ui(dtmp1,0)!=0){
                fprintf(stderr,"ERR: remainder1 in NUDUPL\n");
            }
        }
        else {
            mpz_divexact(g, g, d);
        }

        mpz_set(b2, f.b);
        mpz_set(v2, d1);
        mpz_set(f.a, a2);
    }  
    else {
        if(count&1){
            mpz_neg(v,v);
            mpz_neg(d,d);
        }

        mpz_mul(e, f.c, v);
        mpz_mul(dtmp1, f.b, d);
        mpz_add(e, e, dtmp1);
        if(PARANOIA_CHECKS){
            mpz_tdiv_qr(e, dtmp1, e, f.a);
            if(mpz_cmp_ui(dtmp1,0)!=0){
                fprintf(stderr,"ERR: remainder2 in NUDUPL\n");
            }
        }
        else {
            mpz_divexact(e, e, f.a);
        }


        mpz_mul(g,e,v2);
        mpz_sub(g,g,f.b);
        if(PARANOIA_CHECKS){
            mpz_tdiv_qr(g,dtmp1,g,v);
            if(mpz_cmp_ui(dtmp1,0)!=0){
                fprintf(stderr,"ERR: remainder3 in NUDUPL\n");
            }
        }
        else {
            mpz_divexact(g,g,v);
        }

        mpz_mul(dtmp1, v, g);
        mpz_mul(b2, e, v2);
        mpz_add(b2, b2, dtmp1);

        if(mpz_cmp_ui(d1,1) != 0){
            mpz_mul(b2,b2,d1);
            mpz_mul(v,v,d1);
            mpz_mul(v2,v2,d1);
        }

        mpz_mul(dtmp1, e, v);
        mpz_add(f.a, a2, dtmp1);
    }

    // out.b = (d+v3)^2 - a2 - c2 + b2
    mpz_add(dtmp1, d, v3);
    mpz_mul(dtmp1, dtmp1, dtmp1);
    mpz_sub(dtmp1, dtmp1, a2);
    mpz_sub(dtmp1, dtmp1, c2);
    mpz_add(f.b, dtmp1, b2);

    // out.c = g*v2 + c2
    mpz_mul(dtmp1, g, v2);
    mpz_add(f.c, dtmp1, c2);

    // finish reduction
    reduce();
}

void initialize_globals()
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
    mpz_init2(u, 4096);
    mpz_init2(v, 4096);

    mpz_init2(L, 4096);
    mpz_init2(etmp_q, 4096);
    mpz_init2(etmp_r, 4096);
    mpz_init2(etmp, 4096);

    mpz_init2(c, 4096);
    mpz_init2(v2, 4096);
    mpz_init2(v3, 4096);
    mpz_init2(a2, 4096);
    mpz_init2(b2, 4096);
    mpz_init2(c2, 4096);
    mpz_init2(d1, 4096);
    mpz_init2(g, 4096);
    mpz_init2(d, 4096);
    mpz_init2(e, 4096);
    mpz_init2(dtmp1, 4096);


    // initialize "generator"
    mpz_init_set_ui(f.a, 2);
    mpz_init_set_ui(f.b, 1);
    // f.c = (1-d)/8
    mpz_init(f.c);
    mpz_ui_sub(f.c,1,discriminant);
    mpz_div_2exp(f.c, f.c, 3);

    // L = ((1-d)/4)^(1/4) + 1
    mpz_ui_sub(L,1,discriminant);
    mpz_div_2exp(L,L,2);
    mpz_root(L,L,4);
    mpz_add_ui(L,L,1);
}


int main(int argc, char* argv[]) {

    int ret = mpz_init_set_str(discriminant, argv[1], 0);
    if(PARANOIA_CHECKS && ret){
        fprintf(stderr,"ERR: read discriminant failed\n");
        _exit(1);
    }
    uint64_t iterations = atoll(argv[2]);

    initialize_globals();

    // Performs the VDF squaring iterations
    for (uint64_t i=0; i < iterations; i++) {
        nudupl();
    }

    // Output a and b of resulting form
    mpz_out_str(stdout, 10, f.a);
    fputc('\n',stdout);
    mpz_out_str(stdout, 10, f.b);
    fflush(stdout);

    _exit(0);
}

