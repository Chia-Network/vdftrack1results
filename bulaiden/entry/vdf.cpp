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

/**
compile:
g++ -O3 vdf.cpp -lgmpxx -lgmp

run:
time sh ./run.sh -0xdc2a335cd2b355c99d3d8d92850122b3d8fe20d0f5360e7aaaecb448960d57bcddfee12a229bbd8d370feda5a17466fc725158ebb78a2a7d37d0a226d89b54434db9c3be9a9bb6ba2c2cd079221d873a17933ceb81a37b0665b9b7e247e8df66bdd45eb15ada12326db01e26c861adf0233666c01dec92bbb547df7369aed3b1fbdff867cfc670511cc270964fbd98e5c55fbe0947ac2b9803acbfd935f3abb8d9be6f938aa4b4cc6203f53c928a979a2f18a1ff501b2587a93e95a428a107545e451f0ac6c7f520a7e99bf77336b1659a2cb3dd1b60e0c6fcfffc05f74cfa763a1d0af7de9994b6e35a9682c4543ae991b3a39839230ef84dae63e88d90f457 2097152

expected output:
41925311306490234810413138012137571066190117321086581650264103905212443247266638302300216176280962581197588248081055726676846585755174597863581644696718682501212285820220079934416655823371811162467870568376504352001314014329129123905108990627789940142020588151887088845810858316508373198920329374052520999187
34852345576682667643688406428836586142306150733999584872507440542733753870999851189854421199460417700803159204522889615033067098373432317909891523566449002340814898056506838421999772358720363985014768637377270846456207130003363936723439348271781377460453326857670664463211470390811155611448465292430048563927
**/

#include <iostream>
#include <gmpxx.h>
#include <cassert>

#include <flint/fmpz.h>
//#include "flint.h"
//#include <mpir.h>

#include <time.h>
#include <vector>

using namespace std;

struct form {
	// y = ax^2 + bxy + y^2
	mpz_t a;
	mpz_t b;
	mpz_t c;
	// mpz_t d; // discriminant
};

mpz_t RT;
mpz_t D, L, XQ, XR;
mpz_t negative_a, denom, r, old_b, ra, s, p, old_a;
mpz_t G, dx, dy, By, Dy, x, y, t1, t2, bx, by, ax, ay, q, t, Q1;
form F;

fmpz_t fy, fx, fby, fbx, fL;

ostream& operator<<(ostream& os, const form& f) {
	return os << "a: " <<  f.a << endl << "b: " << f.b << endl << "c: " << f.c << endl;
}

inline void normalize(form& f) {
    mpz_neg(r, f.a);
    if (mpz_cmp(f.b, r) > 0 && mpz_cmp(f.b, f.a) <= 0) {
        // Already normalized
        return;
    }
    mpz_sub(r, f.a, f.b);
    mpz_mul_2exp(ra, f.a, 1);
    mpz_fdiv_q(r, r, ra);
    mpz_mul(ra, r, f.a);
    mpz_addmul(f.c, ra, r);
    mpz_addmul(f.c, r, f.b);
    mpz_mul_2exp(ra, ra, 1);
    mpz_add(f.b, f.b, ra);
}

inline void reduce(form& f) {
    normalize(f);
    int cmp;
    while (((cmp = mpz_cmp(f.a, f.c)) > 0) ||
           (cmp == 0 && mpz_sgn(f.b) < 0)) {
        mpz_add(s, f.c, f.b);
        mpz_mul_2exp(p, f.c, 1);
        mpz_fdiv_q(s, s, p);
        mpz_set(old_a, f.a);
        mpz_set(old_b, f.b);
        mpz_set(f.a, f.c);
        mpz_neg(f.b, f.b);
        mpz_mul(p, s, f.c);
        mpz_mul_2exp(p, p, 1);
        mpz_add(f.b, f.b, p);
        mpz_mul(p, old_b, s);
        mpz_mul(s, s, s);
        mpz_mul(f.c, f.c, s);
        mpz_sub(f.c, f.c, p);
        mpz_add(f.c, f.c, old_a);
    }
    normalize(f);
}

long mpz_bits(mpz_t x)
{
	if (x->_mp_size == 0) return 0;
	return mpz_sizeinbase(x, 2);
}

void mpz_addmul_si(mpz_t r, mpz_t x, long u)
{
	if (u >= 0)
		mpz_addmul_ui(r, x, u);
	else
		mpz_submul_ui(r, x, -u);
}

// https://www.researchgate.net/publication/221451638_Computational_aspects_of_NUCOMP
void gmp_nudupl(form& f) {

	mpz_gcdext(G, y, NULL, f.b, f.a);

	mpz_divexact(By, f.a, G);
	mpz_divexact(Dy, f.b, G);

	mpz_mul(bx, y, f.c);
	mpz_mod(bx, bx, By);

	mpz_set(by, By);

	if (mpz_cmpabs(by, L) <= 0) {
		mpz_mul(dx, bx, Dy);
		mpz_sub(dx, dx, f.c);
		mpz_divexact(dx, dx, By);
		mpz_mul(f.a, by, by);
		mpz_mul(f.c, bx, bx);
		mpz_add(t, bx, by);
		mpz_mul(t, t, t);
		mpz_sub(f.b, f.b, t);
		mpz_add(f.b, f.b, f.a);
		mpz_add(f.b, f.b, f.c);
		mpz_mul(t, G, dx);
		mpz_sub(f.c, f.c, t);
		return;
	}

	fmpz_set_mpz(fy, y);
	fmpz_set_mpz(fx, x);
	fmpz_set_mpz(fby, by);
	fmpz_set_mpz(fbx, bx);
	fmpz_set_mpz(fL, L);

	fmpz_xgcd_partial(fy, fx, fby, fbx, fL);

	fmpz_get_mpz(y, fy);
	fmpz_get_mpz(x, fx);
	fmpz_get_mpz(by, fby);
	fmpz_get_mpz(bx, fbx);
	fmpz_get_mpz(L, fL);


	//gmp_xgcd_partial(y, x, by, bx, L);

	mpz_neg(x, x);
	if (mpz_sgn(x) > 0) {
		mpz_neg(y, y);
	} else {
		mpz_neg(by, by);
	}

	mpz_mul(ax, G, x);
	mpz_mul(ay, G, y);

	mpz_mul(t, Dy, bx);
	mpz_submul(t, f.c, x);
	mpz_divexact(dx, t, By);
	mpz_mul(Q1, y, dx);
	mpz_add(dy, Q1, Dy);
	mpz_add(f.b, dy, Q1);
	mpz_mul(f.b, f.b, G);
	mpz_divexact(dy, dy, x);
	mpz_mul(f.a, by, by);
	mpz_mul(f.c, bx, bx);
	mpz_add(t, bx, by);
	mpz_submul(f.b, t, t);
	mpz_add(f.b, f.b, f.a);
	mpz_add(f.b, f.b, f.c);
	mpz_submul(f.a, ay, dy);
	mpz_submul(f.c, ax, dx);
}

// DOESN'T NEED TO BE OPTIMIZED
inline void generator_for_discriminant(form& x, mpz_t& d) {
	mpz_t denom;
	mpz_init(denom);
	mpz_set_ui(x.a, 2);
	mpz_set_ui(x.b, 1);
	mpz_mul(x.c, x.b, x.b);
	mpz_sub(x.c, x.c, d);
	mpz_mul_ui(denom, x.a, 4);
	mpz_fdiv_q(x.c, x.c, denom);
	reduce(x);
	mpz_clear(denom);
}

int main(int argc, char* argv[]) {

	mpz_init(RT);
	mpz_inits(D, L, XQ, XR, NULL);
	mpz_inits(negative_a, denom, r, old_b, ra, s, p, old_a, NULL);
	mpz_inits(G, dx, dy, By, Dy, x, y, t1, t2, bx, by, ax, ay, q, t, Q1, NULL);
	mpz_inits(F.a, F.b, F.c, NULL);


	fmpz_init(fy);
	fmpz_init(fx);
	fmpz_init(fby);
	fmpz_init(fbx);
	fmpz_init(fL);

	mpz_set_str(D, argv[1], 0);
	generator_for_discriminant(F, D);
	uint64_t n = stoi(argv[2]);

	mpz_abs(L, D);
	mpz_root(L, L, 4);

	for (int i=0; i<n; i++) {
		gmp_nudupl(F);
		reduce(F);
	}

	cout << F.a << endl << F.b;
	cout << flush;
}
