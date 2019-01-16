#include <iostream>
#include <cassert>
#include "gmp.h"
#include <chrono>

using namespace std;
//__________________________________________

struct form 
{
  // y = ax^2 + bxy + y^2
  mpz_t a;
  mpz_t b;
  mpz_t c;
  mpz_t d; // discriminant
};
//__________________________________________

ostream& operator<<(ostream& os, const form& f) 
{
  return os << "a: " <<  f.a << endl << "b: " << f.b << endl << "c: " << f.c << endl;
}
//__________________________________________

mpz_t denom, denom1, old_a, old_b, s, x, x1,  g, d, e, a, b, q, m, mu, negative_a, r, ra, ra1;
//__________________________________________

inline void normalize(form& f) {
    mpz_neg(negative_a, f.a);
    if (mpz_cmp(f.b, negative_a) > 0 && mpz_cmp(f.b, f.a) <= 0) {
        return;
    }
    
    mpz_sub(r, f.a, f.b); //r=a-b
    mpz_mul_ui(denom, f.a, 2); //denom=2*a
    mpz_fdiv_q(r, r, denom); // r = r / denom

    mpz_mul(ra, r, f.a); //ra=a*r

    //new C
    mpz_add(ra1, ra, f.b); //ra1=a*r+b
    mpz_mul(ra1,ra1,r); //ra1=(a*r+b)*r
    mpz_add(f.c, f.c, ra1); // c=c+(a*r+b)*r

    //new B
    mpz_add(f.b, f.b, ra); //b=b+a*r
    mpz_add(f.b, f.b, ra); //b=b+2*a*r
}
//__________________________________________

inline void reduce(form& f) 
{
  normalize(f);
  while ((mpz_cmp(f.a, f.c) > 0) || (mpz_cmp(f.a, f.c) == 0 && mpz_cmp_si(f.b, 0) < 0))
  {
    mpz_add(s, f.c, f.b); // s=b+c
    mpz_mul_ui(denom1, f.c, 2); // x=2*c
    mpz_fdiv_q(s, s, denom1); // s=s/x

    mpz_set(old_a, f.a); // old_a=a
    mpz_set(old_b, f.b); // old_b=b

    // new A
    mpz_set(f.a, f.c); //a=c
   
    // new B
    mpz_neg(f.b, f.b); // b = -b
    mpz_mul(x, s, f.c); // x=s*c
    mpz_add(f.b, f.b, x); // b=b+s*c
    mpz_add(f.b, f.b, x); // b=b+2*s*c

    // new C
    mpz_sub(f.c, x, old_b); //c=cs-b
    mpz_mul(f.c, f.c, s); //c=(cs-b)*s
    mpz_add(f.c, f.c, old_a); //c=(cs-b)*s+a
  }
  normalize(f);

}
//__________________________________________

inline form generator_for_discriminant(mpz_t* d) 
{
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

  reduce(x);
  return x;
}
//__________________________________________

// Returns mu, solving for x:  ax = b mod m
// such that x = u + vn (n are all integers). Assumes that mu and v are initialized.
// Faster version without check, and without returning v
inline void solve_linear_congruence(mpz_t& mu, mpz_t& a, mpz_t& b, mpz_t& m) 
{
  mpz_gcdext(g, d, e, a, m);
  mpz_fdiv_q(q, b, g);
  mpz_mul(mu, q, d);
  mpz_mod(mu, mu, m);
}
//__________________________________________

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
inline void square(form &f1) 
{
  solve_linear_congruence(mu, f1.b, f1.c, f1.a);
  mpz_mul(x1, f1.b, mu); // m=b*mu
  mpz_sub(x1, x1, f1.c); // m=b*mu-c
  mpz_fdiv_q(x1, x1, f1.a); // m=(b*mu-c)/a

  // New c
  mpz_mul(f1.c, mu, mu); //c=mu^2 
  mpz_sub(f1.c, f1.c, x1); //c=mu^2-m

  // New b
  mpz_mul(a, mu, f1.a); // a=a*mu
  mpz_sub(f1.b, f1.b, a); // b=b-a*mu
  mpz_sub(f1.b, f1.b, a); // b=b-2*a*mu

  // New a
  mpz_mul(f1.a, f1.a, f1.a); //a=a^2 
  reduce(f1);
}
//__________________________________________

int main(int argc, char* argv[])
{
  mpz_inits(denom, denom1, old_a, old_b, s, x, x1, negative_a, r, ra, ra1, a, g, d, e, b, q, m, mu, NULL);

  mpz_t discriminant;
  mpz_init_set_str(discriminant, argv[1], 0);
  uint64_t iterations = stoi(argv[2]);

  form x = generator_for_discriminant(&discriminant);
  for (uint64_t i=0; i < iterations; i++)  
    square(x);
 
  // Outputs a and b of final element
  cout << x.a << endl << x.b;
  return 0;
}
