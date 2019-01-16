/**
Copyright (C) 2018 Markku Pulkkinen

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

#ifndef VDFFORMGMP_H
#define VDFFORMGMP_H

#include <cassert>
#include <functional>
#include <iostream>
#include <x86intrin.h>

#include "gmp.h"

/**
 * @brief The ClassGroup data struct for VDF variables a, b, c and discriminant.
 * Optimal size because it fits into single entry of 64 byte wide cache line.
 */
struct alignas(64) ClassGroup {
  mpz_t a;
  mpz_t b;
  mpz_t c;
  mpz_t d;
};

/**
 * @brief ClassGroupContext struct - placeholder for variables
 * in classgroup arithmetic operations. Uses two cache
 * line entries, 128 bytes.
 */
struct alignas(64) ClassGroupContext {
  mpz_t a;
  mpz_t b;
  mpz_t c;
  mpz_t mu;
  /* ----- cache line boundary ----- */
  mpz_t m;
  mpz_t r;
  mpz_t s;
  mpz_t v;

  ClassGroupContext(uint32_t numBits = 4096) {
    mpz_init2(a, numBits);
    mpz_init2(b, numBits);
    mpz_init2(c, numBits);
    mpz_init2(mu, numBits);
    mpz_init2(m, numBits);
    mpz_init2(r, numBits);
    mpz_init2(s, numBits);
    mpz_init2(v, numBits);
  }

  ~ClassGroupContext() { mpz_clears(a, b, c, mu, m, r, s, v, NULL); }
};

/**
 * @brief The VdfFormGmp class that use GNU MP library for VDF.
 * This is also aligned to default cache line width. Uses three
 * cache line entries, 192 bytes.
 */
class alignas(64) VdfFormGmp {

public:
  /**
   * @brief VdfFormGmp - default constructor.
   */
  VdfFormGmp() {
    numBits = 4096;
    init();
  }

  ~VdfFormGmp() { clear(); }

  /**
   * @brief VdfFormGmp - copy constructor.
   */
  VdfFormGmp(const VdfFormGmp &rhs) {
    numBits = rhs.numBits;
    init();
    assign(rhs);
  }

  /**
   * @brief VdfFormGmp - generates classgroup values using given input params.
   * Sets also initial mem alloc size for mpz_t type variables
   */
  VdfFormGmp(uint64_t a_, uint64_t b_, const char *discriminator,
             uint16_t numBits_ = 4096) {
    numBits = numBits_;
    init();
    generate(a_, b_, discriminator);
  }

  /**
   * @brief operator = overload.
   */
  VdfFormGmp &operator=(const VdfFormGmp &rhs) {
    if (&rhs == this)
      return *this;

    numBits = rhs.numBits;
    clear();
    init();
    assign(rhs);
    return *this;
  }

  /**
   * @brief A, B, C, D - getters for classgroup values.
   */
  inline const mpz_t &A() { return cg.a; }
  inline const mpz_t &B() { return cg.b; }
  inline const mpz_t &C() { return cg.c; }
  inline const mpz_t &D() { return cg.d; }

  /**
   * @brief square() - does repeated squaring iterations cnt times.
   * @param cnt - iteration counter
   *
   * This algorithm is the same as the composition/multiply algorithm,
   * but simplified to where both inputs are equal (squaring). It also
   * assumes that the discriminant is a negative prime. Algorithm:
   *
   * 1. solve for mu: b(mu) = c mod a
   * 2. compose A = a^2, B = B - 2a * mu, C = mu^2 - (b * mu - c)/a
   * 3. reduce f(A, B, C)
   **/
  inline void square(const uint64_t &cnt) {

    // prefetch to ensure that data including functors is cached
    _mm_prefetch(this, _MM_HINT_T0);

    // create and init local static struct for fast access,
    // this should play well with cache, all variables needed
    // by algorithm are close to each other.
    static ClassGroupContext t;
    mpz_set(t.a, cg.a);
    mpz_set(t.b, cg.b);
    mpz_set(t.c, cg.c);

    unsigned long long i(cnt);
    while (i--) {
      // solve linear congruence
      /////////////////////////////////////////////////////////////////////////
      mpz_gcdext(t.r, t.s, nullptr, t.b, t.a);
      // mu = (c/r) * s)%a ( == ((c/r)%a * s%a)%a )
      if (mpz_cmp_ui(t.r, 1) == 0) {
        mpz_mul(t.m, t.c, t.s);
        mpz_mod(t.mu, t.m, t.a);
      } else {
        mpz_mdiv(t.mu, t.c, t.r);
        mpz_mul(t.m, t.mu, t.s);
        mpz_mod(t.mu, t.m, t.a);
      }
      // compose new A, B, C
      /////////////////////////////////////////////////////////////////////////
      // ((b * mu) - c)/a -> s
      mpz_mul(t.m, t.b, t.mu);
      mpz_sub(t.m, t.m, t.c);
      mpz_mdiv(t.s, t.m, t.a);
      // (mu * a) -> r
      mpz_mul(t.r, t.mu, t.a);

      // new b = b - (2*(mu * a))
      mpz_submul_ui(t.b, t.r, 2);
      // new c = (mu * mu) - s
      mpz_pow_ui(t.c, t.mu, 2);
      mpz_sub(t.c, t.c, t.s);
      // new a = a * a
      mpz_pow_ui(t.a, t.a, 2);

      // reduction
      /////////////////////////////////////////////////////////////////////////
      {
        normalize(t);
        // in-place inner loop for reduce(a bit faster than functor call):
        // Instrumentation shows that this loop is the most critical part
        // for performance. Almost 90% of time is spend here for long
        // iterations, so every GMP API function call counts, less is more.
        // Avoid multiply and divide, use left/right shifts if possible !
        int res;
        while (1) {
          // reduced: (a < c) || (a == c && b >= 0)
          // not reduced: a > c || (a == c && b < 0)
          // branch check, looks scary because of many brackets but compiler
          // can optimize it so that mpz_sgn is called only if necessary
          res = mpz_cmp(t.a, t.c);
          if (!((res > 0) || !((res != 0) || (mpz_sgn(t.b) >= 0))))
            break;

          // (c + b)/2c == (1 + (b/c))/2 -> s
          mpz_mdiv(t.r, t.b, t.c);
          mpz_add_ui(t.r, t.r, 1);
          mpz_div_2exp(t.s, t.r, 1);
          // cs -> m
          mpz_mul(t.m, t.c, t.s);
          // 2cs -> r
          mpz_mul_2exp(t.r, t.m, 1);
          // (cs - b) -> m
          mpz_sub(t.m, t.m, t.b);

          // new b = -b + 2cs
          mpz_sub(t.b, t.r, t.b);
          // new a = c, c = a
          mpz_swap(t.a, t.c);
          // new c = c + cs^2 - bs ( == c + (s * ( cs - b)))
          mpz_addmul(t.c, t.s, t.m);
        };
      }
    }
    // do last normalize
    normalize(t);

    // copy new values into class group
    mpz_set(cg.a, t.a);
    mpz_set(cg.b, t.b);
    mpz_set(cg.c, t.c);
  }

private:
  using cgcontext_fn_t = std::function<void(ClassGroupContext &)>;
  using cggenerate_fn_t = std::function<void(uint64_t, uint64_t, const char *)>;

  // Internal variables - most often used 1st for cache efficiency
  /////////////////////////////////////////////////////////////////////////////
  // functors for classgroup arithmetic
  cgcontext_fn_t normalize;
  cgcontext_fn_t reduce;
  /* ----- cache line boundary ------------------------------- */
  // classgroup data
  ClassGroup cg;
  /* ----- cache line boundary ------------------------------- */
  // functor for initial classgroup generation
  cggenerate_fn_t generate;
  // fill to align with cache line boundary
  uint8_t fill[30];
  // initial allocation size for mpz_t variables
  uint16_t numBits{4096};

  // Internal methods - clear, assign, init
  /////////////////////////////////////////////////////////////////////////////
  /**
   * @brief clear - frees mem allocated for classgroup data.
   */
  inline void clear() { mpz_clears(cg.a, cg.b, cg.c, cg.d, NULL); }

  /**
   * @brief assign - assigning new values for classgroup.
   */
  inline void assign(const VdfFormGmp &rhs) {
    mpz_set(cg.a, rhs.cg.a);
    mpz_set(cg.b, rhs.cg.b);
    mpz_set(cg.c, rhs.cg.c);
    mpz_set(cg.d, rhs.cg.d);
  }

  /**
   * @brief init - allocates mem for variables and initializes functors
   * for repeated squaring algorithm.
   */
  inline void init() {
    mpz_init2(cg.a, numBits);
    mpz_init2(cg.b, numBits);
    mpz_init2(cg.c, numBits);
    mpz_init2(cg.d, numBits);

    // lambda expressions to init generate, normalize and reduce functors
    // bound lambda is faster than inlined method

    generate = [&](uint64_t a_, uint64_t b_, const char *discriminator) {
      static ClassGroupContext t(numBits);
      mpz_set_str(cg.d, discriminator, 0);
      mpz_set_ui(t.a, a_);
      mpz_set_ui(t.b, b_);

      // c = (b * b - d)/4a
      mpz_pow_ui(t.m, t.b, 2);
      mpz_sub(t.m, t.m, cg.d);
      mpz_mul_2exp(t.r, t.a, 2);
      mpz_mdiv(t.c, t.m, t.r);

      // do initial reduction
      normalize(t);
      reduce(t);
      normalize(t);

      mpz_set(cg.a, t.a);
      mpz_set(cg.b, t.b);
      mpz_set(cg.c, t.c);
    };

    normalize = [](ClassGroupContext &t) {
      // branch check for normalization
      // normalized: b > -a && b <= a
      // not normalized: b <= -a || b > a
      if (!(mpz_cmp(t.b, t.a) > 0 || ([&t]() -> bool {
              mpz_neg(t.m, t.a);
              return (mpz_cmp(t.b, t.m) <= 0);
            }())))
        return;

      // (a - b)/2a == (1 - (b/a))/2 -> r
      mpz_mdiv(t.s, t.b, t.a);
      mpz_neg(t.s, t.s);
      mpz_add_ui(t.s, t.s, 1);
      mpz_div_2exp(t.r, t.s, 1);
      // ar -> m
      mpz_mul(t.m, t.r, t.a);
      // ar + b -> v
      mpz_add(t.v, t.m, t.b);

      // new a = a
      // new b = b + 2ar
      mpz_addmul_ui(t.b, t.m, 2);
      // new c = c + ar^2 + br (c == c + (ar + b)*r)
      mpz_addmul(t.c, t.v, t.r);
    };

    reduce = [](ClassGroupContext &t) {
      int res;
      while (1) {
        // reduced: (a < c) || (a == c && b >= 0)
        // not reduced: a > c || (a == c && b < 0)
        res = mpz_cmp(t.a, t.c);
        if (!((res > 0) || !((res != 0) || (mpz_sgn(t.b) >= 0))))
          break;

        // (c + b)/2c == (1 + (b/c))/2 -> s
        mpz_mdiv(t.r, t.b, t.c);
        mpz_add_ui(t.r, t.r, 1);
        mpz_div_2exp(t.s, t.r, 1);
        // cs -> m
        mpz_mul(t.m, t.c, t.s);
        // 2cs -> r
        mpz_mul_2exp(t.r, t.m, 1);
        // (cs - b) -> m
        mpz_sub(t.m, t.m, t.b);

        // new b = -b + 2cs
        mpz_sub(t.b, t.r, t.b);
        // new a = c, c = a
        mpz_swap(t.a, t.c);
        // new c = c + cs^2 - bs ( == c + (s * ( cs - b)))
        mpz_addmul(t.c, t.s, t.m);
      }
    };
  }
};
#endif // VDFFORMGMP_H
