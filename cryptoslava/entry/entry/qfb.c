#include "qfb.h"

#define LOG(msg, ...) (fprintf(stderr, msg "\n", ##__VA_ARGS__))

/* Find such g, x that a*x == g (mod n). Optimize for g == 1, so that x =
 * a^(-1) (mod n) . */
void fast_gcdinv(fmpz_t g, fmpz_t x, const fmpz_t a, const fmpz_t n)
{
	int ret = fmpz_invmod(x, a, n);
	if (ret) {
		fmpz_one(g);
		return;
	}

	fmpz_gcdinv(g, x, a, n);
}

void qfb_nudupl2(qfb_t r, qfb_t f, fmpz_t D, fmpz_t L)
{
	fmpz_t a1, c1, cb, k, s, t, u2, v1, v2;

	fmpz_init(a1); fmpz_init(c1);
	fmpz_init(cb);
	fmpz_init(k);
	fmpz_init(s);
	fmpz_init(t); fmpz_init(u2); fmpz_init(v1); fmpz_init(v2);

	/* nucomp calculation */

	/* a1 = a */
	fmpz_set(a1, f->a);
	/* c1 = c */
	fmpz_set(c1, f->c);

	fmpz_zero(k);

	/* b < 0 */
	if (fmpz_sgn(f->b) < 0) {
		fmpz_neg(f->b, f->b);
		/* s = gcd(abs(b), a); v2 = inv(b) (mod a) */
		fast_gcdinv(s, v2, f->b, a1);
		fmpz_neg(f->b, f->b);
		fmpz_neg(v2, v2);
	} else {
		fast_gcdinv(s, v2, f->b, a1);
	}

	fmpz_mul(k, v2, c1);
	fmpz_neg(k, k);

	if (!fmpz_is_one(s)) {
		fmpz_fdiv_q(a1, a1, s);
		fmpz_mul(c1, c1, s);
		LOG("s is not one!");
	}

	/* k = -(c*inv(b)) (mod a) */
	fmpz_fdiv_r(k, k, a1);

	if (fmpz_cmp(a1, L) < 0) {
		fmpz_mul(t, a1, k);

		fmpz_mul(r->a, a1, a1);

		fmpz_mul_2exp(cb, t, 1);
		fmpz_add(cb, cb, f->b);

		fmpz_add(r->c, f->b, t);
		fmpz_mul(r->c, r->c, k);
		fmpz_add(r->c, r->c, c1);

		fmpz_fdiv_q(r->c, r->c, a1);
	} else {
		fmpz_t m2, r1, r2, co1, co2, temp;

		fmpz_init(m2); fmpz_init(r1); fmpz_init(r2);
		fmpz_init(co1); fmpz_init(co2); fmpz_init(temp);

		fmpz_set(r2, a1);
		/* r1 = k */
		fmpz_swap(r1, k);

		/* Satisfies co2*r1 - co1*r2 == +/- r2_orig */
		fmpz_xgcd_partial(co2, co1, r2, r1, L);

		/* m2 = b * r1 */
		fmpz_mul(m2, f->b, r1);
		fmpz_submul(m2, c1, co1);
		/* m2 = (b*r1 - c1*co1) / a1 */
		fmpz_divexact(m2, m2, a1);

		/* new_a = r1^2 */
		fmpz_mul(r->a, r1, r1);
		/* new_a = new_a - co1 * m2 */
		fmpz_submul(r->a, co1, m2);
		if (fmpz_sgn(co1) >= 0)
			fmpz_neg(r->a, r->a);

		fmpz_mul(cb, r->a, co2);
		fmpz_submul(cb, a1, r1);
		/* cb = a1*r1 - new_a*co2 */
		fmpz_neg(cb, cb);
		/* cb = 2 * (a1*r1 - new_a*co2) */
		fmpz_mul_2exp(cb, cb, 1);
		fmpz_divexact(cb, cb, co1);
		fmpz_sub(cb, cb, f->b);
		fmpz_mul_2exp(temp, r->a, 1);
		fmpz_fdiv_r(cb, cb, temp);

		fmpz_mul(r->c, cb, cb);
		fmpz_sub(r->c, r->c, D);
		fmpz_divexact(r->c, r->c, r->a);
		fmpz_tdiv_q_2exp(r->c, r->c, 2);

		if (fmpz_sgn(r->a) < 0) {
			fmpz_neg(r->a, r->a);
			fmpz_neg(r->c, r->c);
		}

		fmpz_clear(m2); fmpz_clear(r1); fmpz_clear(r2);
		fmpz_clear(co1); fmpz_clear(co2); fmpz_clear(temp);
	}

	fmpz_set(r->b, cb);

	fmpz_clear(cb);
	fmpz_clear(k);
	fmpz_clear(s);
	fmpz_clear(t); fmpz_clear(u2); fmpz_clear(v2);
	fmpz_clear(a1); fmpz_clear(c1);
}
