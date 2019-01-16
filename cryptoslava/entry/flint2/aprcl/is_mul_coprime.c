/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

/*
    Returns 1 if gcd(x * y, n) == 1; otherwise returns 0.
*/
int
is_mul_coprime_ui_ui(ulong x, ulong y, const fmpz_t n)
{
    ulong a, rem, result;
    fmpz_t m;

    fmpz_init(m);

    rem = fmpz_tdiv_ui(n, x);
    a = n_gcd(x, rem);              /* a = gcd(x, n % x) */
    fmpz_tdiv_q_ui(m, n, a);        /* m = n / a */
    rem = fmpz_tdiv_ui(m, y);
    result = a * n_gcd(y, rem);     /* result = gcd(y, m % y) */

    fmpz_clear(m);

    /* q and r must be small, so q*r must get into ulong. */
    if (result == 1 || fmpz_equal_ui(n, result))
        return 1;
    return 0;
}


/*
    Returns 1 if gcd(x * y, n) == 1; otherwise returns 0.
*/
int
is_mul_coprime_ui_fmpz(ulong x, const fmpz_t y, const fmpz_t n)
{
    int is_coprime;
    ulong a, rem;
    fmpz_t m, result;

    fmpz_init(m);
    fmpz_init(result);

    rem = fmpz_tdiv_ui(n, x);
    a = n_gcd(x, rem);              /* a = gcd(x, n % x) */
    fmpz_tdiv_q_ui(m, n, a);        /* m = n / a */
    fmpz_fdiv_r(result, m, y);      /* result = m % y */
    fmpz_gcd(result, result, n);    /* result = gcd(y, m % y) */

    is_coprime = 0;
    if (fmpz_is_one(result) || fmpz_equal(n, result))
        is_coprime = 1;

    fmpz_clear(m);
    fmpz_clear(result);

    return is_coprime;
}
