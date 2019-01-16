/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "profiler.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("reconstruct_fmpz_2....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        int result;
        int modresult;
        int special_case;
        fmpq_t x, y;
        fmpz_t mod, res, N, D, t;
        mpz_t tmp;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(mod);
        fmpz_init(res);
        fmpz_init(N);
        fmpz_init(D);
        fmpz_init(t);
        mpz_init(tmp);

        fmpq_randtest(x, state, 100);

        fmpz_abs(N, fmpq_numref(x));
        fmpz_set(D, fmpq_denref(x));

        /* Randomly generate larger bounds */
        if (n_randint(state, 2))
        {
            fmpz_randtest_not_zero(t, state, 100);
            fmpz_abs(t, t);
            fmpz_mul(N, N, t);
        }
        if (n_randint(state, 2))
        {
            fmpz_randtest_not_zero(t, state, 100);
            fmpz_abs(t, t);
            fmpz_mul(D, D, t);
        }

        fmpz_mul(mod, N, D);
        fmpz_mul_ui(mod, mod, UWORD(2));
        /* Next prime greater than or equal */
        fmpz_get_mpz(tmp, mod);
        flint_mpz_sub_ui(tmp, tmp, UWORD(1));
        mpz_nextprime(tmp, tmp);
        fmpz_set_mpz(mod, tmp);

        modresult = fmpq_mod_fmpz(res, x, mod);
        result = fmpq_reconstruct_fmpz_2(y, res, mod, N, D);

        /* Special case: both 1 and -1 have residue 1 mod 2.
           There's probably no particular reason to disallow this. */
        special_case = (fmpz_cmp_ui(mod, UWORD(2)) == 0 &&
                        fmpz_get_si(&x->num) == WORD(-1) &&
                        fmpz_cmp_ui(&x->den, UWORD(1)) == 0);

        if (special_case)
        {
            if (!modresult || !result ||
                !fmpz_is_one(&y->num) || !fmpz_is_one(&y->den))
            {
                flint_printf("FAIL: special case: -1 mod 2\n");
                abort();
            }
        }
        else if (!modresult || !result || !fmpq_equal(x, y))
        {
            flint_printf("FAIL: reconstruction failed\n");
            flint_printf("input = ");
            fmpq_print(x);
            flint_printf("\nmodulus = ");
            fmpz_print(mod);
            flint_printf("\nresidue = ");
            fmpz_print(res);
            flint_printf("\nreconstructed = ");
            fmpq_print(y);
            flint_printf("\nfmpq_mod_fmpz return value = %d", modresult);
            flint_printf("\nfmpq_reconstruct_fmpz return value = %d", result);
            flint_printf("\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(mod);
        fmpz_clear(res);
        fmpz_clear(N);
        fmpz_clear(D);
        fmpz_clear(t);
        mpz_clear(tmp);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
