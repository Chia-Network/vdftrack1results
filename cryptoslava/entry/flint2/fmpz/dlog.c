/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

double
fmpz_dlog(const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return log(*x);
    }
    else
    {
        double s;

#if defined(__MPIR_VERSION)
        slong e;
#else
        long e;
#endif        

        s = mpz_get_d_2exp(&e, COEFF_TO_PTR(*x));
       
        return log(s) + e * 0.69314718055994530942;
    }
}
