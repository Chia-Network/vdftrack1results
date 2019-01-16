/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, const fmpz_poly_t g)
{
    slong i;

    fmpz_mod_poly_fit_length(f, g->length);
    _fmpz_mod_poly_set_length(f, g->length);

    for (i = 0; i < g->length; i++)
    {
        fmpz_mod(f->coeffs + i, g->coeffs + i, &(f->p));
    }

    _fmpz_mod_poly_normalise(f);
}
