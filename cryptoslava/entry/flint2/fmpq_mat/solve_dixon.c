/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

int
fmpq_mat_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    fmpz_mat_t Anum;
    fmpz_mat_t Bnum;
    fmpz_mat_t Xnum;
    fmpz_t mod;
    int success;

    fmpz_mat_init(Anum, A->r, A->c);
    fmpz_mat_init(Bnum, B->r, B->c);
    fmpz_mat_init(Xnum, B->r, B->c);
    fmpz_init(mod);

    fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
    success = fmpz_mat_solve_dixon(Xnum, mod, Anum, Bnum);
    if (success)
        success = fmpq_mat_set_fmpz_mat_mod_fmpz(X, Xnum, mod);

    fmpz_mat_clear(Anum);
    fmpz_mat_clear(Bnum);
    fmpz_mat_clear(Xnum);
    fmpz_clear(mod);

    return success;
}
