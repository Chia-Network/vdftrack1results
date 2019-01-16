/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    Try to set G to the gcd of A and B using Brown's alogrithm M.
    This function switches to a big primes version if needed.
    It should only really fail if the dense size of the inputs is too large.
*/
int nmod_mpoly_gcd_brown(nmod_mpoly_t G,
                            const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    nmod_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    nmod_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (nmod_mpoly_is_zero(A, ctx)) {
        nmod_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (nmod_mpoly_is_zero(B, ctx)) {
        nmod_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    nmod_mpolyd_ctx_init(dctx, nvars);
    success = nmod_mpolyd_ctx_set_for_gcd(dctx, A, B, ctx);
    if (!success)
    {
        nmod_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    nmod_mpolyd_init(Ad, nvars);
    nmod_mpolyd_init(Bd, nvars);
    nmod_mpolyd_init(Gd, nvars);
    nmod_mpolyd_init(Abar, nvars);
    nmod_mpolyd_init(Bbar, nvars);

    nmod_mpoly_convert_to_nmod_mpolyd(Ad, dctx, A, ctx);
    nmod_mpoly_convert_to_nmod_mpolyd(Bd, dctx, B, ctx);
    success = nmod_mpolyd_gcd_brown_smprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
    if (!success)
    {
        nmod_mpoly_convert_to_nmod_mpolyd(Ad, dctx, A, ctx);
        nmod_mpoly_convert_to_nmod_mpolyd(Bd, dctx, B, ctx);
        success = nmod_mpolyd_gcd_brown_lgprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
        if (!success) {
            nmod_mpoly_zero(G, ctx);
        } else {
            nmod_mpoly_convert_from_nmod_mpolyd(G, ctx, Gd, dctx);
        }
    } else
    {
        nmod_mpoly_convert_from_nmod_mpolyd(G, ctx, Gd, dctx);
    }

    nmod_mpolyd_clear(Bbar);
    nmod_mpolyd_clear(Abar);
    nmod_mpolyd_clear(Gd);
    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Ad);

cleanup_stage1:

    nmod_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (!nmod_mpoly_is_zero(G, ctx))
        nmod_mpoly_make_monic(G, G, ctx);

    return success;
}
