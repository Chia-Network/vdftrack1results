/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "thread_pool.h"
#ifdef _OPENMP
#include <omp.h>
#endif

FLINT_TLS_PREFIX int _flint_num_threads = 1;
#pragma omp threadprivate(_flint_num_threads)

int flint_get_num_threads()
{
    return _flint_num_threads;
}

void flint_set_num_threads(int num_threads)
{
    _flint_num_threads = num_threads;
    if (global_thread_pool_initialized)
    {
        if (!thread_pool_set_size(global_thread_pool, num_threads - 1))
        {
            flint_throw(FLINT_ERROR,
               "flint_set_num_threads called while global thread pool in use");
        }
    }
    else
    {
        thread_pool_init(global_thread_pool, num_threads - 1);
        global_thread_pool_initialized = 1;
    }
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
}

void flint_parallel_cleanup()
{
    int needs_cleanup = 1;
#pragma omp master
    needs_cleanup = 0;

    if (needs_cleanup)
        flint_cleanup();
}
