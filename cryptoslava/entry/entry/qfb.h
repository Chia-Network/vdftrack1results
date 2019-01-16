#include <flint/fmpz.h>

typedef struct qfb
{
    fmpz_t a;
    fmpz_t b;
    fmpz_t c;
} qfb;

typedef qfb qfb_t[1];

static inline void qfb_init(qfb_t q)
{
   fmpz_init(q->a);
   fmpz_init(q->b);
   fmpz_init(q->c);
}

static inline void qfb_clear(qfb_t q)
{
   fmpz_clear(q->a);
   fmpz_clear(q->b);
   fmpz_clear(q->c);
}

static inline void qfb_set(qfb_t f, qfb_t g)
{
   fmpz_set(f->a, g->a);
   fmpz_set(f->b, g->b);
   fmpz_set(f->c, g->c);
}

/* qfb.c */
void qfb_nudupl2(qfb_t r, qfb_t f, fmpz_t D, fmpz_t L);

/* reduce.c */
void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D);
