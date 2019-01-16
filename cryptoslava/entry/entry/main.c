#include "qfb.h"
#include "params.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

#define LOG(msg, ...) (fprintf(stderr, msg "\n", ##__VA_ARGS__))
#define FATAL(msg, ...) (LOG(msg, ##__VA_ARGS__), exit(1))

static char verify_buf[1024];

static int verify_result(qfb_t f, size_t n)
{
	int ret = 0;
	size_t i;
	struct verify_params *params = NULL;
	for (i = 0; i < sizeof(vp) / sizeof(*vp); i++) {
		if (vp[i].n == n) {
			params = &vp[i];
			break;
		}
	}
	if (!params) {
		LOG("Cannot verify for n=%zu", n);
		return 0;
	}

	fmpz_get_str(verify_buf, 10, f->a);
	if (strcmp(verify_buf, params->a)) {
		LOG("Wrong a, expected:\n%s\n\nfound:\n%s\n",
				params->a, verify_buf);
		ret = 1;
	}
	fmpz_get_str(verify_buf, 10, f->b);
	if (strcmp(verify_buf, params->b)) {
		LOG("Wrong b, expected:\n%s\n\nfound:\n%s\n",
				params->b, verify_buf);
		ret = 1;
	}
	if (!ret)
		LOG("Verification successful!");
	return ret;
}

static void log_sizes(qfb_t f, size_t n)
{
	LOG("Sizes of coeffs at iter %zu: a=%d b=%d c=%d", n,
			(int) fmpz_bits(f->a), (int) fmpz_bits(f->b),
			(int) fmpz_bits(f->c));
}

static int do_squarings(qfb_t f, fmpz_t D, size_t n)
{
	size_t i;
	fmpz_t L;

	fmpz_init(L);
	/* L = -D */
	fmpz_neg(L, D);
	/* L = 4th_root(L) = 4th_root(|D|) */
	fmpz_root(L, L, 4);

	/*log_sizes(f, 0);*/
	for (i = 0; i < n; i++) {
		qfb_nudupl2(f, f, D, L);
		/*log_sizes(f, i + 1);*/
	}
	qfb_reduce(f, f, D);
	/*log_sizes(f, i + 1);*/
	fmpz_clear(L);
	return 0;
}

static int do_work(qfb_t f, const char *disc, size_t n)
{
	int ret = 0;
	fmpz_t D;

	fmpz_init(D);
	if (fmpz_set_str(D, disc, 0)) {
		ret = 1;
		goto out;
	}

	/*fmpz_init(t);*/

	fmpz_set_ui(f->a, 2);
	fmpz_set_ui(f->b, 1);
	/* c = -D */
	fmpz_neg(f->c, D);
	/* c = c + 1 = 1 - D */
	fmpz_add_ui(f->c, f->c, 1);
	/* c = c / 8 = (1 - D) / 8 = (b^2 - D) / (4a) */
	fmpz_divexact_ui(f->c, f->c, 8);

	if (!(ret = do_squarings(f, D, n))) {
		fmpz_print(f->a);
		puts("");
		fmpz_print(f->b);
		puts("");
	}


out:
	fmpz_clear(D);
	return ret;
}

static void ts_diff(struct timespec *d, struct timespec *a, struct timespec *b)
{
	d->tv_sec = a->tv_sec - b->tv_sec;
	d->tv_nsec = a->tv_nsec - b->tv_nsec;
	if (d->tv_nsec < 0) {
		d->tv_sec -= 1;
		d->tv_nsec += 1000000000;
	}
}

int main(int argc, char **argv)
{
	int ret, is_bench;
	size_t n_iter;
	char *disc = BENCH_DISC;
	qfb_t f;
	struct timespec start_time, end_time, diff;

	if (argc < 3)
		FATAL("2 args required");

	if (!(n_iter = strtoul(argv[2], NULL, 0)))
		FATAL("Bad n_iter");

	is_bench = !strcmp(argv[1], "bench");

	if (!is_bench) {
		disc = argv[1];
	}

	clock_gettime(CLOCK_MONOTONIC, &start_time);
	qfb_init(f);
	ret = do_work(f, disc, n_iter);
	clock_gettime(CLOCK_MONOTONIC, &end_time);

	ts_diff(&diff, &end_time, &start_time);
	LOG("Time: %ld.%09ld", diff.tv_sec, diff.tv_nsec);

	if (is_bench && !ret)
		ret = verify_result(f, n_iter);
	qfb_clear(f);
	return ret;
}
