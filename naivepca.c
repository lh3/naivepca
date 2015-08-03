#include <zlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include "ksort.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

typedef struct {
	double ev;
	int i;
} evsrt_t;

#define evsrt_lt(a, b) ((a).ev > (b).ev)
KSORT_INIT(ev, evsrt_t, evsrt_lt)

int n_eigen_symm(double *_a, int n, double *eval);

int main(int argc, char *argv[])
{
	int c, dret, lineno = 0, n_rows = 0, m_rows = 0, n_cols = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int8_t **C = 0;
	double **M, *X;
	char **names = 0;

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (argc - optind == 0) {
		fprintf(stderr, "Usage: naivepca <in.txt>\n");
		return 1;
	}

	fp = gzopen(argv[optind], "r");
	ks = ks_init(fp);

	// read the matrix into C
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		int8_t *q;
		char *p, *name = str.s;
		int i;
		++lineno;
		for (p = str.s; *p && *p != '\t' && *p != ' '; ++p);
		if (*p) {
			*p++ = 0;
			for (; *p && (*p == '\t' || *p == ' '); ++p);
		}
		if (*p == 0) {
			fprintf(stderr, "WARNING: line %d has one field; skipped.\n", lineno);
			continue;
		}
		if (n_cols != 0) {
			if (n_cols != str.s + str.l - p) {
				fprintf(stderr, "WARNING: line %d has a different number of columns; skipped.\n", lineno);
				continue;
			}
		} else n_cols = str.s + str.l - p;
		if (n_rows == m_rows) {
			m_rows = m_rows? m_rows<<1 : 16;
			C = (int8_t**)realloc(C, m_rows * sizeof(int8_t*));
			names = (char**)realloc(names, m_rows * sizeof(char*));
		}
		names[n_rows] = strdup(name);
		q = C[n_rows++] = (int8_t*)calloc(n_cols, sizeof(double));
		for (i = 0; i < n_cols; ++i) {
			if (p[i] >= '0' && p[i] <= '9') q[i] = p[i] - '0';
			else q[i] = -1;
		}
	}
	fprintf(stderr, "MESSAGE: read %d samples and %d sites\n", n_rows, n_cols);

	{ // normalize the matrix into M
		int i, j, *sum, *cnt;
		double *mu;
		sum = (int*)calloc(n_cols, sizeof(int));
		cnt = (int*)calloc(n_cols, sizeof(int));
		mu = (double*)calloc(n_cols, sizeof(double));
		for (i = 0; i < n_rows; ++i) {
			int8_t *q = C[i];
			for (j = 0; j < n_cols; ++j)
				if (q[j] >= 0) sum[j] += q[j], ++cnt[j];
		}
		for (j = 0; j < n_cols; ++j)
			if (cnt[j] > 0) mu[j] = (double)sum[j] / cnt[j];
		M = (double**)calloc(n_rows, sizeof(double*));
		for (i = 0; i < n_rows; ++i) {
			int8_t *q = C[i];
			double *r;
			r = M[i] = (double*)calloc(n_cols, sizeof(double));
			for (j = 0; j < n_cols; ++j)
				r[j] = q[j] < 0? 0. : (q[j] - mu[j]) / sqrt(.5 * mu[j] * (1. - .5 * mu[j]));
		}
		free(sum); free(cnt); free(mu);
		for (i = 0; i < n_rows; ++i) free(C[i]);
		free(C);
	}

	{ // multiplication
		int i, j, k;
		X = (double*)calloc(n_rows * n_rows, sizeof(double));
		for (i = 0; i < n_rows; ++i) {
			double *zi = M[i];
			for (j = 0; j <= i; ++j) {
				double t = 0., *zj = M[j];
				for (k = 0; k < n_cols; ++k)
					t += zi[k] * zj[k];
				X[i*n_rows + j] = X[j*n_rows + i] = t / n_cols;
			}
		}
		for (i = 0; i < n_rows; ++i) free(M[i]);
		free(M);
	}

	{ // print eigan vectors
		double *ev;
		int i, j;
		evsrt_t *evsrt;
		ev = (double*)calloc(n_rows, sizeof(double));
		evsrt = (evsrt_t*)calloc(n_rows, sizeof(evsrt_t));
		n_eigen_symm(X, n_rows, ev);
		for (i = 0; i < n_rows; ++i)
			evsrt[i].ev = ev[i], evsrt[i].i = i;
		ks_introsort(ev, n_rows, evsrt);
		for (i = 0; i < n_rows; ++i) {
			printf("%s", names[i]);
			for (j = 0; j < n_rows; ++j)
				printf("\t%.6f", X[i*n_rows + evsrt[j].i]);
			putchar('\n');
			free(names[i]);
		}
		free(ev);
		free(X); free(names);
	}
	
	ks_destroy(ks);
	gzclose(fp);
	return 0;
}
