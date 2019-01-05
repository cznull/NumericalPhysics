#include "pch.h"
#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>

#define ERRD 1e-100
#define ERRF 1e-10

template <class T>
int hc(T *m, int n);

template <class T>
int givensqr(T *m, T *q, int n) {
	T a, b, c;
	int i, j, k;
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			a = m[i*n + i];
			b = m[j*n + i];
			c = a * a + b * b;
			if (c < ERRD) {
				continue;
			}
			c = 1.0 / sqrt(c);
			a *= c;
			b *= c;
			for (k = i; k < n; k++) {
				c = m[i*n + k] * a + m[j*n + k] * b;
				m[j*n + k] = m[j*n + k] * a - m[i*n + k] * b;
				m[i*n + k] = c;
			}
			for (k = 0; k < n; k++) {
				c = q[k*n + i] * a + q[k*n + j] * b;
				q[k*n + j] = q[k*n + j] * a - q[k*n + i] * b;
				q[k*n + i] = c;
			}
		}
	}
	return 0;
}

template <class T>
int givensqr_m(T *m, T *q, int n) {
	T a, b, c;
	int i, j, k;
	hc(q, n);
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			a = m[i*n + i];
			b = m[j*n + i];
			c = a * a + b * b;
			if (c < ERRD) {
				continue;
			}
			c = 1.0 / sqrt(c);
			a *= c;
			b *= c;
			for (k = i; k < n; k++) {
				c = m[i*n + k] * a + m[j*n + k] * b;
				m[j*n + k] = m[j*n + k] * a - m[i*n + k] * b;
				m[i*n + k] = c;
			}
			for (k = 0; k < n; k++) {
				c = q[i*n + k] * a + q[j*n + k] * b;
				q[j*n + k] = q[j*n + k] * a - q[i*n + k] * b;
				q[i*n + k] = c;
			}
		}
	}
	hc(q, n);
	return 0;
}

int givensqr_m(double *m, double *q, int n) {
	double a, b, c;
	__m256d a4, b4, x1, x2;
	int i, j, k;
	hc(q, n);
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			a = m[i*n + i];
			b = m[j*n + i];
			c = a * a + b * b;
			if (c < ERRD) {
				continue;
			}
			c = 1.0 / sqrt(c);
			a *= c;
			b *= c;
			a4 = _mm256_set1_pd(a);
			b4 = _mm256_set1_pd(b);
			for (k = i; k < (n - 3); k += 4) {
				x1 = _mm256_loadu_pd(m + i * n + k);
				x2 = _mm256_loadu_pd(m + j * n + k);
				_mm256_storeu_pd(m + j * n + k, _mm256_sub_pd(_mm256_mul_pd(x2, a4), _mm256_mul_pd(x1, b4)));
				_mm256_storeu_pd(m + i * n + k, _mm256_add_pd(_mm256_mul_pd(x1, a4), _mm256_mul_pd(x2, b4)));
			}
			for (; k < n; k++) {
				c = m[i*n + k] * a + m[j*n + k] * b;
				m[j*n + k] = m[j*n + k] * a - m[i*n + k] * b;
				m[i*n + k] = c;
			}
			for (k = 0; k < (n - 3); k += 4) {
				x1 = _mm256_loadu_pd(q + i * n + k);
				x2 = _mm256_loadu_pd(q + j * n + k);
				_mm256_storeu_pd(q + j * n + k, _mm256_sub_pd(_mm256_mul_pd(x2, a4), _mm256_mul_pd(x1, b4)));
				_mm256_storeu_pd(q + i * n + k, _mm256_add_pd(_mm256_mul_pd(x1, a4), _mm256_mul_pd(x2, b4)));
			}
			for (; k < n; k++) {
				c = q[i*n + k] * a + q[j*n + k] * b;
				q[j*n + k] = q[j*n + k] * a - q[i*n + k] * b;
				q[i*n + k] = c;
			}
		}
	}
	hc(q, n);
	return 0;
}

int givensqr_m(float *m, float *q, int n) {
	float a, b, c;
	__m256 a4, b4, x1, x2;
	int i, j, k;
	hc(q, n);
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			a = m[i*n + i];
			b = m[j*n + i];
			c = a * a + b * b;
			if (c < ERRF) {
				continue;
			}
			c = 1.0 / sqrt(c);
			a *= c;
			b *= c;
			a4 = _mm256_set1_ps(a);
			b4 = _mm256_set1_ps(b);
			for (k = i; k < (n - 7); k += 8) {
				x1 = _mm256_loadu_ps(m + i * n + k);
				x2 = _mm256_loadu_ps(m + j * n + k);
				_mm256_storeu_ps(m + j * n + k, _mm256_sub_ps(_mm256_mul_ps(x2, a4), _mm256_mul_ps(x1, b4)));
				_mm256_storeu_ps(m + i * n + k, _mm256_add_ps(_mm256_mul_ps(x1, a4), _mm256_mul_ps(x2, b4)));
			}
			for (; k < n; k++) {
				c = m[i*n + k] * a + m[j*n + k] * b;
				m[j*n + k] = m[j*n + k] * a - m[i*n + k] * b;
				m[i*n + k] = c;
			}
			for (k = 0; k < (n - 7); k += 8) {
				x1 = _mm256_loadu_ps(q + i * n + k);
				x2 = _mm256_loadu_ps(q + j * n + k);
				_mm256_storeu_ps(q + j * n + k, _mm256_sub_ps(_mm256_mul_ps(x2, a4), _mm256_mul_ps(x1, b4)));
				_mm256_storeu_ps(q + i * n + k, _mm256_add_ps(_mm256_mul_ps(x1, a4), _mm256_mul_ps(x2, b4)));
			}
			for (; k < n; k++) {
				c = q[i*n + k] * a + q[j*n + k] * b;
				q[j*n + k] = q[j*n + k] * a - q[i*n + k] * b;
				q[i*n + k] = c;
			}
		}
	}
	hc(q, n);
	return 0;
}

template <class T>
int hc(T *s, T *d, int n) {
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			d[i*n + j] = s[j*n + i];
		}
	}
	return 0;
}

template <class T>
int hc(T *m, int n) {
	int i, j;
	T x;
	for (i = 1; i < n; i++) {
		for (j = 0; j < i; j++) {
			x = m[i*n + j];
			m[i*n + j] = m[j*n + i];
			m[j*n + i] = x;
		}
	}
	return 0;
}

template <class T>
int identitymatrix(T *x, int n) {
	int i;
	memset(x, 0, n*n * sizeof(T));
	for (i = 0; i < n; i++) {
		x[i*n + i] = 1.0;
	}
	return 0;
}

template <class T>
int pv(T *v, int n) {
	int i;
	for (i = 0; i < n - 1; i++) {
		printf("%f,", v[i]);
	}
	printf("%f\n", v[i]);
	return 0;
}

template <class T>
int pm(T *m, int a, int b) {
	int i;
	for (i = 0; i < a; i++) {
		pv(m + b * i, b);
	}
	return 0;
}

template <class T>
int copym(T *d, T *s, int n) {
	memcpy(d, s, n*n * sizeof(T));
	return 0;
}

template <class T>
int mult(T *l, T *r, T *a, int m, int n, int o, int p) {
	int i, j, k;
	T t;
	for (i = 0; i < m; i++) {
		for (j = 0; j < p; j++) {
			t = 0;
			for (k = 0; k < n; k++) {
				t += l[i*n + k] * r[k*p + j];
			}
			a[i*p + j] = t;
		}
	}
	return 0;
}

typedef float current;

int main()
{
	int size, i;
	current *m, *q, *ws;

	if (scanf("%d", &size) == 1) {

		m = (current*)malloc(size*size * sizeof(current));
		q = (current*)malloc(size*size * sizeof(current));
		ws = (current*)malloc(size*size * sizeof(current));

		identitymatrix(q, size);
		for (i = 0; i < size*size; i++) {
			if (scanf("%f", m + i) != 1) {
				return 0;
			}
		}

		givensqr_m(m, q, size);
		printf("q:\n");
		pm(q, size, size);
		printf("r:\n");
		pm(m, size, size);
		mult(q, m, ws, size, size, size, size);
		//printf("qr:\n");
		//pm(ws, size, size);

		free(m);
		free(q);
	}
}
