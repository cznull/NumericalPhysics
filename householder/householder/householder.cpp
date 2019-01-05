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
int householderqr(T *m, T *ws, T *q, int n) {//ws:n
	T s, a, aq;
	int i, j, k;
	for (i = 0; i < n - 1; i++) {
		s = 0;
		for (j = i; j < n; j++) {
			s += m[j*n + i] * m[j*n + i];
		}
		ws[i] = m[i*n + i] - sqrt(s);
		for (j = i + 1; j < n; j++) {
			ws[j] = m[j*n + i];
		}
		s = s - m[i*n + i] * m[i*n + i] + ws[i] * ws[i];
		if (s < ERRD) {
			continue;
		}
		for (j = i; j < n; j++) {
			a = 0;
			for (k = i; k < n; k++) {
				a += ws[k] * m[k*n + j];
			}
			a *= (2.0 / s);
			for (k = i; k < n; k++) {
				m[k*n + j] -= ws[k] * a;
			}
		}
		for (j = 0; j < n; j++) {
			a = 0;
			for (k = i; k < n; k++) {
				a += ws[k] * q[j*n + k];
			}
			a *= (2.0 / s);
			for (k = i; k < n; k++) {
				q[j*n + k] -= ws[k] * a;
			}
		}
	}
	return 0;
}

template <class T>
int householderqr_m(T *m, T *ws, T *q, int n) {//ws:n
	T s, a, aq;
	int i, j, k;
	hc(m, n);
	for (i = 0; i < n - 1; i++) {
		s = 0;
		for (j = i; j < n; j++) {
			s += m[i*n + j] * m[i*n + j];
		}
		ws[i] = m[i*n + i] - sqrt(s);
		for (j = i + 1; j < n; j++) {
			ws[j] = m[i*n + j];
		}
		s = s - m[i*n + i] * m[i*n + i] + ws[i] * ws[i];
		if (s < ERRD) {
			continue;
		}
		for (j = i; j < n; j++) {
			a = 0;
			for (k = i; k < n; k++) {
				a += ws[k] * m[j*n + k];
			}
			a *= (2.0 / s);
			for (k = i; k < n; k++) {
				m[j*n + k] -= ws[k] * a;
			}
		}
		for (j = 0; j < n; j++) {
			a = 0;
			for (k = i; k < n; k++) {
				a += ws[k] * q[j*n + k];
			}
			a *= (2.0 / s);
			for (k = i; k < n; k++) {
				q[j*n + k] -= ws[k] * a;
			}
		}
	}
	hc(m, n);
	return 0;
}

int householderqr_m(double *m, double *ws, double *q, int n) {//ws:n
	double s, a, aq;
	__m256d a4;
	int i, j, k;
	hc(m, n);
	for (i = 0; i < n - 1; i++) {
		s = 0;
		for (j = i; j < n; j++) {
			s += m[i*n + j] * m[i*n + j];
		}
		ws[i] = m[i*n + i] - sqrt(s);
		for (j = i + 1; j < n; j++) {
			ws[j] = m[i*n + j];
		}
		s = s - m[i*n + i] * m[i*n + i] + ws[i] * ws[i];
		if (s < ERRD) {
			continue;
		}
		for (j = i; j < n; j++) {
			a4 = _mm256_set1_pd(0.0);
			for (k = i; k < (n - 3); k += 4) {
				a4 = _mm256_add_pd(a4, _mm256_mul_pd(_mm256_loadu_pd(ws + k), _mm256_loadu_pd(m + j * n + k)));
			}
			a = a4.m256d_f64[0] + a4.m256d_f64[1] + a4.m256d_f64[2] + a4.m256d_f64[3];
			for (; k < n; k++) {
				a += ws[k] * m[j*n + k];
			}
			a *= (-2.0 / s);
			a4 = _mm256_set1_pd(a);
			for (k = i; k < (n - 3); k += 4) {
				_mm256_storeu_pd(m + j * n + k, _mm256_add_pd(_mm256_loadu_pd(m + j * n + k), _mm256_mul_pd(_mm256_loadu_pd(ws + k), a4)));
			}
			for (; k < n; k++) {
				m[j*n + k] += ws[k] * a;
			}
		}
		for (j = 0; j < n; j++) {
			a4 = _mm256_set1_pd(0.0);
			for (k = i; k < (n - 3); k += 4) {
				a4 = _mm256_add_pd(a4, _mm256_mul_pd(_mm256_loadu_pd(ws + k), _mm256_loadu_pd(q + j * n + k)));
			}
			a = a4.m256d_f64[0] + a4.m256d_f64[1] + a4.m256d_f64[2] + a4.m256d_f64[3];
			for (; k < n; k++) {
				a += ws[k] * q[j*n + k];
			}
			a *= (-2.0 / s);
			a4 = _mm256_set1_pd(a);
			for (k = i; k < (n - 3); k += 4) {
				_mm256_storeu_pd(q + j * n + k, _mm256_add_pd(_mm256_loadu_pd(q + j * n + k), _mm256_mul_pd(_mm256_loadu_pd(ws + k), a4)));
			}
			for (; k < n; k++) {
				q[j*n + k] += ws[k] * a;
			}
		}
	}
	hc(m, n);
	return 0;
}

int householderqr_m(float *m, float *ws, float *q, int n) {//ws:n
	float s, a, aq;
	__m256 a4;
	int i, j, k;
	hc(m, n);
	for (i = 0; i < n - 1; i++) {
		s = 0;
		for (j = i; j < n; j++) {
			s += m[i*n + j] * m[i*n + j];
		}
		ws[i] = m[i*n + i] - sqrt(s);
		for (j = i + 1; j < n; j++) {
			ws[j] = m[i*n + j];
		}
		s = s - m[i*n + i] * m[i*n + i] + ws[i] * ws[i];
		if (s < ERRF) {
			continue;
		}
		for (j = i; j < n; j++) {
			a4 = _mm256_set1_ps(0.0);
			for (k = i; k < (n - 7); k += 8) {
				a4 = _mm256_add_ps(a4, _mm256_mul_ps(_mm256_loadu_ps(ws + k), _mm256_loadu_ps(m + j * n + k)));
			}
			a = a4.m256_f32[0] + a4.m256_f32[1] + a4.m256_f32[2] + a4.m256_f32[3] + a4.m256_f32[4] + a4.m256_f32[5] + a4.m256_f32[6] + a4.m256_f32[7];
			for (; k < n; k++) {
				a += ws[k] * m[j*n + k];
			}
			a *= (-2.0 / s);
			a4 = _mm256_set1_ps(a);
			for (k = i; k < (n - 7); k += 8) {
				_mm256_storeu_ps(m + j * n + k, _mm256_add_ps(_mm256_loadu_ps(m + j * n + k), _mm256_mul_ps(_mm256_loadu_ps(ws + k), a4)));
			}
			for (; k < n; k++) {
				m[j*n + k] += ws[k] * a;
			}
		}
		for (j = 0; j < n; j++) {
			a4 = _mm256_set1_ps(0.0);
			for (k = i; k < (n - 7); k += 8) {
				a4 = _mm256_add_ps(a4, _mm256_mul_ps(_mm256_loadu_ps(ws + k), _mm256_loadu_ps(q + j * n + k)));
			}
			a = a4.m256_f32[0] + a4.m256_f32[1] + a4.m256_f32[2] + a4.m256_f32[3] + a4.m256_f32[4] + a4.m256_f32[5] + a4.m256_f32[6] + a4.m256_f32[7];
			for (; k < n; k++) {
				a += ws[k] * q[j*n + k];
			}
			a *= (-2.0 / s);
			a4 = _mm256_set1_ps(a);
			for (k = i; k < (n - 7); k += 8) {
				_mm256_storeu_ps(q + j * n + k, _mm256_add_ps(_mm256_loadu_ps(q + j * n + k), _mm256_mul_ps(_mm256_loadu_ps(ws + k), a4)));
			}
			for (; k < n; k++) {
				q[j*n + k] += ws[k] * a;
			}
		}
	}
	hc(m, n);
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

		householderqr_m(m, ws, q, size);
		printf("q:\n");
		pm(q, size, size);
		printf("r:\n");
		pm(m, size, size);
		mult(q, m, ws, size, size, size, size);
		//printf("qr:\n");
		//pm(ws, size, size);

		free(m);
		free(q);
		free(ws);
	}
}
