#include <stdio.h>
#include <math.h>
#include <string.h>

template <class T>
int pm(T *m, int a, int b);

template <class T>
int solve(T *m, T *y, T *x, int n) {
	int i, j, k;
	T a;
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			a = -m[j*n + i] / m[i*n + i];
			m[j*n + i] = 0;
			for (k = i + 1; k < n; k++) {
				m[j*n + k] += a * m[i*n + k];
			}
			y[j] += a * y[i];
		}
	}
	for (i = 0; i < n; i++) {
		a = 0;
		for (j = 0; j < i; j++) {
			a += x[n - 1 - j] * m[(n - 1 - i)*n + (n - 1 - j)];
		}
		x[n - 1 - i] = (y[n - 1 - i] - a) / m[(n - 1 - i)*n + (n - 1 - i)];
	}
	return 0;
}


template <class T>
int solvep(T *m, T *y, T *x, int n) {
	int i, j, k, xm, ym;
	T a;
	int *change = new int[n];
	T *xtemp = new T[n];
	for (i = 0; i < n; i++) {
		change[i] = i;
	}
	for (i = 0; i < n - 1; i++) {
		a = abs(m[i*n + i]);
		xm = i;
		ym = i;
		for (j = i ; j < n; j++) {
			for (k = i; k < n; k++) {
				if (fabs(m[j*n + k])>a) {
					xm = k;
					ym = j;
					a = fabs(m[j*n + k]);
				}
			}
		}
		for (j = 0; j < n; j++) {
			a = m[j*n + i];
			m[j*n + i] = m[j*n + xm];
			m[j*n + xm] = a;
		}
		j = change[i];
		change[i] = change[xm];
		change[xm] = j;
		for (k = 0; k < n; k++) {
			a = m[i*n + k];
			m[i*n + k] = m[ym*n + k];
			m[ym*n + k] = a;
		}
		a = y[i];
		y[i] = y[ym];
		y[ym] = a;
		for (j = i + 1; j < n; j++) {
			a = -m[j*n + i] / m[i*n + i];
			m[j*n + i] = 0;
			for (k = i + 1; k < n; k++) {
				m[j*n + k] += a * m[i*n + k];
			}
			y[j] += a * y[i];
		}
	}
	for (i = 0; i < n; i++) {
		a = 0;
		for (j = 0; j < i; j++) {
			a += xtemp[n - 1 - j] * m[(n - 1 - i)*n + (n - 1 - j)];
		}
		xtemp[n - 1 - i] = (y[n - 1 - i] - a) / m[(n - 1 - i)*n + (n - 1 - i)];
	}
	for (i = 0; i < n; i++) {
		x[change[i]] = xtemp[i];
	}
	delete change;
	delete xtemp;
	return 0;
}

template <class T>
int cholesky(T *hm, T *u, int n) {
	int i, j, k;
	T a, b;
	memset(u, 0, n*n * sizeof(T));
	u[0] = sqrt(hm[0]);
	for (i = 1; i < n; i++) {
		b = 0;
		for (j = 0; j < i; j++) {
			a = 0;
			for (k = 0; k < j; k++) {
				a += u[k*n + i] * u[k*n + j];
			}
			u[j*n + i] = (hm[j*n + i] - a) / u[j*n + j];
			b += u[j*n + i] * u[j*n + i];
		}
		if (hm[i*n + i] < b) {
			return -1;
		}
		u[i*n + i] = sqrt(hm[i*n + i] - b);
	}
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

template <class T>
int solve3(void) {
	return 0;
}

template <class T>
int solveut(T *u, T *y, T *x, int n) {
	int i, j;
	T a;
	for (i = 0; i < n; i++) {
		a = 0;
		for (j = 0; j < i; j++) {
			a += x[j] * u[j*n + i];
		}
		x[i] = (y[i] - a) / u[i*n + i];
	}
	return 0;
}

template <class T>
int identityvector(T *x, int n) {
	int i;
	for (i = 0; i < n; i++) {
		x[i] = 1.0;
	}
	return 0;
}

template <class T>
int solveu(T *u, T *y, T *x, int n) {
	int i, j;
	T a;
	for (i = 0; i < n; i++) {
		a = 0;
		for (j = 0; j < i; j++) {
			a += x[n - 1 - j] * u[(n - 1 - i)*n + n - 1 - j];
		}
		x[n - 1 - i] = (y[n - 1 - i] - a) / u[(n - 1 - i)*n + n - 1 - i];
	}
	return 0;
}

template <class T>
int hilbert(T *m, int n) {
	int i, j;
	for (i = 0; i < n; i++) {
		m[i] = 1.0 / (i+1);
	}
	for (i = 1; i < n; i++) {
		for (j = 0; j < n - 1; j++) {
			m[i*n + j] = m[i*n + j - n + 1];
		}
		m[i*n + j] = 1.0 / (i + j + 1);
	}
	return 0;
}

template <class T>
int pv(T *v, int n) {
	int i;
	for (i = 0; i < n-1; i++) {
		printf("%f,",v[i]);
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

typedef double current;

int main()
{
	int size = 7;
	current m[100] = { 1,1,1,0,1,1,0,0,1 };
	current hm[9] = { 2,1,0,1,2,1,0,1,2 };
	current u[100],l[9],a[9];
	current y[10] = { 1,1,1,1,1,1,1,1,1,1 };
	current x[10];

	for (size = 1; size <= 10; size++) {

		identityvector(y, size);
		hilbert(m, size);
		pm(m, size, size);
		printf("\n");
		solve(m, y, x, size);
		pv(x, size);
		printf("\n");

		identityvector(y, size);
		hilbert(m, size);
		solvep(m, y, x, size);
		pv(x, size);
		printf("\n");

		identityvector(y, size);
		hilbert(m, size);
		cholesky(m, u, size);
		pm(u, size, size);
		printf("\n");
		solveut(u, y, x, size);
		solveu(u, x, x, size);
		pv(x, size);
		printf("\n");
	}
}