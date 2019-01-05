#include <stdio.h>
#include <math.h>

#define PI 3.1415926535897932
#define SUM1 6
#define SUM2 2
#define TN 4096
#define INTE 0.0625

template <class T,class F>
T integrate(F f, T l, T yl, T r, T yr, int depth, T err) {
	T m, ym, d;
	m = 0.5*(l + r);
	ym = f(m);
	d = ym - 0.5*(yl + yr);
	if (depth && (d > err || d < -err || (l - r)>INTE)) {
		return integrate(f, l, yl, m, ym, depth - 1, err) + integrate(f, m, ym, r, yr, depth - 1, err);
	}
	else {
		return (ym*4.0 + yl + yr)*(1.0 / 6.0)*(r - l);
	}
}

template <class T>
T func1(T q, int x, int y, int z) {
	T r = x * x + y * y + z * z;
	return exp(q - r) / (r - q);
}

template <class T>
T func2(T q,T t) {
	T r = sqrt(t);
	r = r * r*r;
	return (exp(t*q) - 1) / r;
}

template <class T>
T func3(T q, int x, int y, int z,T t) {
	T r = x * x + y * y + z * z;
	return exp(t*q)*exp(-PI * PI*r / t);
}

template <class T>
T sum1(T q) {
	int x, y, z;
	T s = 0;
	for (x = -SUM1; x <= SUM1; x++) {
		for (y = -SUM1; y <= SUM1; y++) {
			for (z = -SUM1; z <= SUM1; z++) {
				s += func1(q, x, y, z);
			}
		}
	}
	return s;
}

template <typename T>
T sum2(T q, T t) {
	int x, y, z;
	T s = 0;
	for (x = -SUM2; x <= SUM2; x++) {
		for (y = -SUM2; y <= SUM2; y++) {
			for (z = -SUM2; z <= SUM2; z++) {
				s += func3(q, x, y, z, t);
			}	
		}
	}
	s -= func3(q, 0, 0, 0, t);
	return s;
}

template <typename T>
struct fx{
	double(*f)(double);
};

template <typename T>
T func(T q) {
	T yl1, yr1;
	T yl2, yr2;
	T i1, i2;
	T t0;
	t0 = 1.0 / TN;
	yl1 = func2(q, t0);
	yr1 = func2(q, 1.0);
	i1 = integrate([q](double t)->double {return func2(q, t); }, t0, yl1, 1.0, yr1, 20, 1e-10);
	i1 += 2.0 * sqrt(t0)*q + (1.0 / 3.0)*t0*sqrt(t0)*q*q + (1.0 / 15.0)*t0*t0*sqrt(t0)*q*q*q;

	yl2 = 0.0;
	yr2 = sum2(q, 1.0);
	i2 = integrate(
		[q](double t)->double {
		T a;
		a = sqrt(t);
		a = a * a*a;
		return sum2(q, t) / a;
	},
		0.0, yl2, 1.0, yr2, 20, 1e-10);

	return 0.5 / sqrt(PI)*sum1(q) - PI + 0.5*PI*i1 +0.5*sqrt(PI)*i2;
}

int main()
{
	int i;
	double q,y;
	double l, r;
	l = 0.01, r = 0.99;
	for (i = 0; i < 53; i++) {
		q = 0.5*(l + r);
		y = func(q) - PI * sqrt(PI)*(1.0 + 0.25*q);
		if (y > 0) {
			r = q;
		}
		else {
			l = q;
		}
	}
	printf("%.17f", q);
}