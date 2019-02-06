#include "vec.h"

template <typename T, typename F>
T rk4(T x, double t, double h, F f) {
	T k1 = h * f(x, t);
	T k2 = h * f(x + 0.5*k1, t + 0.5*h);
	T k3 = h * f(x + 0.5*k2, t + 0.5*h);
	T k4 = h * f(x + k3, t + h);
	return x + (1.0 / 6.0)*(k1 + k2 + k2 + k3 + k3 + k4);
}

template <typename T, typename F>
void sd(T x, double t, double h, int n, F f, std::vector<double2> &re) {
	for (int i = 0; i < n; i++) {
		re.push_back(x);
		x = rk4(x, t, h, f);
		t += h;
		re.push_back(x);
	}
}

int main(){
		double2 x0[5] = { {0.8, 0.8}, {1.0, 1.0}, {1.2, 1.2}, {1.4, 1.4}, {1.6, 1.6} };
		std::vector<double2> re;

		for (int i = 0; i < 5; i++) {
			sd(x0[i], 0, (1.0 / 1024), 65536, [](double2 x, double t)->double2 {return{ (2.0 / 3.0)*x.x - (4.0 / 3.0)*x.x*x.y,x.x*x.y - x.y }; }, re);
		}
    return 0;
}
