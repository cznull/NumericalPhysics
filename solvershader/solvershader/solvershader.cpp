// solvershader.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <stdio.h>
#include <Windows.h>
#include <gl/glew.h>
#include <math.h>

HDC hdc1;
HGLRC m_hrc;
GLuint v1, v2, m;
GLuint shadermult, programmult, shadernorm, programnorm;
GLuint mn, nn;

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);

const char *ssmult ="#version 430 core\nlayout(local_size_x = 32) in;layout(r32f, binding = 0) uniform image1D v1;layout(r32f, binding = 1) uniform image1D v2;layout(r32f, binding = 2) uniform image2D m;uniform int n;void main() {int i;float a;ivec2 pos;a = 0.0;pos.x = int(gl_GlobalInvocationID.x);for (i = 0; i < n; i++) {pos.y = i;a += imageLoad(m, pos.xy).x*imageLoad(v1, i).x;}imageStore(v2, pos.x, vec4(a));}";

const char *ssnorm ="#version 430 core\nlayout(local_size_x = 1) in;layout(r32f, binding = 0) uniform image1D v1;layout(r32f, binding = 1) uniform image1D v2;uniform int n;void main() {int i;	float a, b;	a = 0.0;for (i = 0; i < n; i++) {b = imageLoad(v1, i).x;a += b * b;	}a = 1.0 / sqrt(a);for (i = 0; i < n; i++) {b = imageLoad(v1, i).x;imageStore(v2, i, vec4(b*a));}}";

ATOM MyRegisterClass(HINSTANCE hInstance)
{
	WNDCLASSEX wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style = CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc = WndProc;
	wcex.cbClsExtra = 0;
	wcex.cbWndExtra = 0;
	wcex.hInstance = hInstance;
	wcex.hIcon = LoadIcon(hInstance, NULL);
	wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
	wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	wcex.lpszMenuName = NULL;
	wcex.lpszClassName ="test";
	wcex.hIconSm = LoadIcon(wcex.hInstance, NULL);

	return RegisterClassEx(&wcex);
}

BOOL InitInstance(HINSTANCE hInstance, int nCmdShow,HWND &hWnd)
{
	hWnd = CreateWindow("test", "test", WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

	if (!hWnd)
	{
		return FALSE;
	}

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	return TRUE;
}

int init(int n) {
	int i;
	char s[1024];
	PIXELFORMATDESCRIPTOR pfd = {
			sizeof(PIXELFORMATDESCRIPTOR),
			1,
			PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER | PFD_STEREO,
			PFD_TYPE_RGBA,
			24,
			0,0,0,0,0,0,0,0,
			0,
			0,0,0,0,
			32,
			0,0,
			PFD_MAIN_PLANE,
			0,0,0,0
	};
	HWND hwnd; 
	HINSTANCE hInstance = ::GetModuleHandle(NULL);

	MyRegisterClass(hInstance);
	InitInstance(hInstance, SW_SHOW, hwnd);

	hdc1 = GetDC(hwnd);
	int uds = ::ChoosePixelFormat(hdc1, &pfd);
	::SetPixelFormat(hdc1, uds, &pfd);
	m_hrc = ::wglCreateContext(hdc1);
	::wglMakeCurrent(hdc1, m_hrc);	
	glewInit();

	i = glGetError();
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_TEXTURE_1D);

	i = glGetError();
	glGenTextures(1, &v1);
	glGenTextures(1, &v2);
	glGenTextures(1, &m); 

	i = glGetError();

	glBindTexture(GL_TEXTURE_1D,v1);
	glTexStorage1D(GL_TEXTURE_1D, 1, GL_R32F, n);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glBindTexture(GL_TEXTURE_1D, v2);
	glTexStorage1D(GL_TEXTURE_1D, 1, GL_R32F, n);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glBindTexture(GL_TEXTURE_2D, m);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32F, n,n);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	shadermult = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(shadermult, 1, &ssmult, NULL);
	glCompileShader(shadermult);
	glGetInfoLogARB(shadermult, 1024, &i, s);
	printf("%s\n", s);

	programmult = glCreateProgram();
	glAttachShader(programmult, shadermult);
	glLinkProgram(programmult);
	glGetInfoLogARB(programmult, 1024, &i, s);
	printf("%s\n", s);

	glUseProgram(programmult);
	mn = glGetUniformLocation(programmult, "n");
	glUniform1i(mn, n);
	glBindImageTexture(2, m, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32F);

	shadernorm = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(shadernorm, 1, &ssnorm, NULL);
	glCompileShader(shadernorm);
	glGetInfoLogARB(shadernorm, 1024, &i, s);
	printf("%s\n", s);

	programnorm = glCreateProgram();
	glAttachShader(programnorm, shadernorm);
	glLinkProgram(programnorm);
	glGetInfoLogARB(programnorm, 1024, &i, s);
	printf("%s\n", s);

	glUseProgram(programnorm);
	nn = glGetUniformLocation(programnorm, "n");
	glUniform1i(nn, n);

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
int identityvector(T *x, int n) {
	int i;
	for (i = 0; i < n; i++) {
		x[i] = 1.0;
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
int rand(T *m, int n) {
	int i;
	for (i = 0; i < n*n; i++) {
		m[i] = rand()*(1.0 / RAND_MAX);
	}
	return 0;
}

template <class T>
int ocsh(T *m, int n) {
	int i;
	memset(m, 0, n*n * sizeof(T));
	m[n - 1] = -1.0;
	m[0] = 2.0;
	m[1] = -1.0;
	for (i = 1; i < n - 1; i++) {
		m[i*n + i - 1] = -1.0;
		m[i*n + i] = 2.0;
		m[i*n + i + 1] = -1.0;
	}
	m[(n - 1)*n + n - 2] = -1.0;
	m[(n - 1)*n + n - 1] = 2.0;
	m[(n - 1)*n] = -1.0;
	return 0;
}

template <class T>
int norm2(T *v, T &x, int n) {
	x = 0.0;
	int i;
	for (i = 0; i < n; i++) {
		x += v[i] * v[i];
	}
	return 0;
}

int hc(float *m, int n) {
	int i, j;
	float x;
	for (i = 1; i < n; i++) {
		for (j = 0; j < i; j++) {
			x = m[i*n + j];
			m[i*n + j] = m[j*n + i];
			m[j*n + i] = x;
		}
	}
	return 0;
}

typedef float current;

int compute(int n,int r,float *v) {
	int i;
	int err;
	for (i = 0; i < r; i++) {
		glUseProgram(programmult);
		glBindImageTexture(0, v1, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32F);
		glBindImageTexture(1, v2, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_R32F);
		glDispatchCompute((n + 31) / 32, 1, 1);

		glFinish();
		glUseProgram(programnorm);
		glBindImageTexture(0, v2, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32F);
		glBindImageTexture(1, v1, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_R32F);
		glDispatchCompute(1, 1, 1);
		glFinish();
	}
	glBindTexture(GL_TEXTURE_1D, v1);
	glGetTexImage(GL_TEXTURE_1D, 0,GL_RED, GL_FLOAT, v);
	return 0;
}
	
int main()
{
	int n, r, i;
	current *hm;
	current *x, *x1;
	current l;

	scanf("%d%d", &n, &r);

	hm = (current*)malloc(n*n * sizeof(current));
	x = (current*)malloc(n * sizeof(current));
	x1 = (current*)malloc(n * sizeof(current));

	ocsh(hm, n);
	for (i = 0; i < n; i++) {
		x[i] = 0;
	}
	x[0] = 1;

	if (n > 512) {
		init(n);

		glBindTexture(GL_TEXTURE_1D, v1);
		glTexSubImage1D(GL_TEXTURE_1D, 0, 0, n, GL_RED, GL_FLOAT, x);

		glBindTexture(GL_TEXTURE_2D, m);
		hc(hm, n);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, n, n, GL_RED, GL_FLOAT, hm);
		hc(hm, n);
		compute(n, r, x);
	}
	else {
		for (i = 0; i < r; i++) {
			mult(hm, x, x1, n, n, n, 1);
			norm2(x1, l, n);
			l = 1.0 / sqrt(l);
			mult(x1, &l, x, n, 1, 1, 1);
		}
	}
	mult(hm, x, x1, n, n, n, 1);
	mult(x, x1, &l, 1, n, n, 1);
	pv(x, n);
	printf("%f", l);
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{
	case WM_PAINT:
	{
		PAINTSTRUCT ps;
		HDC hdc = BeginPaint(hWnd, &ps);
		EndPaint(hWnd, &ps);
	}
	break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}
