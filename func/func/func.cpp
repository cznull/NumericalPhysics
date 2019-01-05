#include "stdafx.h"
#include "func.h"

#include <gl/glew.h>
#include <math.h>
#include <vector>
#include "vec.h"
#include <sstream>
#include <COmmdlg.h>

#define MAX_LOADSTRING 100
#define PI 3.14159265358979324

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name

HDC hdc1, hdc2;
HGLRC m_hrc;
int mx, my, cx, cy;
double ang1, ang2, len, cenx, ceny, cenz;
GLuint vbo[10];
int size;
float3 color[10] = { {0.7,0.0,0.0},{0.5,0.5,1.0},{0.5,1.0,0.5},{1.0,1.0,0.5},{0.5,0.0,0.5},{0.0,0.5,0.0},{0.0,0.0,0.5},{0.5,1.0,1.0},{1.0,0.7,0.0},{1.0,0.5,1.0} };
std::vector<float3> keypointsrecord[4];

std::stringstream s;

// Forward declarations of functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

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
		/*pm(m, n, n);
		printf("\n");*/
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

int solve3(double *m, double *y, double *x, int n, double2 *workspace) {
	int i;
	workspace[0].x = m[1];
	for (i = 1; i < n; i++) {
		workspace[i].y = m[i * 3] / workspace[i - 1].x;
		workspace[i].x = m[i * 3 + 1] - workspace[i].y*m[i * 3 - 1];
	}
	x[0] = y[0];
	for (i = 1; i < n; i++) {
		x[i] = y[i] - workspace[i].y*x[i - 1];
	}
	x[n - 1] = x[n - 1] / workspace[n - 1].x;
	for (i = n - 2; i >= 0; i--) {
		x[i] = (x[i] - m[i * 3 + 2] * x[i + 1]) / workspace[i].x;
	}
	return 0;
}

int solve3_(double *m, double *y, double *x, int n, double2 *workspace) {
	double *m1;
	int i;
	m1 = (double*)malloc(n*n * sizeof(double));
	memset(m1, 0, n*n * sizeof(double));
	for (i = 1; i < n - 1; i++) {
		m1[i*n + i - 1] = m[i * 3];
		m1[i*n + i] = m[i * 3 + 1];
		m1[i*n + i + 1] = m[i * 3 + 2];
	}
	m1[n - 1] = m[0];
	m1[0] = m[1];
	m1[1] = m[2];
	m1[(n - 1)*n + n - 2] = m[(n - 1) * 3];
	m1[(n - 1)*n + n - 1] = m[(n - 1) * 3 + 1];
	m1[(n - 1)*n] = m[(n - 1) * 3 + 2];
	solve(m1, y, x, n);
	free(m1);
	return 0;
}

int getspline3para(double2 *keypoints,int n, double *workspace, int type, double3 *para) {
	double *m, *y, *x, hl, hr;
	int i,j;
	m = workspace;
	y = workspace + 3 * n;
	x = workspace + 4 * n;
	hl = keypoints[1].x - keypoints[0].x;
	for (i = 1; i < n - 1; i++) {
		hr = keypoints[i + 1].x - keypoints[i].x;
		m[i * 3 + 0] = hl;
		m[i * 3 + 1] = 2.0*(hl+hr);
		m[i * 3 + 2] = hr;
		y[i] = 6.0*((keypoints[i + 1].y - keypoints[i].y) / hr - (keypoints[i].y - keypoints[i - 1].y) / hl);
		hl = hr;
	}
	switch (type) {
	case 0://y'''==0
		m[0] = 0;
		m[1] = 1;
		m[2] = -1;
		m[3 * (n - 1)] = -1;
		m[3 * (n - 1) + 1] = 1;
		m[3 * (n - 1) + 2] = 0;
		y[0] = 0;
		y[n - 1] = 0;
		solve3(m, y, x, n, (double2*)(workspace + 5 * n));
		break;
	case 1://y''
		m[0] = 0;
		m[1] = m[3]-m[5];
		m[2] = m[4]+2.0*m[5];
		m[3 * (n - 1)] = m[3 * (n - 2) + 1] + 2.0 * m[3 * (n - 2)];
		m[3 * (n - 1) + 1] = m[3 * (n - 2) + 2] - m[3 * (n - 2)];
		m[3 * (n - 1) + 2] = 0;
		y[0] = y[1];
		y[n - 1] = y[n-2];
		solve3(m, y, x, n, (double2*)(workspace + 5 * n));
		break;
	case  2:
		hl = keypoints[n - 1].x - keypoints[n - 2].x;
		hr = keypoints[1].x - keypoints[0].x;
		m[0] = hl;
		m[1] = 2.0*(hl+hr);
		m[2] = hr;
		y[0] = 6.0*((keypoints[1].y - keypoints[0].y) / hr - (keypoints[n - 1].y - keypoints[n - 2].y) / hl);
		solve3_(m, y, x, n-1, (double2*)(workspace + 5 * n));
		x[n - 1] = x[0];
		break;
	}
	for (i = 0; i < n - 1; i++) {
		hr = keypoints[i + 1].x - keypoints[i].x;
		para[i].x = (x[i + 1] - x[i]) / (6 * hr);
		para[i].y = 0.5*x[i];
		para[i].z = (keypoints[i + 1].y - keypoints[i].y) / hr - (para[i].x*hr + para[i].y)*hr;
	}
	return 0;
}

int getfuncspline3(double2 *keypoints, int n, float3 *lines, int count, double l, double r,int type) {
	double *workspace;
	double3 *para;
	float3 v;
	int i, cur;
	double dx;
	workspace = (double *)malloc(7 * n * sizeof(double));
	para = (double3 *)malloc((n - 1) * sizeof(double3));
	getspline3para(keypoints, n, workspace, type, para);

	cur = 0;
	v.x = l;
	while (cur<n - 1 && v.x > keypoints[cur + 1].x) {
		cur++;
	}
	dx = v.x - keypoints[cur].x;
	v.y = ((para[cur].x*dx + para[cur].y)*dx + para[cur].z)*dx + keypoints[cur].y;
	v.z = 0;
	lines[0] = v;
	for (i = 1; i < count; i++) {
		v.x = l + (r-l)*i / count;
		while (cur<n - 1 && v.x > keypoints[cur + 1].x) {
			cur++;
		}
		dx = v.x - keypoints[cur].x;
		v.y = ((para[cur].x*dx + para[cur].y)*dx + para[cur].z)*dx + keypoints[cur].y;
		v.z = 0;
		lines[i * 2 - 1] = v;
		lines[i * 2] = v;
	}
	v.x = l + (r - l)*i / count;
	while (cur<n - 1 && v.x > keypoints[cur + 1].x) {
		cur++;
	}
	dx = v.x - keypoints[cur].x;
	v.y = ((para[cur].x*dx + para[cur].y)*dx + para[cur].z)*dx + keypoints[cur].y;
	v.z = 0;
	lines[count * 2 - 1] = v;

	free(workspace);
	free(para);
	return 0;
}

int getparaspline3(double2 *keypointsx, double2 *keypointsy, int n, float3 *lines, int count, double l, double r,int type) {
	double *workspace;
	double3 *parax,*paray;
	float3 v;
	int i, cur;
	double dx,t;
	workspace = (double *)malloc(7 * n * sizeof(double));
	parax = (double3 *)malloc((n - 1) * sizeof(double3));
	paray = (double3 *)malloc((n - 1) * sizeof(double3));
	getspline3para(keypointsx, n, workspace, type, parax);
	getspline3para(keypointsy, n, workspace, type, paray);

	cur = 0;
	t = l;
	while (cur<n - 1 && t > keypointsx[cur + 1].x) {
		cur++;
	}
	dx = t - keypointsx[cur].x;
	v.x = ((parax[cur].x*dx + parax[cur].y)*dx + parax[cur].z)*dx + keypointsx[cur].y;
	v.y = ((paray[cur].x*dx + paray[cur].y)*dx + paray[cur].z)*dx + keypointsy[cur].y;
	v.z = 0;
	lines[0] = v;
	for (i = 1; i < count; i++) {
		t = l + (r - l)*i / count;
		while (cur<n - 1 && t > keypointsx[cur + 1].x) {
			cur++;
		}
		dx = t - keypointsx[cur].x;
		v.x = ((parax[cur].x*dx + parax[cur].y)*dx + parax[cur].z)*dx + keypointsx[cur].y;
		v.y = ((paray[cur].x*dx + paray[cur].y)*dx + paray[cur].z)*dx + keypointsy[cur].y;
		v.z = 0;
		lines[i * 2 - 1] = v;
		lines[i * 2] = v;
	}
	t = l + (r - l)*i / count;
	while (cur<n - 1 && t > keypointsx[cur + 1].x) {
		cur++;
	}
	dx = t - keypointsx[cur].x;
	v.x = ((parax[cur].x*dx + parax[cur].y)*dx + parax[cur].z)*dx + keypointsx[cur].y;
	v.y = ((paray[cur].x*dx + paray[cur].y)*dx + paray[cur].z)*dx + keypointsy[cur].y;
	v.z = 0;
	lines[count * 2 - 1] = v;

	free(workspace);
	free(parax);
	free(paray);
	return 0;
}

double getfuncvalue1(double2 *keypoints, int n, double *workspace, double x) {
	int i, j;
	for (i = 0; i < n; i++) {
		workspace[i] = keypoints[i].y;
	}
	for (i = 1; i < n; i++) {
		for (j = 0; j < n - i; j++) {
			workspace[j] = ((x - keypoints[j].x)*workspace[j+1] + (keypoints[j + i].x - x)*workspace[j]) / (keypoints[j + i].x - keypoints[j].x);
		}
	}
	return workspace[0];
}

double getfunccvalue1(double *ck, int n, double x) {
	int i, j;
	double a, b, c;
	a = ck[n - 1];
	b = ck[n - 2] + a * 2.0*x;
	for (i = n - 3; i > 0; i--) {
		c = ck[i] + 2.0*x*b - a;
		a = b;
		b = c;
	}
	return 0.5*ck[0] + x * b - a;
}


template <class F>
int getkeypoint(double2 *keypoints, int count, F f, double l, double r, std::vector<float3> &re) {
	int i;
	double x;
	double *workspace;
	workspace = (double *)malloc(count * sizeof(double));
	for (i = 0; i < count; i++) {
		x = l + (r - l) *i / (count - 1);
		keypoints[i] = { x,f(x) };
		re.push_back({ x,f(x),0 });
	}
	free(workspace);
	return 0;
}

template <class F>
int pts(double2 *keypoints, int count, F f, double l, double r, std::stringstream &s) {
	int i;
	double x;
	double *workspace;
	workspace = (double *)malloc(count * sizeof(double));
	s << "x\tf(x)\tP20(x)\n";
	for (i = 0; i < count - 1; i++) {
		x = keypoints[i].x;
		s << x << "\t" << f(x) << "\t" << getfuncvalue1(keypoints, count, workspace, x) << "\t" << abs(getfuncvalue1(keypoints, count, workspace, x) - f(x)) << "\n";
		x = 0.5*(x + keypoints[i + 1].x);
		s << x << "\t" << f(x) << "\t" << getfuncvalue1(keypoints, count, workspace, x) << "\t" << abs(getfuncvalue1(keypoints, count, workspace, x) - f(x)) << "\n";
	}
	x = keypoints[i].x;
	s << x << "\t" << f(x) << "\t" << getfuncvalue1(keypoints, count, workspace, x) << "\t" << abs(getfuncvalue1(keypoints, count, workspace, x) - f(x)) << "\n\n";
	free(workspace);
	return 0;
}

template <class F>
int ptsc(double2 *keypoints, double *ck, int count, F f, double l, double r, std::stringstream &s) {
	int i;
	double x;
	double *workspace;
	workspace = (double *)malloc(count * sizeof(double));
	s << "x\tf(x)\tC(x)\n";
	for (i = 0; i < count - 1; i++) {
		x = keypoints[i].x;
		s << x << "\t" << f(x) << "\t" << getfunccvalue1(ck,count,x) << "\t" << abs(getfunccvalue1(ck, count, x) - f(x)) << "\n";
		x = 0.5*(x + keypoints[i + 1].x);
		s << x << "\t" << f(x) << "\t" << getfunccvalue1(ck, count, x) << "\t" << abs(getfunccvalue1(ck, count, x) - f(x)) << "\n";
	}
	x = keypoints[i].x;
	s << x << "\t" << f(x) << "\t" << getfunccvalue1(ck, count, x) << "\t" << abs(getfunccvalue1(ck, count, x) - f(x)) << "\n\n";
	free(workspace);
	return 0;
}

int ptspara(double2 *keypointsx, double2 *keypointsy, int count, std::stringstream &s) {
	int i;
	s << "t\tx(t)\ty(t)\n";
	for (i = 0; i < count; i++) {
		s << keypointsx[i].x << "\t" << keypointsx[i].y << "\t" << keypointsy[i].y << "\n";
	}
	s << "\n";
	return 0;
}

template <class F>
int pts3(double2 *keypoints, int count, F f, double l, double r, std::stringstream &s, int type) {
	double *workspace;
	double3 *para;
	float3 v;
	int i, cur;
	double dx;
	workspace = (double *)malloc(7 * count * sizeof(double));
	para = (double3 *)malloc((count - 1) * sizeof(double3));
	getspline3para(keypoints, count, workspace, type, para);

	s << "x\tf(x)\tspline(x)\n";
	cur = 0;
	for (i = 0; i < count - 1; i++) {
		v.x = keypoints[i].x;
		while (cur<count - 1 && v.x > keypoints[cur + 1].x) {
			cur++;
		}
		dx = v.x - keypoints[cur].x;
		v.y = ((para[cur].x*dx + para[cur].y)*dx + para[cur].z)*dx + keypoints[cur].y;
		s << v.x << "\t" << f(v.x) << "\t" << v.y << "\t" << abs(f(v.x) - v.y) << "\n";
		v.x = 0.5*(v.x + keypoints[i + 1].x);
		while (cur<count - 1 && v.x > keypoints[cur + 1].x) {
			cur++;
		}
		dx = v.x - keypoints[cur].x;
		v.y = ((para[cur].x*dx + para[cur].y)*dx + para[cur].z)*dx + keypoints[cur].y;
		s << v.x << "\t" << f(v.x) << "\t" << v.y << "\t" << abs(f(v.x) - v.y) << "\n";
	}
	v.x = keypoints[i].x;
	while (cur<count - 1 && v.x > keypoints[cur + 1].x) {
		cur++;
	}
	dx = v.x - keypoints[cur].x;
	v.y = ((para[cur].x*dx + para[cur].y)*dx + para[cur].z)*dx + keypoints[cur].y;
	s << v.x << "\t" << f(v.x) << "\t" << v.y << "\t" << abs(f(v.x) - v.y) << "\n\n";
	free(workspace);
	free(para);
	return 0;
}


int getkeypointc(double2 *keypoints,double *ck, int count, double(*f)(double x), double l, double r, std::vector<float3> &re) {
	int i,j;
	double x;
	for (i = 0; i < count; i++) {
		x = cos(PI *(i + 0.5) / count)*(r - l)*0.5 + (r + l)*0.5;
		keypoints[i] = { x,f(x) };
		re.push_back({ x,f(x),0 });
	}
	for (i = 0; i < count; i++) {
		ck[i] = 0;
		for (j = 0; j < count; j++) {
			ck[i] += keypoints[j].y*cos(PI*i*(j + 0.5) / count);
		}
		ck[i] *= (2.0 / count);
	}
	return 0;
}

int getparakeypoint(double2 *keypointsx, double2 *keypointsy, int count, double(*fx)(double x), double(*fy)(double x), double l, double r, std::vector<float3> &re) {
	int i;
	double x;
	for (i = 0; i < count; i++) {
		x = l + (r - l) *i / (count - 1);
		keypointsx[i] = { x,fx(x) };
		keypointsy[i] = { x,fy(x) };
		re.push_back({ fx(x),fy(x),0 });
	}
	return 0;
}

int getfuncc(double2 *keypoints, double *ck, int n, float3 *lines, int count, double l, double r) {
	int i;
	float3 v;
	double *workspace;
	workspace = (double *)malloc(n * sizeof(double));
	v.x = l;
	v.y = getfunccvalue1(ck, n, v.x);
	v.z = 0;
	lines[0] = v;
	for (i = 1; i < count; i++) {
		v.x = l + (r - l)*i / count;
		v.y = getfunccvalue1(ck, n, v.x);
		v.z = 0;
		lines[i * 2 - 1] = v;
		lines[i * 2] = v;
	}
	v.x = l + (r - l)*i / count;
	v.y = getfunccvalue1(ck, n, v.x);
	v.z = 0;
	lines[count * 2 - 1] = v;
	free(workspace);
	return 0;
}

int getfunc(double2 *keypoints, int n, float3 *lines,int count, double l, double r) {
	int i;
	float3 v;
	double *workspace;
	workspace = (double *)malloc(n * sizeof(double));
	v.x = l;
	v.y = getfuncvalue1(keypoints,n, workspace, v.x);
	v.z = 0;
	lines[0] = v;
	for (i = 1; i < count; i++) {
		v.x = l + (r-l)*i / count;
		v.y = getfuncvalue1(keypoints,n, workspace, v.x);
		v.z = 0;
		lines[i * 2 - 1] = v;
		lines[i * 2] = v;
	}
	v.x = l + (r - l)*i / count;
	v.y = getfuncvalue1(keypoints, n, workspace, v.x);
	v.z = 0;
	lines[count * 2 - 1] = v;
	free(workspace);
	return 0;
}

template<typename Lx, typename Ly>
int getoran(float3 *lines, int count, Lx fx, Ly fy, double l, double r) {
	int i;
	float3 v;
	float t;
	t = l;
	v.x = fx(t);
	v.y = fy(t);
	v.z = 0;
	lines[0] = v;
	for (i = 1; i < count; i++) {
		t = l + (r-l)*i / count;
		v.x = fx(t);
		v.y = fy(t);
		lines[i * 2 - 1] = v;
		lines[i * 2] = v;
	}
	t = l + (r - l)*i / count;
	v.x = fx(t);
	v.y = fy(t);
	lines[count * 2 - 1] = v;
	return 0;
}

void draw(void) {
	int i;
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(0x00004100);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	double m[16] = { 1.0 / len,0,0,0,0,1.0 / len,0,0,0,0,1.0 / len,0,-cenx / len,-ceny / len,-cenz / len - 1,1 };

	gluLookAt(len * cos(ang1)*cos(ang2) + cenx, len * sin(ang2) + ceny, len * sin(ang1)*cos(ang2) + cenz, cenx, ceny, cenz, 0, cos(ang2), 0);
	//gluLookAt(0, 0, 1, 0, 0, 0, 0, cos(ang2), 0);
	//glScalef(1.0 / len, 1.0 / len, 1.0 / len);
	//glTranslatef(-cenx, -ceny, -cenz);
	//glLoadMatrixd(m);
	glBegin(GL_LINES);
	glColor3f(0.3, 0.3, 0.3);
	glVertex3f(1.0f, 0.0f, -0.0f);
	glVertex3f(-1.0f, 0.0f, -0.0f);
	glVertex3f(0.0f, 1.0f, -0.0f);
	glVertex3f(0.0f, -1.0f, -0.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, -1.0f);
	glEnd();

	for (i = 0; i < 10; i++) {
		glBindBuffer(GL_ARRAY_BUFFER, vbo[i]);
		glVertexPointer(3, GL_FLOAT, 0, NULL);
		glColor3f(color[i].x, color[i].y, color[i].z);
		glDrawArrays(GL_LINES, 0, size * 2);
	}
	glBegin(GL_LINES);
	glColor3f(0.5, 0.5, 0.5);
	for (float3 x : keypointsrecord[3]) {
		glVertex3f(x.x + len * 0.01f, x.y, x.z);
		glVertex3f(x.x - len * 0.01f, x.y, x.z);
		glVertex3f(x.x, x.y + len * 0.01f, x.z);
		glVertex3f(x.x, x.y - len * 0.01f, x.z);
		glVertex3f(x.x, x.y, x.z + len * 0.01f);
		glVertex3f(x.x, x.y, x.z - len * 0.01f);
	}
	glEnd();
	SwapBuffers(wglGetCurrentDC());
	glFinish();
}

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Place code here.

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_FUNC, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_FUNC));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_FUNC));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_FUNC);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // Store instance handle in our global variable

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{
	case WM_COMMAND:
	{
		int wmId = LOWORD(wParam);
		// Parse the menu selections:
		switch (wmId)
		{
		case IDM_ABOUT:
			DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
			break;
		case IDM_EXIT:
			DestroyWindow(hWnd);
			break;
		case ID_FILE_SAVE: {
			TCHAR szBuffer[MAX_PATH] = { 0 };
			OPENFILENAME ofn = { 0 };
			FILE *fi;
			ofn.lStructSize = sizeof(ofn);
			ofn.hwndOwner = hWnd;
			ofn.lpstrFilter = NULL;
			ofn.lpstrInitialDir = NULL;
			ofn.lpstrFile = szBuffer;
			ofn.nMaxFile = sizeof(szBuffer) / sizeof(*szBuffer);
			ofn.nFilterIndex = 0;
			ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_EXPLORER;
			if (GetSaveFileName(&ofn)) {
				if (!fopen_s(&fi, szBuffer, "wb")) {
					fwrite(s.str().c_str(), 1, s.str().length(), fi);
					fclose(fi);
				}
			}
			break;
		}
		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
		}
	}
	break;
	case WM_PAINT:
	{
		PAINTSTRUCT ps;
		HDC hdc = BeginPaint(hWnd, &ps);
		draw();
		// TODO: Add any drawing code that uses hdc here...
		EndPaint(hWnd, &ps);
	}
	break; 
	case WM_CREATE: {
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
		hdc1 = GetDC(hWnd);
		hdc2 = GetDC(NULL);
		int uds = ::ChoosePixelFormat(hdc1, &pfd);
		::SetPixelFormat(hdc1, uds, &pfd);
		m_hrc = ::wglCreateContext(hdc1);
		::wglMakeCurrent(hdc1, m_hrc);
		glewInit();
		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_COLOR, GL_DST_COLOR);
		glBlendEquation(GL_FUNC_ADD);
		glBlendEquation(GL_MAX);
		glEnableClientState(GL_VERTEX_ARRAY);
		glGenBuffers(10, vbo);
		((bool(_stdcall*)(int))wglGetProcAddress("wglSwapIntervalEXT"))(0);
		size = 131072;
		//size = 1024;
		float3 *buffer;
		double2 *keypoints, *keypoints1;
		double *ck;
		keypoints = (double2*)malloc(1024 * sizeof(double2));
		keypoints1 = (double2*)malloc(1024 * sizeof(double2));
		ck = (double*)malloc(1024 * sizeof(double2));
		buffer = (float3*)malloc(size * 2 * sizeof(float3));

		getoran(buffer, size, [](double x)-> double {return x; }, [](double x)-> double {return 1.0 / (1.0 + 25.0*x*x); }, -1, 1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);

		getkeypoint(keypoints, 21, [](double x)-> double {return 1.0 / (1.0 + 25.0*x*x); }, -1, 1, keypointsrecord[0]);
		pts(keypoints, 21, [](double x)-> double {return 1.0 / (1.0 + 25.0*x*x); }, -1, 1, s);
		getfunc(keypoints, 21, buffer, size, -1, 1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);



		getkeypointc(keypoints,ck, 20, [](double x)-> double {return 1.0 / (1.0 + 25.0*x*x); }, -1, 1, keypointsrecord[1]);
		ptsc(keypoints, ck,20, [](double x)-> double {return 1.0 / (1.0 + 25.0*x*x); }, -1, 1, s);
		getfunc(keypoints, 20, buffer, size, -1, 1);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);

		getkeypoint(keypoints, 21, [](double x)-> double {return 1.0 / (1.0 + 25.0*x*x); }, -1, 1, keypointsrecord[2]);
		pts3(keypoints, 21, [](double x)-> double {return 1.0 / (1.0 + 25.0*x*x); }, -1, 1, s, 0);
		getfuncspline3(keypoints, 21, buffer, size, -1, 1, 0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);

		getoran(buffer, size, [](double x)-> double {return 1.0 *(1.0 - cos(x))*cos(x); }, [](double x)-> double {return 1.0 *(1.0 - cos(x))*sin(x); }, 0, 2 * PI);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[4]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);
		getoran(buffer, size, [](double x)-> double {return x; }, [](double x)-> double {return 1.0 *(1.0 - cos(x))*cos(x); }, 0, 2 * PI);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[5]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);
		getoran(buffer, size, [](double x)-> double {return x; }, [](double x)-> double {return 1.0 *(1.0 - cos(x))*sin(x); }, 0, 2 * PI);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[6]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);

		getparakeypoint(keypoints, keypoints1, 9, [](double x)-> double {return 1.0 *(1.0 - cos(x))*cos(x); }, [](double x)-> double {return 1.0 *(1.0 - cos(x))*sin(x); }, 0, 2 * PI, keypointsrecord[3]);
		ptspara(keypoints, keypoints1, 9, s);
		getparaspline3(keypoints, keypoints1, 9, buffer, size, 0, 2 * PI, 2);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[7]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);
		getfuncspline3(keypoints, 9, buffer, size, 0, 2 * PI, 2);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[8]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);
		getfuncspline3(keypoints1, 9, buffer, size, 0, 2 * PI, 2);
		glBindBuffer(GL_ARRAY_BUFFER, vbo[9]);
		glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);
		/*int i=0;
		for (i = 0; i < 10; i++) {
			auto f = [i](double x)-> double {return integrate1(x, i, 0, 2 * PI); };
			auto f1= [i](double x)-> double {return -pow(x,i)*2 * PI/i; };
			getoran(buffer, size, [](double x)-> double {return x; },f, 0, 1);
			glBindBuffer(GL_ARRAY_BUFFER, vbo[i]);
			glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float3), buffer, GL_STATIC_DRAW);
		}*/

		ang1 = PI / 2;
		ang2 = 0;
		len = 2.5;
		cenx = 0.0;
		ceny = 0.0;
		cenz = 0.0;
		free(buffer);
		free(keypoints);
		free(keypoints1);
		free(ck);
		break;
	}
	case WM_SIZE: {
		cx = lParam & 0xffff;
		cy = (lParam & 0xffff0000) >> 16;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-(float)cx / (cx + cy) *len*0.0078125, (float)cx / (cx + cy) *len*0.0078125, -(float)cy / (cx + cy)*len*0.0078125, (float)cy / (cx + cy) *len*0.0078125, len*0.00390625, len*100.0);
		glViewport(0, 0, cx, cy);
		break;
	}
	case WM_MOUSEMOVE: {
		int x, y, f;
		f = 0;
		x = (lParam & 0xffff);
		y = ((lParam & 0xffff0000) >> 16);
		if (MK_LBUTTON&wParam) {
			f = 1;
			//ang1 += (x - mx)*0.002;
			//ang2 += (y - my)*0.002;
		}
		if (MK_RBUTTON&wParam) {
			double l;
			f = 1;
			l = len * 4.0 / (cx + cy);
			cenx += l * (-(x - mx)*sin(ang1) - (y - my)*sin(ang2)*cos(ang1));
			ceny += l * ((y - my)* cos(ang2));
			cenz += l * ((x - mx)*cos(ang1) - (y - my)*sin(ang2)*sin(ang1));
		}
		mx = x;
		my = y;
		if (f) {
			draw();
		}
		break;
	}
	case WM_MOUSEWHEEL: {
		short m;
		double a,x,y,z,l;
		m = (wParam & 0xffff0000) >> 16;
		a = exp(-m * 0.001);
		l = len * 4.0 / (cx + cy);
		x= l * ((mx - cx*0.5+0.5)*sin(ang1) + (my - cy*0.5+0.5)*sin(ang2)*cos(ang1));
		y= l * (-(my - cy * 0.5 + 0.5)* cos(ang2));
		z= l * (-(mx - cx * 0.5 + 0.5)*cos(ang1) + (my - cy * 0.5 + 0.5)*sin(ang2)*sin(ang1));
		len *= a;
		cenx += x * (1 - a);
		ceny += y * (1 - a);
		cenz += z * (1 - a);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-(float)cx / (cx + cy) *len*0.0078125, (float)cx / (cx + cy) *len*0.0078125, -(float)cy / (cx + cy)*len*0.0078125, (float)cy / (cx + cy) *len*0.0078125, len*0.00390625, len*100.0);
		draw();
		break;
	}
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
