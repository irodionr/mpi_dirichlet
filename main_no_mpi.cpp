#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 100 // array size
#define PI 3.14159265359

// f1, f2, g1, g2 - boundary conditions
// x = 0..a, y = 0
double f1(double x) {
	return x*x;
}

// x = 0..a, y = b
double f2(double x) {
	return 4*cos(2*PI*x);
}

// x = 0, y = 0..b
double g1(double y) {
	return 4*y;
}

// x = a, y = 0..b
double g2(double y) {
	return 4*cos(2*PI*y);
}

int main(int argc, char *argv[]) {
	// .csv file for storing solution u(x, y) 
	FILE* fp;
	
	// 0 <= x <= a, 0 <= y <= b
	// h - step, e - epsilon, max - maximal improvement
	// u0, u - solution u(x, y)
	// x, y - axes
	double a, b, h, e, max, u0[N][N], u[N][N], x[N], y[N];

	// i, j - loop variables
	// count - amount of loops until max <= e
	int i, j, count;

	clock_t t0, t1;

	fp = fopen("u.csv", "w");
	if (!fp) {
		perror("File opening failed");

		return 1;
	}

	a = 2.0;
	b = 1.0;
	h = 1.0 / N;
	e = 0.001;
	max = e + 1; // initial value must be > e
	count = 0;

	t0 = clock();

	// evenly fill x, y from 0 to a, 0 to b with N values
	for (i = 0; i < N; i++) {
		x[i] = i*h*a;
		y[i] = i*h*b;
	}

	// put boundary conditions into u0 and copy to u
	for (i = 0; i < N; i++) {
		u0[i][0] = f1(x[i]);
		u0[0][i] = g1(y[i]);
		u0[i][N-1] = f2(x[i]);
		u0[N-1][i] = g2(y[i]);

		u[i][0] = u0[i][0];
		u[0][i] = u0[0][i];
		u[i][N-1] = u0[i][N-1];
		u[N-1][i] = u0[N-1][i];
	}

	// initial approximation
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			u0[i][j] = 0.0;
		}
	}

	// main iterations
	while (max > e) {
		// calculate u(x, y) as average of adjacent values
		for (i = 1; i < N - 1; i++) {
			for (j = 1; j < N - 1; j++) {
				u[i][j] = 0.25*(u0[i+1][j]+u0[i-1][j]+u0[i][j+1]+u0[i][j-1]);
			}
		}

		// find maximal improvement from u0 to u
		max = fabs(u[0][0]-u0[0][0]);
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				if (fabs(u[i][j]-u0[i][j]) > max) {
					max = fabs(u[i][j]-u0[i][j]);
				}
			}
		}

		// copy u to u0 for next iteration
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				u0[i][j] = u[i][j];
			}
		}

		count++;
	}

	t1 = clock();

	// print solution u(x, y) to .csv file
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fprintf(fp, "%.3f,", u[j][i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// print stats to console
	printf("%dx%d array\nCount = %d\nTime = %f", N, N, count, (double)(t1-t0)/CLOCKS_PER_SEC);

	fclose(fp);

	return 0;
}
