#include <iostream>
#include <fstream>
#include <cstdlib>

#include "matrix.hh"


Matrix hilbert(unsigned n) {
	Matrix M = Matrix(n, n);
	for (unsigned i = 0; i < n; i++){
		for(unsigned j = 0; j < n; j++){
			M.set(i,j,1.0/(i+j+1));
		}
	}
	return M;
}

float absv(float x) {
	float y = 0;
	if (x > 0) {
		y = x;
	}
	else {
		y = -x;
	}
	return y;
}

double time_classic(unsigned n) {
	clock_t timer;
	double calctime;
	Matrix M = hilbert(n);
	timer = clock();
	Matrix M0 = inverse(M);
	timer = clock() - timer;
	calctime = ((float)timer)/CLOCKS_PER_SEC;
	return calctime;
}

int precision_classic(unsigned n) {
	unsigned precision = 0;
	unsigned current;
	Matrix M = hilbert(n);
	Matrix I = Id(n);
	Matrix M0 = inverse(M);
	Matrix M1 = M*M0;
	Matrix M2 = M1 - I;
	for (unsigned i=0; i<n; i++) {
		for (unsigned j=0; j<n; j++) {
			current = absv(M2.get(i, j));
			if (precision < current) {
				precision = current;
			}
		}
	}
	return precision;
}

double time_gauss(unsigned n) {
	clock_t timer;
	double calctime;
	Matrix M = hilbert(n);
	timer = clock();
	Matrix M0 = inverse_gauss(M);
	timer = clock() - timer;
	calctime = ((float)timer)/CLOCKS_PER_SEC;
	return calctime;
}

int precision_gauss(unsigned n) {
	unsigned precision = 0;
	unsigned current;
	Matrix M = hilbert(n);
	Matrix I = Id(n);
	Matrix M0 = inverse_gauss(M);
	Matrix M1 = M*M0;
	Matrix M2 = M1 - I;
	for (unsigned i=0; i<n; i++) {
		for (unsigned j=0; j<n; j++) {
			current = absv(M2.get(i, j));
			if (precision < current) {
				precision = current;
			}
		}
	}
	return precision;
}

int main () {
	ofstream f;
	f.open("data_classic.dat");
	for (unsigned n = 2; n < 10; n++)  {
		f << n << " " << time_classic(n) << " " << precision_classic(n)
			<< "\n";
	}
	f.close();
	f.open("data_gauss.dat");
	for (unsigned n = 2; n < 21; n++)  {
		f << n << " " << time_gauss(n) << " " << precision_gauss(n)
			<< "\n";
	}
	f.close();
	return 0;
}
