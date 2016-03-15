#include <iostream>
#include <iomanip>

#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>

#include <vector>

#include "matrix.hh"

using namespace std;

Matrix::Matrix(unsigned n, unsigned p) {
  size_i = n;
  size_j = p;
  assert(1 <= size_i);
  assert(1 <= size_j);

  contents = vector<vector<scalar_t>/**/>(size_i);
  for (unsigned i = 0; i < size_i; i++)
    contents.at(i) = vector<scalar_t>(size_j, 0.0);	// Initialization to 0.0
}

unsigned Matrix::get_size_i() const {
  return size_i;
}

unsigned Matrix::get_size_j() const {
  return size_j;
}

void Matrix::set(unsigned i, unsigned j, scalar_t x) {
  assert(0 <= i && i < size_i);
  assert(0 <= j && j < size_j);
  contents.at(i).at(j) = x;
}

Matrix::scalar_t Matrix::get(unsigned i, unsigned j) const {
  assert(0 <= i && i < size_i);
  assert(0 <= j && j < size_j);
  return contents.at(i).at(j);
}

void Matrix::print() const {
  for (unsigned i = 0; i < size_i; i++) {
    for (unsigned j = 0; j < size_j; j++) {
      cout << setprecision(2) << setw(8) << contents.at(i).at(j);
    }
    cout << endl;
  }
  cout << "___________________________" << endl;
}

void swapline(Matrix& M1, unsigned i1, unsigned i2){
  unsigned size_j = M1.get_size_j();
  for (unsigned j = 0;j<size_j;j++){
    Matrix::scalar_t temp = M1.get(i2,j);
    M1.set(i2,j,M1.get(i1,j));
    M1.set(i1,j,temp);
  }
}

void multalphaline(Matrix& M1, unsigned i, Matrix::scalar_t alpha){
  unsigned size_j = M1.get_size_j();
  for (unsigned j = 0;j<size_j;j++)
    M1.set(i,j,M1.get(i,j)*alpha);
}

void addline(Matrix& M1, unsigned i1, unsigned i2, Matrix::scalar_t alpha){
  unsigned size_j = M1.get_size_j();
  for (unsigned j = 0;j<size_j;j++)
    M1.set(i1,j,M1.get(i1,j) + M1.get(i2,j) * alpha);
}

/*****************************************************/

Matrix operator+(const Matrix& M1, const Matrix& M2) {
  assert(M1.get_size_i() == M2.get_size_i());
  assert(M1.get_size_j() == M2.get_size_j());
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix M(size_i, size_j);

  for (unsigned i = 0; i < size_i; i++) {
    for (unsigned j = 0; j < size_j; j++) {
      Matrix::scalar_t x =  M1.get(i, j) + M2.get(i, j);
      M.set(i, j, x);
    }
  }

  return M;
}

Matrix operator-(const Matrix& M1, const Matrix& M2) {
  assert(M1.get_size_i() == M2.get_size_i());
  assert(M1.get_size_j() == M2.get_size_j());
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix M(size_i, size_j);

  for (unsigned i = 0; i < size_i; i++) {
    for (unsigned j = 0; j < size_j; j++) {
      Matrix::scalar_t x =  M1.get(i, j) - M2.get(i, j);
      M.set(i, j, x);
    }
  }

  return M;
}

Matrix operator*(Matrix::scalar_t a, const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix M(size_i, size_j);

  for (unsigned i = 0; i < size_i; i++) {
    for (unsigned j = 0; j < size_j; j++) {
      Matrix::scalar_t x =  M1.get(i, j) *a;
      M.set(i, j, x);
    }
  }

  return M;
}

Matrix operator*(const Matrix& M1, const Matrix& M2) {
  assert(M1.get_size_j() == M2.get_size_i());
  unsigned size_i1 = M1.get_size_i();
  unsigned size_j2 = M2.get_size_j();
  unsigned size_i2 = M2.get_size_i();
  Matrix M = Matrix(size_i1,size_j2);
  Matrix::scalar_t m;
  for(unsigned i = 0;i<size_i1;i++){
    for(unsigned j = 0;j<size_j2;j++){
      m = 0;
      for(unsigned k = 0;k<size_i2;k++){
          m += M1.get(i,k) * M2.get(k,j);
      }
      M.set(i,j,m);
    }
  }
  return M;
}

Matrix Id(unsigned n) {
  unsigned size = n;
  Matrix M(size, size);
  for (unsigned i = 0; i < size; i++) {
    M.set(i, i, 1.0);
  }
  return M;
}

Matrix::scalar_t norm(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix::scalar_t m = 0;
  for (unsigned i = 0; i < size_i; i++)
    for (unsigned j = 0; j < size_j; j++)
      if (abs(M1.get(i,j) > m))
       m = M1.get(i,j);
  return m;
}

static Matrix submatrix(const Matrix& M1, unsigned a, unsigned b) {	//Note it is static!
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert (0 <= a && a < size_i);
  assert (0 <= b && b < size_j);
  assert(size_i >= 2);
  assert(size_j >= 2);
  Matrix M = Matrix(size_i-1, size_j-1);
  for (unsigned i=0; i<size_i-1; i++) {
    for (unsigned j=0; j<size_j-1; j++) {
      if (i<a) {
        if (j<b)
          M.set(i, j, M1.get(i, j));
	else
	  M.set(i, j, M1.get(i, j+1));
      }
      else {
        if (j<b)
	  M.set(i, j, M1.get(i+1, j));
        else
	  M.set(i, j, M1.get(i+1, j+1));
      }
    }
  }
  return M;
}

static int toggle(unsigned i) {		//Note it is static!
  if (i%2)
    return 1;
  else
    return (-1);
}

Matrix::scalar_t determinant(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert(size_i == size_j);
  if (size_i == 1)
    return M1.get(0,0);
  Matrix::scalar_t det = 0;
  for (unsigned i = 0 ; i< size_i;i++){
    det += toggle(i) * M1.get(i,0) * determinant(submatrix(M1,i,0));
  }
  return det;
}

Matrix transpose(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  Matrix M = Matrix(size_j,size_i);
  for (unsigned i=0;i<size_i;i++){
    for (unsigned j=0;j<size_j;j++){
      M.set(j,i,M1.get(i,j));
    }
  }
  return M;
}

Matrix inverse(const Matrix& M1) {
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert(size_i == size_j);

  Matrix::scalar_t det = determinant(M1);
  if (det == 0) {
    cerr << "Cannot invert a matrix with null determinant" << endl;
    M1.print();
    exit(1);
  }

  unsigned size = size_i;
  Matrix cofactors(size, size);
  Matrix M2 = transpose(M1);	// Do not forget it!

  for (unsigned i = 0; i < size; i++) {
    for (unsigned j = 0; j < size; j++) {
      Matrix M = submatrix(M2, i, j);

      Matrix::scalar_t x = determinant(M);
      Matrix::scalar_t y =  toggle(i+j) * (x / det);
      cofactors.set(i, j, y);
    }
  }

  return cofactors;
}

Matrix inverse_gauss(const Matrix& M1){
  unsigned size_i = M1.get_size_i();
  unsigned size_j = M1.get_size_j();
  assert(size_i == size_j);
  Matrix M1bis = Matrix(size_i,size_j);
  Matrix inverse = Id(size_i);


  for (unsigned i = 0;i<size_i;i++)
    for (unsigned j = 0;j<size_j;j++)
      M1bis.set(i,j,M1.get(i,j)); // création d'une copie de M1 pour opération en place


  unsigned r = 0;
  for (unsigned j = 0;j<size_j;j++){
    unsigned imax = r;
    for (unsigned i = r;i<size_i;i++)//Récuperation du maximum
      if (M1bis.get(i,j) > M1bis.get(imax,j))
        imax = i;


    if (M1bis.get(imax,j) != 0){
      Matrix::scalar_t buff = M1bis.get(imax,j);
      multalphaline(M1bis,imax,1/buff);
      multalphaline(inverse,imax,1/buff);
      swapline(M1bis,imax,r);
      swapline(inverse,imax,r);
      for (unsigned i = 0;i < size_i;i++){
        if (i != r){
          buff = M1bis.get(i,j);
          addline(M1bis,i,r,-1*buff);
          addline(inverse,i,r,-1*buff);
        }
      }
      r++;
    }
  }
  return inverse;
}
