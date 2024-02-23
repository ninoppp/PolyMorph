#include "lib/Eigen/Dense"
#include "lib/Eigen/Sparse"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

double k = 0.1; // degradation coefficient
double D = 1; // diffusion coefficient
int N = 10; // grid size
int M = N*N; // number of grid points
double dx = 0.1;
double neighbour_contribution = D/(dx*dx);

VectorXd U;

int main()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;

  U = VectorXd(N*N); // solution vector
  MatrixXd A = MatrixXd(M, M); // systems matrix

  // fill A
  for(int i = 0; i < M; i++) {
    A(i, i) += -4*D/(dx*dx) - k;
    if(i != 0) {
      A(i-1, i) += neighbour_contribution;
      A(i, i-1) += neighbour_contribution;
    }
    if (i != M-1) {
      A(i+1, i) += neighbour_contribution;
      A(i, i+1) += neighbour_contribution;
    }
  }
  
  //std::cout << A << std::endl;

  VectorXd B = VectorXd(M); // right hand side
  for(int i = 0; i < N; i++) {  // first row
    B(i) = 1;
  }

  //std::cout << B << std::endl;
  
  // solve
  U = A.colPivHouseholderQr().solve(B);
  MatrixXd U_matrix = MatrixXd(N, N);
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      U_matrix(i, j) = U(i*N + j);
    }
  } 

  std::cout << U_matrix << std::endl;
}

double R(double u) {
  return -k*u;
}

double &u(int i, int j) {
  return U[i*N + j];
}

