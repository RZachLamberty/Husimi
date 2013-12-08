#ifndef MATRIX_FUNCTIONS_H
#define MATRIX_FUNCTIONS_H

#include"global.h"
#include"matrix.h"

using namespace std;

// Diagonalize A and report eigenvalues and orthogonal eigenvecs

void general_real_diagonalize
                (Matrix &a, 
                 std::vector< complex<double> > &eigs, 
                 Matrix &eigenvecs);

// Diagonalize symmetric matrix A and report eigenvalues and orthogonal eigenvecs

void symmetric_diagonalize
                (Matrix &a, 
                 std::vector<double> &eigs, 
                 Matrix &eigenvecs);

// Orthogonalize eigenvectors

void gram_schmidt(Matrix &eigenvecs, std::vector< complex<double> > &eigs);

void eta_orthogonalize( std::vector< std::vector<double> > &eigenvecs,
			std::vector<double> const &eta);

// C=A*B

void real_matrix_multiply( Matrix &a, 
                           Matrix &b, 
                           Matrix &c);

void perform_lu(Matrix &mat, Matrix &lu, std::vector<int> &ipiv);
void invert_real_matrix(Matrix &mat, Matrix &inverse);

// Multiply matrix times a vector

void real_matrix_times_vector( Matrix &a,
                              std::vector<double> &b,
                              std::vector<double> &c);

// Multiply vector times vector
void real_vector_times_vector(std::vector<double> &a,
                              std::vector<double> &b,
                              Matrix &c);

// Singular Value decomposition

void svd( Matrix &a, 
          std::vector<double> &eigs, 
          Matrix &u, Matrix &vt);

void matrix_lanczos(int 			      		iterations,
		    int 			      		how_many_eigenvecs, 
                    Matrix 					&a,
		    std::vector<double> 			&eigs,
		    Matrix 					&eigenvecs,
	            bool 					ipr);
#endif
