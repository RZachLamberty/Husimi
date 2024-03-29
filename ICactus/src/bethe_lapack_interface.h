#ifndef BETHE_LAPACK_INTERFACE_HEADER
#define BETHE_LAPACK_INTERFACE_HEADER

extern "C" void dscal_(const int *, const double *, double *, const int *);
extern "C" double ddot_(const int *, const double *, const int *, const double *, const int *);
extern "C" void daxpy_(const int *, const double *, const double *, const int *, double *, const int *);
extern "C" void dgemm_(const char *, const char *, const int *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);
extern "C" void dgeev_(char *, char *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
extern "C" void dggev_(char *, char *, int *, double *, int *, double *, int *, double *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
extern "C" void dgesvd_(char *, char *, int *, int *, double *, int *, double *, double *, int *, double *, int *, double *, int *, int *);
extern "C" void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
extern "C" void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *);
extern "C" void dgetrf_(int *, int *, double *, int *, int * , int *);
extern "C" void dgetri_(int *, double *, int *, int * , double *, int *, int *);

inline void dscal(const int n, const double a, double *x, const int incx)
{
    dscal_(&n, &a, x, &incx);
}


inline double ddot(const int n, const double *x, const int incx, const double *y, const int incy)
{
    return ddot_(&n, x, &incx, y, &incy);
}


inline void daxpy(const int n, const double a, const double *x, const int incx, double *y, const int incy)
{
    daxpy_(&n, &a, x, &incx, y, &incy);
}

// Matrix multiplication  AB = C
inline void dgemm(const char transa, const char transb, const int m, 
                    const int n, const int k,
                    const double alpha, const double *a, 
                    const int lda, const double *b, 
                    const int ldb,
                    const double beta, double *c, const int ldc)
{
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void dgeev(char jobvl, char jobvr, int n, double *a, int lda,
                double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr,
                double *work, int lwork, int & info)
{
    dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
}

// Added by HJC for solving generalized eigen equation
// A x= lambda B x
// (alphaR(j)+alphaI(j)*i)/Beta(j) are the eigenvalues
// http://www.netlib.org/lapack/double/dggev.f for more detailed definitions

inline void dggev(char jobvl, char jobvr, int n, double *a, int lda,
                double *b, int ldb, double *alphar, double *alphai, double *beta, 
                double *vl, int ldvl, double *vr, int ldvr,
                double *work, int lwork, int & info)
{
    dggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
}


inline void dgesvd(char jobu, char jobvt, int m, int n, double *a, int lda,
                 double *s, double *u, int ldu, double *vt, int ldvt,
                 double *work, int lwork, int & info)
{
    dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
}

inline void dgesv(int n, int nrhs, double *a, int lda,
                int *ipiv, double *b, int ldb, int & info)
{
    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
}

inline void dsyev(char jobz, char uplo, int n, double *a, int lda,
                double * w, double * work, int lwork, int & info)
{
    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
}

inline void dgetrf(int m, int n, double *a, int lda, int * ipiv, int & info)
{
    dgetrf_(&m, &n, a, &lda, ipiv, &info);
}

// Finding inverse from LU factorized matrix  
inline void dgetri(int n, double *a, int lda, int * ipiv, double *work, int lwork, int & info)
{
    dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
}

#endif
