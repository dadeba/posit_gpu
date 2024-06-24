/* LICENSE HEADER MANAGED BY add-license-header
 *
 *    Copyright 2024  Naohito Nakasato
 *    All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AN
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
 */

#include <algorithm>
#include <cmath>
#include <cfloat>
#include "softposit_cpp.h"

double e_time(void);
double N1_norm_diff(float *cc, float *dd, const int n);
double N1_norm_diff(double *cc, float *dd, const int n);
double N1_norm_diff_L(float *cc, float *dd, const int n);
double N1_norm_diff_L(double *cc, float *dd, const int n);
void Dsymm(double *cc, const int n);
void Ssymm(float *cc, const int n);
void Psymm(posit32 *cc, const int n);
void dump(double *cc, const int n);
double L2_norm(double  *cc, const int n);
double L2_norm(float   *cc, const int n);
double L2_norm(posit32 *cc, const int n);
double gen_po2(double *AA_d, int n, double factor, bool normal);
double gen_po(double *AA_d, int n, double factor, bool normal);
double gen_ge(double *AA_d, int n, double factor, bool normal);
double gen_po(double *AA_d, int n, double factor);
double gen_po2(double *AA_d, int n, double factor);
double gen_ge(double *AA_d, int n, double factor);

void gen_rand(double *AA_d, int n, int ld, double factor);
void gen_rand(double *AA_d, int n, int ld, double factor, bool normal);

void transpose_naive(double *, const int);
void transpose_cache(double *, const int);


void initOpenCL(void);

extern "C" {
  void spotrf_ (char * UPLO, int* N, float*  A, int* lda, int* INFO );
  void spotrs_ (char * UPLO, int* N, int* NRHS, float* A, int* lda, float* B, int* ldb, int* INFO );
  void sgemv_  (char * TRANS, int* M, int* N, float* alpha, float* A, int* lda,
		float* X, int* INCX, float* beta, float* B, int* INCY);
  void sgetrf_ (int* M, int* N, float* A, int* lda, int* IPIV, int* INFO );
  void sgetrs_ (char * TRANS, int* N, int* NRHS, float * A, int* lda, int * IPIV, float * B, int* ldb, int* INFO );
  
  void dpotrf_ (char * UPLO, int* N, double* A, int* lda, int* INFO );
  void dpotrs_ (char * UPLO, int* N, int* NRHS, double* A, int* lda, double* B, int* ldb, int* INFO );
  void dgemv_  (char * TRANS, int* M, int* N, double* alpha, double* A, int* lda,
		double* X, int* INCX, double* beta, double* B, int* INCY);
  void dgetrf_ (int* M, int* N, double* A, int* lda, int* IPIV, int* INFO );
  void dgetrs_ (char * TRANS, int* N, int* NRHS, double* A, int* lda, int * IPIV, double* B, int* ldb, int* INFO );

  void dgemm_(char * TRANSA, char * TRANSB, int* M, int* N, int *K, double *alpha, double* A, int* lda,
	      double* b, int* ldb, double *beta, double *c, int *ldc);
  void dgecon_(char * NORM, int *N, double * A, int *lda, double *ANORM, double *RCOND, double *work, int *iwork, int *INFO);
  double dlange_(char * NORM, int *M, int *N, double * A, int *lda, double *work);

  void dpocon_(char * UPLO, int *N, double * A, int *lda, double *ANORM, double *RCOND, double *work, int *iwork, int *INFO);
};

#ifndef ___MPLAPACK_BUILD_WITH__POSIT__
#define ___MPLAPACK_BUILD_WITH__POSIT__
#endif
 
