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

#include <iostream>
#include <random>
#include <chrono>
#include <cstring>
#include <malloc.h>
#include <unistd.h>

#include "def.h"
#include "mpblas.h"
#include "mplapack.h"

int main(int narg, char **argv)
{
  double factor = 1.0;

  using Clock = std::chrono::high_resolution_clock;
  using std::chrono::nanoseconds;
  using std::chrono::duration_cast;
  
  int n;

  if (narg == 2) factor = atof(argv[1]);

  int nb = 32;
  if (narg == 3) nb = atoi(argv[2]);
  
#ifdef BENCH_CPU
#else
  initOpenCL();
#endif

  //  int n_test[] = {100, 200, 300, 400, 500,600,700,800,900,1000, 1300, 1600, 2100, 2559};
  //  int n_test[] = {300, 400, 500, 700, 900, 1119, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 8000};
  //  int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1119, 1400, 1700, 2200, 2500};
  //  int n_test[] = {16*100};
  //  int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1200, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};

  //  int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1000, 1200, 1400, 1700, 2000, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};

  //int n_test[] = {900, 1200, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};
  //  int n_test[] = {300, 600, 900, 1200, 1500, 1800, 2100, 2500, 3000, 4000, 5000, 6000, 7000, 8000};

  //  int n_test[] = {1500, 1800, 2100, 2500, 3000, 4000, 5000, 6000, 7000, 8000};

  //  int n_test[] = {10000,15000,20000};
    
  //  int n_test[] = {300, 600, 900, 1200, 1500, 1800, 2100, 2500, 3000, 4000};

  int n_test[] = {100, 1000, 2000, 4000, 8000};
  
  for(long unsigned int tt = 0; tt < sizeof(n_test)/sizeof(int); tt++) {
    n = n_test[tt];  
    size_t size_in_byte = n*n*sizeof(float)*2; // need to allocate a large buffer
  
    REAL *aa, *bb, *AA;
    aa = (REAL *)memalign(64, size_in_byte);
    bb = (REAL *)memalign(64, size_in_byte);
    AA = (REAL *)memalign(64, size_in_byte);
    
    float *aa_f, *AA_f;
    aa_f = (float *)memalign(64, size_in_byte);
    AA_f = (float *)memalign(64, size_in_byte);

    double *aa_d, *AA_d, *xx_d;
    aa_d = (double *)memalign(64, size_in_byte*2);
    AA_d = (double *)memalign(64, size_in_byte*2);
    xx_d = (double *)memalign(64, sizeof(double)*n);

    double err_ppp, err_fff;
    double *bb_p, *bb_f;
    bb_p = (double *)memalign(64, size_in_byte*2);
    bb_f = (double *)memalign(64, size_in_byte*2);

    double rcond_a = gen_ge(aa_d, n, factor, true);
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	aa  [j*n+i] = aa_d[j*n + i];
	bb  [j*n+i] = aa_d[j*n + i];
	aa_f[j*n+i] = aa_d[j*n + i];
      }
    }
    // copy original matrix
    memcpy(AA,   aa,   sizeof(REAL)*n*n);
    memcpy(AA_f, aa_f, sizeof(float)*n*n);
    memcpy(AA_d, aa_d, sizeof(double)*n*n);
    
    int ipiv[10000];
    int ipiv_p[10000];
    int ipiv_f[10000];
    int ipiv_d[10000];
    int err;

    auto st = Clock::now();
#ifdef BENCH_CPU
    Rgetrf(n, n, aa, n, ipiv, err);
#endif
    auto en = Clock::now();
    auto time_cpu = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    
#ifndef BENCH_CPU
    st = Clock::now();
    Rgetrf_opencl(n, n, bb, n, ipiv_p, err, nb);
    en = Clock::now();
    auto time_ocl = duration_cast<nanoseconds>(en-st).count()/1.0e9;
#endif
    
    st = Clock::now();
    sgetrf_(&n, &n, aa_f, &n, ipiv_f, &err);
    en = Clock::now();
    auto time_b32 = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    st = Clock::now();
    dgetrf_(&n, &n, aa_d, &n, ipiv_d, &err);
    en = Clock::now();
    auto time_b64 = duration_cast<nanoseconds>(en-st).count()/1.0e9;

#ifndef BENCH_CPU
    for(int i = 0; i < n; i++) xx_d[i] = 1.0/sqrt((double)n);
	  
    // err_p for POSIT using Rgemm_opencl : bb
    int nrhs = 1;
    int ldb = n;
    int inc = 1;
    double err_p = 0.0;
    {
      REAL *b  = (REAL *)memalign(64, sizeof(REAL)*n);
      REAL *b_ = (REAL *)memalign(64, sizeof(REAL)*n);
      REAL *x  = (REAL *)memalign(64, sizeof(REAL)*n);
      for(int i = 0; i < n; i++) b[i] = 0.0;
      for(int i = 0; i < n; i++) x[i] = xx_d[i];

      REAL alpha = 1.0;
      REAL beta  = 0.0;
      // compute vector b
      Rgemv("N", n, n, alpha, AA, n, x, inc, beta, b, inc);
      // copy vector
      memcpy(x, b, sizeof(REAL)*n);

      // solve Ax = b
      Rgetrs("N", n, nrhs, bb, n, ipiv_p, x, ldb, err);
      // compute b_' = Ax
      Rgemv("N", n, n, alpha, AA, n, x, inc, beta, b_, inc);

      // b - b'
      for(int i = 0; i < n; i++) {
	bb_p[i] = b_[i].toDouble();
	b_[i] = b[i] - b_[i];
      }
      // compute |b - Ax'|/|b| as L2 norm
      err_p = L2_norm(b_,n)/L2_norm(b,n);
      //      std::cout << "POSIT " << L2_norm(b_,n)/L2_norm(b,n) << "\n";
      free(b);
      free(b_);
      free(x);
    }

    double err_f = 0.0;
    {
      float *b  = (float *)memalign(64, sizeof(float)*n);
      float *b_ = (float *)memalign(64, sizeof(float)*n);
      float *x  = (float *)memalign(64, sizeof(float)*n);
      for(int i = 0; i < n; i++) b[i] = 0.0;
      for(int i = 0; i < n; i++) x[i] = xx_d[i];

      float alpha = 1.0;
      float beta  = 0.0;
      sgemv_((char *)"N", &n, &n, &alpha, AA_f, &n, x, &inc, &beta, b, &inc);
      memcpy(x, b, sizeof(float)*n);
      sgetrs_((char *)"N", &n, &nrhs, aa_f, &n, ipiv_f, x, &ldb, &err);
      sgemv_((char *)"N", &n, &n, &alpha, AA_f, &n, x, &inc, &beta, b_, &inc);
      for(int i = 0; i < n; i++) {
	bb_f[i] = (double)b_[i];
	b_[i] = b[i] - b_[i];
      }
      err_f = L2_norm(b_,n)/L2_norm(b,n);
      //      std::cout << "float " << L2_norm(b_,n)/L2_norm(b,n) << "\n";      
      free(b);
      free(b_);
      free(x);
    }

    double err_d = 0.0;
    {
      double *b  = (double *)memalign(64, sizeof(double)*n);
      double *b_ = (double *)memalign(64, sizeof(double)*n);
      double *x  = (double *)memalign(64, sizeof(double)*n);
      for(int i = 0; i < n; i++) b[i] = 0.0;
      for(int i = 0; i < n; i++) x[i] = xx_d[i];

      double alpha = 1.0;
      double beta  = 0.0;
      dgemv_((char *)"N", &n, &n, &alpha, AA_d, &n, x, &inc, &beta, b, &inc);
      memcpy(x, b, sizeof(double)*n);
      dgetrs_((char *)"N", &n, &nrhs, aa_d, &n, ipiv_d, x, &ldb, &err);
      dgemv_((char *)"N", &n, &n, &alpha, AA_d, &n, x, &inc, &beta, b_, &inc);
      for(int i = 0; i < n; i++) {
	b_[i] = b[i] - b_[i];
      }
      err_d = L2_norm(b_,n)/L2_norm(b,n);
      //      std::cout << "double " << L2_norm(b_,n)/L2_norm(b,n) << "\n";      

      for(int i = 0; i < n; i++) {
	b_[i] = b[i] - bb_p[i];
      }
      err_ppp = L2_norm(b_,n)/L2_norm(b,n);
      for(int i = 0; i < n; i++) {
	b_[i] = b[i] - bb_f[i];
      }
      err_fff = L2_norm(b_,n)/L2_norm(b,n);

      free(b);
      free(b_);
      free(x);
    }
#endif
    
    std::cout << n << " ";
#ifdef BENCH_CPU
    std::cout << time_cpu << "\t" << time_b32 << "\t" << time_b64 << "\t";
    //    std::cout << time_cpu << "\t" << time_ocl << "\t" << time_b32 << "\t" << time_b64 << "\t";
    std::cout << std::endl;
#else
    std::cout << time_ocl << "\t" << time_b32 << "\t" << time_b64 << "\t";
    //    std::cout << err_p << "\t" << err_f << "\t" << err_d << "\t" << log10(err_f/err_p) << "\n";
    std::cout << "\t" << err_ppp << "\t" << err_fff << "\t" << log10(err_fff/err_ppp)<< "\n";
#endif
    //    std::cout << err_p << "\t" << err_f << "\t" << "(" << log10(err_f/err_p) << ")" << "\t" << err_d << "\n";

    /*
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	rr_f[j*n+i] = (float)aa[j*n + i].toDouble();
      }
    }

    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	rr_f[j*n+i] = (float)aa[j*n + i].toDouble();
      }
    }

    std::cout << N1_norm_diff(aa_f, rr_f, n) << "\t" << N1_norm_diff(aa_d, rr_f, n) << "\t";
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	rr_f[j*n+i] = (float)bb[j*n + i].toDouble();
      }
    }
    std::cout << N1_norm_diff(aa_f, rr_f, n) << "\t" << N1_norm_diff(aa_d, rr_f, n) << "\n";
    */
    //    std::cout << std::endl;
    std::cout << std::flush;
    
    free(aa);
    free(bb);
    free(AA);
    free(aa_f);
    free(AA_f);
    free(aa_d);
    free(AA_d);
    free(xx_d);

    sleep(5);
  }

  return 1;
}
