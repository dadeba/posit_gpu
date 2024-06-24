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
  
  //  int n_test[] = {100, 200, 300, 400, 500,600,700,800,900,1000, 1300, 1600, 2100, 2559};
  //  int n_test[] = {300, 400, 500, 700, 900, 1119, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 8000};
  //  int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1119, 1400, 1700, 2200, 2500};
  //  int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1200, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};
  //int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1200, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};
  //  int n_test[] = {16*100};

#ifdef BENCH_CPU
#else
  initOpenCL();
#endif
  
  //  int n_test[] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,
  //		  1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,3000,3500,4000,4500,5000,6000,7000,8000};

  int n_test[] = {100,200,300,400,500,750,1000,1500,2000,2500,3000,3500,4000,6000,8000};

  
  for(long unsigned int tt = 0; tt < sizeof(n_test)/sizeof(int); tt++) {
    n = n_test[tt];  
    size_t size_in_byte = n*n*sizeof(float)*2; // need to allocate a large buffer
    
    REAL *aa, *bb, *cc;
    REAL *cc2, *cc3, *cc4;
    aa = (REAL *)memalign(64, size_in_byte);
    bb = (REAL *)memalign(64, size_in_byte);
    cc  = (REAL *)memalign(64, size_in_byte);
    cc2 = (REAL *)memalign(64, size_in_byte);
    cc3 = (REAL *)memalign(64, size_in_byte);
    cc4 = (REAL *)memalign(64, size_in_byte);
    
    double *aa_d, *bb_d, *cc_d;
    aa_d = (double *)memalign(64, size_in_byte*2);
    bb_d = (double *)memalign(64, size_in_byte*2);
    cc_d = (double *)memalign(64, size_in_byte*2);
    
    double rcond_a = gen_ge(aa_d, n, factor, true);
    double rcond_b = gen_ge(bb_d, n, factor, true);
    double rcond_c = gen_ge(cc_d, n, factor, true);

    double alpha_d = aa_d[0];
    double beta_d  = bb_d[0];
    //    alpha_d = beta_d = 1.0;
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	aa[j*n+i] = aa_d[j*n + i];
	bb[j*n+i] = bb_d[j*n + i];
	cc [j*n+i] = cc_d[j*n + i];
	cc2[j*n+i] = cc_d[j*n + i];
	cc3[j*n+i] = cc_d[j*n + i];
	cc4[j*n+i] = cc_d[j*n + i];
      }
    }

    char str_a[4] = "", str_b[4]= "";
    char str_aa[5] = "nntt", str_bb[5]= "ntnt";
    //    int w = tt % 4;
    int w = 0;
    sprintf(str_a, "%c",  str_aa[w]);
    sprintf(str_b, "%c",  str_bb[w]);
    //    printf("%s %s\t", str_a, str_b);

    REAL alpha_p, beta_p;
    alpha_p = alpha_d;
    beta_p  = beta_d;
    
    auto st = Clock::now();
#ifdef BENCH_CPU
    Rgemm((char *)str_a, (char *)str_b, n, n, n, alpha_p, aa, n, bb, n, beta_p, cc, n);
#endif
    auto en = Clock::now();
    auto time_cpu = duration_cast<nanoseconds>(en-st).count()/1.0e9;

#ifndef BENCH_CPU
    st = Clock::now();
    RgemmG((char *)str_a, (char *)str_b, n, n, n, alpha_p, aa, n, bb, n, beta_p, cc2, n);
    en = Clock::now();
    auto time_ocl = duration_cast<nanoseconds>(en-st).count()/1.0e9;
    st = Clock::now();
    RgemmG((char *)str_a, (char *)str_b, n, n, n, alpha_p, aa, n, bb, n, beta_p, cc3, n);
    en = Clock::now();
    time_ocl += duration_cast<nanoseconds>(en-st).count()/1.0e9;
    st = Clock::now();
    RgemmG((char *)str_a, (char *)str_b, n, n, n, alpha_p, aa, n, bb, n, beta_p, cc4, n);
    en = Clock::now();
    time_ocl += duration_cast<nanoseconds>(en-st).count()/1.0e9;
    time_ocl /= 3.0;
#endif
    
    st = Clock::now();
    dgemm_((char *)str_a, (char *)str_b, &n, &n, &n, &alpha_d, aa_d, &n, bb_d, &n, &beta_d, cc_d, &n);
    en = Clock::now();
    auto time_b64 = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    double ops = 2.0*pow((double)n, 3.0);
    
    std::cout << n << " ";
    //    std::cout << time_cpu << "\t" << time_ocl << "\t" << time_b64 << "\t";
    //    std::cout << ops/time_cpu/1.0e9 << "\t" << ops/time_ocl/1.0e9 << "\t" << ops/time_b64/1.0e9 << "\t";
    //    std::cout << "\t" << ops/time_ocl/1.0e9 << "\t" << ops/time_b64/1.0e9 << "\t";

#ifdef BENCH_CPU
    std::cout << "\t" << ops/time_cpu/1.0e9 << "\t";
#else
    std::cout << "\t" << ops/time_ocl/1.0e9 << "\t";
#endif


    double norm_b64 = L2_norm(cc_d, n);
      
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	aa_d[j*n+i] = cc_d[j*n + i] - cc [j*n + i].toDouble();
	bb_d[j*n+i] = cc_d[j*n + i] - cc2[j*n + i].toDouble();
      }
    }
    double norm_p_cpu = L2_norm(aa_d, n);
    double norm_p_ocl = L2_norm(bb_d, n);    

    //    std::cout << norm_p_cpu/norm_b64 << "\t" << norm_p_ocl/norm_b64;
    std::cout << norm_p_ocl/norm_b64;
    //    std::cout << "\t" << norm_p_ocl;
    
    std::cout << std::endl;
    std::cout << std::flush;
    
    free(aa);
    free(bb);
    free(cc);
    free(cc2);
    free(cc3);
    free(cc4);    
    free(aa_d);
    free(bb_d);
    free(cc_d);
  }

  return 1;
}
