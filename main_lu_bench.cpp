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

  //int n_test[] = {900, 1200, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};
  //  int n_test[] = {300, 600, 900, 1200, 1500, 1800, 2100, 2500, 3000, 4000, 5000,
  //		  6000, 7000, 8000, 10000, 12000, 14000, 16000};

  //  int n_test[] = {10000,15000,20000};
  int n_test[] = {8000};
    
  //  int n_test[] = {300, 600, 900, 1200, 1500, 1800, 2100, 2500, 3000, 4000};

  for(long unsigned int tt = 0; tt < sizeof(n_test)/sizeof(int); tt++) {
    n = n_test[tt];  
    size_t size_in_byte = n*n*sizeof(float)*2; // need to allocate a large buffer
  
    REAL *aa, *bb, *AA;
    bb = (REAL *)memalign(64, size_in_byte);
    
    double *aa_d, *AA_d, *xx_d;
    aa_d = (double *)memalign(64, size_in_byte*2);

    double rcond_a = gen_ge(aa_d, n, factor, true);
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	bb[j*n+i] = aa_d[j*n + i];
      }
    }
    
    int ipiv[10000];
    int err;

    auto st = Clock::now();
    Rgetrf_opencl(n, n, bb, n, ipiv, err, nb);
    auto en = Clock::now();
    auto time_ocl = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    double ops = (2.0/3.0)*pow((double)n, 3.0);
    std::cout << n << "\t" << time_ocl << "\t" << ops/time_ocl/1.0e9 << "\t";

    std::cout << std::endl;
    std::cout << std::flush;
    
    free(bb);
    free(aa_d);
  }

  return 1;
}
