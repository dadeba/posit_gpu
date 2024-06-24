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
  using Clock = std::chrono::high_resolution_clock;
  using std::chrono::nanoseconds;
  using std::chrono::duration_cast;
  
  //  int n_test[] = {100, 200, 300, 400, 500,600,700,800,900,1000, 1300, 1600, 2100, 2559};
  //  int n_test[] = {300, 400, 500, 700, 900, 1119, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 8000};
  //  int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1119, 1400, 1700, 2200, 2500};
  //  int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1200, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};
  //int n_test[] = {100, 200, 300, 400, 500, 700, 900, 1200, 1400, 1700, 2200, 2500, 3000, 4000, 5000, 6000, 7000, 8000};
  //  int n_test[] = {8};
  /*
  int n_test[] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,
		  1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,3000,3500,4000,4500,5000,6000,7000,8000,
		  15000, 20000};
  */
  int n_test[] = {32, 64, 128, 1024, 2048, 4096, 8192, 8192*2};

  for(long unsigned int tt = 0; tt < sizeof(n_test)/sizeof(int); tt++) {
    int n = n_test[tt];  
    size_t size_in_byte = n*n*sizeof(float)*2; // need to allocate a large buffer
    
    double *aa_d, *bb_d, *cc_d;
    aa_d = (double *)memalign(64, size_in_byte*2);
    bb_d = (double *)memalign(64, size_in_byte*2);
    cc_d = (double *)memalign(64, size_in_byte*2);
    
    double rcond_a = gen_ge(aa_d, n, 1.0);

    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	aa_d[j*n+i] = i;
      }
    }
    
    memcpy(bb_d, aa_d, sizeof(double)*n*n);

    auto st = Clock::now();
    auto en = Clock::now();
    
    st = Clock::now();
    transpose_cache(bb_d, n);
    transpose_cache(bb_d, n);
    transpose_cache(bb_d, n);
    en = Clock::now();
    auto time_cache = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    st = Clock::now();
    transpose_naive(aa_d, n);
    transpose_naive(aa_d, n);
    transpose_naive(aa_d, n);    
    en = Clock::now();
    auto time_naive = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    /*
    dump(aa_d, n);
    puts("");
    dump(bb_d, n);
    */
    
    for(int j = 0; j < n; j++) {
      for(int i = 0; i < n; i++) {
	cc_d[j*n+i] = aa_d[j*n + i] - bb_d[j*n + i];
      }
    }
    double norm_trans = L2_norm(cc_d, n);

    std::cout << n << "\t" << time_naive << "\t" << time_cache << "\t" << norm_trans;

    
    std::cout << std::endl;
    std::cout << std::flush;
    
    free(aa_d);
    free(bb_d);
    free(cc_d);
  }

  return 1;
}
