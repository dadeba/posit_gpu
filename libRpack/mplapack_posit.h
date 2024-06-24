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

#ifndef _MPLAPACK_POSIT_H_
#define _MPLAPACK_POSIT_H_

#include "softposit_cpp.h"
#include <cfloat> // dummy for FLT_MIN

void Rpotrf(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info);
void Rpotrf2(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info);
bool Risnan(REAL const din);
void Rpotrs(const char *uplo, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, INTEGER &info);

void Rpotrf_opencl(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info);
void Rpotrf2_opencl(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info);

void Rgetrs(const char *trans, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, INTEGER *ipiv, REAL *b, INTEGER const ldb, INTEGER &info);

void Rgetrf(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info);
void Rgetrf_opencl(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info, INTEGER);

void Rgetrf2(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info);
void Rgetrf2_opencl(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info);

void Rlaswp(INTEGER const n, REAL *a, INTEGER const lda, INTEGER const k1, INTEGER const k2, INTEGER *ipiv, INTEGER const incx);
INTEGER iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4);
INTEGER iMieeeck(INTEGER const &ispec, REAL const &zero, REAL const &one);
INTEGER iMparmq(INTEGER const ispec, const char *name, const char *opts, INTEGER const n, INTEGER const ilo, INTEGER const ihi, INTEGER const lwork);
#endif
