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

#ifndef _MPBLAS_POSIT_H_
#define _MPBLAS_POSIT_H_

#include "softposit_cpp.h"

void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, posit32 const alpha, posit32 *a, mplapackint const lda, posit32 *b, mplapackint const ldb, posit32 const beta, posit32 *c, mplapackint const ldc);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, posit32 const alpha, posit32 *a, mplapackint const lda, posit32 *b, mplapackint const ldb);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, posit32 const alpha, posit32 *a, mplapackint const lda, posit32 const beta, posit32 *c, mplapackint const ldc);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, posit32 const alpha, posit32 *a, mplapackint const lda, posit32 *x, mplapackint const incx, posit32 const beta, posit32 *y, mplapackint const incy);

void Mxerbla(const char *srname, int info);
void Rscal(mplapackint const n, posit32 const da, posit32 *dx, mplapackint const incx);
mplapackint iRamax(mplapackint const n, posit32 *dx, mplapackint const incx);
bool Mlsame(const char *a, const char *b);

posit32 abs(posit32 x);

void RgemmG(const char *transa, const char *transb,
	    mplapackint const m, mplapackint const n, mplapackint const k,
	    posit32 const alpha,
	    posit32* A, mplapackint const lda,
	    posit32* B, mplapackint const ldb,
	    posit32 const beta,
	    posit32* C, mplapackint const ldc);
#endif
