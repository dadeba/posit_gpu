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

/*============================================================================

This C source file is part of the SoftPosit Posit Arithmetic Package
by S. H. Leong (Cerlane).

Copyright 2017, 2018 A*STAR.  All rights reserved.

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3d, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions, and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions, and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the University nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

/*
   ushort -- 2
   uint   -- 4
   ulong  -- 8
*/

#define signP32UI( a ) ((bool) ((uint) (a)>>31))
#define signregP32UI( a ) ((bool) (((uint) (a)>>30) & 0x1))
#define packToP32UI(regime, expA, fracA) ( (uint) regime + (uint) expA + ((uint)(fracA)) )

typedef struct { uint v; } posit32_t;
union ui32_p32 { uint ui; posit32_t p; };

posit32_t p32_add(posit32_t, posit32_t);
posit32_t p32_mul(posit32_t, posit32_t);

posit32_t softposit_addMagsP32( uint uiA, uint uiB ) {
	ushort regA, regB;
        ulong frac64A=0, frac64B=0;
        uint fracA=0, regime, tmp;
        bool sign, regSA, regSB, rcarry=0, bitNPlusOne=0, bitsMore=0;
        char kA=0;
        int expA;
        short shiftRight;
        union ui32_p32 uZ;

        sign = signP32UI( uiA );
        if (sign){
                uiA = -uiA & 0xFFFFFFFF;
                uiB = -uiB & 0xFFFFFFFF;
        }

        if ((int)uiA < (int)uiB){
                uiA ^= uiB;
                uiB ^= uiA;
                uiA ^= uiB;
        }
        regSA = signregP32UI( uiA );
	regSB = signregP32UI( uiB );

	tmp = (uiA<<2)&0xFFFFFFFF;
        if (regSA){
                while (tmp>>31){
                        kA++;
                        tmp= (tmp<<1) & 0xFFFFFFFF;
                }
        }
        else{
                kA=-1;
                while (!(tmp>>31)){
                        kA--;
                        tmp= (tmp<<1) & 0xFFFFFFFF;
                }
                tmp&=0x7FFFFFFF;
        }

        expA = tmp>>29; //to get 2 bits
      //frac64A = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	frac64A = ((0x40000000UL | tmp<<1) & 0x7FFFFFFFUL) <<32;
        shiftRight = kA;

        tmp = (uiB<<2) & 0xFFFFFFFF;
        if (regSB){
                while (tmp>>31){
                        shiftRight--;
                        tmp= (tmp<<1) & 0xFFFFFFFF;
                }
        }
        else{
                shiftRight++;
                while (!(tmp>>31)){
                        shiftRight++;
                        tmp= (tmp<<1) & 0xFFFFFFFF;
                }
                tmp&=0x7FFFFFFF;
        }
        //frac64B = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	frac64B = ((0x40000000UL | tmp<<1) & 0x7FFFFFFFUL) <<32;
        //This is 4kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
        shiftRight = (shiftRight<<2) + expA - (tmp>>29);

        //Manage CLANG (LLVM) compiler when shifting right more than number of bits
        (shiftRight>63) ? (frac64B=0): (frac64B >>= shiftRight); //frac64B >>= shiftRight

        frac64A += frac64B;

        rcarry = 0x8000000000000000 & frac64A; //first left bit
        if (rcarry){
                expA++;
                if (expA>3){
                        kA ++;
                        expA&=0x3;
                }
                frac64A>>=1;
        }
        if(kA<0){
                regA = -kA;
                regSA = 0;
                regime = 0x40000000>>regA;
        }
        else{
                regA = kA+1;
                regSA=1;
                regime = 0x7FFFFFFF - (0x7FFFFFFF>>regA);
        }

        if(regA>30){
                //max or min pos. exp and frac does not matter.
                (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
        }
        else{
                //remove hidden bits
                frac64A = (frac64A & 0x3FFFFFFFFFFFFFFF) >>(regA + 2) ; // 2 bits exp

                fracA = frac64A>>32;

                if (regA<=28){
                        bitNPlusOne |= (0x80000000 & frac64A) ;
                        expA <<= (28-regA);
                }
                else {
                        if (regA==30){
                                bitNPlusOne = expA&0x2;
                                bitsMore = (expA&0x1);
                                expA = 0;
                        }
                        else if (regA==29){
                                bitNPlusOne = expA&0x1;
                                expA>>=1;
                        }
                        if (fracA>0){
                                fracA=0;
                                bitsMore =1;
                        }
                }

                uZ.ui = packToP32UI(regime, expA, fracA);
                //n+1 frac bit is 1. Need to check if another bit is 1 too if not round to even
                if (bitNPlusOne){
                        (0x7FFFFFFF & frac64A) ? (bitsMore=1) : (bitsMore=0);
                        uZ.ui += (uZ.ui&1) | bitsMore;
                }
        }
        if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
        return uZ.p;
}

posit32_t softposit_subMagsP32( ulong uiA, ulong uiB ) {
	ushort regA, regB;
	ulong frac64A=0, frac64B=0;
	uint fracA=0, regime, tmp;
	bool sign, regSA, regSB, ecarry=0, bitNPlusOne=0, bitsMore=0;
	char kA=0;
	int expA=0;
	short shiftRight;
	union ui32_p32 uZ;

	sign = signP32UI( uiA );
	if (sign)
		uiA = -uiA & 0xFFFFFFFF;
	else
		uiB = -uiB & 0xFFFFFFFF;

	if (uiA==uiB){ //essential, if not need special handling
		uZ.ui = 0;
		return uZ.p;
	}
	if ((int)uiA < (int)uiB){
		uiA ^= uiB;
		uiB ^= uiA;
		uiA ^= uiB;
		(sign) ? (sign = 0 ) : (sign=1); //A becomes B
	}
	regSA = signregP32UI( uiA );
	regSB = signregP32UI( uiB );

	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}

	expA = tmp>>29; //to get 2 bits
	//frac64A = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	frac64A = ((0x40000000UL | tmp<<1) & 0x7FFFFFFFUL) <<32;
	shiftRight = kA;


	tmp = (uiB<<2) & 0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}

	}
	else{
		shiftRight++;
		while (!(tmp>>31)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;

	}
	//frac64B = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	frac64B = ((0x40000000UL | tmp<<1) & 0x7FFFFFFFUL) <<32;

	//This is 4kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
	shiftRight = (shiftRight<<2) + expA - (tmp>>29);
	if (shiftRight>63){
		uZ.ui = uiA;
		if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
		return uZ.p;
	}
	else
		(frac64B >>= shiftRight);

	frac64A -= frac64B;

	while((frac64A>>59)==0){
		kA--;
		frac64A<<=4;
	}
	ecarry = (0x4000000000000000 & frac64A);//(0x4000000000000000 & frac64A)>>62;
	while (!ecarry){
		if (expA==0){
			kA--;
			expA=3;
		}
		else
			expA--;
		frac64A<<=1;
		ecarry = (0x4000000000000000 & frac64A);
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFF - (0x7FFFFFFF>>regA);
	}
	if(regA>30){
		//max or min pos. exp and frac does not matter.
		(regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
	}
	else{
		//remove hidden bits
		frac64A = (frac64A & 0x3FFFFFFFFFFFFFFF) >>(regA + 2) ; // 2 bits exp

		fracA = frac64A>>32;

		if (regA<=28){
			bitNPlusOne |= (0x80000000 & frac64A) ;
			expA <<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1;
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}

		}

		uZ.ui = packToP32UI(regime, expA, fracA);
		//n+1 frac bit is 1. Need to check if another bit is 1 too if not round to even
		if (bitNPlusOne){
			(0x7FFFFFFF & frac64A) ? (bitsMore=1) : (bitsMore=0);
			uZ.ui += (uZ.ui&1) | bitsMore;
		}
	}
	if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}

posit32_t p32_add( posit32_t a, posit32_t b ){
    union ui32_p32 uA, uB, uZ;
    ulong uiA, uiB;

    uA.p = a;
    uiA = uA.ui;
    uB.p = b;
    uiB = uB.ui;

    //Zero or infinity
    if (uiA==0 || uiB==0){ // Not required but put here for speed
       uZ.ui = uiA | uiB;

       return uZ.p;
    }
    else if ( uiA==0x80000000 || uiB==0x80000000 ){
	//printf("in infinity\n");
	uZ.ui = 0x80000000;

	return uZ.p;
    }

    //different signs
    if ((uiA^uiB)>>31)
    	return softposit_subMagsP32(uiA, uiB);
    else
	return softposit_addMagsP32(uiA, uiB);
}

posit32_t p32_mul( posit32_t pA, posit32_t pB ){
	union ui32_p32 uA, uB, uZ;
	uint uiA, uiB;
	uint regA, fracA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int expA;
	char kA=0;
	ulong frac64Z;

	uA.p = pA;
	uiA = uA.ui;
	uB.p = pB;
	uiB = uB.ui;

	//NaR or Zero
	if ( uiA==0x80000000 || uiB==0x80000000 ){
		uZ.ui = 0x80000000;
		return uZ.p;
	}
	else if (uiA==0 || uiB==0){
		uZ.ui = 0;
		return uZ.p;
	}

	signA = signP32UI( uiA );
	signB = signP32UI( uiB );
	signZ = signA ^ signB;

	if(signA) uiA = (-uiA & 0xFFFFFFFF);
	if(signB) uiB = (-uiB & 0xFFFFFFFF);

	regSA = signregP32UI(uiA);
	regSB = signregP32UI(uiB);


	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){

		while (tmp>>31){

			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA = tmp>>29; //to get 2 bits
	fracA = ((tmp<<1) | 0x40000000) & 0x7FFFFFFF;

	tmp = (uiB<<2)&0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA--;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA += tmp>>29;
	frac64Z = (ulong) fracA * (((tmp<<1) | 0x40000000) & 0x7FFFFFFF);

	if (expA>3){
		kA++;
		expA&=0x3; // -=4
	}

	rcarry = frac64Z>>61;//3rd bit of frac64Z
	if (rcarry){
		expA++;
		if (expA>3){
			kA ++;
			expA&=0x3;
		}
		frac64Z>>=1;
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFF - (0x7FFFFFFF>>regA);
	}

	if(regA>30){
		//max or min pos. exp and frac does not matter.
		(regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
	}
	else{
		//remove carry and rcarry bits and shift to correct position (2 bits exp, so + 1 than 16 bits)
		frac64Z = (frac64Z&0xFFFFFFFFFFFFFFF) >> regA;
		fracA = (uint) (frac64Z>>32);

		if (regA<=28){
			bitNPlusOne |= (0x80000000 & frac64Z);
			expA<<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1; //taken care of by the pack algo
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}

		}
		//sign is always zero
		uZ.ui = packToP32UI(regime, expA, fracA);
		//n+1 frac bit is 1. Need to check if another bit is 1 too if not round to even
		if (bitNPlusOne){
			(0x7FFFFFFF & frac64Z) ? (bitsMore=1) : (bitsMore=0);
			uZ.ui += (uZ.ui&1) | bitsMore;
		}
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}

// blocking https://cnugteren.github.io/tutorial/pages/page4.html

#define BLOCKA As[TS][TS]
#define BLOCKB Bs[TS][TS]

#define STOREBLOCK { As[l_i][l_j] = aa; Bs[l_j][l_i] = bb; }
#define LOADBLOCK  { aa = As[l_i][p];  bb = Bs[l_j][p]; }

//#define STOREBLOCK { As[l_j][l_i] = aa; Bs[l_j][l_i] = bb; }
//#define LOADBLOCK  { aa = As[p][l_i];  bb = Bs[l_j][p]; }

//            i   j
// Matrix A : m * k  
// Matrix B : k * n
// Matrix C : m * n

// Column major layout
//
// if A is not transposed : shape m-by-k, buffer lda*k
// if A is transposed     : shape k-by-m, buffer lda*m
//
// if B is not transposed : shape k-by-n, buffer ldb*n
// if B is transposed     : shape n-by-k, buffer ldb*k
//
// C : shape m-by-n, buffer ldc*n
//

// m x k
#define LOAD_A_N { aa.v = ((i   < m) && (t_j < k)) ? A[i   + t_j*lda] : 0x0; }
// k x m
#define LOAD_A_T { aa.v = ((i   < m) && (t_j < k)) ? A[t_j + i*lda] : 0x0; }

// k x n
#define LOAD_B_N { bb.v = ((t_i < k) && (  j < n)) ? B[t_i + j*ldb] : 0x0; }
// n x k
#define LOAD_B_T { bb.v = ((t_i < k) && (  j < n)) ? B[j + t_i*ldb] : 0x0; }

//#define STORERESULT { if (i < m && j < n) { C[i + j * ldc] = temp.v; } }

#define STORERESULT { if (i < m && j < n) { C[i + j * ldc] = store_al_be(C[i + j * ldc], alpha, beta, temp); } }

uint store_al_be(const uint C, const uint alpha, const uint beta, posit32_t res)
{
  posit32_t Ctmp, altmp, betmp;

  Ctmp.v = C;
  altmp.v = alpha;
  betmp.v = beta;

  posit32_t tmp1, tmp2;

  tmp1 = p32_mul(res,  altmp);
  tmp2 = p32_mul(Ctmp, betmp);
  res = p32_add(tmp1, tmp2);

  return res.v;
}

__kernel void gemm_NN_X(const int m, const int n, const int k,
        	        const uint alpha, 
			__global uint *restrict A, const int lda,
			__global uint *restrict B, const int ldb,
			const uint beta, 
			__global uint *restrict C, const int ldc
		      )
{
  const int l_i = get_local_id(0);  // row
  const int l_j = get_local_id(1);  // col
  const int i = TS*get_group_id(0) + l_i; //const int i = get_global_id(0);
  const int j = TS*get_group_id(1) + l_j; //const int j = get_global_id(1);

  __local posit32_t BLOCKA;
  __local posit32_t BLOCKB;

  posit32_t temp;
  temp.v = 0x0;
  
  int numtiles = k/TS;
  if (k % TS != 0) numtiles++;
  for (int t = 0; t < numtiles; t++) {
    posit32_t aa, bb;

    const int t_i = TS*t + l_i;
    const int t_j = TS*t + l_j;

    LOAD_A_N;
    LOAD_B_N;

    STOREBLOCK;
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int p = 0; p < TS; p++) {
      posit32_t aa, bb, tmp;

      LOADBLOCK;

      tmp = p32_mul(aa, bb);
      temp = p32_add(tmp, temp);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  STORERESULT;
  return;
}

__kernel void gemm_TN_X(const int m, const int n, const int k,
        	        const uint alpha, 
			__global uint *restrict A, const int lda,
			__global uint *restrict B, const int ldb,
			const uint beta, 
			__global uint *restrict C, const int ldc
		      )
{
  const int l_i = get_local_id(0);  // row
  const int l_j = get_local_id(1);  // col
  const int i = TS*get_group_id(0) + l_i; //const int i = get_global_id(0);
  const int j = TS*get_group_id(1) + l_j; //const int j = get_global_id(1);

  __local posit32_t BLOCKA;
  __local posit32_t BLOCKB;

  posit32_t temp;
  temp.v = 0x0;

  int numtiles = k/TS;
  if (k % TS != 0) numtiles++;
  for (int t = 0; t < numtiles; t++) {
    posit32_t aa, bb;

    const int t_i = TS*t + l_i;
    const int t_j = TS*t + l_j;

    LOAD_A_T;
    LOAD_B_N;

    STOREBLOCK;
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int p = 0; p < TS; p++) {
      posit32_t aa, bb, tmp;

      LOADBLOCK;

      tmp = p32_mul(aa, bb);
      temp = p32_add(tmp, temp);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  STORERESULT;
  return;
}

__kernel void gemm_NT_X(const int m, const int n, const int k,
        	        const uint alpha, 
			__global uint *restrict A, const int lda,
			__global uint *restrict B, const int ldb,
			const uint beta, 
			__global uint *restrict C, const int ldc
		      )
{
  const int l_i = get_local_id(0);  // row
  const int l_j = get_local_id(1);  // col
  const int i = TS*get_group_id(0) + l_i; //const int i = get_global_id(0);
  const int j = TS*get_group_id(1) + l_j; //const int j = get_global_id(1);

  __local posit32_t BLOCKA;
  __local posit32_t BLOCKB;

  posit32_t temp;
  temp.v = 0x0;

  int numtiles = k/TS;
  if (k % TS != 0) numtiles++;
  for (int t = 0; t < numtiles; t++) {
    posit32_t aa, bb;

    const int t_i = TS*t + l_i;
    const int t_j = TS*t + l_j;

    LOAD_A_N;
    LOAD_B_T;

    STOREBLOCK;  
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int p = 0; p < TS; p++) {
      posit32_t aa, bb, tmp;

      LOADBLOCK;
      
      tmp = p32_mul(aa, bb);
      temp = p32_add(tmp, temp);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  STORERESULT;
  return;
}

__kernel void gemm_TT_X(const int m, const int n, const int k,
        	        const uint alpha, 
			__global uint *restrict A, const int lda,
			__global uint *restrict B, const int ldb,
			const uint beta, 
			__global uint *restrict C, const int ldc
		      )
{
  const int l_i = get_local_id(0);  // row
  const int l_j = get_local_id(1);  // col
  const int i = TS*get_group_id(0) + l_i; //const int i = get_global_id(0);
  const int j = TS*get_group_id(1) + l_j; //const int j = get_global_id(1);

  __local posit32_t BLOCKA;
  __local posit32_t BLOCKB;

  posit32_t temp;
  temp.v = 0x0;

  int numtiles = k/TS;
  if (k % TS != 0) numtiles++;
  for (int t = 0; t < numtiles; t++) {
    posit32_t aa, bb;

    const int t_i = TS*t + l_i;
    const int t_j = TS*t + l_j;

    LOAD_A_T;
    LOAD_B_T;

    STOREBLOCK; 
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int p = 0; p < TS; p++) {
      posit32_t aa, bb, tmp;

      LOADBLOCK;

      tmp = p32_mul(aa, bb);
      temp = p32_add(tmp, temp);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  STORERESULT;
  return;
}

__kernel void gemm_TN(const int m, const int n, const int k,
		      const uint alpha, 
		      __global uint *restrict A, const int lda,
		      __global uint *restrict B, const int ldb,
		      const uint beta, 
		      __global uint *restrict C, const int ldc
		      )
{
  int i, j, l;

  i = get_global_id(0);
  j = get_global_id(1);

  if (i > m-1) return;
  if (j > n-1) return;
  
  // Rgemm_TN_omp.cpp from mpackapck
  // C = A^t B, ignore ALPHA, BETA
  posit32_t temp;
  temp.v = 0x0;

  for (l = 0; l < k; l++) {
    posit32_t aa, bb, tmp;

    aa.v = A[l + i * lda];
    bb.v = B[l + j * ldb];

    tmp = p32_mul(aa, bb);
    temp = p32_add(tmp, temp);
  }

  C[i + j * ldc] = temp.v;
}
