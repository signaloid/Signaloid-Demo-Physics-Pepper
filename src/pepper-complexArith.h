/*
 *	Copyright (c) 2026, Signaloid.
 *
 *	Permission is hereby granted, free of charge, to any person obtaining a copy
 *	of this software and associated documentation files (the "Software"), to deal
 *	in the Software without restriction, including without limitation the rights
 *	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the Software is
 *	furnished to do so, subject to the following conditions:
 *
 *	The above copyright notice and this permission notice shall be included in all
 *	copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *	SOFTWARE.
 */

#pragma once

#include "pepper-vec4.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	Complex number as a pair of distributional doubles.
 */
typedef struct
{
	double  re;
	double  im;
} Complex;

/**
 *	@brief	Complex four-vector (Lorentz vector with complex components).
 *		Used for Berends-Giele gluon currents and polarization vectors.
 */
typedef struct
{
	Complex c[kFourVectorIndexMax];
} CVec4;

/**
 *	@brief	Weyl spinor with two complex components.
 */
typedef struct
{
	Complex u1;
	Complex u2;
} WeylSpinor;

/*
 *	Complex arithmetic
 */
Complex
complexMake(double re, double im);

Complex
complexFromReal(double realPart);

Complex
complexZero(void);

Complex
complexI(void);

Complex
complexAdd(Complex left, Complex right);

Complex
complexSub(Complex left, Complex right);

Complex
complexMul(Complex left, Complex right);

Complex
complexConj(Complex value);

Complex
complexNeg(Complex value);

Complex
complexScale(double scale, Complex value);

Complex
complexDiv(Complex left, Complex right);

double
complexNorm2(Complex value);

/*
 *	Complex four-vector arithmetic
 */
CVec4
cvec4Make(Complex energy, Complex momentumX, Complex momentumY, Complex momentumZ);

CVec4
cvec4Zero(void);

CVec4
cvec4FromVec4(Vec4 realVector);

CVec4
cvec4Add(CVec4 left, CVec4 right);

CVec4
cvec4Sub(CVec4 left, CVec4 right);

CVec4
cvec4Scale(Complex scale, CVec4 vector);

CVec4
cvec4ScaleReal(double scale, CVec4 vector);

CVec4
cvec4Neg(CVec4 vector);

void
cvec4AddAssign(CVec4 * dest, CVec4 src);

/**
 *	@brief	Lorentz dot product of two complex four-vectors.
 *		Uses metric (+,-,-,-).
 */
Complex
cvec4Dot(CVec4 left, CVec4 right);

/*
 *	Weyl spinor operations
 */

/**
 *	@brief	Construct a Weyl spinor from a lightlike momentum.
 *
 *	@param	helicity	: +1 or -1.
 *	@param	momentum	: Lightlike four-momentum.
 *	@return			: Weyl spinor.
 */
WeylSpinor
weylFromMomentum(int helicity, Vec4 momentum);

/**
 *	@brief	Spinor product: a.u1 * b.u2 - a.u2 * b.u1.
 */
Complex
weylProduct(WeylSpinor left, WeylSpinor right);

/**
 *	@brief	Sigma-mu contraction: maps two Weyl spinors to a complex four-vector.
 *		Implements the sigma^mu_{alpha dot-beta} contraction from
 *		the Weyl-van-der-Waerden spinor formalism.
 */
CVec4
vt(WeylSpinor left, WeylSpinor right);

#ifdef __cplusplus
}
#endif
