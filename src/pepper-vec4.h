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

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	Component indices of a Lorentz four-vector, in the order the
 *		contravariant components are stored: p^mu = (E, px, py, pz).
 *
 *	Used by both `Vec4` below and `CVec4` in `pepper-complexArith.h`.
 */
typedef enum
{
	kFourVectorIndexE = 0,
	kFourVectorIndexPx,
	kFourVectorIndexPy,
	kFourVectorIndexPz,
	kFourVectorIndexMax,
} FourVectorIndex;

/**
 *	@brief	Real four-vector with Lorentz metric (+,-,-,-).
 */
typedef struct
{
	double p[kFourVectorIndexMax];
} Vec4;

Vec4
vec4Make(double energy, double momentumX, double momentumY, double momentumZ);

Vec4
vec4Zero(void);

Vec4
vec4Add(Vec4 left, Vec4 right);

Vec4
vec4Sub(Vec4 left, Vec4 right);

Vec4
vec4Scale(double scale, Vec4 vector);

Vec4
vec4Neg(Vec4 vector);

/**
 *	@brief	Lorentz dot product with metric (+,-,-,-).
 */
double
vec4Dot(Vec4 left, Vec4 right);

/**
 *	@brief	Invariant mass squared: p^mu p_mu.
 */
double
vec4M2(Vec4 vector);

/**
 *	@brief	Transverse momentum squared: px^2 + py^2.
 */
double
vec4PT2(Vec4 vector);

/**
 *	@brief	Rapidity: 0.5 * ln((E+pz)/(E-pz)).
 */
double
vec4Y(Vec4 vector);

/**
 *	@brief	Light-cone plus component: E + pz.
 */
double
vec4PPlus(Vec4 vector);

/**
 *	@brief	Light-cone minus component: E - pz.
 */
double
vec4PMinus(Vec4 vector);

#ifdef __cplusplus
}
#endif
