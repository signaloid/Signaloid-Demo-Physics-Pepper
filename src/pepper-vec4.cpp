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

#include <math.h>
#include "pepper-vec4.h"


Vec4
vec4Make(double energy, double momentumX, double momentumY, double momentumZ)
{
	Vec4 result;

	result.p[kFourVectorIndexE]     = energy;
	result.p[kFourVectorIndexPx]    = momentumX;
	result.p[kFourVectorIndexPy]    = momentumY;
	result.p[kFourVectorIndexPz]    = momentumZ;

	return result;
}

Vec4
vec4Zero(void)
{
	Vec4 result;

	result.p[kFourVectorIndexE]     = 0.0;
	result.p[kFourVectorIndexPx]    = 0.0;
	result.p[kFourVectorIndexPy]    = 0.0;
	result.p[kFourVectorIndexPz]    = 0.0;

	return result;
}

Vec4
vec4Add(Vec4 left, Vec4 right)
{
	Vec4 result;

	result.p[kFourVectorIndexE]     = left.p[kFourVectorIndexE] + right.p[kFourVectorIndexE];
	result.p[kFourVectorIndexPx]    = left.p[kFourVectorIndexPx] + right.p[kFourVectorIndexPx];
	result.p[kFourVectorIndexPy]    = left.p[kFourVectorIndexPy] + right.p[kFourVectorIndexPy];
	result.p[kFourVectorIndexPz]    = left.p[kFourVectorIndexPz] + right.p[kFourVectorIndexPz];

	return result;
}

Vec4
vec4Sub(Vec4 left, Vec4 right)
{
	Vec4 result;

	result.p[kFourVectorIndexE]     = left.p[kFourVectorIndexE] - right.p[kFourVectorIndexE];
	result.p[kFourVectorIndexPx]    = left.p[kFourVectorIndexPx] - right.p[kFourVectorIndexPx];
	result.p[kFourVectorIndexPy]    = left.p[kFourVectorIndexPy] - right.p[kFourVectorIndexPy];
	result.p[kFourVectorIndexPz]    = left.p[kFourVectorIndexPz] - right.p[kFourVectorIndexPz];

	return result;
}

Vec4
vec4Scale(double scale, Vec4 vector)
{
	Vec4 result;

	result.p[kFourVectorIndexE]     = scale * vector.p[kFourVectorIndexE];
	result.p[kFourVectorIndexPx]    = scale * vector.p[kFourVectorIndexPx];
	result.p[kFourVectorIndexPy]    = scale * vector.p[kFourVectorIndexPy];
	result.p[kFourVectorIndexPz]    = scale * vector.p[kFourVectorIndexPz];

	return result;
}

Vec4
vec4Neg(Vec4 vector)
{
	Vec4 result;

	result.p[kFourVectorIndexE]     = -vector.p[kFourVectorIndexE];
	result.p[kFourVectorIndexPx]    = -vector.p[kFourVectorIndexPx];
	result.p[kFourVectorIndexPy]    = -vector.p[kFourVectorIndexPy];
	result.p[kFourVectorIndexPz]    = -vector.p[kFourVectorIndexPz];

	return result;
}

double
vec4Dot(Vec4 left, Vec4 right)
{
	return left.p[kFourVectorIndexE] * right.p[kFourVectorIndexE]
	       - left.p[kFourVectorIndexPx] * right.p[kFourVectorIndexPx]
	       - left.p[kFourVectorIndexPy] * right.p[kFourVectorIndexPy]
	       - left.p[kFourVectorIndexPz] * right.p[kFourVectorIndexPz];
}

double
vec4M2(Vec4 vector)
{
	return vec4Dot(vector, vector);
}

double
vec4PT2(Vec4 vector)
{
	return vector.p[kFourVectorIndexPx] * vector.p[kFourVectorIndexPx]
	       + vector.p[kFourVectorIndexPy] * vector.p[kFourVectorIndexPy];
}

double
vec4Y(Vec4 vector)
{
	double  lightConePlus   = vec4PPlus(vector);
	double  lightConeMinus  = vec4PMinus(vector);

	return 0.5 * log(lightConePlus / lightConeMinus);
}

double
vec4PPlus(Vec4 vector)
{
	return vector.p[kFourVectorIndexE] + vector.p[kFourVectorIndexPz];
}

double
vec4PMinus(Vec4 vector)
{
	return vector.p[kFourVectorIndexE] - vector.p[kFourVectorIndexPz];
}
