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
#include <uxhw.h>
#include "pepper-complexArith.h"

/*
 *	Complex number operations
 */
Complex
complexMake(double re, double im)
{
	Complex result;

	result.re   = re;
	result.im   = im;

	return result;
}

Complex
complexFromReal(double realPart)
{
	Complex result;

	result.re   = realPart;
	result.im   = 0.0;

	return result;
}

Complex
complexZero(void)
{
	Complex result;

	result.re   = 0.0;
	result.im   = 0.0;

	return result;
}

Complex
complexI(void)
{
	Complex result;

	result.re   = 0.0;
	result.im   = 1.0;

	return result;
}

Complex
complexAdd(Complex left, Complex right)
{
	Complex result;

	result.re   = left.re + right.re;
	result.im   = left.im + right.im;

	return result;
}

Complex
complexSub(Complex left, Complex right)
{
	Complex result;

	result.re   = left.re - right.re;
	result.im   = left.im - right.im;

	return result;
}

Complex
complexMul(Complex left, Complex right)
{
	Complex result;

	result.re   = left.re * right.re - left.im * right.im;
	result.im   = left.re * right.im + left.im * right.re;

	return result;
}

Complex
complexConj(Complex value)
{
	Complex result;

	result.re   = value.re;
	result.im   = -value.im;

	return result;
}

Complex
complexNeg(Complex value)
{
	Complex result;

	result.re   = -value.re;
	result.im   = -value.im;

	return result;
}

Complex
complexScale(double scale, Complex value)
{
	Complex result;

	result.re   = scale * value.re;
	result.im   = scale * value.im;

	return result;
}

Complex
complexDiv(Complex left, Complex right)
{
	double  denominator = right.re * right.re + right.im * right.im;
	Complex result;

	result.re   = (left.re * right.re + left.im * right.im) / denominator;
	result.im   = (left.im * right.re - left.re * right.im) / denominator;

	return result;
}

double
complexNorm2(Complex value)
{
	return value.re * value.re + value.im * value.im;
}

/*
 *	Complex four-vector operations
 */
CVec4
cvec4Make(Complex energy, Complex momentumX, Complex momentumY, Complex momentumZ)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = energy;
	result.c[kFourVectorIndexPx]    = momentumX;
	result.c[kFourVectorIndexPy]    = momentumY;
	result.c[kFourVectorIndexPz]    = momentumZ;

	return result;
}

CVec4
cvec4Zero(void)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = complexZero();
	result.c[kFourVectorIndexPx]    = complexZero();
	result.c[kFourVectorIndexPy]    = complexZero();
	result.c[kFourVectorIndexPz]    = complexZero();

	return result;
}

CVec4
cvec4FromVec4(Vec4 realVector)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = complexFromReal(realVector.p[kFourVectorIndexE]);
	result.c[kFourVectorIndexPx]    = complexFromReal(realVector.p[kFourVectorIndexPx]);
	result.c[kFourVectorIndexPy]    = complexFromReal(realVector.p[kFourVectorIndexPy]);
	result.c[kFourVectorIndexPz]    = complexFromReal(realVector.p[kFourVectorIndexPz]);

	return result;
}

CVec4
cvec4Add(CVec4 left, CVec4 right)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = complexAdd(left.c[kFourVectorIndexE], right.c[kFourVectorIndexE]);
	result.c[kFourVectorIndexPx]    = complexAdd(left.c[kFourVectorIndexPx], right.c[kFourVectorIndexPx]);
	result.c[kFourVectorIndexPy]    = complexAdd(left.c[kFourVectorIndexPy], right.c[kFourVectorIndexPy]);
	result.c[kFourVectorIndexPz]    = complexAdd(left.c[kFourVectorIndexPz], right.c[kFourVectorIndexPz]);

	return result;
}

CVec4
cvec4Sub(CVec4 left, CVec4 right)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = complexSub(left.c[kFourVectorIndexE], right.c[kFourVectorIndexE]);
	result.c[kFourVectorIndexPx]    = complexSub(left.c[kFourVectorIndexPx], right.c[kFourVectorIndexPx]);
	result.c[kFourVectorIndexPy]    = complexSub(left.c[kFourVectorIndexPy], right.c[kFourVectorIndexPy]);
	result.c[kFourVectorIndexPz]    = complexSub(left.c[kFourVectorIndexPz], right.c[kFourVectorIndexPz]);

	return result;
}

CVec4
cvec4Scale(Complex scale, CVec4 vector)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = complexMul(scale, vector.c[kFourVectorIndexE]);
	result.c[kFourVectorIndexPx]    = complexMul(scale, vector.c[kFourVectorIndexPx]);
	result.c[kFourVectorIndexPy]    = complexMul(scale, vector.c[kFourVectorIndexPy]);
	result.c[kFourVectorIndexPz]    = complexMul(scale, vector.c[kFourVectorIndexPz]);

	return result;
}

CVec4
cvec4ScaleReal(double scale, CVec4 vector)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = complexScale(scale, vector.c[kFourVectorIndexE]);
	result.c[kFourVectorIndexPx]    = complexScale(scale, vector.c[kFourVectorIndexPx]);
	result.c[kFourVectorIndexPy]    = complexScale(scale, vector.c[kFourVectorIndexPy]);
	result.c[kFourVectorIndexPz]    = complexScale(scale, vector.c[kFourVectorIndexPz]);

	return result;
}

CVec4
cvec4Neg(CVec4 vector)
{
	CVec4 result;

	result.c[kFourVectorIndexE]     = complexNeg(vector.c[kFourVectorIndexE]);
	result.c[kFourVectorIndexPx]    = complexNeg(vector.c[kFourVectorIndexPx]);
	result.c[kFourVectorIndexPy]    = complexNeg(vector.c[kFourVectorIndexPy]);
	result.c[kFourVectorIndexPz]    = complexNeg(vector.c[kFourVectorIndexPz]);

	return result;
}

void
cvec4AddAssign(CVec4 * dest, CVec4 src)
{
	dest->c[kFourVectorIndexE]  = complexAdd(dest->c[kFourVectorIndexE], src.c[kFourVectorIndexE]);
	dest->c[kFourVectorIndexPx] = complexAdd(dest->c[kFourVectorIndexPx], src.c[kFourVectorIndexPx]);
	dest->c[kFourVectorIndexPy] = complexAdd(dest->c[kFourVectorIndexPy], src.c[kFourVectorIndexPy]);
	dest->c[kFourVectorIndexPz] = complexAdd(dest->c[kFourVectorIndexPz], src.c[kFourVectorIndexPz]);
}

Complex
cvec4Dot(CVec4 left, CVec4 right)
{
	/*
	 *	Lorentz metric (+,-,-,-):
	 *	left^mu right_mu = l0*r0 - l1*r1 - l2*r2 - l3*r3
	 */
	Complex result = complexMul(left.c[kFourVectorIndexE], right.c[kFourVectorIndexE]);

	result  = complexSub(result, complexMul(left.c[kFourVectorIndexPx], right.c[kFourVectorIndexPx]));
	result  = complexSub(result, complexMul(left.c[kFourVectorIndexPy], right.c[kFourVectorIndexPy]));
	result  = complexSub(result, complexMul(left.c[kFourVectorIndexPz], right.c[kFourVectorIndexPz]));

	return result;
}

/*
 *	Weyl spinor operations
 */
WeylSpinor
weylFromMomentum(int helicity, Vec4 momentum)
{
	WeylSpinor  spinor;
	double      energy = momentum.p[kFourVectorIndexE];

	/*
	 *	This Weyl spinor basis uses the px component as the light-cone "z"
	 *	axis:
	 *	  p+  = E  + px
	 *	  p-  = E  - px
	 *	  pT  = py + i*pz
	 *
	 *	px is perpendicular to the physical beam axis (stored in pz
	 *	by the phase-space generator), so boosts along the beam do not touch
	 *	the spinor's light-cone components and the Weyl construction is
	 *	manifestly invariant under beam-axis boosts.
	 */
	double  pp      = energy + momentum.p[kFourVectorIndexPx];
	double  pm      = energy - momentum.p[kFourVectorIndexPx];
	double  ptRe    = momentum.p[kFourVectorIndexPy];
	double  ptIm    = (helicity > 0) ? momentum.p[kFourVectorIndexPz] : -momentum.p[kFourVectorIndexPz];

	/*
	 *	Complex square root: sqrt(x) = sqrt(|x|) for x >= 0, i*sqrt(|x|) for x < 0.
	 *	Incoming-particle momenta (negative energy) produce pp, pm < 0 here.
	 *
	 *	On UxHw the sqrt argument is passed through LimitDistributionSupport
	 *	to trim the slightly-negative tail that distributional arithmetic
	 *	produces at the (y=0, phi=pi) corner of outgoing phase space, where
	 *	pp = pT*(cosh(y)+cos(phi)) touches zero. Without the clamp, sqrt()
	 *	on that tail emits NaN that contaminates the whole distribution.
	 *
	 *	The clamps live inside the branch that matches the sign of pp / pm
	 *	so the Laplace core never sees a distribution with zero support in
	 *	the limits (which would otherwise fall back to a NaN particle).
	 */
	const double kSqrtFloor = 1e-30;

	Complex u1;
	Complex u2;

	if (pp >= 0.0)
	{
		double ppSafe = UxHwDoubleLimitDistributionSupport(pp, kSqrtFloor, INFINITY);
		u1 = complexMake(sqrt(ppSafe), 0.0);
	}
	else
	{
		double negPpSafe = UxHwDoubleLimitDistributionSupport(-pp, kSqrtFloor, INFINITY);
		u1 = complexMake(0.0, sqrt(negPpSafe));
	}

	if (pm >= 0.0)
	{
		double pmSafe = UxHwDoubleLimitDistributionSupport(pm, kSqrtFloor, INFINITY);
		u2 = complexMake(sqrt(pmSafe), 0.0);
	}
	else
	{
		double negPmSafe = UxHwDoubleLimitDistributionSupport(-pm, kSqrtFloor, INFINITY);
		u2 = complexMake(0.0, sqrt(negPmSafe));
	}

	/*
	 *	When the transverse momentum is nonzero, re-derive u2 from u1 and pT
	 *	so that the spinor product satisfies <p|p> = 0 to machine precision.
	 */
	double  scaleTolerance  = fabs(energy) * 1e-12;
	int     ptNonzero       = (fabs(ptRe) > scaleTolerance) || (fabs(ptIm) > scaleTolerance);
	int     u1Nonzero       = (fabs(u1.re) > scaleTolerance) || (fabs(u1.im) > scaleTolerance);

	if (ptNonzero && u1Nonzero)
	{
		Complex pt = complexMake(ptRe, ptIm);

		u2 = complexDiv(pt, u1);
	}

	/*
	 *	Phase convention for negative-norm light-cone components (incoming
	 *	particles): positive-helicity spinors get rotated by -i, negative-
	 *	helicity spinors by +i. This matches the Python Weyl convention.
	 */
	if (pp < 0.0 || pm < 0.0)
	{
		Complex t;

		if (helicity < 0)
		{
			/* u -> i*u : (a + b*i) -> -b + a*i */
			t   = complexMake(-u1.im, u1.re); u1 = t;
			t   = complexMake(-u2.im, u2.re); u2 = t;
		}
		else
		{
			/* u -> -i*u : (a + b*i) -> b - a*i */
			t   = complexMake(u1.im, -u1.re); u1 = t;
			t   = complexMake(u2.im, -u2.re); u2 = t;
		}
	}

	spinor.u1   = u1;
	spinor.u2   = u2;

	return spinor;
}

Complex
weylProduct(WeylSpinor left, WeylSpinor right)
{
	/*
	 *	Spinor product: <left|right> = left.u1 * right.u2 - left.u2 * right.u1
	 */
	return complexSub(complexMul(left.u1, right.u2), complexMul(left.u2, right.u1));
}

CVec4
vt(WeylSpinor left, WeylSpinor right)
{
	/*
	 *	Sigma-mu contraction in the Weyl basis where index 1 is the
	 *	light-cone "z" axis and (index 2, index 3) forms the transverse
	 *	(real, imaginary) plane:
	 *	  [kFourVectorIndexE]  = l.u1*r.u1 + l.u2*r.u2       (energy)
	 *	  [kFourVectorIndexPx] = l.u1*r.u1 - l.u2*r.u2       (Weyl-longitudinal)
	 *	  [kFourVectorIndexPy] = l.u1*r.u2 + l.u2*r.u1       (Weyl-transverse Re)
	 *	  [kFourVectorIndexPz] = i*(l.u1*r.u2 - l.u2*r.u1)   (Weyl-transverse Im)
	 */
	Complex u1u1    = complexMul(left.u1, right.u1);
	Complex u2u2    = complexMul(left.u2, right.u2);
	Complex u1u2    = complexMul(left.u1, right.u2);
	Complex u2u1    = complexMul(left.u2, right.u1);

	CVec4 result;

	result.c[kFourVectorIndexE]     = complexAdd(u1u1, u2u2);
	result.c[kFourVectorIndexPx]    = complexSub(u1u1, u2u2);
	result.c[kFourVectorIndexPy]    = complexAdd(u1u2, u2u1);
	result.c[kFourVectorIndexPz]    = complexMul(complexI(), complexSub(u1u2, u2u1));

	return result;
}
