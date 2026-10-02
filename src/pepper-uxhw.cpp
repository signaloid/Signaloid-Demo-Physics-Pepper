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

/*
 *
 *	UxHw cross-section kernel. Computes the gg -> gg
 *	differential cross section dsigma in a single distributional pass.
 *
 */

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#ifndef uncertain
 #define uncertain
#endif

#include <uxhw.h>
#ifdef max
 #undef max
#endif
#ifdef min
 #undef min
#endif

extern "C"
{
#include "pepper-kernel.h"
#include "pepper-uxhw.h"
#include "pepper-lhapdfLite.h"
#include "pepper-utilities.h"
}

/*
 *	Args for the per-particle PDF lookup, threaded through the
 *	UxHwDoublePropagateFunction `args` parameter so each particle's scalar
 *	x reaches the lhapdfLite reader. Without that, UxHw collapses the
 *	distributional x to its mean before the reader's bin-selection
 *	conditionals run.
 */
typedef struct
{
	PepperPDF * pdf;
	double      q2;
	int         pid;
} PdfPropagateArgs;

extern "C" void
pdfPropagateCallback(
	void *      args,
	double *    input,
	size_t      nIn,
	double *    output,
	size_t      nOut)
{
	(void) nIn;
	(void) nOut;

	PdfPropagateArgs * lookupArgs = (PdfPropagateArgs *) args;

	output[0] = pepperPDFxfxQ2(lookupArgs->pdf, lookupArgs->pid, input[0], lookupArgs->q2);
}

/**
 *	@brief	Per-particle PDF lookup, wrapped in UxHwDoublePropagateFunction.
 *
 *	Calls pepperPDFxfxQ2 with a scalar x for each particle of the input
 *	distribution, so the output is a joint distribution of xf values
 *	correlated with x. That is what makes the cube integral
 *	E[... * xf(x1) * xf(x2) * ...] correct.
 */
static double
pdfLookupPropagated(double x, double q2)
{
	PdfPropagateArgs lookupArgs;

	lookupArgs.pdf  = pepperKernelPDF();
	lookupArgs.q2   = q2;
	lookupArgs.pid  = 21;   /* gluon */

	uncertain double    in[1] = { x };
	uncertain double    out[1];

	UxHwDoublePropagateFunction(
		&pdfPropagateCallback,
		(void *) &lookupArgs,
		in,
		1,
		out,
		1
	);

	return out[0];
}

/*
 *	Exact per-particle yminV/ymaxV = (-log|ymnRaw|, log|ymxRaw|), wrapped in
 *	a single UxHwDoublePropagateFunction call. Used to feed log() a
 *	strictly-non-negative recoil-rapidity bound without the sqrt(x^2 + eps^2)
 *	surrogate: the scalar fabs runs once per Dirac, so it is exact (no eps
 *	floor) and never branches on a distribution.
 */
extern "C" void
yvBoundsPropagateCallback(void * args, double * input, size_t nIn, double * output, size_t nOut)
{
	(void) args;
	(void) nIn;
	(void) nOut;

	output[0]   = -log(fabs(input[0]));
	output[1]   = log(fabs(input[1]));
}

static void
yvBoundsPropagated(double ymnRaw, double ymxRaw, double * yminVOut, double * ymaxVOut)
{
	uncertain double    in[2] = { ymnRaw, ymxRaw };
	uncertain double    out[2];

	UxHwDoublePropagateFunction(&yvBoundsPropagateCallback, NULL, in, 2, out, 2);

	*yminVOut   = out[0];
	*ymaxVOut   = out[1];
}

/*
 *	Exact Heaviside as plain arithmetic.
 *
 *	copysign(1.0, x) is +1 for x > 0 and for x == +0.0, and -1 for x < 0 and
 *	for x == -0.0, so this returns exactly 1.0 and 0.0 on those branches --
 *	bit-identical to the commented-out callback for every x != 0, with no eps
 *	softening and no approximation.
 */
static double
step(double value)
{
	return 0.5 * (1.0 + copysign(1.0, value));
}

static double
acceptanceGate(
	double  ymnRaw,
	double  ymxRaw,
	double  pp,
	double  pm,
	double  ecmsMinusPp,
	double  ecmsMinusPm)
{
	double  g0  = step(ymnRaw) * step(ymxRaw);
	double  g1  = step(pp) * step(pm);
	double  g2  = step(ecmsMinusPp) * step(ecmsMinusPm);

	return (g0 * g1) * g2;
}

double
generateXsPointUxHw(void)
{
	double  randomPt2               = UxHwDoubleUniformDist(0.0, 1.0);
	double  randomRapidity          = UxHwDoubleUniformDist(0.0, 1.0);
	double  randomRecoilRapidity    = UxHwDoubleUniformDist(0.0, 1.0);

	double  ecms    = kPepperConstantEcms;
	double  ptmin   = kPepperConstantPtMin;
	double  g_s2    = 4.0 * M_PI * pepperKernelAlphasMZ2();
	double  pt2min  = ptmin * ptmin;
	double  sHat    = ecms * ecms;
	double  pt2max  = sHat / 4.0;

	/*
	 *	pT^2: PeakedDist inverse-CDF, importance weight ~ 1/(pT^2)^alpha
	 *	with alpha = 2.
	 */
	const double    alphaImp            = kPepperConstantAlphaImp;
	double          oneMinusAlphaImp    = 1.0 - alphaImp;
	double          cdfLow              = pow(pt2min, oneMinusAlphaImp);
	double          cdfHigh             = pow(pt2max, oneMinusAlphaImp);
	double          pt2                 = pow(cdfLow + randomPt2 * (cdfHigh - cdfLow), 1.0 / oneMinusAlphaImp);
	double          pt2Weight           = (cdfHigh - cdfLow) / (oneMinusAlphaImp * pow(pt2, -alphaImp));

	double weight = (2.0 * M_PI / sHat) * (pt2Weight / (32.0 * M_PI * M_PI * M_PI));

	/*
	 *	y (auxiliary rapidity): uniform on the kinematically-allowed
	 *	range [-ymax_eff, +ymax_eff], where ymax_eff is the smaller
	 *	of the user-specified |y| < 5 cut and the boundary set by
	 *	pT^2 (cosh^2 y_b = sHat/(4 pT^2)).
	 */
	double  innerR  = sHat / (4.0 * pt2);
	double  yb      = log(sqrt(innerR) + sqrt(innerR - 1.0));
	double  ymaxY   = fmin(yb, kPepperConstantYMax);
	double  yminY   = -ymaxY;
	double  y       = yminY + randomRapidity * (ymaxY - yminY);

	weight *= (ymaxY - yminY) * 2.0 * M_PI; /* 2 pi = analytic int d phi */

	double  sqrtPt2 = sqrt(pt2);
	double  Qt      = sqrtPt2;
	double  ey      = exp(y);
	double  emy     = exp(-y);

	/*
	 *	yv (recoil rapidity): the kinematic bound depends on Qt and y,
	 *	and can be negative (rejection region) when the recoil cannot
	 *	balance. Feed log() the exact |bound| (per-particle, via
	 *	yvBoundsPropagated) so it never sees a non-positive input.
	 */
	double  ymnRaw  = ecms / Qt - emy;
	double  ymxRaw  = ecms / Qt - ey;
	double  yminV, ymaxV;
	yvBoundsPropagated(ymnRaw, ymxRaw, &yminV, &ymaxV);

	double yv = yminV + randomRecoilRapidity * (ymaxV - yminV);

	weight *= (ymaxV - yminV);

	double  eyv     = exp(yv);
	double  emyv    = 1.0 / eyv;
	double  pp      = sqrtPt2 * (eyv + ey);         /* light-cone+ of beam */
	double  pm      = sqrtPt2 * (emyv + emy);       /* light-cone- of beam */

	/*
	 *	Closed-form gg -> gg matrix element. The textbook expression
	 *	(9/4) g_s^4 (3 - tu/s^2 - us/t^2 - st/u^2) reduces to a function
	 *	of the single distribution r = exp(yv - y) (the pp^2 and Qt^2
	 *	factors algebraically cancel in every term):
	 *
	 *	    |M|^2 = (9/4) g_s^4 (3 + r + 1/r + r^2 + 1/r^2 - r/(r+1)^2)
	 *
	 */
	double  r           = emy / emyv;
	double  invR        = 1.0 / r;
	double  r2          = r * r;
	double  invR2       = invR * invR;
	double  rp1         = r + 1.0;
	double  rp1sq       = rp1 * rp1;
	double  rOverRp1sq  = r / rp1sq;
	double  bracket     = 3.0 + r;
	bracket += invR;
	bracket += r2;
	bracket += invR2;
	bracket -= rOverRp1sq;
	double me2 = (9.0 / 4.0) * g_s2 * g_s2 * bracket;

	double  dxsPreConst = me2 * weight;
	double  dxs         = dxsPreConst * kPepperConstantGeV2ToPb;

	/*
	 *	Parton x's.
	 */
	double  x1  = pp / ecms;
	double  x2  = pm / ecms;

	double  xf1 = pdfLookupPropagated(x1, kPepperConstantQ2);
	double  xf2 = pdfLookupPropagated(x2, kPepperConstantQ2);

	/*
	 *	Combined PDF luminosity + partonic flux factor:
	 *	    (xf1/x1) (xf2/x2) / (2 x1 x2 ecms^2)
	 *	  = xf1 * xf2 * ecms^2 / (2 pp^2 pm^2)
	 *	using x1 = pp/ecms, x2 = pm/ecms.
	 */
	double xfProduct = xf1 * xf2;
	dxs *= xfProduct;
	dxs *= ecms * ecms;
	double  ppSq    = pp * pp;
	double  pmSq    = pm * pm;
	dxs /= 2.0 * ppSq * pmSq;

	/*
	 *	Acceptance gates: exact per-particle Heaviside and acceptanceGate,
	 *	multiplied in as a single combined factor. It is 1 in the physical
	 *	region and 0 in the rejection region, preserving per-particle
	 *	correlation with dxs.
	 *
	 *	    ymnRaw > 0 / ymxRaw > 0 -> recoil rapidity range valid
	 *	    pp > 0 / pm > 0         -> light-cone fractions positive
	 *	    ecms - pp / ecms - pm   -> light-cone fractions <= ecms
	 */
	dxs *= acceptanceGate(ymnRaw, ymxRaw, pp, pm, ecms - pp, ecms - pm);

	/*
	 *	Self-gate on dxs to zero negative contributions from PDF
	 *	interpolator overshoot or quantization artefacts.
	 */
	dxs *= step(dxs);

	return dxs;
}
