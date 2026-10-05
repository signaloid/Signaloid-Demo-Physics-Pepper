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
 *	========================================================================
 *	pepper-monteCarlo.cpp
 *
 *	Monte Carlo (native sampling) cross-section kernel. Computes the
 *	gg -> gg differential cross section dsigma with conventional Monte
 *	Carlo: each call to generateXsPointMonteCarlo() draws three Uniform[0,1]
 *	random samples via the GSL-backed compat layer and evaluates the
 *	integrand as PLAIN SCALAR arithmetic. No UxHw distributional API is
 *	used. The kinematic gates, |.|, and parton-x clamp reduce to ordinary
 *	fmin / fmax / fabs / step operations, and the PDF lookup is a direct
 *	scalar table read. calculateOutputMonteCarlo (pepper-kernel.cpp)
 *	averages the samples to estimate the integrated cross section.
 *
 *	This mirrors the distributional pepper-uxhw.cpp op-for-op, so the two
 *	kernels compute the same physics: pepper-uxhw.cpp propagates the whole
 *	Uniform[0,1]^3 cube in one pass and contracts to the first moment,
 *	while this file samples the cube point-by-point and contracts to the
 *	sample mean.
 *	========================================================================
 */

#include <cmath>
#include <cstdlib>

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
#include "pepper-monteCarlo.h"
#include "pepper-lhapdfLite.h"
}

double
generateXsPointMonteCarlo(void)
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
	 *	with alpha = 2 (matches the dominant pT^2 dependence of the
	 *	matrix element, so the post-Jacobian integrand is nearly flat).
	 */
	const double    alphaImp            = kPepperConstantAlphaImp;
	double          oneMinusAlphaImp    = 1.0 - alphaImp;
	double          cdfLow              = pow(pt2min, oneMinusAlphaImp);
	double          cdfHigh             = pow(pt2max, oneMinusAlphaImp);
	double          pt2                 = pow(cdfLow + randomPt2 * (cdfHigh - cdfLow), 1.0 / oneMinusAlphaImp);
	double          pt2Weight           = (cdfHigh - cdfLow) / (oneMinusAlphaImp * pow(pt2, -alphaImp));

	double weight = (2.0 * M_PI / sHat) * (pt2Weight / (32.0 * M_PI * M_PI * M_PI));

	/*
	 *	y (auxiliary rapidity): uniform on [-ymax_eff, +ymax_eff], where
	 *	ymax_eff is the smaller of the user |y| < 5 cut and the pT^2
	 *	boundary (cosh^2 y_b = sHat/(4 pT^2)). On a scalar sample the smaller
	 *	of is a plain fmin.
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
	double  emy     = 1.0 / ey;

	/*
	 *	yv (recoil rapidity): the kinematic bound depends on Qt and y and
	 *	goes negative in the rejection region. Feed log() the exact
	 *	|bound| via fabs so it never sees a non-positive input. The step
	 *	gates below zero the rejected contribution.
	 */
	double  ymnRaw  = ecms / Qt - emy;
	double  ymxRaw  = ecms / Qt - ey;
	double  ymnSafe = fabs(ymnRaw);
	double  ymxSafe = fabs(ymxRaw);

	double  yminV   = -log(ymnSafe);
	double  ymaxV   = log(ymxSafe);
	double  yv      = yminV + randomRecoilRapidity * (ymaxV - yminV);

	weight *= (ymaxV - yminV);

	double  eyv     = exp(yv);
	double  emyv    = 1.0 / eyv;
	double  pp      = sqrtPt2 * (eyv + ey);         /* light-cone+ of beam */
	double  pm      = sqrtPt2 * (emyv + emy);       /* light-cone- of beam */

	/*
	 *	Closed-form gg -> gg matrix element, a function of the single
	 *	ratio r = exp(yv - y):
	 *
	 *	    |M|^2 = (9/4) g_s^4 (3 + r + 1/r + r^2 + 1/r^2 - r/(r+1)^2)
	 */
	double  r       = emy / emyv;
	double  invR    = 1.0 / r;
	double  r2      = r * r;
	double  invR2   = invR * invR;
	double  rp1     = r + 1.0;
	double  bracket = 3.0 + r + invR + r2 + invR2 - r / (rp1 * rp1);
	double  me2     = (9.0 / 4.0) * g_s2 * g_s2 * bracket;

	double dxs = me2 * weight * kPepperConstantGeV2ToPb;

	/*
	 *	Parton x's, fed to the PDF lookup unclamped. The kinematics keep
	 *	x inside the CT18NLO grid, so the interpolator never
	 *	extrapolates. See the matching comment in pepper-uxhw.cpp.
	 */
	double  x1  = pp / ecms;
	double  x2  = pm / ecms;

	double  xf1 = pepperPDFxfxQ2(pepperKernelPDF(), 21, x1, kPepperConstantQ2);
	double  xf2 = pepperPDFxfxQ2(pepperKernelPDF(), 21, x2, kPepperConstantQ2);

	/*
	 *	Combined PDF luminosity + partonic flux factor:
	 *	    (xf1/x1) (xf2/x2) / (2 x1 x2 ecms^2)
	 *	  = xf1 * xf2 * ecms^2 / (2 pp^2 pm^2).
	 */
	dxs *= xf1 * xf2;
	dxs *= ecms * ecms;
	dxs /= 2.0 * pp * pp * pm * pm;

	/*
	 *	Acceptance gates: scalar Heaviside indicators multiplied in. Each
	 *	factor is 1 in the physical region and 0 in its rejection region.
	 *
	 *	    ymnRaw > 0 / ymxRaw > 0 -> recoil rapidity range valid
	 *	    pp > 0 / pm > 0         -> light-cone fractions positive
	 *	    ecms - pp / ecms - pm   -> light-cone fractions <= ecms
	 */
	dxs *= (ymnRaw > 0.0) ? 1.0 : 0.0;
	dxs *= (ymxRaw > 0.0) ? 1.0 : 0.0;
	dxs *= (pp > 0.0) ? 1.0 : 0.0;
	dxs *= (pm > 0.0) ? 1.0 : 0.0;
	dxs *= (ecms - pp > 0.0) ? 1.0 : 0.0;
	dxs *= (ecms - pm > 0.0) ? 1.0 : 0.0;

	/*
	 *	Self-gate: zero negative contributions from PDF interpolator
	 *	overshoot (relu).
	 */
	dxs *= (dxs > 0.0) ? 1.0 : 0.0;

	return dxs;
}

double
generateXsMonteCarlo(
	size_t      numberOfMonteCarloIterations,
	double *    monteCarloOutputSamples)
{
	double sample = 0.0;

	for (size_t iteration = 0; iteration < numberOfMonteCarloIterations; iteration++)
	{
		sample = generateXsPointMonteCarlo();

		/*
		 *	Scalar NaN guard: a rare sample deep in the rejection region
		 *	can make log()/sqrt() emit NaN. This is a comparison on a
		 *	single scalar sample, not a branch on a distribution.
		 */
		if (sample != sample)
		{
			sample = 0.0;
		}

		monteCarloOutputSamples[iteration] = sample;
	}

	return sample;
}
