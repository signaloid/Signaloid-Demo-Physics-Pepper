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
#include "pepper-phaseSpace.h"

/**
 *	@brief	Importance sampling exponent for the pT^2 distribution.
 */
static const double kAlphaImportance = 2.0;

/**
 *	@brief	Importance sampling via the inverse-CDF method for a distribution
 *		proportional to (a + k*x)^(-cn).
 *
 *		Based on TChannel.PeakedDist from the Python implementation.
 *
 *	@param	offset	: Offset parameter.
 *	@param	exponent: Power-law exponent.
 *	@param	xLow	: Lower bound of x.
 *	@param	xHigh	: Upper bound of x.
 *	@param	scale	: Scale parameter.
 *	@param	uniform	: Uniform random number in [0, 1].
 *	@param	x	: Output sampled value.
 *	@param	weight	: Output Jacobian weight.
 */
static void
peakedDist(
	double      offset,
	double      exponent,
	double      xLow,
	double      xHigh,
	double      scale,
	double      uniform,
	double *    x,
	double *    weight)
{
	double oneMinusExponent = 1.0 - exponent;

	if (oneMinusExponent != 0.0)
	{
		double  aminVal = pow(offset + scale * xLow, oneMinusExponent);
		double  amaxVal = pow(offset + scale * xHigh, oneMinusExponent);
		double  mapped  = aminVal + uniform * (amaxVal - aminVal);

		*x      = scale * (pow(mapped, 1.0 / oneMinusExponent) - offset);
		*weight = (amaxVal - aminVal) / (scale * oneMinusExponent * pow(offset + scale * (*x), -exponent));
	}
	else
	{
		double  aminVal = offset + scale * xLow;
		double  amaxVal = offset + scale * xHigh;

		*x      = scale * (aminVal * pow(amaxVal / aminVal, uniform) - offset);
		*weight = log(amaxVal / aminVal) / (scale * pow(offset + scale * (*x), -exponent));
	}
}

void
generatePhaseSpacePoint(
	Particle    particles[kPartonIndexMax],
	double      ecms,
	double      pt2min,
	double      ymax,
	double *    weight)
{
	double  S       = ecms * ecms;
	double  pt2max  = S / 4.0;
	double  psw     = 2.0 * M_PI / S;

	/*
	 *	Generate pT^2, rapidity, and azimuthal angle for parton 3
	 *	(the loop runs only for i=3 in the 2->2 case).
	 */
	double  rn0 = UxHwDoubleUniformDist(0.0, 1.0);
	double  rn1 = UxHwDoubleUniformDist(0.0, 1.0);
	double  rn2 = UxHwDoubleUniformDist(0.0, 1.0);

	double pt2, wpt2;

	peakedDist(0.0, kAlphaImportance, pt2min, pt2max, 1.0, rn0, &pt2, &wpt2);
	psw *= wpt2 / (32.0 * M_PI * M_PI * M_PI);

	/*
	 *	Rapidity bounds from kinematics.
	 */
	double  ratio   = sqrt(S / (4.0 * pt2));
	double  yb      = log(ratio + sqrt(ratio * ratio - 1.0));
	double  ymaxEff = (ymax < yb) ? ymax : yb;
	double  yminEff = (-ymax > -yb) ? -ymax : -yb;

	double  y   = yminEff + rn1 * (ymaxEff - yminEff);
	double  phi = 2.0 * M_PI * rn2;

	psw *= (ymaxEff - yminEff) * 2.0 * M_PI;

	/*
	 *	Construct four-momentum for parton 3.
	 *
	 *	Derive cosh(y) and sinh(y) from a single exp(y) evaluation (and
	 *	its reciprocal) rather than calling cosh()/sinh() independently.
	 *	Downstream Lorentz-invariant cross products evaluate light-cone
	 *	sums of the form cosh(y) + sinh(y), which analytically collapse
	 *	to exp(y). UxHw COMMON_ANCESTOR tracking can only reconstruct
	 *	that collapse when both terms derive from one shared exp(y)
	 *	distribution. Otherwise cosh(y) + sinh(y) loses precision at
	 *	extreme |y|.
	 */
	double  expy    = exp(y);
	double  expmy   = exp(-y);
	double  coshy   = 0.5 * (expy + expmy);
	double  sinhy   = 0.5 * (expy - expmy);
	double  sqrtPt2 = sqrt(pt2);

	Vec4 p3 = vec4Make(
		sqrtPt2 * coshy,
		sqrtPt2 * cos(phi),
		sqrtPt2 * sin(phi),
		sqrtPt2 * sinhy
	);

	Vec4 psum = p3;

	/*
	 *	Kinematic variables for the recoil system. In the 2->2 case psum is
	 *	exactly p3, which we just built from (pt2, y, phi) and which is
	 *	lightlike by construction. Using that construction directly (rather
	 *	than vec4M2 / vec4PT2 / vec4Y on psum) avoids correlation loss on
	 *	UxHw: a distributional v.p[i] * v.p[i] is evaluated as a product of
	 *	two random variables with the same distribution rather than as a
	 *	square, so ptj2 can pick up negative support and cosh^2 y - sinh^2 y
	 *	can drift off 1. sqrt(ptj2) then emits NaN. Skipping those calls
	 *	keeps the recoil invariant mass at 0 and its pT^2 at pt2 (neither
	 *	is needed below, so neither is materialised), and yj = y,
	 *	Qt = sqrtPt2, mt = sqrtPt2 strictly by construction.
	 */
	double  yj  = y;
	double  Qt  = sqrtPt2;
	double  mt  = sqrtPt2;

	/*
	 *	Rapidity bounds for the last parton (parton 2),
	 *	determined by momentum conservation.
	 */
	double  ymn = ecms / Qt * (1.0 - mt / ecms * exp(-yj));
	double  ymx = ecms / Qt * (1.0 - mt / ecms * exp(yj));

	/*
	 *	Truncate ymn and ymx to positive values with
	 *	UxHwDoubleLimitDistributionSupport. These are distributional on
	 *	UxHw, and hard branching on a distributional value evaluates on
	 *	the mean and discards valid parts of the distribution. Limiting
	 *	the support to the physically valid positive region instead lets
	 *	invalid configurations produce zero weight naturally.
	 */
	ymn = UxHwDoubleLimitDistributionSupport(ymn, 1e-30, INFINITY);
	ymx = UxHwDoubleLimitDistributionSupport(ymx, 1e-30, INFINITY);

	ymn = -log(ymn);
	ymx = log(ymx);

	double  rn3 = UxHwDoubleUniformDist(0.0, 1.0);
	double  yv  = ymn + rn3 * (ymx - ymn);

	psw *= ymx - ymn;

	/*
	 *	Construct parton 2 momentum from momentum conservation.
	 *	Use a single exp(yv)/exp(-yv) pair as the shared ancestor for
	 *	cosh(yv) and sinh(yv), so the downstream light-cone sum
	 *	cosh(yv) + sinh(yv) collapses back to exp(yv) under UxHw
	 *	COMMON_ANCESTOR tracking.
	 */
	double  expyv   = exp(yv);
	double  expmyv  = exp(-yv);
	double  coshyv  = 0.5 * (expyv + expmyv);
	double  sinhyv  = 0.5 * (expyv - expmyv);

	Vec4 p2 = vec4Make(
		Qt * coshyv,
		-psum.p[kFourVectorIndexPx],
		-psum.p[kFourVectorIndexPy],
		Qt * sinhyv
	);

	/*
	 *	Light-cone components of each outgoing parton, computed directly
	 *	from exp(+/-) to avoid the cosh + sinh cancellation. These are
	 *	the values that the Berends-Giele propagator cache reads via
	 *	particles[i].plus / .minus below.
	 */
	double  p3Plus  = sqrtPt2 * expy;
	double  p3Minus = sqrtPt2 * expmy;
	double  p2Plus  = Qt * expyv;
	double  p2Minus = Qt * expmyv;

	/*
	 *	Reconstruct initial-state momenta from energy/momentum conservation.
	 *	Computing pp, pm via particle.plus / particle.minus keeps the
	 *	calculation in light-cone form (sums of positive exp values) rather
	 *	than going through vec4PPlus / vec4PMinus of a Cartesian sum, which
	 *	would again need cosh + sinh cancellations on UxHw.
	 */
	double  pp  = p2Plus + p3Plus;
	double  pm  = p2Minus + p3Minus;

	/*
	 *	Energy conservation constraints: pp and pm must be in (0, ecms].
	 *	Rather than branching (which evaluates on the mean for distributional
	 *	values), we limit pp and pm to the valid range. Configurations outside
	 *	this range will have their contribution suppressed.
	 */
	pp  = UxHwDoubleLimitDistributionSupport(pp, 1e-30, ecms);
	pm  = UxHwDoubleLimitDistributionSupport(pm, 1e-30, ecms);

	Vec4    p0  = vec4Make(-pp / 2.0, 0.0, 0.0, -pp / 2.0);
	Vec4    p1  = vec4Make(-pm / 2.0, 0.0, 0.0, pm / 2.0);

	/*
	 *	Populate particles[].plus / .minus directly. Beam particles have
	 *	one light-cone component zero by construction. Outgoing particles
	 *	take their exp-derived values from above.
	 *
	 *	r = emy / emyv = exp(yv - y).
	 *
	 *	This is the SINGLE distributional quantity the closed-form gg->gg
	 *	|M|^2 depends on (everything else algebraically cancels). Computing
	 *	r here, directly from the two exp results we already have, makes it
	 *	a clean function of rn1 (through expmy) and rn3 (through expmyv)
	 *	-- no sqrtPt2 prefactor that UxHw would have to cancel between
	 *	numerator and denominator at the ME^2 site. Stashed on particles[kPartonIndexIncomingA]
	 *	as a shared event-level scalar.
	 */
	double rKine = expmy / expmyv;

	particles[kPartonIndexIncomingA].mom    = p0;
	particles[kPartonIndexIncomingA].plus   = -pp;  /* p0.E + p0.pz = -pp */
	particles[kPartonIndexIncomingA].minus  = 0.0;  /* p0.E - p0.pz = 0   */
	particles[kPartonIndexIncomingA].rKine  = rKine;
	particles[kPartonIndexIncomingA].pol    = 0;

	particles[kPartonIndexIncomingB].mom    = p1;
	particles[kPartonIndexIncomingB].plus   = 0.0;  /* p1.E + p1.pz = 0  */
	particles[kPartonIndexIncomingB].minus  = -pm;  /* p1.E - p1.pz = -pm */
	particles[kPartonIndexIncomingB].rKine  = 0.0;
	particles[kPartonIndexIncomingB].pol    = 0;

	particles[kPartonIndexOutgoingA].mom    = p2;
	particles[kPartonIndexOutgoingA].plus   = p2Plus;   /* Qt * exp(yv)  */
	particles[kPartonIndexOutgoingA].minus  = p2Minus;  /* Qt * exp(-yv) */
	particles[kPartonIndexOutgoingA].rKine  = 0.0;
	particles[kPartonIndexOutgoingA].pol    = 0;

	particles[kPartonIndexOutgoingB].mom    = p3;
	particles[kPartonIndexOutgoingB].plus   = p3Plus;   /* sqrtPt2 * exp(y)  */
	particles[kPartonIndexOutgoingB].minus  = p3Minus;  /* sqrtPt2 * exp(-y) */
	particles[kPartonIndexOutgoingB].rKine  = 0.0;
	particles[kPartonIndexOutgoingB].pol    = 0;

	*weight = psw;
}

int
passesCuts(
	Particle    particles[kPartonIndexMax],
	double      ptmin,
	double      ymax)
{
	double pt2min = ptmin * ptmin;

	for (int parton = kPartonIndexOutgoingA; parton < kPartonIndexMax; parton++)
	{
		if (vec4PT2(particles[parton].mom) < pt2min)
		{
			return 0;
		}

		if (fabs(vec4Y(particles[parton].mom)) > ymax)
		{
			return 0;
		}
	}

	return 1;
}
