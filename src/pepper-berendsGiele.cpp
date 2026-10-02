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
#include "pepper-berendsGiele.h"

/**
 *	@brief	The two color orderings of the outgoing pair summed over by the
 *		Berends-Giele amplitude, indexing both the permutation table and
 *		the 2x2 color-factor matrix that contracts the two amplitudes.
 */
typedef enum
{
	kColorOrderingIndexIdentity = 0,
	kColorOrderingIndexSwapped,
	kColorOrderingIndexMax,
} ColorOrderingIndex;

/**
 *	@brief	The independent helicity configurations of gg -> gg. The
 *		remaining configurations are related to these by parity and are
 *		accounted for by the color factors.
 */
typedef enum
{
	kHelicityConfigIndexPPMM = 0,
	kHelicityConfigIndexPMMP,
	kHelicityConfigIndexPMPM,
	kHelicityConfigIndexMax,
} HelicityConfigIndex;

/**
 *	@brief	Number of outgoing partons, i.e. those the color orderings
 *		permute. The two incoming partons keep their positions.
 */
#define kOutgoingPartonCount (kPartonIndexMax - kPartonIndexOutgoingA)

/**
 *	@brief	Gauge vector fixed for all external gluons: k = (1, 0, 1, 0).
 */
static const Vec4 kGaugeVector = { { 1.0, 0.0, 1.0, 0.0 } };

/**
 *	@brief	Compute positive-helicity polarization vector EP(p) given gauge spinors.
 *		EP = VT(kp, p_minus) / (sqrt(2) * conj(<km | p_minus>))
 *		where p_minus = Weyl(-1, p).
 */
static CVec4
polarizationPlus(Vec4 p, WeylSpinor kp, WeylSpinor km)
{
	WeylSpinor  pm              = weylFromMomentum(-1, p);
	CVec4       polarization    = vt(kp, pm);
	Complex     denom           = complexScale(sqrt(2.0), complexConj(weylProduct(km, pm)));

	return cvec4Scale(complexDiv(complexFromReal(1.0), denom), polarization);
}

/**
 *	@brief	Compute negative-helicity polarization vector EM(p) given gauge spinors.
 *		EM = VT(p_plus, km) / (sqrt(2) * conj(<kp | p_plus>))
 *		where p_plus = Weyl(+1, p).
 */
static CVec4
polarizationMinus(Vec4 p, WeylSpinor kp, WeylSpinor km)
{
	WeylSpinor  pp              = weylFromMomentum(+1, p);
	CVec4       polarization    = vt(pp, km);
	Complex     denom           = complexScale(sqrt(2.0), complexConj(weylProduct(kp, pp)));

	return cvec4Scale(complexDiv(complexFromReal(1.0), denom), polarization);
}

/**
 *	@brief	Initialize an external gluon current from a particle momentum and helicity.
 *		Sets the gauge vector to k = (1, 0, 1, 0) and computes gauge spinors.
 *
 *		For incoming particles (negative energy), the positive-energy momentum |p|
 *		is used to construct the polarization vector so that it is transverse:
 *		epsilon · |p| = 0. The original (negative-energy) momentum is stored
 *		in current->p
 *		for the recursion momentum bookkeeping.
 */
static void
constructExternalCurrent(GluonCurrent * current, const Particle * particle, int helicity)
{
	current->p      = particle->mom;
	current->plus   = particle->plus;
	current->minus  = particle->minus;
	current->p2     = 0.0;  /* external gluons are lightlike */
	current->kp     = weylFromMomentum(+1, kGaugeVector);
	current->km     = weylFromMomentum(-1, kGaugeVector);

	if (helicity >= 0)
	{
		current->j = polarizationPlus(particle->mom, current->kp, current->km);
	}

	if (helicity <= 0)
	{
		current->j = polarizationMinus(particle->mom, current->kp, current->km);
	}
}

/**
 *	@brief	Multiply the current by the Feynman-gauge scalar propagator: -i / p^2.
 *
 *	On UxHw the distributional representation of p^2 typically has a small
 *	tail that crosses zero even when the analytic value is single-signed.
 *	A near-zero p^2 sample makes 1/p^2 explode, contaminating amplitudes
 *	by many orders of magnitude. We apply a Regge-style regulator that
 *	clamps |p^2| to be at least kP2MinAbs GeV^2, on whichever side of zero
 *	the distribution's mean lives. This biases the final cross-section
 *	mildly downward (cutting the t -> 0 tail that perturbation theory does
 *	not trust anyway) but keeps the UxHw mean numerically tractable.
 *	On native this branch picks the same side every sample, so the clamp
 *	is a no-op for any physical p^2 value typical phase-space produces.
 */
static const double kP2MinAbs = 1.0e2;  /* 100 GeV^2 ~ (1/9) * (pT cut)^2 */

static void
addPropagator(GluonCurrent * current)
{
	double p2 = current->p2;

	if (p2 >= 0.0)
	{
		p2 = UxHwDoubleLimitDistributionSupport(p2, kP2MinAbs, INFINITY);
	}
	else
	{
		p2 = UxHwDoubleLimitDistributionSupport(p2, -INFINITY, -kP2MinAbs);
	}

	/*
	 *	Multiply j by -i/p2:
	 *	(-i/p2) * (re + i*im) = (im/p2) + i*(-re/p2)
	 */
	Complex factor = complexMake(0.0, -1.0 / p2);

	current->j = cvec4Scale(factor, current->j);
}

/**
 *	@brief	Triple gluon vertex V3.
 *
 *	Implements:
 *	  jc = dot(a.j, b.j) * (a.p - b.p) + dot(a.j, 2*b.p + a.p) * b.j
 *	       - dot(b.j, 2*a.p + b.p) * a.j
 *	  c.j += -g3 * jc
 *	  c.p  = a.p + b.p
 */
static void
tripleGluonVertex(
	GluonCurrent *  leftCurrent,
	GluonCurrent *  rightCurrent,
	GluonCurrent *  target,
	double          g3)
{
	/*
	 *	Term 1: dot(a.j, b.j) * (a.p - b.p)
	 */
	Complex dAB     = cvec4Dot(leftCurrent->j, rightCurrent->j);
	Vec4    pDiff   = vec4Sub(leftCurrent->p, rightCurrent->p);
	CVec4   term1   = cvec4Scale(dAB, cvec4FromVec4(pDiff));

	/*
	 *	Term 2: dot(a.j, 2*b.p + a.p) * b.j
	 */
	Vec4    mom2    = vec4Add(vec4Scale(2.0, rightCurrent->p), leftCurrent->p);
	Complex dA2BpA  = cvec4Dot(leftCurrent->j, cvec4FromVec4(mom2));
	CVec4   term2   = cvec4Scale(dA2BpA, rightCurrent->j);

	/*
	 *	Term 3: dot(b.j, 2*a.p + b.p) * a.j
	 */
	Vec4    mom3    = vec4Add(vec4Scale(2.0, leftCurrent->p), rightCurrent->p);
	Complex dB2ApB  = cvec4Dot(rightCurrent->j, cvec4FromVec4(mom3));
	CVec4   term3   = cvec4Scale(dB2ApB, leftCurrent->j);

	/*
	 *	jc = term1 + term2 - term3
	 *	c.j += -g3 * jc
	 */
	CVec4 vertexCurrent = cvec4Sub(cvec4Add(term1, term2), term3);

	cvec4AddAssign(&target->j, cvec4ScaleReal(-g3, vertexCurrent));

	/*
	 *	Composite momentum (Cartesian) and its cached light-cone pair.
	 *	All four are simple sums that preserve structure on UxHw.
	 */
	target->p       = vec4Add(leftCurrent->p, rightCurrent->p);
	target->plus    = leftCurrent->plus + rightCurrent->plus;
	target->minus   = leftCurrent->minus + rightCurrent->minus;

	/*
	 *	Cached (p)^2 via the light-cone identity
	 *		2 * a.b = a.plus * b.minus + a.minus * b.plus
	 *			- 2 * a.px * b.px - 2 * a.py * b.py
	 *	(a+b)^2 = a^2 + b^2 + 2 * a.b.
	 *
	 *	The light-cone term a.plus * b.minus + a.minus * b.plus is a sum
	 *	of two cross-products of *independent* distributional quantities
	 *	(no a/a self-multiplication, no cosh-sinh cancellation), so UxHw
	 *	preserves its sign exactly. For every n=2 composite that receives
	 *	a propagator in the 4-gluon case at least one constituent is along
	 *	the beam axis (px = py = 0), which zeroes the transverse cross
	 *	terms and keeps the result free of cos^2 + sin^2 cancellations too.
	 */
	double twoAB = leftCurrent->plus * rightCurrent->minus
	               + leftCurrent->minus * rightCurrent->plus
	               - 2.0 * leftCurrent->p.p[kFourVectorIndexPx] * rightCurrent->p.p[kFourVectorIndexPx]
	               - 2.0 * leftCurrent->p.p[kFourVectorIndexPy] * rightCurrent->p.p[kFourVectorIndexPy];

	target->p2 = leftCurrent->p2 + rightCurrent->p2 + twoAB;
}

/**
 *	@brief	Quartic gluon vertex V4.
 *
 *	Implements:
 *	  jc = dot(e.j, b.j) * a.j + dot(e.j, a.j) * b.j - 2 * dot(a.j, b.j) * e.j
 *	  c.j += i * g3^2 * jc
 */
static void
quarticGluonVertex(
	GluonCurrent *  leftCurrent,
	GluonCurrent *  middleCurrent,
	GluonCurrent *  rightCurrent,
	GluonCurrent *  target,
	double          g3)
{
	Complex dEB = cvec4Dot(middleCurrent->j, rightCurrent->j);
	Complex dEA = cvec4Dot(middleCurrent->j, leftCurrent->j);
	Complex dAB = cvec4Dot(leftCurrent->j, rightCurrent->j);

	CVec4   term1   = cvec4Scale(dEB, leftCurrent->j);
	CVec4   term2   = cvec4Scale(dEA, rightCurrent->j);
	CVec4   term3   = cvec4Scale(complexScale(2.0, dAB), middleCurrent->j);

	CVec4 vertexCurrent = cvec4Sub(cvec4Add(term1, term2), term3);

	/*
	 *	Multiply by i * g3^2.
	 *	Note: the quartic vertex does not have the extra sqrt(2) per coupling
	 *	that the triple vertex has.
	 */
	Complex factor = complexMake(0.0, g3 * g3);

	cvec4AddAssign(&target->j, cvec4Scale(factor, vertexCurrent));
}

/**
 *	@brief	Run the full Berends-Giele recursion for 4 gluons and return
 *		the color-ordered amplitude.
 *
 *	The recursion fills a 4x4 grid currents[begin][end]:
 *	- Diagonal currents[i][i]: external currents, set from particle data.
 *	- Off-diagonal currents[begin][begin+span-1]: span-particle currents
 *	  built from the V3 and V4 vertices.
 *
 *	@param	event	: Array of 4 particles with momenta and helicities.
 *	@param	g3	: Strong coupling sqrt(4*pi*alpha_s).
 *	@return		: Color-ordered amplitude A (complex).
 */
static Complex
berendsGieleAmplitude(Particle event[kPartonIndexMax], double g3)
{
	GluonCurrent    currents[kPartonIndexMax][kPartonIndexMax];
	int             gluonCount = kPartonIndexMax;

	/*
	 *	Initialize diagonal: external currents.
	 */
	for (int gluon = 0; gluon < gluonCount; gluon++)
	{
		constructExternalCurrent(&currents[gluon][gluon], &event[gluon], event[gluon].pol);
	}

	/*
	 *	Bottom-up recursion: build span-particle currents for
	 *	span = 2, 3, ..., gluonCount-1.
	 */
	for (int span = 2; span < gluonCount; span++)
	{
		for (int begin = 0; begin < gluonCount - span; begin++)
		{
			currents[begin][begin + span - 1].j = cvec4Zero();

			for (int split = begin; split < begin + span - 1; split++)
			{
				tripleGluonVertex(
					&currents[begin][split],
					&currents[split + 1][begin + span - 1],
					&currents[begin][begin + span - 1],
					g3
				);

				for (int secondSplit = split + 1; secondSplit < begin + span - 1; secondSplit++)
				{
					quarticGluonVertex(
						&currents[begin][split],
						&currents[split + 1][secondSplit],
						&currents[secondSplit + 1][begin + span - 1],
						&currents[begin][begin + span - 1],
						g3
					);
				}
			}

			/*
			 *	Add propagator for intermediate currents (not the final one).
			 */
			if (span < gluonCount - 1)
			{
				addPropagator(&currents[begin][begin + span - 1]);
			}
		}
	}

	/*
	 *	The amplitude is the dot product of the two (gluonCount-1)-particle
	 *	currents.
	 */
	return cvec4Dot(currents[gluonCount - 1][gluonCount - 1].j, currents[0][gluonCount - 2].j);
}

/*
 *	Closed-form gg -> gg matrix element. Color- and spin-averaged with the
 *	identical-particle 1/2! included, so it can be used as a drop-in
 *	replacement for the original Berends-Giele computeME2 below.
 *
 *	    |M|^2 = (9/4) g_s^4 (3 - tu/s^2 - us/t^2 - st/u^2)
 *
 *	The Mandelstam invariants are read off pepper's all-outgoing convention
 *	as products of light-cone variables stored in the Particle struct by
 *	the phase-space generator:
 *	    s = (p0+p1)^2 = pp*pm                = particles[kPartonIndexIncomingA].plus * particles[kPartonIndexIncomingB].minus
 *	    t = (p0+p2)^2 = -pp*sqrtPt2*exp(-yv) = particles[kPartonIndexIncomingA].plus * particles[kPartonIndexOutgoingA].minus
 *	    u = (p0+p3)^2 = -pp*sqrtPt2*exp(-y)  = particles[kPartonIndexIncomingA].plus * particles[kPartonIndexOutgoingB].minus
 *	Each is a product of two *independent* clean light-cone quantities,
 *	with no spinor products, cosh/sinh identities, intermediate complex
 *	four-vectors, or propagator 1/p^2, which the algebraic reduction has
 *	folded into the matrix element.
 *
 *	Verified to match the original Berends-Giele computeME2 (now
 *	`computeME2_bg`) to ~6 decimals across the integrand.
 *
 */
double
computeME2(Particle particles[kPartonIndexMax], double g3)
{
	/*
	 *	Closed-form gg -> gg |M|^2, written as a function of the SINGLE
	 *	event-level distribution r = exp(yv - y).
	 *
	 *	The textbook expression
	 *		(9/4) g_s^4 (3 - tu/s^2 - us/t^2 - st/u^2)
	 *	algebraically collapses to
	 *		(9/4) g_s^4 (3 + r + 1/r + r^2 + 1/r^2 - r/(r+1)^2)
	 *	once you substitute
	 *		s =  pp * Qt * (emy + emyv)
	 *		t = -pp * Qt *        emyv
	 *		u = -pp * Qt *        emy
	 *	and observe that pp^2 and Qt^2 cancel in every ratio.
	 *
	 */
	double r = particles[kPartonIndexIncomingA].rKine;

	/*
	 *	Clamp r to its physical range. r = exp(yv - y) > 0 by
	 *	construction, but UxHw's representation can place
	 *	a small fraction of its support at slightly negative or
	 *	near-zero values. 1/r and 1/r^2 then spike with the wrong sign
	 *	and skew the bracket mean, which came out 0.44 x truth on
	 *	Jupiter without the clamp. With |y| < 5 and |yv| < 5, r is
	 *	bounded above by exp(10) ~ 22 000, so [1e-5, 1e5] is
	 *	conservative.
	 */
	r = UxHwDoubleLimitDistributionSupport(r, 1.0e-5, 1.0e5);

	double  r2      = r * r;
	double  invR    = 1.0 / r;
	double  invR2   = invR * invR;
	double  rp1     = r + 1.0;
	double  rp1sq   = rp1 * rp1;

	double bracket = 3.0 + r + invR + r2 + invR2 - r / rp1sq;

	double  g_s2    = g3 * g3;
	double  g_s4    = g_s2 * g_s2;
	double  me2     = (9.0 / 4.0) * g_s4 * bracket;

	return me2;
}

/*
 *	Original Berends-Giele path, kept for reference and debugging.
 *	generateXsPoint does not call it. It stays in the translation unit so a
 *	future build can #ifdef-switch between the two to validate the
 *	closed-form replacement. The formula above is the algebraic reduction
 *	of the loop below for the 4-gluon color- and spin-averaged sum,
 *	including the identical-particle 1/2! factor.
 */
__attribute__((unused))
static double
computeME2_bg(Particle particles[kPartonIndexMax], double g3)
{
	int     orderings[kColorOrderingIndexMax][kOutgoingPartonCount] = {
		{ kPartonIndexOutgoingA, kPartonIndexOutgoingB },
		{ kPartonIndexOutgoingB, kPartonIndexOutgoingA }
	};
	double  colorFactors[kColorOrderingIndexMax][kColorOrderingIndexMax] = {
		{ 72.0, 36.0 },
		{ 36.0, 72.0 }
	};
	int     helicities[kHelicityConfigIndexMax][kPartonIndexMax] = {
		{ +1, +1, -1, -1 },
		{ +1, -1, -1, +1 },
		{ +1, -1, +1, -1 }
	};

	double me2 = 0.0;

	for (int helicityConfig = 0; helicityConfig < kHelicityConfigIndexMax; helicityConfig++)
	{
		Complex amplitudes[kColorOrderingIndexMax];

		for (int ordering = 0; ordering < kColorOrderingIndexMax; ordering++)
		{
			Particle permuted[kPartonIndexMax];

			permuted[kPartonIndexIncomingA]     = particles[kPartonIndexIncomingA];
			permuted[kPartonIndexIncomingA].pol = helicities[helicityConfig][kPartonIndexIncomingA];
			permuted[kPartonIndexIncomingB]     = particles[kPartonIndexIncomingB];
			permuted[kPartonIndexIncomingB].pol = helicities[helicityConfig][kPartonIndexIncomingB];

			for (int outgoing = 0; outgoing < kOutgoingPartonCount; outgoing++)
			{
				int partonIndex = orderings[ordering][outgoing];

				permuted[kPartonIndexOutgoingA + outgoing]      = particles[partonIndex];
				permuted[kPartonIndexOutgoingA + outgoing].pol  = helicities[helicityConfig][partonIndex];
			}

			amplitudes[ordering] = berendsGieleAmplitude(permuted, g3);
		}

		for (int ordering = 0; ordering < kColorOrderingIndexMax; ordering++)
		{
			for (int otherOrdering = 0; otherOrdering < kColorOrderingIndexMax; otherOrdering++)
			{
				Complex prod = complexMul(amplitudes[ordering], complexConj(amplitudes[otherOrdering]));

				me2 += colorFactors[ordering][otherOrdering] * prod.re;
			}
		}
	}

	me2 *= 2.0 / (256.0 * 2.0);

	return me2;
}
