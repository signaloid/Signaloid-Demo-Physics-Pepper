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

#include "pepper-berendsGiele.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	Generate a 2 -> 2 phase-space point using the t-channel
 *		parameterization from arXiv:2302.10449.
 *
 *		Draws its randomness from 4 calls to UxHwDoubleUniformDist(0, 1),
 *		which return full Uniform[0,1] distributions on Signaloid
 *		hardware and single samples on the compatibility layer.
 *
 *	@param	particles	: Output array of 4 partons, momenta filled in.
 *	@param	ecms		: Center-of-mass energy (GeV).
 *	@param	pt2min		: Minimum pT^2 (GeV^2).
 *	@param	ymax		: Maximum absolute rapidity.
 *	@param	weight		: Output phase-space weight.
 */
void
generatePhaseSpacePoint(
	Particle    particles[kPartonIndexMax],
	double      ecms,
	double      pt2min,
	double      ymax,
	double *    weight);

/**
 *	@brief	Check whether all final-state partons pass pT and rapidity cuts.
 *
 *	@param	particles	: The four partons.
 *	@param	ptmin		: Minimum transverse momentum (GeV).
 *	@param	ymax		: Maximum absolute rapidity.
 *	@return			: 1 if all pass, 0 if any fail.
 */
int
passesCuts(
	Particle    particles[kPartonIndexMax],
	double      ptmin,
	double      ymax);

#ifdef __cplusplus
}
#endif
