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

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	Draw one Monte Carlo sample of the gg -> gg differential cross
 *		section dsigma with plain scalar arithmetic (no UxHw
 *		distributional API). Three Uniform[0,1] samples drive the phase
 *		space, and the kinematic gates, clamps, and PDF lookups are
 *		ordinary scalar operations (fmin, fabs, step, direct lookup).
 *
 *	@return	double	: One per-event differential cross section sample (pb).
 */
double
generateXsPointMonteCarlo(void);

/**
 *	@brief	Fill `monteCarloOutputSamples` with `numberOfMonteCarloIterations`
 *		independent per-event differential cross section samples.
 *
 *	@param	numberOfMonteCarloIterations	: Number of samples to draw.
 *	@param	monteCarloOutputSamples		: Output array of that length.
 *	@return	double				: The last sample drawn.
 */
double
generateXsMonteCarlo(
	size_t      numberOfMonteCarloIterations,
	double *    monteCarloOutputSamples);

#ifdef __cplusplus
}
#endif
