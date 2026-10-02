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
#include "pepper-complexArith.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	Indices of the four partons of the gg -> gg process, in the
 *		order the phase-space generator fills them.
 *
 *	Momenta follow pepper's all-outgoing convention, so the two incoming
 *	partons carry negative energy. `kPartonIndexMax` doubles as the parton
 *	multiplicity of the process.
 */
typedef enum
{
	kPartonIndexIncomingA = 0,
	kPartonIndexIncomingB,
	kPartonIndexOutgoingA,
	kPartonIndexOutgoingB,
	kPartonIndexMax,
} PartonIndex;

/**
 *	@brief	Particle with momentum, light-cone pair, and helicity polarization.
 *
 *	The light-cone pair (plus, minus) = (E + pz, E - pz) sits alongside the
 *	Cartesian Vec4 so dot products with beam-axis momenta can skip the
 *	cosh(y) + sinh(y) reconstruction, which loses precision on UxHw's
 *	distributional representation at large |y|. The phase-space generator
 *	fills plus/minus straight from exp(y) and exp(-y), never from sums of
 *	cosh and sinh.
 */
typedef struct
{
	Vec4    mom;
	double  plus;   /* E + pz, light-cone forward component	*/
	double  minus;  /* E - pz, light-cone backward component */
	double  rKine;  /* event-level r = exp(yv - y) = emy/emyv,
	                 * the SINGLE distribution the closed-form
	                 * gg->gg ME^2 depends on. Filled once per
	                 * event by the phase-space generator. Only
	                 * particles[0].rKine is read.		*/
	int pol;        /* +1 or -1				*/
} Particle;

/**
 *	@brief	Gluon current for the Berends-Giele recursion.
 *
 *	p2 caches (p)^2, built incrementally from the light-cone identity
 *		2 * a.b = a.plus * b.minus + a.minus * b.plus
 *			- 2 * a.px * b.px - 2 * a.py * b.py
 *	so no pair of distributional factors has to cancel via
 *	cosh^2 y - sinh^2 y = 1. On UxHw that keeps the sign of the
 *	propagator's p^2 correct even in the far tails of the rapidity
 *	distribution, where cosh(y) and sinh(y) are each ~ e^|y| / 2 and their
 *	raw sum can flip sign from representation error.
 */
typedef struct
{
	CVec4       j;      /* Current vector (complex four-vector)	*/
	Vec4        p;      /* Associated momentum (always real)	*/
	double      plus;   /* E + pz, light-cone forward		*/
	double      minus;  /* E - pz, light-cone backward		*/
	double      p2;     /* Cached (p)^2, light-cone arithmetic	*/
	WeylSpinor  kp;     /* Positive-helicity gauge spinor	*/
	WeylSpinor  km;     /* Negative-helicity gauge spinor	*/
} GluonCurrent;

/**
 *	@brief	Color- and helicity-summed |M|^2 for gg -> gg, via the
 *		Berends-Giele recursion.
 *
 *	@param	particles	: The four partons, two incoming then two outgoing.
 *	@param	g3		: Strong coupling constant sqrt(4*pi*alpha_s).
 *	@return			: The squared matrix element (real, positive).
 */
double
computeME2(Particle particles[kPartonIndexMax], double g3);

#ifdef __cplusplus
}
#endif
