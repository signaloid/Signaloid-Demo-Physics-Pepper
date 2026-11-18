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

#include <stddef.h>
#include "pepper-lhapdfLite.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	Output variable indices, selected with the `-S` flag.
 *
 *	    -S 0  kOutputVariableIndexXs : the integrated gg -> gg cross section
 *	                                   (a Scalar). In UxHw mode, the first
 *	                                   moment of the propagated
 *	                                   differential-cross-section
 *	                                   distribution. In Monte Carlo mode,
 *	                                   the mean of the per-event
 *	                                   differential cross sections.
 */
typedef enum
{
	kOutputVariableIndexXs = 0,
	kOutputVariableIndexMax,
} OutputVariableIndex;

/*
 *	Physics constants, shared by the UxHw and Monte Carlo cross-section
 *	kernels (pepper-uxhw.cpp and pepper-monteCarlo.cpp).
 */
#define kPepperConstantEcms     13000.0
#define kPepperConstantPtMin    30.0
#define kPepperConstantYMax     5.0
#define kPepperConstantMZ       91.1876
#define kPepperConstantQ2       (kPepperConstantMZ * kPepperConstantMZ)
#define kPepperConstantGeV2ToPb 3.89379656e8

/*
 *	Importance-sampling exponent for the pT^2 PeakedDist (matches the
 *	dominant pT^2 dependence of the matrix element).
 */
#define kPepperConstantAlphaImp 2.0

/**
 *	@brief	Load the CT18NLO PDF set and compute g3 = sqrt(4 * pi * alpha_s(MZ^2)).
 *
 *		Must be called once before calculateOutputUxHw() or
 *		calculateOutputMonteCarlo().
 *
 *	@param	dataPath	: Directory holding the PDF grid data
 *				  ("inputs" locally, "." on SCCE).
 */
void
initializeLHAPDF(const char * dataPath);

/**
 *	@brief	Clean up LHAPDF resources.
 */
void
cleanupLHAPDF(void);

/**
 *	@brief	The loaded PDF handle, valid after initializeLHAPDF(). Shared by
 *		both cross-section kernels.
 */
PepperPDF *
pepperKernelPDF(void);

/**
 *	@brief	alpha_s(MZ^2), valid after initializeLHAPDF(). Shared by both
 *		cross-section kernels.
 */
double
pepperKernelAlphasMZ2(void);

/**
 *	@brief	UxHw-mode kernel. Propagates three Uniform[0,1] phase-space axes
 *		through the differential cross section op-by-op in a single pass,
 *		then reads the integrated cross section off the first moment of
 *		the resulting joint. Writes the scalar into
 *		`outputVariables[kOutputVariableIndexXs]` and
 *		`monteCarloOutputSamples[0]`.
 *
 *	@param	outputVariables		: Array of size `kOutputVariableIndexMax` to fill.
 *	@param	monteCarloOutputSamples	: Single-element array for the result.
 *	@return	double			: Returns the integrated cross section (pb).
 */
double
calculateOutputUxHw(
	double *    outputVariables,
	double *    monteCarloOutputSamples);

/**
 *	@brief	Monte Carlo kernel. Draws
 *		`numberOfMonteCarloIterations` independent
 *		phase-space samples into `monteCarloOutputSamples` with plain
 *		scalar arithmetic (no UxHw distributional API) and reports their
 *		mean as the integrated cross section. Writes the scalar into
 *		`outputVariables[kOutputVariableIndexXs]`.
 *
 *	@param	numberOfMonteCarloIterations	: Number of Monte Carlo iterations to run.
 *	@param	outputVariables		: Array of size `kOutputVariableIndexMax` to fill.
 *	@param	monteCarloOutputSamples	: `numberOfMonteCarloIterations` doubles, filled with per-event dsigma.
 *	@return	double			: Returns the integrated cross section (pb).
 */
double
calculateOutputMonteCarlo(
	size_t      numberOfMonteCarloIterations,
	double *    outputVariables,
	double *    monteCarloOutputSamples);

#ifdef __cplusplus
}
#endif
