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
 *	pepper-kernel.cpp
 *
 *	Shared LHAPDF lifecycle and the mode dispatch for the gg -> gg cross
 *	section. Two computation kernels sit behind a common interface:
 *
 *	  calculateOutputUxHw          single distributional pass, contracted
 *	                               to the first moment (pepper-uxhw.cpp).
 *	  calculateOutputMonteCarlo    conventional Monte Carlo sampling,
 *	                               contracted to the sample mean
 *	                               (pepper-monteCarlo.cpp).
 *
 *	main() (pepper.cpp) selects between them based on `isMonteCarloMode`.
 *	========================================================================
 */

#include <cmath>
#include <cstdlib>
#include <cstdio>

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
#include "pepper-monteCarlo.h"
#include "pepper-lhapdfLite.h"
#include "common.h"
}

/*
 *	Module-level state initialised by initializeLHAPDF(). The lhapdfLite
 *	reader (pure C, no iostream / YAML / exceptions) loads the CT18NLO
 *	grid and provides xf(x, Q^2) and alpha_s(Q^2) via log-bicubic and
 *	cubic interpolation respectively. Shared by both cross-section kernels
 *	via the pepperKernelPDF() / pepperKernelAlphasMZ2() accessors.
 */
static PepperPDF *  gPdf        = NULL;
static double       gG3         = 0.0;
static double       gAlphasMZ2  = 0.0;

void
initializeLHAPDF(const char * dataPath)
{
	const char *    pdfSetName = "CT18NLO";
	char            pdfDir[1024];

	snprintf(pdfDir, sizeof(pdfDir), "%s/%s", dataPath, pdfSetName);

	gPdf = pepperPDFLoad(pdfDir, pdfSetName, 0);

	if (gPdf == NULL)
	{
		fprintf(
			stderr, "Error: Failed to load PDF set \"%s\" from \"%s\".\n",
			pdfSetName, pdfDir
		);
		exit(EXIT_FAILURE);
	}

	gAlphasMZ2  = pepperPDFAlphasQ2(gPdf, kPepperConstantQ2);
	gG3         = sqrt(4.0 * M_PI * gAlphasMZ2);
}

void
cleanupLHAPDF(void)
{
	if (gPdf != NULL)
	{
		pepperPDFFree(gPdf);
		gPdf = NULL;
	}
}

PepperPDF *
pepperKernelPDF(void)
{
	return gPdf;
}

double
pepperKernelAlphasMZ2(void)
{
	return gAlphasMZ2;
}

double
calculateOutputUxHw(
	double *    outputVariables,
	double *    monteCarloOutputSamples)
{
	/*
	 *	dxs is the differential-cross-section distribution over the
	 *	Uniform[0,1]^3 input cube. The full propagated joint on SCCE, a
	 *	single sample on the compat (native) layer.
	 */
	double dxs = generateXsPointUxHw();

	/*
	 *	Scalar output: the integrated cross section is the first moment
	 *	of dxs. On SCCE this is the in-hardware cube integral. The NaN
	 *	check is a scalar test on the contracted moment (never a branch on
	 *	a distribution): it catches the rare compat-layer sample that fell
	 *	deep in the rejection region where log()/sqrt() produced NaN.
	 */
	double integral = UxHwDoubleNthMoment(dxs, 1);

	if (integral != integral)
	{
		integral = 0.0;
	}

	outputVariables[kOutputVariableIndexXs] = integral;
	monteCarloOutputSamples[0]              = integral;

	return integral;
}

double
calculateOutputMonteCarlo(
	size_t      numberOfMonteCarloIterations,
	double *    outputVariables,
	double *    monteCarloOutputSamples)
{
	MeanAndVariance meanAndVariance;

	/*
	 *	Draw the per-event differential cross section samples, then report
	 *	the integrated cross section as their sample mean (which converges
	 *	to the cube integral).
	 */
	generateXsMonteCarlo(numberOfMonteCarloIterations, monteCarloOutputSamples);

	meanAndVariance = calculateMeanAndVarianceOfDoubleSamples(
		monteCarloOutputSamples,
		numberOfMonteCarloIterations
	);

	outputVariables[kOutputVariableIndexXs] = meanAndVariance.mean;

	return meanAndVariance.mean;
}
