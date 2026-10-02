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

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	Opaque handle to a loaded PDF grid.
 */
typedef struct PepperPDF PepperPDF;

/**
 *	@brief	Load a PDF grid from an LHAPDF lhagrid1 .dat file and .info file.
 *
 *	@param	dataDir		: Directory containing the PDF set (e.g. "inputs/CT18NLO-minimal").
 *	@param	setName		: PDF set name (e.g. "CT18NLO-minimal").
 *	@param	member		: Member number (0 for central).
 *	@return			: Handle to the loaded PDF, or NULL on failure.
 */
PepperPDF *
pepperPDFLoad(const char * dataDir, const char * setName, int member);

/**
 *	@brief	Free a loaded PDF grid.
 */
void
pepperPDFFree(PepperPDF * pdf);

/**
 *	@brief	Evaluate x * f(x, Q^2) for a given parton flavor.
 *
 *	Uses log-bicubic interpolation matching LHAPDF's "logcubic" interpolator.
 *
 *	@param	pdf	: Loaded PDF handle.
 *	@param	pid	: Parton ID (21 for gluon).
 *	@param	x	: Bjorken x.
 *	@param	q2	: Factorization scale Q^2 in GeV^2.
 *	@return		: x * f(x, Q^2).
 */
double
pepperPDFxfxQ2(const PepperPDF * pdf, int pid, double x, double q2);

/**
 *	@brief	Evaluate alpha_s(Q^2) using interpolation from the .info file.
 *
 *	@param	pdf	: Loaded PDF handle.
 *	@param	q2	: Scale Q^2 in GeV^2.
 *	@return		: alpha_s(Q^2).
 */
double
pepperPDFAlphasQ2(const PepperPDF * pdf, double q2);

#ifdef __cplusplus
}
#endif
