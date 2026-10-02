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
 *	Lightweight reimplementation of the subset of LHAPDF6 needed for Pepper:
 *	- Load a single lhagrid1 .dat file (one subgrid, one member)
 *	- Log-bicubic interpolation of xf(x, Q^2) for a single parton flavor
 *	- Cubic interpolation of alpha_s(Q^2) from .info metadata
 *
 *	Pure C. No C++ exceptions, no iostream, no YAML parser, no stat/opendir.
 *	Parses .info as plain text (key: value lines) instead of YAML.
 */

#include "pepper-lhapdfLite.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define kMaxLine        4096
#define kMaxKnots       256
#define kMaxAlphaSKnots 64
#define kMaxFlavors     13

/*
 *	Value of alpha_s(M_Z) used when a PDF set ships fewer than two
 *	alpha_s knots, leaving nothing to interpolate. This is the PDG world
 *	average, which is what LHAPDF falls back to in the same situation.
 */
#define kAlphaSMZFallback 0.118

/*
 *	Cap on the magnitude of an interpolated alpha_s. Any |alpha_s| at or
 *	above this is unphysical (the coupling is O(0.1) across the whole
 *	grid) and indicates the interpolation has run away in an
 *	extrapolation region, so the result is clamped rather than returned.
 */
#define kAlphaSMaxMagnitude 2.0

struct PepperPDF
{
	/*
	 *	Grid knots.
	 */
	double  xs[kMaxKnots];
	double  logxs[kMaxKnots];
	double  q2s[kMaxKnots];
	double  logq2s[kMaxKnots];
	size_t  nX;
	size_t  nQ2;

	/*
	 *	Flavor mapping: pids[i] is the PDG ID for column i.
	 */
	int     pids[kMaxFlavors];
	size_t  nPids;

	/*
	 *	Grid data: xf[ix * nQ2 + iq2] for the gluon (pid 21).
	 *	Stored x-major to match LHAPDF's internal layout.
	 */
	double * xf;

	/*
	 *	Alpha_s interpolation knots from .info file.
	 */
	double  asQ[kMaxAlphaSKnots];
	double  asVal[kMaxAlphaSKnots];
	double  asLogQ2[kMaxAlphaSKnots];
	size_t  nAlphaS;
};

/*
 *	Hermite cubic interpolation (same formula as LHAPDF).
 *
 *	The four weights are the standard cubic Hermite basis polynomials on
 *	t in [0, 1], named here as they are in the literature:
 *
 *	    p(t) = h00(t)*valueLow + h10(t)*derivLow
 *	           + h01(t)*valueHigh + h11(t)*derivHigh
 *
 *	h00 and h01 weight the values at the low and high ends of the
 *	interval, h10 and h11 the corresponding derivatives. Their
 *	coefficients are fixed by the interpolation conditions
 *	p(0) = valueLow, p(1) = valueHigh, p'(0) = derivLow,
 *	p'(1) = derivHigh, so they define the
 *	basis rather than being tunable parameters.
 */
static double
interpolateCubic(double t, double valueLow, double derivLow, double valueHigh, double derivHigh)
{
	double  t2  = t * t;
	double  t3  = t2 * t;

	double  h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
	double  h10 = t3 - 2.0 * t2 + t;
	double  h01 = -2.0 * t3 + 3.0 * t2;
	double  h11 = t3 - t2;

	return h00 * valueLow + h10 * derivLow + h01 * valueHigh + h11 * derivHigh;
}

/*
 *	Find index i such that knots[i] <= value < knots[i+1].
 *	Returns 0 below the first knot, knotCount-2 above the last.
 */
static size_t
findBelow(const double * knots, size_t knotCount, double value)
{
	if (value <= knots[0])
	{
		return 0;
	}

	if (value >= knots[knotCount - 1])
	{
		return knotCount - 2;
	}

	size_t  low     = 0;
	size_t  high    = knotCount - 1;

	while (high - low > 1)
	{
		size_t middle = (low + high) / 2;

		if (knots[middle] <= value)
		{
			low = middle;
		}
		else
		{
			high = middle;
		}
	}

	return low;
}

/*
 *	Derivative dxf/dlogx at grid point (xIndex, q2Index), by central,
 *	forward, or backward difference at the interior, low edge, and high edge.
 */
static double
dxfDlogx(const PepperPDF * pdf, size_t xIndex, size_t q2Index)
{
	size_t xKnotCount = pdf->nX;

	if (xIndex > 0 && xIndex < xKnotCount - 1)
	{
		double backwardSlope = (pdf->xf[xIndex * pdf->nQ2 + q2Index] - pdf->xf[(xIndex - 1) * pdf->nQ2 + q2Index])
		                       / (pdf->logxs[xIndex] - pdf->logxs[xIndex - 1]);
		double forwardSlope = (pdf->xf[(xIndex + 1) * pdf->nQ2 + q2Index] - pdf->xf[xIndex * pdf->nQ2 + q2Index])
		                      / (pdf->logxs[xIndex + 1] - pdf->logxs[xIndex]);

		return (backwardSlope + forwardSlope) / 2.0;
	}
	else if (xIndex == 0)
	{
		return (pdf->xf[1 * pdf->nQ2 + q2Index] - pdf->xf[0 * pdf->nQ2 + q2Index])
		       / (pdf->logxs[1] - pdf->logxs[0]);
	}
	else
	{
		return (pdf->xf[xIndex * pdf->nQ2 + q2Index] - pdf->xf[(xIndex - 1) * pdf->nQ2 + q2Index])
		       / (pdf->logxs[xIndex] - pdf->logxs[xIndex - 1]);
	}
}

/*
 *	Parse the alpha_s Q and value arrays from the .info file.
 *	Looks for "AlphaS_Qs:" and "AlphaS_Vals:" lines.
 */
static int
parseAlphaSFromInfo(PepperPDF * pdf, const char * infoPath)
{
	FILE *  f = fopen(infoPath, "r");
	char    line[kMaxLine];
	int     foundQs     = 0;
	int     foundVals   = 0;

	if (f == NULL)
	{
		return -1;
	}

	while (fgets(line, kMaxLine, f) != NULL)
	{
		if (strncmp(line, "AlphaS_Qs:", 10) == 0)
		{
			/*
			 *	Parse array: [ val, val, ... ]
			 */
			char * cursor = strchr(line, '[');

			if (cursor == NULL)
			{
				continue;
			}

			cursor++;
			pdf->nAlphaS = 0;

			while (pdf->nAlphaS < kMaxAlphaSKnots)
			{
				char *  end;
				double  parsedValue = strtod(cursor, &end);

				if (end == cursor)
				{
					break;
				}

				pdf->asQ[pdf->nAlphaS]      = parsedValue;
				pdf->asLogQ2[pdf->nAlphaS]  = log(parsedValue * parsedValue);
				pdf->nAlphaS++;
				cursor = end;

				while (*cursor == ',' || *cursor == ' ' || *cursor == '\t')
				{
					cursor++;
				}

				if (*cursor == ']')
				{
					break;
				}
			}

			foundQs = 1;
		}
		else if (strncmp(line, "AlphaS_Vals:", 12) == 0)
		{
			char * cursor = strchr(line, '[');

			if (cursor == NULL)
			{
				continue;
			}

			cursor++;
			size_t n = 0;

			while (n < kMaxAlphaSKnots)
			{
				char *  end;
				double  parsedValue = strtod(cursor, &end);

				if (end == cursor)
				{
					break;
				}

				pdf->asVal[n] = parsedValue;
				n++;
				cursor = end;

				while (*cursor == ',' || *cursor == ' ' || *cursor == '\t')
				{
					cursor++;
				}

				if (*cursor == ']')
				{
					break;
				}
			}

			foundVals = 1;
		}
	}

	fclose(f);

	return (foundQs && foundVals) ? 0 : -1;
}

PepperPDF *
pepperPDFLoad(const char * dataDir, const char * setName, int member)
{
	PepperPDF * pdf = (PepperPDF *) calloc(1, sizeof(PepperPDF));
	char        path[1024];
	char        line[kMaxLine];
	FILE *      f;
	int         gluonCol;

	if (pdf == NULL)
	{
		return NULL;
	}

	/*
	 *	Load the .dat grid file.
	 */
	snprintf(path, sizeof(path), "%s/%s_%04d.dat", dataDir, setName, member);
	f = fopen(path, "r");

	if (f == NULL)
	{
		fprintf(stderr, "pepper-lhapdfLite: cannot open %s\n", path);
		free(pdf);

		return NULL;
	}

	/*
	 *	Skip header lines until the first "---" separator.
	 */
	while (fgets(line, kMaxLine, f) != NULL)
	{
		if (strncmp(line, "---", 3) == 0)
		{
			break;
		}
	}

	/*
	 *	Read x knots.
	 */
	if (fgets(line, kMaxLine, f) == NULL)
	{
		goto fail;
	}

	{
		char * cursor = line;

		pdf->nX = 0;

		while (pdf->nX < kMaxKnots)
		{
			char *  end;
			double  parsedValue = strtod(cursor, &end);

			if (end == cursor)
			{
				break;
			}

			pdf->xs[pdf->nX]    = parsedValue;
			pdf->logxs[pdf->nX] = log(parsedValue);
			pdf->nX++;
			cursor = end;
		}
	}

	/*
	 *	Read Q knots (file stores Q, we convert to Q^2).
	 */
	if (fgets(line, kMaxLine, f) == NULL)
	{
		goto fail;
	}

	{
		char * cursor = line;

		pdf->nQ2 = 0;

		while (pdf->nQ2 < kMaxKnots)
		{
			char *  end;
			double  parsedValue = strtod(cursor, &end);

			if (end == cursor)
			{
				break;
			}

			pdf->q2s[pdf->nQ2]      = parsedValue * parsedValue;
			pdf->logq2s[pdf->nQ2]   = log(parsedValue * parsedValue);
			pdf->nQ2++;
			cursor = end;
		}
	}

	/*
	 *	Read flavor PIDs.
	 */
	if (fgets(line, kMaxLine, f) == NULL)
	{
		goto fail;
	}

	{
		char * cursor = line;

		pdf->nPids = 0;

		while (pdf->nPids < kMaxFlavors)
		{
			char *  end;
			long    parsedValue = strtol(cursor, &end, 10);

			if (end == cursor)
			{
				break;
			}

			pdf->pids[pdf->nPids] = (int) parsedValue;
			pdf->nPids++;
			cursor = end;
		}
	}

	/*
	 *	Find the gluon column index.
	 *	Declared here (before any goto fail targets) for C++ compatibility.
	 */
	gluonCol = -1;

	for (size_t ii = 0; ii < pdf->nPids; ii++)
	{
		if (pdf->pids[ii] == 21)
		{
			gluonCol = (int) ii;

			break;
		}
	}

	if (gluonCol < 0)
	{
		fprintf(stderr, "pepper-lhapdfLite: gluon (pid 21) not found in flavor list\n");
		goto fail;
	}

	/*
	 *	Allocate xf array (x-major: xf[ix * nQ2 + iq2]).
	 */
	pdf->xf = (double *) malloc(pdf->nX * pdf->nQ2 * sizeof(double));

	if (pdf->xf == NULL)
	{
		goto fail;
	}

	/*
	 *	Read data rows. File is x-major: for each x, nQ2 rows.
	 */
	{
		size_t  row             = 0;
		size_t  expectedRows    = pdf->nX * pdf->nQ2;

		while (row < expectedRows && fgets(line, kMaxLine, f) != NULL)
		{
			if (strncmp(line, "---", 3) == 0)
			{
				break;
			}

			/*
			 *	Parse columns to find the gluon value.
			 */
			char *  cursor  = line;
			int     col     = 0;

			while (col <= gluonCol)
			{
				char *  end;
				double  parsedValue = strtod(cursor, &end);

				if (end == cursor)
				{
					break;
				}

				if (col == gluonCol)
				{
					pdf->xf[row] = parsedValue;
				}

				col++;
				cursor = end;
			}

			row++;
		}

		if (row != expectedRows)
		{
			fprintf(
				stderr, "pepper-lhapdfLite: expected %zu data rows, got %zu\n",
				expectedRows, row
			);
			goto fail;
		}
	}

	fclose(f);

	/*
	 *	Parse alpha_s from the .info file.
	 */
	snprintf(path, sizeof(path), "%s/%s.info", dataDir, setName);

	if (parseAlphaSFromInfo(pdf, path) != 0)
	{
		fprintf(stderr, "pepper-lhapdfLite: warning: could not parse alpha_s from %s\n", path);
	}

	fprintf(
		stderr, "pepper-lhapdfLite: loaded %s member %d (%zu x, %zu Q2, %zu alphaS knots)\n",
		setName, member, pdf->nX, pdf->nQ2, pdf->nAlphaS
	);

	return pdf;

fail:
	fclose(f);

	if (pdf->xf != NULL)
	{
		free(pdf->xf);
	}

	free(pdf);

	return NULL;
}

void
pepperPDFFree(PepperPDF * pdf)
{
	if (pdf != NULL)
	{
		free(pdf->xf);
		free(pdf);
	}

	return;
}

double
pepperPDFxfxQ2(const PepperPDF * pdf, int pid, double x, double q2)
{
	(void) pid; /* We only store the gluon column */

	double  logx    = log(x);
	double  logq2   = log(q2);

	size_t  xIndex  = findBelow(pdf->logxs, pdf->nX, logx);
	size_t  q2Index = findBelow(pdf->logq2s, pdf->nQ2, logq2);

	/*
	 *	Clamp indices so xIndex+1 and q2Index+1 are valid.
	 */
	if (xIndex >= pdf->nX - 1)
	{
		xIndex = pdf->nX - 2;
	}

	if (q2Index >= pdf->nQ2 - 1)
	{
		q2Index = pdf->nQ2 - 2;
	}

	size_t  q2Stride    = pdf->nQ2;
	size_t  q2IndexMax  = pdf->nQ2 - 1;

	/*
	 *	Fall back to bilinear if fewer than 4 Q2 knots.
	 */
	if (pdf->nQ2 < 4)
	{
		double  logx0           = pdf->logxs[xIndex];
		double  logx1           = pdf->logxs[xIndex + 1];
		double  valueAtQ2Low    = pdf->xf[xIndex * q2Stride + q2Index] + (logx - logx0) / (logx1 - logx0)
		                          * (pdf->xf[(xIndex + 1) * q2Stride + q2Index] - pdf->xf[xIndex * q2Stride + q2Index]);
		double valueAtQ2High = pdf->xf[xIndex * q2Stride + q2Index + 1] + (logx - logx0) / (logx1 - logx0)
		                       * (pdf->xf[(xIndex + 1) * q2Stride + q2Index + 1] - pdf->xf[xIndex * q2Stride + q2Index + 1]);

		return valueAtQ2Low + (logq2 - pdf->logq2s[q2Index]) / (pdf->logq2s[q2Index + 1] - pdf->logq2s[q2Index])
		       * (valueAtQ2High - valueAtQ2Low);
	}

	/*
	 *	Log-bicubic interpolation (matches LHAPDF's LogBicubicInterpolator).
	 *	First interpolate cubically in x at each of the 4 Q2 slices,
	 *	then interpolate cubically in Q2.
	 */
	double  dlogx1  = pdf->logxs[xIndex + 1] - pdf->logxs[xIndex];
	double  tlogx   = (logx - pdf->logxs[xIndex]) / dlogx1;

	/*
	 *	Cubic interpolation in x at q2Index and q2Index+1.
	 */
	double sliceLow = interpolateCubic(
		tlogx,
		pdf->xf[xIndex * q2Stride + q2Index], dxfDlogx(pdf, xIndex, q2Index) * dlogx1,
		pdf->xf[(xIndex + 1) * q2Stride + q2Index], dxfDlogx(pdf, xIndex + 1, q2Index) * dlogx1
	);
	double sliceHigh = interpolateCubic(
		tlogx,
		pdf->xf[xIndex * q2Stride + q2Index + 1], dxfDlogx(pdf, xIndex, q2Index + 1) * dlogx1,
		pdf->xf[(xIndex + 1) * q2Stride + q2Index + 1], dxfDlogx(pdf, xIndex + 1, q2Index + 1) * dlogx1
	);

	/*
	 *	Derivatives in Q2 for cubic Q2 interpolation.
	 */
	double  dlogq1  = pdf->logq2s[q2Index + 1] - pdf->logq2s[q2Index];
	double  tlogq   = (logq2 - pdf->logq2s[q2Index]) / dlogq1;
	double  sliceDerivLow, sliceDerivHigh;

	if (q2Index > 0 && q2Index + 1 < q2IndexMax)
	{
		double  dlogq0  = pdf->logq2s[q2Index] - pdf->logq2s[q2Index - 1];
		double  dlogq2  = pdf->logq2s[q2Index + 2] - pdf->logq2s[q2Index + 1];

		double sliceBelow = interpolateCubic(
			tlogx,
			pdf->xf[xIndex * q2Stride + q2Index - 1], dxfDlogx(pdf, xIndex, q2Index - 1) * dlogx1,
			pdf->xf[(xIndex + 1) * q2Stride + q2Index - 1], dxfDlogx(pdf, xIndex + 1, q2Index - 1) * dlogx1
		);
		double sliceAbove = interpolateCubic(
			tlogx,
			pdf->xf[xIndex * q2Stride + q2Index + 2], dxfDlogx(pdf, xIndex, q2Index + 2) * dlogx1,
			pdf->xf[(xIndex + 1) * q2Stride + q2Index + 2], dxfDlogx(pdf, xIndex + 1, q2Index + 2) * dlogx1
		);

		sliceDerivLow   = ((sliceHigh - sliceLow) / dlogq1 + (sliceLow - sliceBelow) / dlogq0) / 2.0;
		sliceDerivHigh  = ((sliceHigh - sliceLow) / dlogq1 + (sliceAbove - sliceHigh) / dlogq2) / 2.0;
	}
	else if (q2Index == 0)
	{
		double dlogq2 = pdf->logq2s[q2Index + 2] - pdf->logq2s[q2Index + 1];

		sliceDerivLow = (sliceHigh - sliceLow) / dlogq1;

		double sliceAbove = interpolateCubic(
			tlogx,
			pdf->xf[xIndex * q2Stride + q2Index + 2], dxfDlogx(pdf, xIndex, q2Index + 2) * dlogx1,
			pdf->xf[(xIndex + 1) * q2Stride + q2Index + 2], dxfDlogx(pdf, xIndex + 1, q2Index + 2) * dlogx1
		);

		sliceDerivHigh = (sliceDerivLow + (sliceAbove - sliceHigh) / dlogq2) / 2.0;
	}
	else
	{
		double dlogq0 = pdf->logq2s[q2Index] - pdf->logq2s[q2Index - 1];

		sliceDerivHigh = (sliceHigh - sliceLow) / dlogq1;

		double sliceBelow = interpolateCubic(
			tlogx,
			pdf->xf[xIndex * q2Stride + q2Index - 1], dxfDlogx(pdf, xIndex, q2Index - 1) * dlogx1,
			pdf->xf[(xIndex + 1) * q2Stride + q2Index - 1], dxfDlogx(pdf, xIndex + 1, q2Index - 1) * dlogx1
		);

		sliceDerivLow = (sliceDerivHigh + (sliceLow - sliceBelow) / dlogq0) / 2.0;
	}

	return interpolateCubic(tlogq, sliceLow, sliceDerivLow * dlogq1, sliceHigh, sliceDerivHigh * dlogq1);
}

double
pepperPDFAlphasQ2(const PepperPDF * pdf, double q2)
{
	if (pdf->nAlphaS < 2)
	{
		return kAlphaSMZFallback;
	}

	/*
	 *	Extrapolation below grid.
	 */
	if (q2 < pdf->asQ[0] * pdf->asQ[0])
	{
		size_t next = 1;

		while (next < pdf->nAlphaS &&
		       pdf->asQ[next] == pdf->asQ[0])
		{
			next++;
		}

		double  dlogq2          = log10(pdf->asQ[next] / pdf->asQ[0]) * 2.0;
		double  deltaLogAlphaS  = log10(pdf->asVal[next] / pdf->asVal[0]);
		double  logSlope        = deltaLogAlphaS / dlogq2;

		return pdf->asVal[0] * pow(q2 / (pdf->asQ[0] * pdf->asQ[0]), logSlope);
	}

	/*
	 *	Extrapolation above grid.
	 */
	if (q2 > pdf->asQ[pdf->nAlphaS - 1] * pdf->asQ[pdf->nAlphaS - 1])
	{
		return pdf->asVal[pdf->nAlphaS - 1];
	}

	/*
	 *	Find index in log(Q^2) array.
	 */
	double  logq2       = log(q2);
	size_t  knotIndex   = findBelow(pdf->asLogQ2, pdf->nAlphaS, logq2);

	if (knotIndex >= pdf->nAlphaS - 1)
	{
		knotIndex = pdf->nAlphaS - 2;
	}

	double  dlogq2      = pdf->asLogQ2[knotIndex + 1] - pdf->asLogQ2[knotIndex];
	double  fraction    = (logq2 - pdf->asLogQ2[knotIndex]) / dlogq2;

	/*
	 *	Derivatives for cubic interpolation.
	 */
	double slopeAtKnot, slopeAtNextKnot;

	if (knotIndex == 0)
	{
		slopeAtKnot = (pdf->asVal[knotIndex + 1] - pdf->asVal[knotIndex]) / dlogq2;
		double  spacingAbove    = pdf->asLogQ2[knotIndex + 2] - pdf->asLogQ2[knotIndex + 1];
		double  slopeAbove      = (pdf->asVal[knotIndex + 2] - pdf->asVal[knotIndex + 1]) / spacingAbove;

		slopeAtNextKnot = (slopeAtKnot + slopeAbove) / 2.0;
	}
	else if (knotIndex == pdf->nAlphaS - 2)
	{
		double  spacingBelow    = pdf->asLogQ2[knotIndex] - pdf->asLogQ2[knotIndex - 1];
		double  slopeBelow      = (pdf->asVal[knotIndex] - pdf->asVal[knotIndex - 1]) / spacingBelow;

		slopeAtKnot     = (slopeBelow + (pdf->asVal[knotIndex + 1] - pdf->asVal[knotIndex]) / dlogq2) / 2.0;
		slopeAtNextKnot = (pdf->asVal[knotIndex + 1] - pdf->asVal[knotIndex]) / dlogq2;
	}
	else
	{
		double  spacingBelow    = pdf->asLogQ2[knotIndex] - pdf->asLogQ2[knotIndex - 1];
		double  slopeBelow      = (pdf->asVal[knotIndex] - pdf->asVal[knotIndex - 1]) / spacingBelow;
		double  spacingAbove    = pdf->asLogQ2[knotIndex + 2] - pdf->asLogQ2[knotIndex + 1];
		double  slopeAbove      = (pdf->asVal[knotIndex + 2] - pdf->asVal[knotIndex + 1]) / spacingAbove;

		slopeAtKnot     = (slopeBelow + (pdf->asVal[knotIndex + 1] - pdf->asVal[knotIndex]) / dlogq2) / 2.0;
		slopeAtNextKnot = ((pdf->asVal[knotIndex + 1] - pdf->asVal[knotIndex]) / dlogq2 + slopeAbove) / 2.0;
	}

	double result = interpolateCubic(
		fraction, pdf->asVal[knotIndex], slopeAtKnot * dlogq2,
		pdf->asVal[knotIndex + 1], slopeAtNextKnot * dlogq2
	);

	return (fabs(result) < kAlphaSMaxMagnitude) ? result : kAlphaSMaxMagnitude;
}
