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

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <uxhw.h>

extern "C"
{
#include "pepper-utilities.h"
#include "pepper-kernel.h"
#include "common.h"
}

#ifdef NO_OS_AVAILABLE

extern "C" void
returnZeroNoOS(void);
#endif

int
main(int argc, char * argv[])
{
	CommandLineArguments        arguments;
	double                      outputVariables[kOutputVariableIndexMax] = { 0 };
	const char *                outputVariableNames[kOutputVariableIndexMax] = {
		"crossSection",
	};
	const char *                outputVariableDescriptions[kOutputVariableIndexMax] = {
		"gg -> gg cross section (pb)",
	};
	kOutputVariableTypeIndex    outputVariableTypes[kOutputVariableIndexMax] = {
		kOutputVariableTypeScalar,
	};
	const char *                applicationDescription  = "Pepper: pp -> gg cross-section";
	double *                    monteCarloOutputSamples = NULL;
	clock_t                     start                   = 0;
	clock_t                     end                     = 0;
	double                      cpuTimeUsedInSeconds    = 0;

	memset(&arguments, 0, sizeof(arguments));

	/*
	 *	Get command line arguments. In the no-OS build there is no command
	 *	line, so this instead applies the hard-coded configuration in
	 *	`setNoOSCommandLineArguments()`.
	 */
	if (getCommandLineArguments(argc, argv, &arguments) != kCommonConstantReturnTypeSuccess)
	{
		return EXIT_FAILURE;
	}

	/*
	 *	Detect LHAPDF data path: on SCCE the cwd is inputs/, so the
	 *	CT18NLO set is found at "./CT18NLO/". For local builds the
	 *	cwd is the repo root, so the data is at "inputs/CT18NLO/".
	 */
	const char * dataPath;

#ifdef NO_OS_AVAILABLE
	/*
	 *	The no-OS build uses the compiled-in grid, so the path is unused.
	 */
	dataPath = ".";
#else

	if (fileExists("CT18NLO/CT18NLO.info"))
	{
		dataPath = ".";
	}
	else
	{
		dataPath = "inputs";
	}
#endif

	initializeLHAPDF(dataPath);

	/*
	 *	MonteCarlo output samples are used even in the UxHw case to store
	 *	the result of the single distributional pass (at index 0).
	 */
	monteCarloOutputSamples = (double *) checkedMalloc(
		arguments.common.numberOfMonteCarloIterations * sizeof(double),
		__FILE__,
		__LINE__
	);

	/*
	 *	Start timing if timing is enabled.
	 */
	if (arguments.common.isTimingEnabled)
	{
		start = clock();
	}

	/*
	 *	Dispatch to the mode-specific kernel.
	 */
	const bool isScalarOutputSelected = isSelectedOutputScalar(
		&arguments.common,
		outputVariableTypes,
		kOutputVariableIndexMax
	);

	if (arguments.common.isMonteCarloMode)
	{
		calculateOutputMonteCarlo(
			arguments.common.numberOfMonteCarloIterations,
			outputVariables,
			monteCarloOutputSamples
		);

		if (isScalarOutputSelected)
		{
			monteCarloOutputSamples[0] = outputVariables[kOutputVariableIndexXs];
		}
	}
	else
	{
		calculateOutputUxHw(outputVariables, monteCarloOutputSamples);
	}

	/*
	 *	Stop timing.
	 */
	if (arguments.common.isTimingEnabled)
	{
		end = clock();
		cpuTimeUsedInSeconds = ((double) (end - start)) / CLOCKS_PER_SEC;
	}

	if (arguments.common.isOutputJSONMode)
	{
		printJSONFormattedOutput(
			&arguments.common,
			monteCarloOutputSamples,
			outputVariables,
			outputVariableNames,
			kOutputVariableIndexMax,
			applicationDescription,
			outputVariableTypes
		);
	}
	else
	{
		printHumanConsumableOutput(
			&arguments.common,
			kOutputVariableIndexMax,
			outputVariables,
			outputVariableNames,
			outputVariableDescriptions,
			monteCarloOutputSamples,
			outputVariableTypes
		);
	}

	/*
	 *	Print timing if timing is enabled.
	 */
	if (arguments.common.isTimingEnabled)
	{
		printf("\nCPU time used: %" SignaloidParticleModifier "lf seconds\n", cpuTimeUsedInSeconds);
	}

	/*
	 *	Save Monte Carlo data to "data.out" if in Monte Carlo mode. The
	 *	cross section is a scalar (the cube integral), so we save only the
	 *	mean (written into monteCarloOutputSamples[0] above).
	 */
	if (arguments.common.isMonteCarloMode)
	{
		size_t samplesToSave = isScalarOutputSelected
		            ? 1
		            : arguments.common.numberOfMonteCarloIterations;

		saveMonteCarloDoubleDataToDataDotOutFile(
			monteCarloOutputSamples,
			(uint64_t) (cpuTimeUsedInSeconds * 1000000),
			samplesToSave
		);
	}
	else
	{
		if (arguments.common.isWriteToFileEnabled)
		{
			if (writeOutputDoubleDistributionsToCSV(
					arguments.common.outputFilePath,
					outputVariables,
					outputVariableNames,
					kOutputVariableIndexMax
			))
			{
				fprintf(
					stderr, "Error: Could not write to output CSV file \"%s\".\n",
					arguments.common.outputFilePath
				);

				return EXIT_FAILURE;
			}
		}
	}

	/*
	 *	Free allocations and clean up.
	 */
	free(monteCarloOutputSamples);

	cleanupLHAPDF();
#ifdef NO_OS_AVAILABLE
	returnZeroNoOS();
#else

	return EXIT_SUCCESS;

#endif
}
