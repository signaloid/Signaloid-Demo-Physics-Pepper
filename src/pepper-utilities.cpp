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
#include <stdlib.h>
#include <string.h>
#include "pepper-utilities.h"
#include "common.h"

bool
fileExists(const char * path)
{
	FILE * f = fopen(path, "r");

	if (f != NULL)
	{
		fclose(f);

		return true;
	}

	return false;
}

void
printUsage(void)
{
	fprintf(stderr, "Pepper: pp -> gg cross-section with Signaloid UxHw distributional arithmetic\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: Valid command-line arguments are:\n");
	fprintf(
		stderr,
		"\t[-i, --input <Path to input CSV file : str>] (Read inputs from file.)\n"
		"\t[-o, --output <Path to output CSV file : str>] (Specify the output file.)\n"
		"\t[-S, --select-output <output : int>] (Compute 0-indexed output, by default 0. "
		"Set to the number of outputs to compute all outputs.)\n"
		"\t[-M, --multiple-executions <Number of executions : int (Default: 1)>] "
		"(Repeatedly execute kernel for benchmarking. Also enables timing.)\n"
		"\t[-T, --time] (Timing mode: Times and prints the timing of the kernel execution.)\n"
		"\t[-v, --verbose] (Verbose mode: Prints extra information about demo execution.)\n"
		"\t[-j, --json] (Print output in JSON format.)\n"
		"\t[-h, --help] (Display this help message.)\n"
	);
	fprintf(stderr, "\n");

	return;
}

CommonConstantReturnType
setDefaultCommandLineArguments(CommandLineArguments * arguments)
{
	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
	*arguments = (CommandLineArguments) {
		.common = (CommonCommandLineArguments) { 0 },
	};
#pragma GCC diagnostic pop

	return kCommonConstantReturnTypeSuccess;
}

#ifdef NO_OS_AVAILABLE
CommonConstantReturnType
setNoOSCommandLineArguments(CommandLineArguments * arguments)
{
	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

	/*
	 *	Start from the defaults so that every demo-specific field is
	 *	initialized, then override the ones the no-OS build fixes. The
	 *	`common` sub-struct is zeroed by this call, so it is set explicitly
	 *	below.
	 */
	if (setDefaultCommandLineArguments(arguments) != kCommonConstantReturnTypeSuccess)
	{
		return kCommonConstantReturnTypeError;
	}

	arguments->common.numberOfMonteCarloIterations  = 1;
	arguments->common.outputSelect                  = kOutputVariableIndexMax;
	arguments->common.isTimingEnabled               = false;
	arguments->common.isMonteCarloMode              = false;
	arguments->common.isOutputJSONMode              = false;
	arguments->common.isWriteToFileEnabled          = false;

	return kCommonConstantReturnTypeSuccess;
}
#endif

CommonConstantReturnType
getCommandLineArguments(int argc, char * argv[], CommandLineArguments * arguments)
{
#ifdef NO_OS_AVAILABLE
	/*
	 *	The no-OS build has no command line to parse, so ignore `argc` and
	 *	`argv` and use the hard-coded configuration instead.
	 */
	(void) argc;
	(void) argv;

	puts("Using hard coded command line arguments");

	return setNoOSCommandLineArguments(arguments);

#else

	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

	if (setDefaultCommandLineArguments(arguments) != kCommonConstantReturnTypeSuccess)
	{
		return kCommonConstantReturnTypeError;
	}

	DemoOption options[] = {
		{ 0 },
	};

	if (parseArgs(argc, argv, &arguments->common, options) != kCommonConstantReturnTypeSuccess)
	{
		fprintf(stderr, "Error: Parsing command line arguments failed.\n");
		printUsage();

		return kCommonConstantReturnTypeError;
	}

	if (arguments->common.isHelpEnabled)
	{
		printUsage();

		exit(EXIT_SUCCESS);
	}

	/*
	 *	If no output is selected, set outputSelect to kOutputVariableIndexMax,
	 *	which triggers computation of all outputs.
	 */
	if (!arguments->common.isOutputSelected)
	{
		arguments->common.outputSelect = kOutputVariableIndexMax;
	}

	/*
	 *	When outputSelect == kOutputVariableIndexMax (all outputs),
	 *	Monte Carlo mode is not supported.
	 */
	if (arguments->common.outputSelect == kOutputVariableIndexMax)
	{
		if (arguments->common.isMonteCarloMode)
		{
			fprintf(stderr, "Error: Please select a single output (-S) when in Monte Carlo mode.\n");

			return kCommonConstantReturnTypeError;
		}
	}
	else if (arguments->common.outputSelect > kOutputVariableIndexMax)
	{
		fprintf(stderr, "Error: Invalid output selection. Supported: -S 0 (cross section, scalar).\n");

		return kCommonConstantReturnTypeError;
	}

	if (arguments->common.isVerbose)
	{
		fprintf(stderr, "Warning: Verbose mode not supported. Continuing in non-verbose mode.\n");
	}

	return kCommonConstantReturnTypeSuccess;

#endif
}
