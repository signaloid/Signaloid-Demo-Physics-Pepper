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
#include <stdbool.h>
#include <inttypes.h>
#include "common.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 *	@brief	No CSV input distributions are used in this application.
 */
typedef enum
{
	kInputDistributionIndexMax = 0,
} InputDistributionIndex;

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

typedef struct
{
	CommonCommandLineArguments common;
} CommandLineArguments;

/**
 *	@brief	Check if a file exists at the given path.
 *
 *	@param	path	: Path to test.
 *	@return		: `true` if the file can be opened for reading, else `false`.
 */
bool
fileExists(const char * path);

/**
 *	@brief	Print out command line usage.
 */
void
printUsage(void);

/**
 *	@brief	Get command line arguments.
 *
 *	@param	argc		: argument count from `main()`.
 *	@param	argv		: argument vector from `main()`.
 *	@param	arguments	: Pointer to struct to store arguments.
 *	@return			: `kCommonConstantReturnTypeSuccess` if successful,
 *				  else `kCommonConstantReturnTypeError`.
 */
CommonConstantReturnType
getCommandLineArguments(int argc, char * argv[], CommandLineArguments * arguments);

#ifdef __cplusplus
}
#endif
