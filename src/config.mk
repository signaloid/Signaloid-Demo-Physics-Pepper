SOURCES =\
	pepper.cpp\
	pepper-utilities.cpp\
	pepper-kernel.cpp\
	pepper-uxhw.cpp\
	pepper-monteCarlo.cpp\
	pepper-lhapdfLite.c\
	common.c

CFLAGS = -std=c11 -Wno-missing-braces -D_GNU_SOURCE

CXXFLAGS = -std=c++11 -D_GNU_SOURCE
