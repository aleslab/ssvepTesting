#ifndef MATLABAUX_H
#define MATLABAUX_H


#include "mex.h"
#include "matrix.h"
#include <memory.h>

// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the MATLABAUX_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// MATLABAUX_API functions as being imported from a DLL, wheras this DLL sees symbols
// defined with this macro as being exported.
#ifndef COMPILE_IN_MATLAB
#ifdef MATLABAUX_EXPORTS
#define MATLABAUX_API __declspec(dllexport)
#else
#define MATLABAUX_API __declspec(dllimport)
#endif
#else
#define	MATLABAUX_API	
#endif



MATLABAUX_API unsigned long maGetUINT32Element(const mxArray* pArray,
unsigned long index);
MATLABAUX_API double maGetDoubleElement(const mxArray* pArray,
unsigned long index);
MATLABAUX_API long maGetINT32Element(const mxArray* pArray, long index);
MATLABAUX_API mxArray* maCreateDoubleElement(double value);
MATLABAUX_API mxArray* maCreateUINT32Element(unsigned long value);
MATLABAUX_API mxArray* maCreateStringElement(const char* temp);
//mxArray *mxCreateCharMatrixFromStrings(int m, const char **str);
MATLABAUX_API mxArray* maCreateBufferElement(char* buff, unsigned long size);


#endif
