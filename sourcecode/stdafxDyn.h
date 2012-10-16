// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//
#ifdef MYMPI
#define DOUBLE_COMPLEX COMPLEX16    //on Georgetown mahcines
#include <mpi.h>
extern MPI::Status istatus;

#endif

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <complex>

typedef std::complex<double> cdouble;
#define MKL_Complex16 cdouble

#ifndef _HEADER_H
#define _HEADER_H

extern const double Pi;

extern int Nsite, Nfermion, Nbzone, K;
extern double * Ham0, * HamT, * EnergT;
extern cdouble *cPT;
extern char BorS;
extern int xtime;
extern int irank, isize;

#endif


// TODO: reference additional headers your program requires here
