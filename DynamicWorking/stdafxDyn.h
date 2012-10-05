// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//
#ifdef MYMPI
#define DOUBLE_COMPLEX COMPLEX16    //on Georgetown machines
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


// variable declaration
extern const double Pi=2.0*acos(0.0);

extern int Nsite;	// number of sites
extern int Nfermion;	// number of particles
extern int Nbzone;	// number of Brillouin zones
extern int K;		//specific momemtum basis K

// 1D single-fermion Hamiltonian (size: Nsite x Nsite) used to store eigenstates
extern double * Ham0;	// initial hamiltonian, with trap centered
extern double * HamT; 	// after the onset of trap displacement

extern double * EnergT;	//Energy level/eigenvalues of HamT

extern cdouble *cPT; 	//time-dependent fermionic slater determinant
extern char BorS;	// marks whether our system is a Boson or Spin
extern int xtime;	// used for time-stepping

// for MPI
extern int irank, isize; 

#endif


// TODO: reference additional headers your program requires here
