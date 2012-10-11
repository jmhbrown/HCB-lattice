/* 
Hamiltonian.cpp

Purpose:	Generates the two parts of our hamiltonian, the non-diagonal 
		time-independent part and the time-dependent diagonal part. 
 
Author: 	Jennifer Brown
 		jmb347@georgetown.edu


NOTES: I'm changing my matrices from double arrays to single arrays. 



*/


#include "stdafxDyn.h" // project-specific header file. 

// not sure if these are needed, they're included in Kai's NoiseCorrelationTrapDyn.cpp
// so I'm including them here as well.
#include <fstream>
#include <sstream>
#include <cstring>
#include <mkl_lapack.h>
#include <mkl_blas.h>

using namespace std;

/*** Variable and Function Declarations ***/

// Lattice constants
int Nsite;	// number of lattice sites
int Nfermion;	// number of particles
int Nbzone;	// number of Brillouin zones
int K=-1;	// specific momentum basis L

// Energy Coefficients
double T;	// hopping
extern double V(int)	// time-dependent lattice potential as of oct. 4, this is not yet defined


// Time-Stepping variables
int xtime;	// intermediate time, used for indexing
int itime;	// initial time. I assume this should be zero.
int ftime;	// final time
int tstep;	// step size

/*** Functions ***/

// State Initi - constructs an idenity matrix of size Nsite * Nsite
void StateInit(int * F) {
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++) {
			F[a*Nsite+b]=0;
		}

		for(int i=0; i<Nsite; i++) {
			F[i*Nsite+i]=1;
		}
	}
}

void KineticHamiltonian(double * H, int * F, char Flag) {

	int * temp;
	temp=(int *)malloc(Nsite*sizeof(int));
	int i, int j;
	double A;	// sign flag

	// fill H with zeros
  // NOT DONE!!! WORK HERE. 	
}
