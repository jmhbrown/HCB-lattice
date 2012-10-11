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

// Function Declarations

void KineticHamiltonian(double *, int **, char);
void PotentialHamiltonian(double *, int **, char);
void StateInit(int **);
void StateCopy(int *, int *, int);
int StateId(int **, int *);

extern void StateEvoWithTime();
extern void MomDen(cdouble *, cdouble *);
extern void NoiseCorr(cdouble *, cdouble *, cdouble *);

void output(cdouble *, cdouble *, cdouble *);
void EnergChk(double *, double*, double *, cdouble *, cdouble *);
void ChkOutput(double, double, double, double);

void MPI_initialize(int, char **);
void MPI_finalize();

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




/*** Function Definitions ***/

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

// Initializes N*N matrix B by copying A to B
void StateCopy(int * A, int * B, int col, int N) {
	for(int i=0; i<N; i++) {
		B[col*N+i]=A[col*N+i];
	}
}

// assigns a state's ID and returns it. This ID becomes the state's column index in the hamiltonian
int StateId(int ** F, int * temp) {
	bool Flag;
	int ID;
	for(ID=0; ID<Nsite; ID++) {
		Flag=true;
		for(int i=0; i<Nsite; i++) {
			if(F[ID*Nsite+i]!=temp[i]) {
				Flag=false;
				break;
			}
		}
		if(Flag) break;
	}
	return ID;
}


// Generates the Kinetic Portion of the Hamiltonian
// 
void KineticHamiltonian(double * H, int * F, char Flag) {

	int * temp;
	temp=(int *)malloc(Nsite*sizeof(int));
	int i;
       	int j;
	double A;	// sign flag

	// fill H with zeros
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++) {
			H[a*Nsite+b]=0;
		}
	}

	// Note that the Kinetic Hamiltonian only has off-diagonal elements
	for(j=0; j<Nsite; j++) {
		for(int k=0; k<Nsite-1; k++) {
			StateCopy(F[
	
	}
		
}
