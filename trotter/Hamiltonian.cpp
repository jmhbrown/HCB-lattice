/* 
Hamiltonian.cpp

Purpose:	Generates the two parts of our hamiltonian, the non-diagonal 
		time-independent part and the time-dependent diagonal part. 
 
Author: 	Jennifer Brown
 		jmb347@georgetown.edu


NOTES: Kai used a combination of double and single arrays to describe matrices.
	I've decided to switch all of these over to 1D arrays, ecausethat's what Lapack works with. 
	I've also decided to switch all functions over to BLAS and LAPACK's built in ones.

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

void KineticHamiltonian(double *, int *, char);
void PotentialHamiltonian(double *, int *, char);

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
double V;	// lattice potential


// Time-Stepping variables
int xtime;	// intermediate time, used for indexing
int itime;	// initial time. I assume this should be zero.
int ftime;	// final time
int tstep;	// step size




/*** Function Definitions ***/


// Generates the Kinetic Portion of the Hamiltonian
// 
void KineticHamiltonian(double * HK, int * F, char Flag) {

	int * temp;
	int i;
	int j;
	double A;	// sign flag

	temp=(int *)malloc(Nsite*sizeof(int));


	// fill HK with zeros
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++) {
			HK[a*Nsite+b]=0;
		}
	}

	// Note that the Kinetic Hamiltonian only has off-diagonal elements
	for(j=0; j<Nsite; j++) {
		StateCopy(F, temp, j, Nsite); // copies jth row of F to temp	
		for(int k=0; k<Nsite-1; k++) {
			if((F[j*Nsite+k]+F[j*Nsite+k+1])==1) {
				temp[k+1]=F[j*Nsite+k];
				temp[k]=F[j*Nsite+k+1];
			}
		}

		// for periodic lattices
		if(Flag=='P') {
			A = -1.0;
			if(Nfermion%2==0) { // anti-periodic case for an even number of particles.
				A = 1.0;
			}
			StateCopy(F, temp, j, Nsite);
			if((F[j*Nsite]+F[j*Nsite+Nsite-1])==1) {
				temp[0]=F[j*Nsite+Nsite-1];
				temp[Nsite-1]=F[j*Nsite];
				i=StateId(F, temp, Nsite);
				HK[i*Nsite+j]=A*T; 
			}
		}           
	}
	free(temp);
	
} // kinetic hamiltonian


// Potential Hamiltonian
void PotentialHamiltonian (double * HP, int * F, char Flag) {
	/*
	Kai's code has a U and a V, for the lattice and hopping potential,
	respectively. The lattice potential doesn't show up in my notes 
	from Professor Rigol's lecture. Here's how U shows up in Kai's 
	potential hamiltonian:
	if(U!=0.0) {
		for(i=0; i<Nsite; i++) {
			for(int k=0; k<Nsite; k++) {
            		H[i*Nsite+i]+=2*cos(2.0*Pi*A*(k+1))*F[i][k]*U+pow(double(k)-xc, 2.0)*F[i][k]*V;
			}
		}
	}
	else { // see uncommented code below:
*/

	for(int i=0; i<Nsite; i++) {
		for(int j = 0; j<Nsite; j++) {
			HP[i*Nsite+j]+=pow((double)j,2.0)*F[i*Nsite+j]*V
		}
	}

} // pontential hamiltonian
