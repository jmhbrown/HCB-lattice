/* 
Hamiltonian.cpp

Purpose:	Generates the two parts of our hamiltonian, the non-diagonal 
		time-independent part and the time-dependent diagonal part. 
 
Author: 	Jennifer Brown
 		jmb347@georgetown.edu


NOTES: Kai used a combination of double and single arrays to describe matrices.
	I've decided to switch all of these over to double matrices. 


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


// Generates the Kinetic Portion of the Hamiltonian
// 
void KineticHamiltonian(double ** H, int ** F, char Flag) {

	int * temp;
	int i;
       	int j;
	double A;	// sign flag

	temp=(int *)malloc(Nsite*sizeof(int));


	// fill H with zeros
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++) {
			H[a][b]=0;
		}
	}

	// Note that the Kinetic Hamiltonian only has off-diagonal elements
	for(j=0; j<Nsite; j++) {
		StateCopy(F, temp, j); // copies jth row of F to temp	
		for(int k=0; k<Nsite-1; k++) {
		if((F[j][k]+F[j][k+1])==1) {
			temp[k+1]=F[j][k];
			temp[k]=F[j][k+1];
		}
		
}
