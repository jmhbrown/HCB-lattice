// Does the calculation P_t = U*e^(-iDt/h)*(U+)*P_0,
// which is the time-dependent fermionic slater determinant
// This was originally found in GreensDyn.cpp

#include "stdafxDyn.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>

void StateEvoWithTime() //calculate the time-dependent fermionic slater determinant
{

	double * P0;  // initial state
	double * Ptemp; // temporary state, used in calculations
	cdouble * cHamT; // complex hamiltonian
	cdouble *cPtemp; // complex version of temporary state

	P0=(double *)malloc(Nfermion*Nsite*sizeof(double));
	Ptemp=(double *)malloc(Nfermion*Nsite*sizeof(double));
	cHamT=(cdouble *)malloc(Nsite*Nsite*sizeof(cdouble));
	cPtemp=(cdouble *)malloc(Nfermion*Nsite*sizeof(cdouble));
	cPT=(cdouble *)malloc(Nfermion*Nsite*sizeof(cdouble));

	//P_0 the initial fermionic slater determinant
	for(int i=0; i<Nfermion*Nsite; i++)
		P0[i]=Ham0[i];

	//Ptemp=(U+)*P_0
	double alpha=1.0, beta=0.0;
	dgemm("C", "N", &Nsite, &Nfermion, &Nsite, &alpha, HamT, &Nsite, P0, &Nsite, &beta, Ptemp, &Nsite);    // This HamT has be diagonalized, it's constructed elsewhere.  -- J 10/2.

	//cPtemp=1e^(-iDt)*Ptemp
	for(int i=0; i<Nfermion; i++)
		for(int j=0; j<Nsite; j++)
			cPtemp[i*Nsite+j]=polar(1.0, -double(xtime)*EnergT[j])*Ptemp[i*Nsite+j];

	//P_t=U*cPtemp
	for(int i=0; i<Nsite*Nsite; i++)
		cHamT[i]=HamT[i];
	cdouble calpha=1.0, cbeta=0.0;
	zgemm("N", "N", &Nsite, &Nfermion, &Nsite, &calpha, cHamT, &Nsite, cPtemp, &Nsite, &cbeta, cPT, &Nsite);

	free(P0);
	free(Ptemp);
	free(cHamT);
	free(cPtemp);
}// determining the state P based on the diagonalized hamiltonian. 
