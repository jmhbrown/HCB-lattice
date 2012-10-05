// Constructs the Momentum density matrix
// Originally code was written by Kai, then put here by Jenny (jmb347@georgetown.edu)
// 
//
//
// These functions were originally in CorrelationsDyn.cpp
//  --J Oct 2012
// CHANGES: replaced corr_2 -> rho, because these refer to the same matrix physically and I was getting confused.
//
//

#include "stdafxDyn.h"
#include <fstream>

using namespace std;

extern void TwoPntGreenMx(cdouble *);
cdouble MomDistFunc(int, cdouble *);
void DenMx(cdouble *);


// Density Matrix -> See eq. (6.13)
void DenMx(cdouble * rho)		//calculate and store <b+b-> in rho[]
{
	//int flag=0;
	cdouble * green_2;
	green_2=(cdouble *)malloc(Nsite*Nsite*sizeof(cdouble));

	TwoPntGreenMx(green_2);

  // (6.13)
	for(int i=0; i<Nsite; i++)
		for(int j=i; j<Nsite; j++)
		{
			if(i==j)
				rho[i*Nsite+j]=1.0-green_2[i*Nsite+j];
			else
            {
                //notice the density matrix is not symmetric anymore
                rho[i*Nsite+j]=green_2[j*Nsite+i];
                rho[j*Nsite+i]=green_2[i*Nsite+j];
            }
		}

	free(green_2);

}

// Momentum Distribution Function, <nk> for a single K.
cdouble MomDistFunc(int k, cdouble * rho)	
{
	cdouble nk=0.0; 
	double theta;

// this is the fourier transform of the density matrix
	for(int i=0; i<Nsite; i++)
		for(int j=0; j<Nsite; j++)
		{
			theta=2*Pi*double(k*(i-j))/double(Nsite);
			//nk+=rho[i*Nsite+j]*cos(theta);    // the real part of nk
			nk+=rho[i*Nsite+j]*polar(1.0, theta);
		}

	nk=nk/double(Nsite); // normalizing the density matrix

	return nk;
}


// Constructs the Momentum Density Matrix rho=<b+b->, Nsite*Nsite; Nk=<nk>, 2*Nbzone*Nsite+1
void MomDen(cdouble * rho, cdouble * Nk){
	DenMx(rho);

	if(irank==0)			//included in stdafxDyn in the line: extern int irank, isize;
		for(int i=0, k=-Nbzone*Nsite; k<=Nbzone*Nsite; k++, i++)
			Nk[i]=MomDistFunc(k, rho);
}

