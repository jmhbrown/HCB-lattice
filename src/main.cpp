#include "stdafxDyn.h"
#include "MyMPI.h"
#include "cdouble.h"
#include "NoiseCorrelationTrapDyn.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>
#include <fstream>
#include <sstream>
#include <cstring>


using namespace std;

int main(int argc, char ** argv)
{
	int ** F;  //single-fermion states
	double * Energ0, GSEnerg=0.0; //Energy levels for single-fermion system, ground state energy for multi-fermion system
	cdouble * Noise, * rho, * Nk; //Noise correlation (k, k'), one-particle density matrix(Nsite*Nsite), momentum distribution function
	double GSEnerg2, KEnerg1, KEnerg2;
	char  PeriodicFlag;
	int Ns, lda, lwork, info0, infoT; //<-----for Lapack use
	double * work;			  //<----/

	MPI_initialize(argc, argv);


	if(irank==0)    //data input to Master
	{
        //input/output setting
        ifstream ip("Parameters");

        ip>>Nfermion>>Nsite>>T>>V>>PeriodicFlag>>Nbzone>>BorS>>K;
		ip>>xdisp>>itime>>ftime>>tstep;
        ip>>ratioA>>ratioB>>U;  // these have to do with the super lattice
        ip.close();
	}


    //Memory allocation
	F=(int **)malloc(Nsite*sizeof(int *));
	for(int i=0; i<Nsite; i++) {
		*(F+i)=(int *)malloc(Nsite*sizeof(int));
	}
	Ham0=(double *)malloc(Nsite*Nsite*sizeof(double));
	HamT=(double *)malloc(Nsite*Nsite*sizeof(double));
	rho=(cdouble *)malloc(Nsite*Nsite*sizeof(cdouble));
	Nk=(cdouble *)malloc((2*Nbzone*Nsite+1)*sizeof(cdouble));
	Noise=(cdouble *)malloc((2*Nbzone*Nsite+1)*(2*Nbzone*Nsite+1)*sizeof(cdouble));
	Energ0=(double *)malloc(Nsite*sizeof(double));
	EnergT=(double *)malloc(Nsite*sizeof(double));
	work=(double *)malloc(3*Nsite*sizeof(double));

	
	StateInit(F);

	HamGen(Ham0, F, 0.0, PeriodicFlag);
	HamGen(HamT, F, xdisp, PeriodicFlag);

	lda=Ns=Nsite;
	lwork=3*Ns;

    // calling Lapack's function dsyev() to solve the single-particle Hamiltonians (init, and final)
	dsyev("V", "U", &Ns, Ham0, &lda, Energ0, work, &lwork, &info0);
	dsyev("V", "U", &Ns, HamT, &lda, EnergT, work, &lwork, &infoT);

	if(info0!=0||infoT!=0) //failed
	{
		cout<<"dsyev failed"<<endl;
		return 0;
	}


    //Otherwise succeeded
	for(int i=0; i<Nfermion; i++) {	
		GSEnerg+=Energ0[i]; //Record the ground state energy for the state of interest (multi-fermion)
	}
	
	//time-evolution in the correlations
	for(xtime=itime; xtime<=ftime; xtime+=tstep) {
		StateEvoWithTime();

		MomDen(rho, Nk);

		NoiseCorr(Nk, rho, Noise);

		if(irank==0) {
			EnergChk(&KEnerg1, &KEnerg2, &GSEnerg2, rho, Nk); // Compares KEnerg1 and KEnerG2. 
			output(rho, Nk, Noise);
			ChkOutput(KEnerg1, KEnerg2, GSEnerg, GSEnerg2);
		}
	}

	MPI_finalize();

//free the memory
	for(int i=0; i<Nsite; i++) {
		free(F[i]);
	}
	free(F);
	free(Energ0);
	free(work);
	free(Ham0);
	free(HamT);
	free(EnergT);
	free(cPT);
	free(rho);
	free(Nk);
	free(Noise);

	return 0;

}
