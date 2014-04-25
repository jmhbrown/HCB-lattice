#include "stdafxDyn.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>
#include <fstream>
#include <sstream>
#include <cstring>


using namespace std;

double T, U, V; //energy coeffecients: hopping, lattice potential, trap potential
char BorS; //indicator of system type: HCB or spin
double * Ham0, * HamT; //1-d single-fermion Hamiltonian matrix (Nsite*Nsite) and later used to store eigenstates;
					 //Ham0: initial(with trap centered); HamT: after onset of trap displacement
double * EnergT;	 //Energy level/eigenvalues of HamT
cdouble * cPT; //time-dependent fermionic slater determinant
int irank, isize; //<-----for MPI use
double ratioA, ratioB;

double xdisp; //displacement of the center of the trap
int xtime, itime, ftime, tstep; //time-variables



int Nsite, Nfermion, Nbzone, K=-1; //# of sites and # of particles and # of Brillouin zones and specific momemtum basis K


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
        ip>>ratioA>>ratioB>>U;
        ip.close();
        
    }

    
    #ifdef MYMPI    //input data sharing through MPI
        
        MPI::COMM_WORLD.Bcast(&Nfermion, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&Nsite, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&T, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&V, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&U, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&ratioA, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&ratioB, 1, MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(&PeriodicFlag, 1, MPI::CHAR, 0);
        MPI::COMM_WORLD.Bcast(&Nbzone, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&BorS, 1, MPI::CHAR, 0);
        MPI::COMM_WORLD.Bcast(&K, 1, MPI::INT, 0);

		MPI::COMM_WORLD.Bcast(&xdisp, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&itime, 1, MPI::INT, 0);
		MPI::COMM_WORLD.Bcast(&ftime, 1, MPI::INT, 0);
		MPI::COMM_WORLD.Bcast(&tstep, 1, MPI::INT, 0);

    #endif

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
    for(int i=0; i<Nfermion; i++)
        GSEnerg+=Energ0[i]; //Record the ground state energy for the state of interest (multi-fermion)

	//time-evolution in the correlations
	for(xtime=itime; xtime<=ftime; xtime+=tstep)
	{
		StateEvoWithTime();

		MomDen(rho, Nk);

		NoiseCorr(Nk, rho, Noise);

		if(irank==0)
        {
            EnergChk(&KEnerg1, &KEnerg2, &GSEnerg2, rho, Nk);
            output(rho, Nk, Noise);
            ChkOutput(KEnerg1, KEnerg2, GSEnerg, GSEnerg2);
        }
	}

    MPI_finalize();

    //free the memory
    for(int i=0; i<Nsite; i++)
        free(F[i]);
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
