/***********************************************************************************************************************************
 *  This program is designed to calculate the noise correlation for multi-HCB and multi-spin1/2 system. The algorithm is based     *
 *  on the slater determinant method for fermions, which constructs the ground state for many-body system with eigenstates for     *
 *  single-particle fermion system. Spin1/2 system can be mapped into the fermion system and HCB system can be mapped into spin1/2 *
 *  system. The Green's functions and four-point correlations can be obtained and finally the noise correlations are calculated.   *
 *  HARMONIC TRAP in presence. Non-equilibrium.                                                                                    *
 ***********************************************************************************************************************************/

#include "stdafxDyn.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>
#include <fstream>
#include <sstream>
#include <cstring>


using namespace std;

/**********************************Declarations of functions and global variables****************************************************/

void MPI_initialize(int, char **);
void StateInit(int **);
void HamGen(double *, int **, double, char);
void StateCopy(int *, int *, int);
int StateId(int **, int *);
extern void StateEvoWithTime();
extern void MomDen(cdouble *, cdouble *, double *, double *, int *);
extern void NoiseCorr(cdouble *, cdouble *, cdouble *, double *, int *);
void output(cdouble *, cdouble *, cdouble *, double, double, double, int, int);
void MPI_finalize();
void EnergChk(double *, double*, double *, cdouble *, cdouble *);
void ChkOutput(double, double, double, double);

/************************General variables*************************/
const double Pi=2.0*acos(0.0);

int Nsite, Nfermion, Nbzone, K=-1; //# of sites and # of particles and # of Brillouin zones and specific momemtum basis K
double T, U, V; //energy coeffecients: hopping, lattice potential, trap potential
char BorS; //indicator of system type: HCB or spin
double * Ham0, * HamT; //1-d single-fermion Hamiltonian matrix (Nsite*Nsite) and later used to store eigenstates;
					 //Ham0: initial(with trap centered); HamT: after onset of trap displacement
double * EnergT;	 //Energy level/eigenvalues of HamT
cdouble * cPT; //time-dependent fermionic slater determinant
int irank, isize; //<-----for MPI use
double ratioA, ratioB;

#ifdef MYMPI
MPI::Status istatus;
#endif

/***********************Dynamic variables**************************/
double xdisp; //displacement of the center of the trap
int xtime, itime, ftime, tstep; //time-variables

/******************************************************The main function**************************************************************/

int main(int argc, char ** argv)
{
    int ** F;  //single-fermion states
    double * Energ0, GSEnerg=0.0; //Energy levels for single-fermion system, ground state energy for multi-fermion system
    cdouble * Noise, * rho, * Nk; //Noise correlation (k, k'), one-particle density matrix(Nsite*Nsite), momentum distribution function
    double GSEnerg2, KEnerg1, KEnerg2;
    char  PeriodicFlag;
    int Ns, lda, lwork, info0, infoT; //<-----for Lapack use
	double * work;			  //<----/
	double xCM, mMax, nMax;
	int kmMax, knMax;

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
    for(int i=0; i<Nsite; i++)
        *(F+i)=(int *)malloc(Nsite*sizeof(int));
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
		
		MomDen(rho, Nk, &xCM, &mMax, &kmMax);

		NoiseCorr(Nk, rho, Noise, &nMax, &knMax);

		if(irank==0)
        {
            EnergChk(&KEnerg1, &KEnerg2, &GSEnerg2, rho, Nk);
            output(rho, Nk, Noise, xCM, mMax, nMax, kmMax, knMax);
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

/****************************************************The end of the main function*****************************************************/
/***************************************************Definition of functions used******************************************************/

void MPI_initialize(int argc, char ** argv)
{
    #ifdef MYMPI
        MPI::Init(argc, argv);
        irank=MPI::COMM_WORLD.Get_rank();
        isize=MPI::COMM_WORLD.Get_size();
    #else
        irank=0;
        isize=1;
    #endif
}

void StateInit(int ** F)
{

    for(int a=0; a<Nsite; a++)
        for(int b=0; b<Nsite; b++)
            F[a][b]=0;

        for (int i=0; i<Nsite; i++)
                F[i][i]=1;          
}

void HamGen(double * H, int ** F, double disp, char Flag)
{
    int * temp;
    temp=(int *)malloc(Nsite*sizeof(int));
    int i, j;
    double A;       //sign flag
    double xc;      //center of the trap

    //fill H with 0s
    for(int a=0; a<Nsite; a++)
        for(int b=0; b<Nsite; b++)
            H[a*Nsite+b]=0;

    //Offdiagonal elements due to hopping
    for(j=0; j<Nsite; j++)
    {
        for(int k=0; k<Nsite-1; k++)
        {
            StateCopy(F[j], temp, Nsite);
            if((F[j][k]+F[j][k+1])==1)
            {
                temp[k+1]=F[j][k];
                temp[k]=F[j][k+1];
                i=StateId(F, temp);
                H[i*Nsite+j]=-T;
            }           
        }

        //Periodic condition applies, which allows the hopping between Site 1 and last site 

        if(Flag=='P') 
        {
            A=-1.0;
            if(Nfermion%2==0)   //antiperiodic case for even number of particles
                A=1.0;

            StateCopy(F[j], temp, Nsite);
            if((F[j][0]+F[j][Nsite-1])==1)
            {
                temp[0]=F[j][Nsite-1];
                temp[Nsite-1]=F[j][0];
                i=StateId(F, temp);
                H[i*Nsite+j]=A*T; 
            }
        }
    }

    //Diangonal elements due to on-site potential
    A=ratioA/ratioB;
    xc=double(Nsite-1)/2.0+disp; //a displaced center of the trap
    
    if(U!=0.0)
    for(i=0; i<Nsite; i++)
        for(int k=0; k<Nsite; k++)
            H[i*Nsite+i]+=2*cos(2.0*Pi*A*(k+1))*F[i][k]*U+pow(double(k)-xc, 2.0)*F[i][k]*V;
    else
    for(i=0; i<Nsite; i++)
        for(int k=0; k<Nsite; k++)
            H[i*Nsite+i]+=pow(double(k)-xc, 2.0)*F[i][k]*V;        

    free(temp);
}


void StateCopy(int * A, int * B, int N) //Initialize B by copying A to B
{
    for(int i=0; i<N; i++)
        B[i]=A[i];
}

int StateId(int ** F, int * temp)
{
    bool Flag;
    int ID;
    for(ID=0; ID<Nsite; ID++)
    {
        Flag=true;
        for(int i=0; i<Nsite;i++)
        {
            if(F[ID][i]!=temp[i])
            {
                Flag=false;
                break;
            }
        }
        if(Flag) break;
    }
    return ID;
}

void EnergChk(double * KEnerg1, double * KEnerg2, double * GSEnerg2, cdouble * rho, cdouble * Nk)
{
	cdouble kin1=0.0, pot1=0.0, tot;
	cdouble kin2=0.0;
    double A=ratioA/ratioB;
    double xc=double(Nsite-1)/2.0+xdisp;

	//following E_k=-t*Sum{rho(i, i+1)+rho(i, i-1)}
	//			E_p= U*Sum{rho(i, i)}
	//			E  = E_k+E_p
	for(int i=0, j1, j2; i<Nsite; i++)
	{
		j1=i+1; j2=i-1;
		if(j1>=Nsite)
			j1-=Nsite;
		if(j2<0)
			j2+=Nsite;
		
		kin1+=-T*(rho[i*Nsite+j1]+rho[i*Nsite+j2]);
        //
        //modified the potential form according to the model    
        // 
		pot1+=U*2.0*cos(2.0*Pi*i*ratioA/ratioB)*rho[i*Nsite+i]+V*pow(double(i)-xc, 2.0)*rho[i*Nsite+i];
	}

	//
    //following E_k=2.0*Sum{cos(ka)*n(k)}
	//
    for(int k=-Nbzone*Nsite, i=0; k<(Nsite+1)/2&&k<Nbzone*Nsite; k++, i++)        //count N[k] within one brilluion zone, [-Pi/2, Pi/2)
        if(k>=-Nsite/2)
			kin2+=-T*2.0*cos(2.0*Pi*k/Nsite)*Nk[i];

	*KEnerg1=real(kin1);
	*KEnerg2=real(kin2);
	*GSEnerg2=real(kin1+pot1);

}


void output(cdouble * rho, cdouble * Nk, cdouble * Noise, double xCM, double mMax, double nMax, int kmMax, int knMax)
{
    double ka, ki, kj; //ka is the momentum in terms of Pi, ka=2*k/Nsite
    ofstream op1, op2, op3, op4, op5, op6;
    ostringstream densfile, momfile, noisefile, smryfile;
    
    densfile<<"Output/DensTrap_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp<<"time"<<xtime;
    momfile<<"Output/MomTrap_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp<<"time"<<xtime;
    noisefile<<"Output/NoiseTrap_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp<<"time"<<xtime;
    smryfile<<"Output/Summary_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp;

    string s0=densfile.str();
    string s1=momfile.str();
    string s2=noisefile.str();
	string s3=smryfile.str();
    
    char dfile[200], mfile[200], nfile[200], sfile[200];
    memset(dfile, '\0', 200); s0.copy(dfile, 100);
    memset(mfile, '\0', 200); s1.copy(mfile, 100);
    memset(nfile, '\0', 200); s2.copy(nfile, 100);
	memset(sfile, '\0', 200); s3.copy(sfile, 100);

//  op1.open("Dens_Mx"); //'os::app' means start from last end
    op2.open(dfile);
    op3.open(mfile);
//  op4.open("Dens_Corr");
    op5.open(nfile);
	op6.open(sfile, ios::app);

/*    for(int i=0; 2*i<(Nsite+1); i++)
        op4<<i<<"\t"<<rho[i]<<endl;

    for(int i=0; i<Nsite; i++)
    {
        for(int j=0; j<Nsite; j++)
            op1<<rho[i*Nsite+j]<<"\t";
        op1<<endl;
    }
*/
    for(int i=0; i<Nsite; i++)
        op2<<xtime<<"\t"<<i<<"\t"<<real(rho[i*Nsite+i])<<endl;

    for(int k=-Nbzone*Nsite, i=0; k<=Nsite/2&&k<=Nbzone*Nsite; k++, i++)        //output N[k] within one brilluion zone, [-Pi/2, Pi/2]
        if(k>=-Nsite/2)
        {
            ka=2*double(k)/double(Nsite);
            op3<<xtime<<"\t"<<ka<<"\t"<<real(Nk[i])<<endl;
        }
    
    //output the noise correlation within pre-set number of brilluion zones
    if(K==-1)       //full matrix of (k, k')
        for(int i=0, k=-Nbzone*Nsite; k<=Nbzone*Nsite; k++, i++)
            for(int j=0, q=-Nbzone*Nsite; q<=Nbzone*Nsite; q++, j++)
            {
                ki=2*double(k)/double(Nsite);
                kj=2*double(q)/double(Nsite);
                op5<<ki<<"\t"<<kj<<"\t"<<real(Noise[i*(2*Nbzone*Nsite+1)+j])<<endl;
            }
    else    //fixed k (or k')
    {
        int i=K+Nbzone*Nsite;
        for(int k=-Nbzone*Nsite, j=0; k<=Nbzone*Nsite; k++, j++)
        {
            ka=2*double(k)/double(Nsite);
            op5<<xtime<<"\t"<<ka<<"\t"<<real(Noise[i*(2*Nbzone*Nsite+1)+j])<<endl;
        }
    }

	ki=2*double(kmMax)/double(Nsite);
	kj=2*double(knMax)/double(Nsite);
	op6<<xtime<<"\t"<<xCM<<"\t"<<ki<<"\t"<<mMax<<"\t"<<kj<<"\t"<<nMax<<endl;

    op1.close();
    op2.close();
    op3.close();
    op4.close();
    op5.close();
	op6.close();
}

void ChkOutput(double ke1, double ke2, double e1, double e2)
{
	ofstream chk;
	ostringstream cfile;

	cfile<<"CheckFile/ChkTrap_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp;
	string s=cfile.str();
	char name[100];
	memset(name, '\0', 100); s.copy(name, 50);

	chk.open(name, ios::app);
//	chk<<"Kinetic1\tKinetic2\tTotal1\t\tTotal2\n";
	chk<<xtime<<"\t"<<ke1<<"\t"<<ke2<<"\t"<<e1<<"\t"<<e2<<endl;

	chk.close();
}

void MPI_finalize()
{
    #ifdef MYMPI
        MPI::Finalize();
    #endif
}

/***********************************************The end of function definitions*******************************************************/
