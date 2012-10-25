/***********************************************************************************************************************************
 *  This program is designed to calculate the noise correlation for multi-HCB and multi-spin1/2 system. The algorithm is based     *
 *  on the slater determinant method for fermions, which constructs the ground state for many-body system with eigenstates for     *
 *  single-particle fermion system. Spin1/2 system can be mapped into the fermion system and HCB system can be mapped into spin1/2 *
 *  system. The Green's functions and four-point correlations can be obtained and finally the noise correlations are calculated.   *
 *  HARMONIC TRAP in presence. Non-equilibrium.                                                                                    * ***********************************************************************************************************************************/

#include "stdafxDyn.h"
#include "cdouble.h"
#include "MyMPI.h"
#include <mkl_lapack.h>
#include <mkl_blas.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <math.h>
#include <stdio.h>

using namespace std;

/************************General variables*************************/
//const double Pi=2.0*acos(0.0);

//Lattice coefficients:
int Nsite; //number of sites
int Nfermion; //number of particles
int Nbzone; //number of Brillouin zones
int K = -1; //specific momemtum basis K
//Energy coefficients:
//Note to self: rename these better, maybe
double T; //hopping potential
double U; //lattice potential
double V; //trap potential
char BorS; //indicates system type: HCB or spin (B or S)
/*1-d single-fermion Hamiltonian matrix (Nsite*Nsite) and later used to store eigenstates*/
double * Ham0; //initial(with trap centered)
double * HamT; //after onset of trap displacement
double * EnergT; //Energy level/eigenvalues of HamT
cdouble * cPT; //time-dependent fermionic slater determinant
double ratioA, ratioB;
//For MPI use:
int irank, isize;

/***********************Dynamic variables**************************/
double xdisp; //displacement of the center of the trap
//Time variables:
int xtime; //
int itime; //
int ftime; //
int tstep; //

//Functions n' stuff
void MPI_initialize(int, char **);
void StateInit(int **);
void HamGen(double *, int **, double, char);
void StateCopy(int *, int *, int);
int StateId(int **, int *);
extern void StateEvoWithTime();
extern void MomDen(cdouble *, cdouble *);
extern void NoiseCorr(cdouble *, cdouble *, cdouble *);
void output(cdouble *, cdouble *, cdouble *);
void MPI_finalize();
void EnergChk(double *, double*, double *, cdouble *, cdouble *);
void ChkOutput(double, double, double, double);


/***Definition of functions used***/

//We think he's trying to initialize an identity matrix of size Nsite*Nsite
void StateInit(int ** F)
{ //you'll need to learn "new" to make dynamically allocated arrays
//you'll have to do it twice to make a 2D one
//(the ** may be for something else, but in general is NOT necessary for 2d dynamic arrays)
//M0AR NOTES: vectors are basically better, standardized dynamic arrays. They live in <vector>. man it!
/* //This was how he did it:
    for(int a=0; a<Nsite; a++)
        for(int b=0; b<Nsite; b++)
            F[a][b]=0;

        for (int i=0; i<Nsite; i++)
                F[i][i]=1;*/
//our [incomplete] shot at it:
//this will likely require us making the dynamic array though
	for (int a=0; a<Nsite; a++) {
		for (int b=0; b<Nsite; b++) {
			if (b==a)
				F[a][b]=1;
			else
				F[a][b]=0;
		}
	}
}//end StateInit

//Probably generating the Hamiltonian
void HamGen(double * H, int ** F, double disp, char Flag) {
	int * temp;
	temp=(int *)malloc(Nsite*sizeof(int));
	int i, j;
	double A; //sign flag
	double xc; //center of the trap

	//fill H with 0s
	for(int a=0; a<Nsite; a++) {
		for(int b=0; b<Nsite; b++)
			H[a*Nsite+b]=0;
	}

	//Off-diagonal elements due to hopping
	for (j=0; j<Nsite; j++) {
		for (int k=0; k<Nsite-1; k++) {
			StateCopy(F[j], temp, Nsite);
			if ((F[j][k]+F[j][k+1])==1) {
				temp[k+1]=F[j][k];
				temp[k]=F[j][k+1];
				i=StateId(F, temp);
				H[i*Nsite+j]=-T;
			}           
		}

	//Periodic condition applies, which allows the hopping between Site 1 and last site 
		if(Flag=='P') {
			A=-1.0;
			if(Nfermion%2==0) //antiperiodic case for even number of particles
				A=1.0;
			StateCopy(F[j], temp, Nsite);
			if((F[j][0]+F[j][Nsite-1])==1) {
				temp[0]=F[j][Nsite-1];
				temp[Nsite-1]=F[j][0];
				i=StateId(F, temp);
				H[i*Nsite+j]=A*T; 
			}//end the second if
		}//end the outer if
	}//end for from 25 lines ago

	//Diangonal elements due to on-site potential
	A=ratioA/ratioB;
	xc=double(Nsite-1)/2.0+disp; //a displaced center of the trap
    
	if(U!=0.0) {
		for(i=0; i<Nsite; i++) {
			for(int k=0; k<Nsite; k++)
				H[i*Nsite+i]+=2*cos(2.0*Pi*A*(k+1))*F[i][k]*U+pow(double(k)-xc, 2.0)*F[i][k]*V;
		}
	}
	else {
		for(i=0; i<Nsite; i++) {
			for(int k=0; k<Nsite; k++)
				H[i*Nsite+i]+=pow(double(k)-xc, 2.0)*F[i][k]*V;
		}
	}     

	free(temp);
}//end HamGen

//This looks ripe for a copy constructor...--K
void StateCopy(int* A, int* B, int N) {//Initialize B by copying A to B
	for(int i=0; i<N; i++)
		B[i]=A[i];
}

int StateId(int ** F, int * temp) {
	bool Flag; //should rename this variable whatever you're looking for and initialize to true --K
	int ID;
	//this seems like an inefficient form of a while loop
	//consider reimplementing as while(!Flag) --K
	for (ID=0; ID<Nsite; ID++) { //while((!Flag)&&(ID<Nsite))
		Flag = true; //won't need if you initialize it as true
		for (int i=0; i<Nsite; i++) {
			if (F[ID][i] != temp[i]) {
				Flag = false;
				break;
			}
		}
		if(Flag) break; //else ID++; (if you use while, delete the line currently here, unnecessary)
	}
	return ID;
}

void EnergChk(double * KEnerg1, double * KEnerg2, double * GSEnerg2, cdouble * rho, cdouble * Nk) {
	cdouble kin1=0.0, pot1=0.0, tot;
	cdouble kin2=0.0;
	double A=ratioA/ratioB;
	double xc=double(Nsite-1)/2.0+xdisp;

	//following E_k=-t*Sum{rho(i, i+1)+rho(i, i-1)}
	//			E_p= U*Sum{rho(i, i)}
	//			E  = E_k+E_p
	for(int i=0, j1, j2; i<Nsite; i++) {
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
	for(int k=-Nbzone*Nsite, i=0; k<(Nsite+1)/2&&k<Nbzone*Nsite; k++, i++) //count N[k] within one brilluion zone, [-Pi/2, Pi/2)
		if(k>=-Nsite/2)
			kin2+=-T*2.0*cos(2.0*Pi*k/Nsite)*Nk[i];

	*KEnerg1=real(kin1);
	*KEnerg2=real(kin2);
	*GSEnerg2=real(kin1+pot1);
}//end EnergChk


void output(cdouble * rho, cdouble * Nk, cdouble * Noise) {
	double ka, ki, kj; //ka is the momentum in terms of Pi, ka=2*k/Nsite
	ofstream op1, op2, op3, op4, op5;
	ostringstream densfile, momfile, noisefile;

	densfile<<"Output/DensTrap_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp<<"time"<<xtime;
	momfile<<"Output/MomTrap_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp<<"time"<<xtime;
	noisefile<<"Output/NoiseTrap_V"<<V<<"L"<<Nsite<<"N"<<Nfermion<<"X"<<xdisp<<"time"<<xtime;

	string s0=densfile.str();
	string s1=momfile.str();
	string s2=noisefile.str();

	char dfile[200], mfile[200], nfile[200];
	memset(dfile, '\0', 200); s0.copy(dfile, 100);
	memset(mfile, '\0', 200); s1.copy(mfile, 100);
	memset(nfile, '\0', 200); s2.copy(nfile, 100);

	//op1.open("Dens_Mx"); //'os::app' means start from last end
	op2.open(dfile);
	op3.open(mfile);
	//op4.open("Dens_Corr");
	op5.open(nfile);

/*	for(int i=0; 2*i<(Nsite+1); i++)
	op4<<i<<"\t"<<rho[i]<<endl;

	for(int i=0; i<Nsite; i++) {
		for(int j=0; j<Nsite; j++)
			op1<<rho[i*Nsite+j]<<"\t";
		op1<<endl;
	}
*/
	for(int i=0; i<Nsite; i++)
		op2<<xtime<<"\t"<<i<<"\t"<<real(rho[i*Nsite+i])<<endl;

	for(int k=-Nbzone*Nsite, i=0; k<=Nsite/2&&k<=Nbzone*Nsite; k++, i++) {//output N[k] within one brilluion zone, [-Pi/2, Pi/2]
		if(k>=-Nsite/2) {
			ka=2*double(k)/double(Nsite);
			op3<<xtime<<"\t"<<ka<<"\t"<<real(Nk[i])<<endl;
		}
	}
    
	//output the noise correlation within pre-set number of brilluion zones
	if(K==-1) {//full matrix of (k, k')
		for(int i=0, k=-Nbzone*Nsite; k<=Nbzone*Nsite; k++, i++) {
			for(int j=0, q=-Nbzone*Nsite; q<=Nbzone*Nsite; q++, j++) {
				ki=2*double(k)/double(Nsite);
				kj=2*double(q)/double(Nsite);
				op5<<ki<<"\t"<<kj<<"\t"<<real(Noise[i*(2*Nbzone*Nsite+1)+j])<<endl;
			}
		}
	}
	else {//fixed k (or k')
		int i=K+Nbzone*Nsite;
		for(int k=-Nbzone*Nsite, j=0; k<=Nbzone*Nsite; k++, j++) {
			ka=2*double(k)/double(Nsite);
			op5<<xtime<<"\t"<<ka<<"\t"<<real(Noise[i*(2*Nbzone*Nsite+1)+j])<<endl;
		}
	}
	//Close all the files we opened
	op1.close();
	op2.close();
	op3.close();
	op4.close();
	op5.close();
}//end output

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
}//end ChkOutput

/***MPI STUFF!***/

//Don't think these are ACTUALLY argc and argv
//you'll have to actually pass those as parameters to this
//when you call this in main)(
void MPI_initialize(int argc, char ** argv)
{
//From http://www.cse.msstate.edu/~ioana/Courses/CS6163/mpi/sld010.htm
//MPI::Init(int &argc, char **& argv)
    #ifdef MYMPI
        MPI::Init(argc, argv);
        irank=MPI::COMM_WORLD.Get_rank();
        isize=MPI::COMM_WORLD.Get_size();
    #else
        irank=0;
        isize=1;
    #endif
}

void MPI_finalize() {
	#ifdef MYMPI
	MPI::Finalize();
	#endif
}
