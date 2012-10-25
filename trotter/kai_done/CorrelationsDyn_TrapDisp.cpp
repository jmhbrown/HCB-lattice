/*	The two-point correlation functions and the four-point functions are implemented here	*/

#include "stdafxDyn.h"
#include <fstream>

using namespace std;

extern void TwoPntGreenMx(cdouble *);
extern cdouble FourPntGreen(int, int, int, int);
cdouble MomDistFunc(int, cdouble *);
void DenMx(cdouble *);
void MomCorr2(cdouble *, cdouble *);
void GreenToMomCorr(int, int, int, int, cdouble, cdouble *, cdouble *);
cdouble GreenToCorr(int, int, int, int, cdouble, cdouble *);
void MomCorrInc(int, int, int, int, cdouble *, cdouble);

/************************************************Function definitions************************************************************/
void DenMx(cdouble * corr_2)		//calculate and store <b+b-> in corr_2[]
{
	int flag=0;
	cdouble * green_2;
	green_2=(cdouble *)malloc(Nsite*Nsite*sizeof(cdouble));

	TwoPntGreenMx(green_2);

	for(int i=0; i<Nsite; i++)
		for(int j=i; j<Nsite; j++)
		{
			if(i==j)
				corr_2[i*Nsite+j]=1.0-green_2[i*Nsite+j];
			else
            {
                //notice the density matrix is not symmetric anymore
                //
                corr_2[i*Nsite+j]=green_2[j*Nsite+i];
                corr_2[j*Nsite+i]=green_2[i*Nsite+j];
            }
		}

	free(green_2);

}

cdouble MomDistFunc(int k, cdouble * rho)	//calulate and return <nk> at k
{
	cdouble nk=0.0;
	double theta;
	
	for(int i=0; i<Nsite; i++)
		for(int j=0; j<Nsite; j++)
		{
			theta=2*Pi*double(k*(i-j))/double(Nsite);
			//nk+=rho[i*Nsite+j]*cos(theta);
			nk+=rho[i*Nsite+j]*polar(1.0, theta);
		}

	nk=nk/double(Nsite);

	return nk;
}

void MomDen(cdouble * rho, cdouble * Nk, double * xCM, double * momMax, int * xMom) //rho=<b+b->, Nsite*Nsite; Nk=<nk>, 2*Nbzone*Nsite+1
{
	DenMx(rho);
		
	if(irank==0)
	{
		*xCM=0.0;
		for(int i=0; i<Nsite; i++)
			*xCM+=real(rho[i*Nsite+i])*i;
		*xCM=*xCM/double(Nfermion)-(double(Nsite)-1.0)/2.0;
//		*xCM=(*xCM-(double(Nsite)-1.0)/2.0)/double(Nfermion);

		*momMax=0.0;
		*xMom=0;
		for(int i=0, k=-Nbzone*Nsite; k<=Nbzone*Nsite; k++, i++)
		{
			Nk[i]=MomDistFunc(k, rho);
			if(k<=(Nsite+1)/2&&k>=-Nsite/2&&real(Nk[i])>=*momMax)
			{
				*momMax=real(Nk[i]);
				*xMom=k;
			}
		}
	}
}

void MomCorr2(cdouble * nk_nq_corr, cdouble * corr_2)
{
	int i, j, l, m;
	cdouble g_4;
	int count=0;
	int flag=0;

	cdouble * nk_nq;
	nk_nq=(cdouble *)malloc((2*Nbzone*Nsite+1)*(2*Nbzone*Nsite+1)*sizeof(cdouble));

	for(int i=0; i<2*Nbzone*Nsite+1; i++)
		for(int j=0; j<2*Nbzone*Nsite+1; j++)
			nk_nq[i*(2*Nbzone*Nsite+1)+j]=0.0;

	for(i=0; i<Nsite-1; i++)		//non-zero elements in the Four-Point Green's Matrix
		for(j=i+1; j<Nsite; j++)
			for(l=i; l<Nsite-1; l++)
			{
				if(l==i)
					m=j;
				else
					m=l+1;

				for(; m<Nsite; m++)
				{
					if(flag++%isize==irank)
					{
						g_4=FourPntGreen(i, j, l, m);

						if(i==l&&j==m)
						{
							GreenToMomCorr(i, j, l, m, g_4, nk_nq, corr_2);
							GreenToMomCorr(i, j, m, l, g_4, nk_nq, corr_2);
							GreenToMomCorr(j, i, l, m, g_4, nk_nq, corr_2);
							GreenToMomCorr(j, i, m, l, g_4, nk_nq, corr_2);
							count+=4;
						}
						else
						{
							GreenToMomCorr(i, j, l, m, g_4, nk_nq, corr_2);
							GreenToMomCorr(i, j, m, l, g_4, nk_nq, corr_2);
							GreenToMomCorr(j, i, l, m, g_4, nk_nq, corr_2);
							GreenToMomCorr(j, i, m, l, g_4, nk_nq, corr_2);
							GreenToMomCorr(m, l, i, j, conj(g_4), nk_nq, corr_2);
							GreenToMomCorr(m, l, j, i, conj(g_4), nk_nq, corr_2);
							GreenToMomCorr(l, m, i, j, conj(g_4), nk_nq, corr_2);
							GreenToMomCorr(l, m, j, i, conj(g_4), nk_nq, corr_2);
							count+=8;
						}
					}

				}
			}

	for(i=0; i<Nsite; i++)		//zero elements in Four-Point Green's matrix
		for(j=0; j<Nsite; j++)
			for(l=0; l<Nsite; l++)
				for(m=0; m<Nsite; m++)
					if(i==j||l==m)
						if(flag++%isize==irank)
						{
							GreenToMomCorr(i, j, l, m, 0.0, nk_nq, corr_2);
							count++;
						}

//	cout<<"rank"<<irank<<": the count of corr matrix is "<<count<<endl;

	#ifdef MYMPI

	if(K==-1)
		MPI::COMM_WORLD.Reduce(nk_nq, nk_nq_corr, (2*Nbzone*Nsite+1)*(2*Nbzone*Nsite+1), MPI::DOUBLE_COMPLEX, MPI::SUM, 0);
	else
	{
		int a=K+Nbzone*Nsite;
		MPI::COMM_WORLD.Reduce((nk_nq+a*(2*Nbzone*Nsite+1)), (nk_nq_corr+a*(2*Nbzone*Nsite+1)), 2*Nbzone*Nsite+1, MPI::DOUBLE_COMPLEX, MPI::SUM, 0);
	}
	#else
	//in non-parallel mode, simply copy nk_nq to nk_nq_corr by reference
	for(int i=0; i<(2*Nbzone*Nsite+1)*(2*Nbzone*Nsite+1); i++)
		nk_nq_corr[i]=nk_nq[i];
	#endif

	free(nk_nq);
}

void GreenToMomCorr(int a, int b, int c, int d, cdouble g_4, cdouble * nq_nk, cdouble * corr_2)
{
	cdouble corr_4;
	int i, j, l, m;
	i=c;
	j=a; 
	l=d;
	m=b;

	corr_4=GreenToCorr(i, j, l, m, g_4, corr_2);
	if(abs(corr_4)>=1.0e-15)
		MomCorrInc(i, j, l, m, nq_nk, corr_4);

}

cdouble GreenToCorr(int i, int j, int l, int m, cdouble green_4, cdouble * corr_2)
{
	cdouble corr_4;

	if(i==j&&l==m)
		if(i==m)
			corr_4=corr_2[i*Nsite+i];
		else
			corr_4=corr_2[i*Nsite+i]+corr_2[m*Nsite+m]+green_4-1.0;
	else if(i==j&&l!=m)
		if(i==m)
			corr_4=0.0;
		else
			corr_4=corr_2[l*Nsite+m]-green_4;
	else if(i!=j&&l==m)
		if(i==m)
			corr_4=0.0;
		else
			corr_4=corr_2[i*Nsite+j]-green_4;
	else
	{
		if(i!=m)
			corr_4=green_4;
		else if(j==l)
			corr_4=1.0-corr_2[l*Nsite+j]-green_4;
		else
			corr_4=corr_2[l*Nsite+j]-green_4;
	}

	if(BorS=='B')
		if(j==l)
			corr_4=corr_2[i*Nsite+m]*2.0-corr_4;

	return corr_4;
}

void MomCorrInc(int i, int j, int l, int m, cdouble * nk_nq, cdouble corr_4)
{
//	double real, imag;
	double theta;
	int a, b, k, q;

	if(K==-1)
		for(a=0, k=-Nbzone*Nsite; k<=Nbzone*Nsite; a++, k++)
			for(b=0,q=-Nbzone*Nsite; q<=Nbzone*Nsite; b++, q++)
			{
				theta=2*Pi*double(k*(i-j)+q*(l-m))/double(Nsite);
				//real=corr_4*cos(theta)/double(Nsite*Nsite);
				//imag=corr_4*sin(theta)/double(Nsite*Nsite);

				nk_nq[a*(2*Nbzone*Nsite+1)+b]+=corr_4*polar(1.0, theta)/double(Nsite*Nsite);
			}

	else
	{
		k=K;
		a=K+Nbzone*Nsite;
		for(b=0,q=-Nbzone*Nsite; q<=Nbzone*Nsite; b++, q++)
		{
			theta=2*Pi*double(k*(i-j)+q*(l-m))/double(Nsite);
			//real=corr_4*cos(theta)/double(Nsite*Nsite);
			//imag=corr_4*sin(theta)/double(Nsite*Nsite);

			nk_nq[a*(2*Nbzone*Nsite+1)+b]+=corr_4*polar(1.0, theta)/double(Nsite*Nsite);
		} 
	}
}

void NoiseCorr(cdouble * nk, cdouble * corr_2, cdouble * noise, double * noiseMax, int * xNoise) //noise(k, k'), 2*Nbzone*Nsite+1
{
	cdouble * nk_nq;
	nk_nq=(cdouble *)malloc((2*Nbzone*Nsite+1)*(2*Nbzone*Nsite+1)*sizeof(cdouble));

	for(int i=0; i<2*Nbzone*Nsite+1; i++)
		for(int j=0; j<2*Nbzone*Nsite+1; j++)
			nk_nq[i*(2*Nbzone*Nsite+1)+j]=0.0;

	MomCorr2(nk_nq, corr_2);

	if(irank==0)
	{
		if(K==-1)
			for(int i=0, k=-Nbzone*Nsite; k<=Nbzone*Nsite; i++, k++)
				for(int j=0, q=-Nbzone*Nsite; q<=Nbzone*Nsite; j++, q++)
					if((q-k)%Nsite==0&&q!=k)
						noise[i*(2*Nbzone*Nsite+1)+j]=nk_nq[i*(2*Nbzone*Nsite+1)+j]-nk[i]*(nk[j]+1.0);
					else
						noise[i*(2*Nbzone*Nsite+1)+j]=nk_nq[i*(2*Nbzone*Nsite+1)+j]-nk[i]*nk[j];
		else
		{
			int k=K;
			int i=K+Nbzone*Nsite;
			*noiseMax=0.0;
			*xNoise=0;
			for(int j=0, q=-Nbzone*Nsite; q<=Nbzone*Nsite; j++, q++)
			{
				if((q-k)%Nsite==0&&q!=k)
					noise[i*(2*Nbzone*Nsite+1)+j]=nk_nq[i*(2*Nbzone*Nsite+1)+j]-nk[i]*(nk[j]+1.0);
				else
					noise[i*(2*Nbzone*Nsite+1)+j]=nk_nq[i*(2*Nbzone*Nsite+1)+j]-nk[i]*nk[j];
				if(real(noise[i*(2*Nbzone*Nsite+1)+j])>=*noiseMax)
				{
					*noiseMax=real(noise[i*(2*Nbzone*Nsite+1)+j]);
					*xNoise=q;
				}
			}
		}
	}

	free(nk_nq);
}

/***********************************************************End******************************************************************/





