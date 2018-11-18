#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>

enum
{
  FLAVS=3
};

typedef struct
{
  int N;
  double d0,d1;
  double (*dens)(double);
  double (*H0)[FLAVS],(*W)[FLAVS];  
} wf_ctx;

double ex(double);
void EigenV(double*, double, double);
void aWF_calc(wf_ctx*, double complex*);

void aWF_calc(wf_ctx *ctx, double complex *Psi)
{
  int i;
  double h,x;
  x=ctx->d0;
  h=(ctx->d1-ctx->d0)/(ctx->N-1);
  double xi_m,xi_p,f_m,f_p;
  double z,p,q;
  double commut[FLAVS][FLAVS],L[2*FLAVS-1];
  double complex A[FLAVS][FLAVS],Eom4[FLAVS][FLAVS],psi_0[FLAVS];
  double complex r0,r1;
  double One[FLAVS][FLAVS]=
  {
    {1., 0., 0.},
    {0., 1., 0.},
    {0., 0., 1.}
  };
  
  for(i=0;i<ctx->N-1;i++)
  {
    xi_m=x+(1.-1./sqrt(3.))*h/2.;
    xi_p=x+(1.+1./sqrt(3.))*h/2.;
    f_m=ctx->dens(xi_m);
    f_p=ctx->dens(xi_p);
    x+=h;
  
    z=0.;
  
    //вычисление коммутатора H0 с W
    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
	for(int j3=0;j3<FLAVS;j3++)
	  commut[j1][j2]=ctx->H0[j1][j3]*ctx->W[j3][j2]-ctx->W[j1][j3]*ctx->H0[j3][j2];
  
    for(int j1=0;j1<FLAVS;j1++)
      {
	for(int j2=0;j2<FLAVS;j2++)
        { 
	   A[j1][j2]=-ctx->H0[j1][j2]-(f_p+f_m)*ctx->W[j1][j2]/2.-
	   I*sqrt(3.)*(f_p-f_m)*commut[j1][j2]*h/12.;
	}
	z+=A[j1][j1]/3.;
      }
    
    //обесшпуривание
    A[0][0]=A[0][0]-z;
    A[1][1]=A[1][1]-z;
    A[2][2]=A[2][2]-z;  
 
    //По Крамеру находит переменную q=Det
    q =A[0][0]*A[1][1]*A[2][2];
    q+=A[0][1]*A[1][2]*A[2][0];
    q+=A[0][2]*A[1][0]*A[2][1];
    q-=A[0][0]*A[1][2]*A[2][1];
    q-=A[0][1]*A[1][0]*A[2][2];
    q-=A[0][2]*A[1][1]*A[2][0];
    
    p=0.;

    for(int j1=0;j1<FLAVS;j1++)
      p+=(A[j1][0]*A[0][j1]+A[j1][1]*A[1][j1]+A[j1][2]*A[2][j1])/2.;

    EigenV(L,p,q);

    r0=-(1.-cexp(I*L[3]*h))/L[3];
    r1=-(-r0-(1.-cexp(I*L[4]*h))/L[4])/(L[3]-L[4]);

    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
	{
	  Eom4[j1][j2]=cexp(I*h*z)*cexp(I*L[0]*h)*
	  ((1.-L[0]*(r0-L[1]*r1))*One[j1][j2]+
	  (r0+L[2]*r1)*A[j1][j2]+
	  r1*(A[j1][0]*A[0][j2]+A[j1][1]*A[1][j2]+A[j1][2]*A[2][j2]));
    }
    
    for(int j1=0;j1<FLAVS;j1++)
    {
      psi_0[j1]=Eom4[j1][0]*Psi[0];
      psi_0[j1]+=Eom4[j1][1]*Psi[1];
      psi_0[j1]+=Eom4[j1][2]*Psi[2];
    }
    
    for(int j1=0;j1<FLAVS;j1++)
      Psi[j1]=psi_0[j1];  
  }
}

double ex(double t)
{
  double g=659560, n=10.54;//g=659560
  return g*exp(-n*t);
}

void EigenV(double *L,double p, double q)
{
  double arg,u;
  arg=acos((3.*q*sqrt(3.))/(2.*p*sqrt(p)))/3.;

  L[0]=cos(arg);
  L[1]=cos(arg-2.*M_PI/3.);

  if(L[1]<L[0])
  {
    u=L[0];
    L[0]=L[1];
    L[1]=u;
  }
  
  L[2]=cos(arg-4.*M_PI/3.);

  if(L[2]<L[0])
  {
    u=L[0];
    L[0]=L[2];
    L[2]=u;

    u=L[1];
    L[1]=L[2];
    L[2]=u;
  }
  else
    if(L[2]<L[1])
      {
        u=L[1];
        L[1]=L[2];
        L[2]=u;
      }
  L[3]=L[1]-L[0]; //a
  L[4]=L[2]-L[0]; //b
  u=2.*sqrt(p/3.);

  for(int i=0;i<2*FLAVS-1;i++)
    L[i]=u*L[i];
}

int main(int argc,char **argv)
{	
  double k=4351960.,b0=0.030554; //k=4351960
  double s12=sqrt(0.308),  
    s13=sqrt(0.0234),
    c12=sqrt(1.-s12*s12),
    c13=sqrt(1.-s13*s13);
  double H0[FLAVS][FLAVS]=
  {
    {0., 0., 0.},
    {0., k*b0, 0.},
    {0., 0., k}
  };
  double W[FLAVS][FLAVS]=
  {
    {c13*c13*c12*c12, c12*s12*c13*c13, c12*c13*s13},
    {c12*s12*c13*c13, s12*s12*c13*c13, s12*c13*s13},
    {c12*s13*c13, s12*c13*s13, s13*s13}
  };
  double complex Psi[FLAVS]=
  {
    c12*c13,
    s12*c13,
    s13
  };
  double d0=0.1,d1=1.;
  double Pee; //средняя вероятность выживания
  int N=10;
  
  if(argc==3)
  {
	d1=atof(argv[1]);
	N=atoi(argv[2]);  
  }
  if((N<0)||(d1<d0))
  {
    printf("Error! N<0 or d1<d0.\n");
    return 1;
  }
  
  wf_ctx ctx;
  ctx.N=N;
  ctx.d0=d0;
  ctx.d1=d1;
  ctx.dens=ex;
  ctx.H0=H0;
  ctx.W=W;
  
  aWF_calc(&ctx,Psi);
  
  Pee=c12*c12*c13*c13*Psi[0]*conj(Psi[0]);
  Pee+=s12*s12*c13*c13*Psi[1]*conj(Psi[1]);
  Pee+=s13*s13*Psi[2]*conj(Psi[2]);  
    
  printf("%lf\t%lf\t%d\n",d1,Pee,N);

return 0;
}
