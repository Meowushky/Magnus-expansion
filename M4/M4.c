#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <stdint.h>

enum
{
  FLAVS=3
};

typedef struct
{
  double d0,d1,tol;
  double (*dens)(double);
  double (*H0)[FLAVS],(*W)[FLAVS],(*cH0W)[FLAVS],(*cH0H0W)[FLAVS],(*cWH0W)[FLAVS];  
} wf_ctx;

typedef struct
{
  double complex Psi[FLAVS];
  uint64_t calls;
  uint64_t step; //double(?) число изменений шагов
} rwf_ctx;

double ex(double);
void EigenV(double*, double, double);
void aWF_calc(wf_ctx*, rwf_ctx*);

void aWF_calc(wf_ctx *ctx, rwf_ctx *res)
{
  double h,x,Er,S1[FLAVS][FLAVS];
  x=ctx->d0;
  h=ctx->tol/2.;
  double xi_m,xi_p,f_m,f_p;
  double z,p,q,F;
  double L[2*FLAVS-1];
  double complex A[FLAVS][FLAVS],Eom4[FLAVS][FLAVS],psi_0[FLAVS],S2[FLAVS][FLAVS];
  double complex r0,r1;
  double One[FLAVS][FLAVS]=
  {
    {1., 0., 0.},
    {0., 1., 0.},
    {0., 0., 1.}
  };
  
  while(x<ctx->d1) 
  {
    if(x+h>ctx->d1) 
      h=ctx->d1-x;
    
    res->calls+=1;
    
    xi_m=x+(1.-1./sqrt(3.))*h/2.;
    xi_p=x+(1.+1./sqrt(3.))*h/2.;
    f_m=ctx->dens(xi_m);
    f_p=ctx->dens(xi_p);
    x+=h;
  
    z=0.;
    
    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
      {
        S1[j1][j2]=-sqrt(3.)*(f_p-f_m)*ctx->cH0W[j1][j2]/12.;
        S2[j1][j2]=I*sqrt(3.)*(f_p-f_m)*(ctx->cH0H0W[j1][j2]+(f_p+f_m)*ctx->cWH0W[j1][j2]/2.)/24.;
      }
    for(int j1=0;j1<FLAVS;j1++)
    {
      for(int j2=0;j2<FLAVS;j2++)
      {
        A[j1][j2]=-ctx->H0[j1][j2]-(f_p+f_m)*ctx->W[j1][j2]/2.+
        I*S1[j1][j2]*h;
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
    
    F=2.*sqrt(p/3.);

    r0=-(1.-cexp(I*L[3]*F*h))/L[3]; // r0/F
    r1=-(-r0-(1.-cexp(I*L[4]*F*h))/L[4])/(L[3]-L[4]);// r1/F^2

    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
      {
        Eom4[j1][j2]=cexp(I*h*z)*cexp(I*L[0]*F*h)*
        ((1.-L[0]*(r0-L[1]*r1))*One[j1][j2]+
        (r0+L[2]*r1)*A[j1][j2]/F+
        r1*(A[j1][0]*A[0][j2]+A[j1][1]*A[1][j2]+A[j1][2]*A[2][j2])/(F*F));
      }
    
    for(int j1=0;j1<FLAVS;j1++)
    {
      psi_0[j1]=Eom4[j1][0]*res->Psi[0];
      psi_0[j1]+=Eom4[j1][1]*res->Psi[1];
      psi_0[j1]+=Eom4[j1][2]*res->Psi[2];
    }

    for(int j1=0;j1<FLAVS;j1++)
      res->Psi[j1]=psi_0[j1];

    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
        for(int j3=0;j3<FLAVS;j3++)
          psi_0[j1]=(S1[j1][j2]+h*S2[j1][j2]+h*h*S1[j1][j3]*S1[j3][j2])*res->Psi[j2];
      
    Er=h*h*sqrt(psi_0[0]*psi_0[0]+psi_0[1]*psi_0[1]+psi_0[2]*psi_0[2]);
    
    if(Er>ctx->tol)
    {
      h=h*0.8*pow(ctx->tol/Er,1./3.);
      res->step+=1;
    }
  }
}

double ex(double t)
{
  double g=659560, n=10.54;
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
}

int main(int argc,char **argv)
{  
  double a=4351960.,b=0.030554,E=1.0,tol=0.001;
  double s12=sqrt(0.308),
    s13=sqrt(0.0234),
    c12=sqrt(1.-s12*s12),
    c13=sqrt(1.-s13*s13);
  double W[FLAVS][FLAVS]=
  {
    {c13*c13*c12*c12, c12*s12*c13*c13, c12*c13*s13},
    {c12*s12*c13*c13, s12*s12*c13*c13, s12*c13*s13},
    {c12*s13*c13, s12*c13*s13, s13*s13}
  };
  double d0=0.1,d1=1.;
  double Pee;
  
  if(argc==4)
  {
    d1=atof(argv[1]);
    tol=atof(argv[2]);
    E=atof(argv[3]);  
  }
  if((d1<d0)||(E<0))
  {
    printf("Error! E<0 or d1<d0.\n");
    return 1;
  }
  
  double H0[FLAVS][FLAVS]=
  {
    {0., 0., 0.},
    {0., a*b/E, 0.},
    {0., 0., a/E}
  };
  double cH0W[FLAVS][FLAVS],cH0H0W[FLAVS][FLAVS],cWH0W[FLAVS][FLAVS];
  
  //вычисление коммутатора H0 с W
  for(int j1=0;j1<FLAVS;j1++)
    for(int j2=0;j2<FLAVS;j2++)
  for(int j3=0;j3<FLAVS;j3++)
    cH0W[j1][j2]=H0[j1][j3]*W[j3][j2]-W[j1][j3]*H0[j3][j2];
      
  //вычисление коммутатора [H0,[H0,W]]
  for(int j1=0;j1<FLAVS;j1++)
    for(int j2=0;j2<FLAVS;j2++)
  for(int j3=0;j3<FLAVS;j3++)
    cH0H0W[j1][j2]=H0[j1][j3]*cH0W[j3][j2]-cH0W[j1][j3]*H0[j3][j2];
      
  //вычисление коммутатора [W,[H0,W]]
  for(int j1=0;j1<FLAVS;j1++)
    for(int j2=0;j2<FLAVS;j2++)
  for(int j3=0;j3<FLAVS;j3++)
    cWH0W[j1][j2]=W[j1][j3]*cH0W[j3][j2]-cH0W[j1][j3]*W[j3][j2];
  
  wf_ctx ctx;
  ctx.d0=d0;
  ctx.d1=d1;
  ctx.dens=ex;
  ctx.H0=H0;
  ctx.W=W;
  ctx.cH0W=cH0W;
  ctx.cH0H0W=cH0H0W;
  ctx.cWH0W=cWH0W;
  ctx.tol=tol;
  
  rwf_ctx res;
  res.Psi[0]=c12*c13;
  res.Psi[1]=s12*c13;
  res.Psi[2]=s13;
  res.calls=0;
  res.step=0;
  
  aWF_calc(&ctx,&res);
  
  Pee=c12*c12*c13*c13*res.Psi[0]*conj(res.Psi[0]);
  Pee+=s12*s12*c13*c13*res.Psi[1]*conj(res.Psi[1]);
  Pee+=s13*s13*res.Psi[2]*conj(res.Psi[2]);  
  
  printf("%lf\t%lf\t%ld\t%ld\n",d1,Pee,res.calls,res.step);
  
  /*
  fprintf(stderr,"%lf+I%lf; ",creal(Psi[0]),cimag(Psi[0]));
  fprintf(stderr,"%lf+I%lf; ",creal(Psi[1]),cimag(Psi[1]));
  fprintf(stderr,"%lf+I%lf; ",creal(Psi[2]),cimag(Psi[2]));
  
  double complex z;
  z=Psi[0]*conj(Psi[0])+Psi[1]*conj(Psi[1])+Psi[2]*conj(Psi[2]);
  fprintf(stderr,"%e+I%e\n",creal(z)-1.0,cimag(z));
  */
  
return 0;
}
