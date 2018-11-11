#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>

enum
{
  FLAVS=3
};

double ex(double);
void EigenV(double*, double, double);

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

int main()
{
  FILE *f1=fopen("M4_2000.dat","w");
	
  double k=4351960.,b0=0.030554; //k=4351960
  double s12=sqrt(0.308),
    s13=sqrt(0.0234),
    c12=sqrt(1.-s12*s12),
    c13=sqrt(1.-s13*s13);
  double One[FLAVS][FLAVS]= //Единичная матрица
  {
    {1., 0., 0.},
    {0., 1., 0.},
    {0., 0., 1.}
  };
  double H0[FLAVS][FLAVS]=
  {
    {0., 0., 0.},
    {0., k*b0, 0.},
    {0., 0., k}
  };

  double commut[FLAVS][FLAVS],W[FLAVS][FLAVS]=
  {
    {c13*c13*c12*c12, c12*s12*c13*c13, c12*c13*s13},
    {c12*s12*c13*c13, s12*s12*c13*c13, s12*c13*s13},
    {c12*s13*c13, s12*c13*s13, s13*s13}
  };
  double L[2*FLAVS-1]; //L Собственные значения обесшпуренной матрицы
  double complex psi_0[FLAVS],Psi[FLAVS]=
  {
    c12*c13,
    s12*c13,
    s13
  };
  //все промежуточные значения и начальное
  double f_m,f_p;
  //d0, d1 - интервал интегрирования, P - сумма квадратов модулей пси
  double d0=0.1,d1=1.0,x,h,P;
  double p,q,z;
  //Для манипуляций с матричной экспонентой, A2 - квадрат матрицы А
  double complex Eom4[FLAVS][FLAVS];
  //Матрица эволюции(?) для 4 порядка (e^omega)
  double xi_m, xi_p; //Точки гаусс-лежандра
  double complex A[FLAVS][FLAVS]; //A2 - квадрат матрицы А
  double (*v)(double);
  double Pee; //средняя вероятность выживания
  
  v=ex;
  
  int N=2000,i;
  double complex r0,r1;
  
  h=(d1-d0)/(N-1);
      
  for(i=0,x=d0;i<N-1;i++)
  {
    xi_m=x+(1.-1./sqrt(3.))*h/2.;
    xi_p=x+(1.+1./sqrt(3.))*h/2.;
    f_m=v(xi_m);
    f_p=v(xi_p);
    x+=h;
  
    z=0.;
  
    //вычисление коммутатора H0 с W
    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
	for(int j3=0;j3<FLAVS;j3++)
	  commut[j1][j2]=H0[j1][j3]*W[j3][j2]-W[j1][j3]*H0[j3][j2];
  
    for(int j1=0;j1<FLAVS;j1++)
      {
	for(int j2=0;j2<FLAVS;j2++)
        { 
	   A[j1][j2]=-H0[j1][j2]-(f_p+f_m)*W[j1][j2]/2.-
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
    
    P=0.;
    
    for(int j1=0;j1<FLAVS;j1++)
    {
      psi_0[j1]=Eom4[j1][0]*Psi[0];
      psi_0[j1]+=Eom4[j1][1]*Psi[1];
      psi_0[j1]+=Eom4[j1][2]*Psi[2];
      P+=creal(psi_0[j1])*creal(psi_0[j1])+cimag(psi_0[j1])*cimag(psi_0[j1]);
    }
    
    for(int j1=0;j1<FLAVS;j1++)
      Psi[j1]=psi_0[j1];
        
   //printf("[%d] %lf \n",i,P); //выводит единицы
  
    Pee=c12*c12*c13*c13*Psi[0]*conj(Psi[0]);
    Pee+=s12*s12*c13*c13*Psi[1]*conj(Psi[1]);
    Pee+=s13*s13*Psi[2]*conj(Psi[2]);
    
    fprintf(f1,"%lf %lf\n",x,Pee);
  }
fclose(f1);

return 0;
}
